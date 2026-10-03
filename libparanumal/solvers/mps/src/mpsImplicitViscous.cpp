/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mps.hpp"

// out = H^{-1} rhs with H = gamma - nu0 Lap (the ins velocity Helmholtz operator), zero initial guess.
// Same steps as ins_t::VelocitySolve, with separate linear solver objects. The boundary data enter as
// in the main solve, so this is linear only for homogeneous velocity boundary data.
void mps_t::PrecVelocitySolve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_out,
                              const dfloat gamma, const dfloat T){

  const dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
  deviceMemory<dfloat> o_UH = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_VH = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_rhsU = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_rhsV = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_WH, o_rhsW;

  platform.linAlg().set(Ntotal*NVfields, (dfloat)0.0, o_out);
  const dfloat nu = flow.nu;
  flow.velocityRhsKernel(mesh.Nelements, mesh.o_wJ, mesh.o_vgeo, mesh.o_sgeo, mesh.o_ggeo, mesh.o_S, mesh.o_D,
                         mesh.o_LIFT, mesh.o_MM, mesh.o_sM, mesh.o_vmapM, mesh.o_EToB, mesh.o_mapB, flow.vTau, T,
                         mesh.o_x, mesh.o_y, mesh.o_z, gamma/nu, nu, o_out, o_rhs,
                         o_UH, o_VH, o_WH, o_rhsU, o_rhsV, o_rhsW);
  platform.linAlg().set(Ntotal, (dfloat)0.0, o_UH);
  platform.linAlg().set(Ntotal, (dfloat)0.0, o_VH);
  flow.uSolver.lambda = gamma/nu;
  flow.vSolver.lambda = gamma/nu;
  if (flow.vDisc_c0) {
    elliptic_t* S[2] = {&flow.uSolver, &flow.vSolver};
    linearSolver_t<dfloat>* L[2] = {&precULinearSolver, &precVLinearSolver};
    deviceMemory<dfloat>* R[2] = {&o_rhsU, &o_rhsV};
    deviceMemory<dfloat>* X[2] = {&o_UH, &o_VH};
    for (int c=0;c<2;++c) {
      deviceMemory<dfloat> o_G  = platform.reserve<dfloat>(S[c]->Ndofs+S[c]->Nhalo);
      deviceMemory<dfloat> o_GX = platform.reserve<dfloat>(S[c]->Ndofs+S[c]->Nhalo);
      S[c]->ogsMasked.Gather(o_G, *R[c], 1, ogs::Add, ogs::Trans);
      platform.linAlg().set(S[c]->Ndofs+S[c]->Nhalo, (dfloat)0.0, o_GX);
      S[c]->Solve(*L[c], o_GX, o_G, flow.velTOL, 5000, 0);
      S[c]->ogsMasked.Scatter(*X[c], o_GX, 1, ogs::NoTrans);
      o_GX.free(); o_G.free();
    }
  } else {
    flow.uSolver.Solve(precULinearSolver, o_UH, o_rhsU, flow.velTOL, 5000, 0);
    flow.vSolver.Solve(precVLinearSolver, o_VH, o_rhsV, flow.velTOL, 5000, 0);
  }
  flow.velocityBCKernel(mesh.Nelements, mesh.o_sgeo, mesh.o_vmapM, mesh.o_mapB, T,
                        mesh.o_x, mesh.o_y, mesh.o_z, nu, flow.vDisc_c0, o_UH, o_VH, o_WH, o_out);
}

// Implicit variable viscosity. The explicit step solved H U0 = f + V(u~), with H = gamma - nu0 Lap and the
// viscous remainder V (mpsViscousQuad2D.okl). The implicit step requires (H - V) U = f, i.e. the correction
// d = U - U0 solves (H - V) d = r0 = V(U0) - V(u~). Right preconditioned flexible GMRES with z = H^{-1} v:
// (H - V) z = v - V(z). At convergence all viscous terms are taken at t^{n+1} (BDF), with the properties of the step.
void mps_t::ImplicitViscousSolve(deviceMemory<dfloat>& o_Uc, const dfloat gamma, const dfloat T){

  const dlong Ntot = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
  const dlong NU = Ntot*NVfields;
  const dlong Nl = mesh.Nelements*mesh.Np*NVfields;   // owned entries (element major layout)
  const int m = viscousImplicitMaxIt;
  linAlg_t& la = platform.linAlg();
  auto V = [&](int i){ return o_gmV + i*NU; };
  auto Z = [&](int i){ return o_gmZ + i*NU; };

  // r0 = V(U0) - V(u~), with V(u~) in o_Visc
  deviceMemory<dfloat> v0 = V(0);
  ViscousRemainder(o_Uc, v0, T);
  la.axpy(NU, (dfloat)-1.0, o_Visc, (dfloat)1.0, v0);
  const dfloat beta = la.norm2(Nl, v0, mesh.comm);
  NiterVisc = 0; viscousImplicitRes = 0.0;
  if (!(beta > 0.0)) return;
  la.scale(NU, (dfloat)(1.0/beta), v0);

  std::vector<dfloat> H((m+1)*m, 0.0), cs(m, 0.0), sn(m, 0.0), g(m+1, 0.0);
  g[0] = beta;
  int k = 0;
  for (int j=0;j<m;++j) {
    deviceMemory<dfloat> zj = Z(j);
    deviceMemory<dfloat> vj = V(j);
    PrecVelocitySolve(vj, zj, gamma, T);
    // w = v_j - V(z_j)
    ViscousRemainder(zj, o_gmW, T);
    la.axpy(NU, (dfloat)1.0, V(j), (dfloat)-1.0, o_gmW);
    // modified Gram-Schmidt
    for (int i=0;i<=j;++i) {
      const dfloat hij = la.innerProd(Nl, o_gmW, V(i), mesh.comm);
      H[i*m+j] = hij;
      la.axpy(NU, -hij, V(i), (dfloat)1.0, o_gmW);
    }
    const dfloat hn = la.norm2(Nl, o_gmW, mesh.comm);
    H[(j+1)*m+j] = hn;
    // Givens rotations
    for (int i=0;i<j;++i) {
      const dfloat a = H[i*m+j], b = H[(i+1)*m+j];
      H[i*m+j] = cs[i]*a + sn[i]*b; H[(i+1)*m+j] = -sn[i]*a + cs[i]*b;
    }
    const dfloat a = H[j*m+j], b = H[(j+1)*m+j], r = sqrt(a*a + b*b);
    cs[j] = a/r; sn[j] = b/r;
    H[j*m+j] = r; H[(j+1)*m+j] = 0.0;
    g[j+1] = -sn[j]*g[j]; g[j] = cs[j]*g[j];
    k = j+1;
    viscousImplicitRes = fabs(g[j+1])/beta;
    if (viscousImplicitRes < viscousImplicitTol || hn <= 0.0) break;
    la.axpy(NU, (dfloat)(1.0/hn), o_gmW, (dfloat)0.0, V(j+1));
  }
  NiterVisc = k;

  // back substitution and U = U0 + sum y_i z_i
  std::vector<dfloat> y(k, 0.0);
  for (int i=k-1;i>=0;--i) {
    dfloat s = g[i];
    for (int l=i+1;l<k;++l) s -= H[i*m+l]*y[l];
    y[i] = s/H[i*m+i];
  }
  for (int i=0;i<k;++i) la.axpy(NU, y[i], Z(i), (dfloat)1.0, o_Uc);
}
