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

#include "ins.hpp"

// Divergence and continuity penalty step (Fehn, Wall & Kronbichler, JCP 2017/2018): given the
// projected velocity uhat, solve
//   (M + dt (A_D + A_C)) u = M uhat
// by conjugate gradients preconditioned with the inverse (lumped) mass matrix. The penalty
// parameters are frozen at uhat. The system is symmetric positive definite and dominated by
// the mass matrix, so few iterations are needed. Acts on the element-local velocity storage.
void ins_t::PenaltyStep(deviceMemory<dfloat>& o_U, const dfloat dt){

  if (!penaltyOn) return;

  const dlong Nlocal = mesh.Nelements*mesh.Np*NVfields;
  const dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
  linAlg_t& la = platform.linAlg();

  // penalty parameters from the projected velocity
  penaltyTauKernel(mesh.Nelements, mesh.o_wJ, penaltyDiv, o_U, o_penSpeed, o_penTauD);
  mesh.halo.Exchange(o_penSpeed, 1);

  deviceMemory<dfloat> o_b  = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_r  = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_z  = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_pp = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_Ap = platform.reserve<dfloat>(Ntotal);

  auto Ax = [&](deviceMemory<dfloat>& o_x, deviceMemory<dfloat>& o_Ax) {
    vTraceHalo.Exchange(o_x, 1);
    penaltyAxKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_sgeo, mesh.o_D, mesh.o_wJ,
                    mesh.o_vmapM, mesh.o_vmapP, mesh.o_EToB, dt, penaltyCont,
                    o_penTauD, o_penSpeed, o_x, o_Ax);
  };

  // b = M uhat, x0 = uhat, r = b - A x0
  penaltyMassKernel(mesh.Nelements, 0, mesh.o_wJ, o_U, o_b);
  Ax(o_U, o_Ap);
  la.axpy(Nlocal, (dfloat)1.0, o_b, (dfloat)0.0, o_r);
  la.axpy(Nlocal, (dfloat)-1.0, o_Ap, (dfloat)1.0, o_r);

  const dfloat bnorm = sqrt(la.innerProd(Nlocal, o_b, o_b, mesh.comm));
  const dfloat tol = 1e-12*std::max(bnorm, (dfloat)1e-30);

  penaltyMassKernel(mesh.Nelements, 1, mesh.o_wJ, o_r, o_z);
  la.axpy(Nlocal, (dfloat)1.0, o_z, (dfloat)0.0, o_pp);
  dfloat rz = la.innerProd(Nlocal, o_r, o_z, mesh.comm);

  NiterPen = 0;
  const int maxIter = 500;
  while (NiterPen < maxIter) {
    const dfloat rnorm = sqrt(la.innerProd(Nlocal, o_r, o_r, mesh.comm));
    if (rnorm <= tol) break;
    Ax(o_pp, o_Ap);
    const dfloat pAp = la.innerProd(Nlocal, o_pp, o_Ap, mesh.comm);
    const dfloat alpha = rz/pAp;
    la.axpy(Nlocal,  alpha, o_pp, (dfloat)1.0, o_U);
    la.axpy(Nlocal, -alpha, o_Ap, (dfloat)1.0, o_r);
    penaltyMassKernel(mesh.Nelements, 1, mesh.o_wJ, o_r, o_z);
    const dfloat rzNew = la.innerProd(Nlocal, o_r, o_z, mesh.comm);
    const dfloat beta = rzNew/rz;
    rz = rzNew;
    la.axpy(Nlocal, (dfloat)1.0, o_z, beta, o_pp);
    NiterPen++;
  }
}
