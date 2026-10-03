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

// smoothed Heaviside, equal to 1 inside (phi<0) and 0 outside
static dfloat Heaviside(const dfloat p, const dfloat eps){
  if (p < -eps) return 1.0;
  if (p >  eps) return 0.0;
  return 0.5*(1.0 - p/eps - sin(M_PI*p/eps)/M_PI);
}

void mps_t::LevelSetErrors(const dfloat T, dfloat& errL2, dfloat& errInf,
                           dfloat& errL1H, dfloat& area, dfloat& areaExact){

  dlong Nlocal = mesh.Nelements*mesh.Np;

  deviceMemory<dfloat> o_phiEx = platform.reserve<dfloat>(Nlocal);
  levelSetExactKernel(mesh.Nelements,
                      T,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      o_phiEx);

  memory<dfloat> phiEx(Nlocal);
  o_phiEx.copyTo(phiEx, Nlocal);
  o_phi.copyTo(phi, Nlocal);

  dfloat thickness = 1.0;
  settings.getSetting("INTERFACE THICKNESS", thickness);
  const dfloat eps = thickness*hmin/mesh.N;

  // phi errors are measured near the interface, on nodes with |phi_exact| <= 2h
  errL2 = 0.0; errInf = 0.0; errL1H = 0.0; area = 0.0; areaExact = 0.0;
  for (dlong n=0;n<Nlocal;++n) {
    const dfloat w  = mesh.wJ[n];
    const dfloat d  = phi[n] - phiEx[n];
    const dfloat H  = Heaviside(phi[n], eps);
    const dfloat He = Heaviside(phiEx[n], eps);

    if (fabs(phiEx[n]) <= 2.0*hmin) {
      errL2  += w*d*d;
      errInf  = std::max(errInf, fabs(d));
    }
    errL1H += w*fabs(H-He);
    area      += w*H;
    areaExact += w*He;
  }

  comm.Allreduce(errL2,  Comm::Sum);
  comm.Allreduce(errInf, Comm::Max);
  comm.Allreduce(errL1H, Comm::Sum);
  comm.Allreduce(area,   Comm::Sum);
  comm.Allreduce(areaExact, Comm::Sum);
  errL2 = sqrt(errL2);
}

void mps_t::Report(dfloat time, int tstep){

  static int frame=0;

  dfloat errL2, errInf, errL1H, area, areaExact;
  LevelSetErrors(time, errL2, errInf, errL1H, area, areaExact);

  hlong Nmarked = NumberOfMarkedElements();
  hlong Ntotal = mesh.Nelements;
  comm.Allreduce(Ntotal, Comm::Sum);

  hlong Nband = NbandElements;
  comm.Allreduce(Nband, Comm::Sum);

  if (twoPhase) {
    // flow diagnostics: max |u|, kinetic energy, pressure inside/outside (means over H)
    const dlong N = mesh.Nelements*mesh.Np;
    memory<dfloat> hu(N*NVfields), hp(N), hH(N), hrho(N);
    flow.o_u.copyTo(hu, N*NVfields);
    flow.o_p.copyTo(hp, N);
    o_Hf.copyTo(hH, N);
    o_rho.copyTo(hrho, N);
    memory<dfloat> hy(N);
    for (dlong n=0;n<N;++n) hy[n] = mesh.y[n];
    dfloat umax=0, ke=0, pin=0, win=0, pout=0, wout=0, Hint=0, yH=0, vH=0;
    dfloat exmin=1e30, exmax=-1e30, eymin=1e30, eymax=-1e30;
    for (dlong e=0;e<mesh.Nelements;++e)
      for (int n=0;n<mesh.Np;++n) {
        const dlong id = e*mesh.Np+n;
        const dfloat u = hu[e*mesh.Np*NVfields+n], v = hu[e*mesh.Np*NVfields+mesh.Np+n];
        const dfloat w = mesh.wJ[id];
        umax = std::max(umax, (dfloat)sqrt(u*u+v*v));
        ke += 0.5*w*hrho[id]*(u*u+v*v);
        Hint += w*hH[id]; yH += w*hH[id]*hy[id]; vH += w*hH[id]*v;
        if (hH[id] > 0.5) {
          exmin = std::min(exmin, mesh.x[id]); exmax = std::max(exmax, mesh.x[id]);
          eymin = std::min(eymin, hy[id]);     eymax = std::max(eymax, hy[id]);
        }
        if (hH[id] > 0.999) { pin += w*hp[id]; win += w; }
        if (hH[id] < 0.001) { pout += w*hp[id]; wout += w; }
      }
    comm.Allreduce(umax, Comm::Max);
    comm.Allreduce(ke, Comm::Sum);
    comm.Allreduce(pin, Comm::Sum);  comm.Allreduce(win, Comm::Sum);
    comm.Allreduce(pout, Comm::Sum); comm.Allreduce(wout, Comm::Sum);
    comm.Allreduce(Hint, Comm::Sum); comm.Allreduce(yH, Comm::Sum); comm.Allreduce(vH, Comm::Sum);
    comm.Allreduce(exmin, Comm::Min); comm.Allreduce(exmax, Comm::Max);
    comm.Allreduce(eymin, Comm::Min); comm.Allreduce(eymax, Comm::Max);
    if (mesh.rank==0 && levelSetMotion==2)
      printf("%10.5f (%6d): inside area %.8e, centroid y %.6f, rise velocity %.6f, extent x [%.5f,%.5f] y [%.5f,%.5f], LS substeps %d, reinits %d, rebuilds %d\n",
             time, tstep, Hint, Hint>0?yH/Hint:0.0, Hint>0?vH/Hint:0.0, exmin, exmax, eymin, eymax,
             NlsSubsteps, Nreinits, Nrebuilds);
    if (surfaceForceSubcell && sigma != 0.0) {
      memory<int> fl(mesh.Nelements); o_stFlag.copyTo(fl);
      hlong nf = 0; for (dlong e=0;e<mesh.Nelements;++e) nf += fl[e];
      comm.Allreduce(nf, Comm::Sum);
      if (mesh.rank==0) printf("%10.5f (%6d): subcell surface force in %lld elements\n", time, tstep, (long long)nf);
    }
    if (mesh.rank==0 && viscousImplicit)
      printf("%10.5f (%6d): implicit viscous FGMRES its %d, residual reduction %.1e\n", time, tstep, NiterVisc, viscousImplicitRes);
    if (mesh.rank==0)
      printf("%10.5f (%6d): max|u| %.3e, KE %.3e, <p>in %.6e, <p>out %.6e, dp %.6e, iters U %d V %d P %d, split passes %d (dp/p %.1e)\n",
             time, tstep, umax, ke, win>0?pin/win:0.0, wout>0?pout/wout:0.0,
             (win>0?pin/win:0.0)-(wout>0?pout/wout:0.0), NiterU, NiterV, NiterP,
             NpressureIters, pressureIterationChange);
  } else if(mesh.rank==0)
    printf("%8.4f (%6d): L2 err %.3e, L1(H) err %.3e, area ratio %.8f, band %lld/%lld, marked %lld, rebuilds %d, reinits %d, mass shift %.2e\n",
           time, tstep, errL2, errL1H, area/areaExact, (long long)Nband, (long long)Ntotal,
           (long long)Nmarked, Nrebuilds, Nreinits, lastMassShift);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_phi.copyTo(phi);
    if (!twoPhase) {
      SetVelocity(time);
      o_U.copyTo(U);
    }

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    o_decayRate.copyTo(decayRate);
    if (twoPhase) {
      const dlong Nt = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
      memory<dfloat> hp(Nt), hrho(Nt), hk(Nt);
      flow.o_u.copyTo(U);
      flow.o_p.copyTo(hp);
      o_rho.copyTo(hrho);
      o_kappa.copyTo(hk);
      PlotFields(phi, U, elementMark, decayRate, std::string(fname),
                 {{"Pressure", hp}, {"Density", hrho}, {"Curvature", hk}});
    } else
      PlotFields(phi, U, elementMark, decayRate, std::string(fname));

    // element data for detector calibration: centroid, decay exponent, mark, min|phi|
    sprintf(fname, "%s_%04d_%04d.elm", name.c_str(), mesh.rank, frame-1);
    FILE *fe = fopen(fname, "w");
    for (dlong e=0;e<mesh.Nelements;++e) {
      dfloat xc=0, yc=0, pmin=1e30;
      for (int n=0;n<mesh.Np;++n) {
        xc += mesh.x[e*mesh.Np+n]/mesh.Np;
        yc += mesh.y[e*mesh.Np+n]/mesh.Np;
        pmin = std::min(pmin, (dfloat)fabs(phi[e*mesh.Np+n]));
      }
      int cut = 0; dfloat pmn=1e30, pmx=-1e30;
      for (int n=0;n<mesh.Np;++n) { pmn=std::min(pmn,phi[e*mesh.Np+n]); pmx=std::max(pmx,phi[e*mesh.Np+n]); }
      cut = (pmn<=0 && pmx>=0);
      fprintf(fe, "%g %g %g %d %g %d\n", xc, yc, decayRate[e], elementMark[e], pmin, cut);
    }
    fclose(fe);
  }
}
