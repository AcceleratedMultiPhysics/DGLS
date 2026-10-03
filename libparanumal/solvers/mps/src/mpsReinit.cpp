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

// Reinitialize phi on the band: evolve u_t + |grad u| = 0 from u = phi and
// v_t + |grad v| = 0 from v = -phi in pseudo time up to reinitDistance, and set phi to
// the signed first arrival time of the zero level (CAMWA 2022).
void mps_t::Reinitialize(){

  const dlong Ntotal2 = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*2;

  // pseudo time step (unit characteristic speed)
  dfloat dtau = reinitCFL*hmin/((mesh.N+1.)*(mesh.N+1.));
  int K = std::max(4, (int)ceil(reinitDistance/dtau));
  dtau = reinitDistance/K;

  o_phi0.copyFrom(o_phi);

  // contact points on the walls, where the contact angle ghost condition applies
  if (o_contactPts.length() == 0) {
    memory<dfloat> big(2*(mesh.Nelements+mesh.totalHaloPairs), (dfloat)1e30);
    o_contactPts = platform.malloc<dfloat>(big);
  }
  if (reinitGhostWall) {
    platform.linAlg().set(2*(mesh.Nelements+mesh.totalHaloPairs), (dfloat)1e30, o_contactPts);   // no stale points
    wallContactPointsKernel(NbandElements, o_bandElements, mesh.o_EToB, mesh.o_vmapM, mesh.o_x, mesh.o_y,
                            o_phi0, o_contactPts);
  }

  reinitInitKernel(NbandElements, o_bandElements, o_phi, o_Q);
  reinitClearFoundKernel((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np, o_found);
  reinitSaveHistoryKernel(NbandElements, o_bandElements, Ntotal2, 0, o_Q, o_Qhist);

  for (int k=1;k<=K;++k) {
    if (reinitMarkMode==2) {
      reinitTraceHalo.Exchange(o_Q, 1);
      levelSetDetectorKernel(NbandElements, o_bandElements, mesh.o_vgeo, o_invV1D,
                             reinitDetectorOn, reinitDetectorOff, (dlong)(2*mesh.Np), (dlong)0,
                             o_Q, o_reinitDecay, o_reinitMark);
      levelSetDetectorKernel(NbandElements, o_bandElements, mesh.o_vgeo, o_invV1D,
                             reinitDetectorOn, reinitDetectorOff, (dlong)(2*mesh.Np), (dlong)mesh.Np,
                             o_Q, o_reinitDecay, o_reinitMarkV);
    }

    // previous step values and element minima for the Hopf-Lax bounds
    o_Qprev.copyFrom(o_Q);
    reinitElementMinKernel(NbandElements, o_bandElements, o_Q, o_Qmin);

    for (int rk=0;rk<5;++rk) {
      reinitTraceHalo.Exchange(o_Q, 1);

      reinitRhsKernel(NbandElements, o_bandElements,
                      mesh.o_vgeo, mesh.o_sgeo, mesh.o_D, o_gllw,
                      mesh.o_vmapP, mesh.o_EToB, o_elementLayer,
                      o_reinitMark, o_reinitMarkV,
                      o_Q, reinitGhostWall, cosContactAngle, (dfloat)(2.0*hmin), o_phi0, mesh.Nelements, o_EToE, mesh.o_x, mesh.o_y, o_contactPts, o_Qrhs);

      reinitUpdateKernel(NbandElements, o_bandElements, dtau, rka[rk], rkb[rk],
                         o_Qrhs, o_Qres, o_Q);
    }

    reinitBoundKernel(NbandElements, o_bandElements, mesh.Nelements, o_EToE, o_elementLayer,
                      o_Qmin, o_Qprev, o_Q);

    reinitSaveHistoryKernel(NbandElements, o_bandElements, Ntotal2, k%Nhist, o_Q, o_Qhist);

    if (k>=2)
      reinitArrivalKernel(NbandElements, o_bandElements, Ntotal2, k, 0, dtau,
                          o_phi0, o_Qhist, o_found, o_arrival);
  }
  reinitArrivalKernel(NbandElements, o_bandElements, Ntotal2, K, 1, dtau,
                      o_phi0, o_Qhist, o_found, o_arrival);

  reinitFinalizeKernel(NbandElements, o_bandElements, reinitDistance,
                       o_found, o_arrival, o_phi);

  // keep the interface where it was
  if (reinitPreserve && bandEnabled)
    reinitPreserveKernel(NbandElements, o_bandElements, reinitPreserve, o_elementLayer,
                         mesh.o_vgeo, mesh.o_D, o_phi0, o_phi);

  Nreinits++;
}
