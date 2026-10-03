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

// fill o_U with the prescribed velocity field at time T, including halo
void mps_t::SetVelocity(const dfloat T){

  velocityFieldKernel(mesh.Nelements,
                      T,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      o_U);

  vTraceHalo.Exchange(o_U, 1);
}

dfloat mps_t::MaxVelocity(){

  dlong N = mesh.Nelements*mesh.Np;
  o_U.copyTo(U);

  dfloat vmax = 0.0;
  for (dlong e=0;e<mesh.Nelements;++e) {
    for (int n=0;n<mesh.Np;++n) {
      dfloat v2 = 0.0;
      for (int fld=0;fld<NVfields;++fld) {
        const dfloat u = U[e*mesh.Np*NVfields + fld*mesh.Np + n];
        v2 += u*u;
      }
      vmax = std::max(vmax, sqrt(v2));
    }
  }
  comm.Allreduce(vmax, Comm::Max);
  (void) N;
  return vmax;
}

//evaluate level set rhs = -u.grad(phi) with upwind surface flux, on the band elements
void mps_t::rhsf(deviceMemory<dfloat>& o_Phi, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // velocity at stage time: prescribed, or the flow velocity extrapolated linearly from
  // U^n and U^{n-1} over the current flow step; optionally tapered near the band edge
  if (twoPhase && levelSetMotion==2) {
    const dlong NU = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
    const dfloat theta = (lsFlowDt > 0.0) ? (T - lsFlowT0)/lsFlowDt : 0.0;
    deviceMemory<dfloat> U0 = o_Uhist + ((flowShift+0)%3)*NU;
    deviceMemory<dfloat> U1 = o_Uhist + ((flowShift+1)%3)*NU;
    const dfloat c1 = (flowStep > 0) ? -theta : 0.0;
    combine3Kernel(NU, 1.0 - c1, U0, c1, U1, 0.0, U1, o_U);
  } else {
    velocityFieldKernel(mesh.Nelements,
                        T,
                        mesh.o_x,
                        mesh.o_y,
                        mesh.o_z,
                        o_U);
  }
  if (bandEnabled && bandTaper)
    bandTaperKernel(NbandElements, o_bandElements, taperStart, taperEnd,
                    bandLayers, o_elementLayer, o_Phi, o_U);
  vTraceHalo.Exchange(o_U, 1);

  // a priori marking at every stage
  MarkElements(o_Phi);

  // phi on neighbor faces is needed for the subcell slopes
  phiTraceHalo.Exchange(o_Phi, 1);

  // reconstructed subcell edge states (nodal copy in unmarked elements)
  levelSetReconstructKernel(NbandElements,
                            o_bandElements,
                            o_gllw,
                            o_elementMark,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            o_elementLayer,
                            o_Phi,
                            o_phiEdge);

  edgeTraceHalo.ExchangeStart(o_phiEdge, 1);

  // DG volume term in unmarked elements, subcell FV interior fluctuations in marked ones.
  // The element face term is the upwind jump of edge states and is added by the surface kernel.
  levelSetVolumeKernel(NbandElements,
                       o_bandElements,
                       mesh.o_vgeo,
                       mesh.o_D,
                       o_gllw,
                       o_elementMark,
                       o_U,
                       o_Phi,
                       o_phiEdge,
                       o_RHS);

  edgeTraceHalo.ExchangeFinish(o_phiEdge, 1);

  levelSetSurfaceKernel(NbandElements,
                        o_bandElements,
                        mesh.o_sgeo,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        mesh.o_EToB,
                        o_EToF,
                        o_elementLayer,
                        T,
                        mesh.o_x,
                        mesh.o_y,
                        mesh.o_z,
                        o_U,
                        o_Phi,
                        o_phiEdge,
                        o_RHS);
}

// one low storage RK4 step, updating only the band elements
void mps_t::Step(const dfloat time, const dfloat _dt){

  for(int rk=0;rk<5;++rk){
    const dfloat stageTime = time + rkc[rk]*_dt;

    rhsf(o_phi, o_rhsPhi, stageTime);

    levelSetUpdateKernel(NbandElements, o_bandElements, _dt, rka[rk], rkb[rk],
                         o_rhsPhi, o_resPhi, o_phi);
  }
}

void mps_t::MarkElements(deviceMemory<dfloat>& o_Phi){

  if (subcellMode!=2) return; // marks are fixed

  levelSetDetectorKernel(NbandElements,
                         o_bandElements,
                         mesh.o_vgeo,
                         o_invV1D,
                         detectorOn,
                         detectorOff,
                         (dlong)mesh.Np,
                         (dlong)0,
                         o_Phi,
                         o_decayRate,
                         o_elementMark);
}

hlong mps_t::NumberOfMarkedElements(){

  o_elementMark.copyTo(elementMark);

  hlong cnt = 0;
  for (dlong l=0;l<NbandElements;++l) cnt += elementMark[bandElements[l]];
  comm.Allreduce(cnt, Comm::Sum);
  return cnt;
}
