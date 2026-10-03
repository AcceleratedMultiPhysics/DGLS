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

dfloat mps_t::GradientDeviation(){

  if (NbandElements==0) return 0.0;

  gradientDeviationKernel(NbandElements, o_bandElements, o_elementLayer,
                          mesh.o_vgeo, mesh.o_wJ, mesh.o_D, o_phi, o_gradDev);

  memory<dfloat> dev(2*NbandElements);
  o_gradDev.copyTo(dev, 2*NbandElements);

  dfloat sd = 0.0, sw = 0.0;
  for (dlong l=0;l<NbandElements;++l) { sd += dev[2*l]; sw += dev[2*l+1]; }
  comm.Allreduce(sd, Comm::Sum);
  comm.Allreduce(sw, Comm::Sum);
  return (sw > 0.0) ? sd/sw : 0.0;
}

void mps_t::MassIntegrals(const dfloat s, dfloat& V, dfloat& Lint){

  dfloat thickness = 1.0;
  settings.getSetting("INTERFACE THICKNESS", thickness);
  const dfloat eps = thickness*hmin/mesh.N;

  massIntegralsKernel(mesh.Nelements, eps, s, o_elementLayer, mesh.o_wJ, o_phi, o_massVL);

  memory<dfloat> VL(2*mesh.Nelements);
  o_massVL.copyTo(VL);

  V = 0.0; Lint = 0.0;
  for (dlong e=0;e<mesh.Nelements;++e) { V += VL[2*e]; Lint += VL[2*e+1]; }
  comm.Allreduce(V, Comm::Sum);
  comm.Allreduce(Lint, Comm::Sum);
}

// Global mass correction (Lee, Dolbow & Mucha 2014): find the constant shift s of phi on
// the band with V(phi + s) = massTarget by Newton iterations, V'(s) = -L(s) where L is
// the smoothed interface length, and apply it.
void mps_t::MassCorrection(){

  dfloat s = 0.0;
  for (int it=0;it<5;++it) {
    dfloat V, Lint;
    MassIntegrals(s, V, Lint);
    if (Lint <= 0.0) break;
    const dfloat ds = (V - massTarget)/Lint;
    s += ds;
    if (fabs(ds) < 1e-12*hmin) break;
  }

  if (s != 0.0)
    massShiftKernel(NbandElements, o_bandElements, s, o_phi);

  lastMassShift = s;
}
