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

// Smooth rigid-body rotation test for convergence rates.
// Domain [-1,1]^2, one revolution per unit time, exact solution known for all t.

#define MPS_OMEGA (2.0*M_PI)
#define MPS_X0 0.4
#define MPS_Y0 0.0
#define MPS_SIGMA 0.3

dfloat mpsRotationPhi0(const dfloat x, const dfloat y){
  const dfloat r2 = (x-MPS_X0)*(x-MPS_X0) + (y-MPS_Y0)*(y-MPS_Y0);
  return 0.5 - exp(-r2/(MPS_SIGMA*MPS_SIGMA));
}

#define mpsPrescribedVelocity2D(t, x, y, u, v) \
{                                         \
  *(u) = -MPS_OMEGA*(y);                  \
  *(v) =  MPS_OMEGA*(x);                  \
}

// exact solution: rotate coordinates back by -omega*t
#define mpsExactLevelSet2D(t, x, y, phi)          \
{                                                 \
  const dfloat ct = cos(MPS_OMEGA*(t));           \
  const dfloat st = sin(MPS_OMEGA*(t));           \
  const dfloat xr =  ct*(x) + st*(y);             \
  const dfloat yr = -st*(x) + ct*(y);             \
  *(phi) = mpsRotationPhi0(xr, yr);               \
}

// inflow data from the exact solution
#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                 \
  mpsExactLevelSet2D(t, x, y, phiB);              \
}
