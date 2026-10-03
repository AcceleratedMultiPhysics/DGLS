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

// Single vortex (LeVeque / Rider-Kothe) with time reversal: domain [-0.5,0.5]^2
// (the classic unit square shifted), circle of radius 0.15 centered at (0,0.25).
// The flow reverses at t = T/2 and the exact solution at t = T is the initial one.

#define MPS_VORTEX_T 8.0
#define MPS_YC 0.25
#define MPS_R  0.15

#define mpsPrescribedVelocity2D(t, x, y, u, v)                                  \
{                                                                          \
  const dfloat X = (x) + 0.5;                                              \
  const dfloat Y = (y) + 0.5;                                              \
  const dfloat ct = cos(M_PI*(t)/MPS_VORTEX_T);                            \
  *(u) = -sin(M_PI*X)*sin(M_PI*X)*sin(2.0*M_PI*Y)*ct;                      \
  *(v) =  sin(M_PI*Y)*sin(M_PI*Y)*sin(2.0*M_PI*X)*ct;                      \
}

// exact only at t = 0 and t = T
#define mpsExactLevelSet2D(t, x, y, phi)                                   \
{                                                                          \
  *(phi) = sqrt((x)*(x) + ((y)-MPS_YC)*((y)-MPS_YC)) - MPS_R;              \
}

// walls have zero normal velocity, so no inflow data is needed
#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                                          \
  *(phiB) = phiM;                                                          \
}
