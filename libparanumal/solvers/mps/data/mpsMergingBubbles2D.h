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

// Merging bubbles: two bubbles of radius MPS_R on the vertical axis x = MPS_XC, the leading one at
// y = MPS_YC1 and the trailing one at y = MPS_YC2, in a liquid column [0,1.5]x[0,4] (phi < 0 inside the
// bubbles; phi = min of the two circle distances). rho 1000/100, mu 10/1, g 0.98, sigma 6.125
// (Eo = 40, Ga = 35 on the diameter 0.5). The trailing bubble is drawn into the wake of the leading one and
// the two coalesce. Walls top and bottom (bc 1), slip sides (bc 4). Initial pressure: hydrostatic in the
// liquid plus the Laplace jump inside each bubble; constants must match the setup file.

#define MPS_XC  0.75
#define MPS_YC1 1.2
#define MPS_YC2 0.6
#define MPS_R   0.25
#define MPS_RHO_LIQUID 1000.0
#define MPS_G 0.98
#define MPS_SIGMA 6.125
#define MPS_HTOP 4.0

// ---- flow (ins) macros ----
#define insInitialConditions2D(nu, t, x, y, u, v, p)                          \
{                                                                             \
  const dfloat r1 = sqrt(((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC1)*((y)-MPS_YC1)); \
  const dfloat r2 = sqrt(((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC2)*((y)-MPS_YC2)); \
  *(u) = 0.0;                                                                 \
  *(v) = 0.0;                                                                 \
  *(p) = MPS_RHO_LIQUID*MPS_G*(MPS_HTOP - (y))                                \
       + ((r1 < MPS_R || r2 < MPS_R) ? MPS_SIGMA/MPS_R : 0.0);                \
}

/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insVelocityDirichletConditions2D(bc, nu, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==3){                 \
    /* directional outflow: zero gradient where fluid leaves, the incoming normal velocity is removed where it enters (backflow) */ \
    const dfloat unB = (uM)*(nx) + (vM)*(ny);                      \
    *(uB) = (unB < 0.0) ? (uM) - unB*(nx) : (uM);                  \
    *(vB) = (unB < 0.0) ? (vM) - unB*(ny) : (vM);                  \
  } else if(bc==4){                 \
    *(uB) = 0.f;                    \
    *(vB) = vM;                     \
  } else if(bc==5){                 \
    *(uB) = uM;                     \
    *(vB) = 0.f;                    \
  }                                 \
}

#define insVelocityNeumannConditions2D(bc, nu, t, x, y, nx, ny, uxM, uyM, vxM, vyM, uxB, uyB, vxB, vyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  } else if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
  } else if(bc==4){                        \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(vxB) = vxM;                          \
    *(vyB) = vyM;                          \
  }                                        \
}

#define insPressureDirichletConditions2D(bc, nu, t, x, y, nx, ny, pM, pB) \
{                                   \
  if(bc==3){                        \
    *(pB) = 0.0;                    \
  } else {                          \
    *(pB) = pM;                     \
  }                                 \
}

#define insPressureNeumannConditions2D(bc, nu, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==3){                               \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
  } else {                                 \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  }                                        \
}

// ---- level set macros ----
#define mpsPrescribedVelocity2D(t, x, y, u, v) \
{                                              \
  *(u) = 0.0;                                  \
  *(v) = 0.0;                                  \
}

#define mpsExactLevelSet2D(t, x, y, phi)                                          \
{                                                                                 \
  const dfloat d1 = sqrt(((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC1)*((y)-MPS_YC1)) - MPS_R; \
  const dfloat d2 = sqrt(((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC2)*((y)-MPS_YC2)) - MPS_R; \
  *(phi) = (d1 < d2) ? d1 : d2;                                                   \
}

#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                                          \
  *(phiB) = phiM;                                                          \
}
