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

// Bubble sliding under an inclined wall: Hysing et al. (2009) test case 1 fluids (rho 1000/100, mu 10/1,
// sigma 24.5, |g| = 0.98, bubble radius 0.25) in a channel [0,5]x[0,1.5] whose top wall is the inclined plate.
// The inclination alpha of the plate to the horizontal is realized by tilting gravity:
//   g = |g| (-sin(alpha), -cos(alpha)),   alpha = 30 degrees: g = (-0.49, -0.848705).
// The bubble (phi < 0 inside) starts at (0.75, 0.6), rises to the plate and slides along it.
// Walls on all sides (bc 1). Initial pressure: hydrostatic for the tilted gravity plus the Laplace jump;
// the constants must match the setup file (GRAVITY X, GRAVITY Y, SURFACE TENSION).

#define MPS_XC 0.75
#define MPS_YC 0.6
#define MPS_R  0.25
#define MPS_RHO_LIQUID 1000.0
#define MPS_GX (-0.49)
#define MPS_GY (-0.848705)
#define MPS_XMAX 5.0
#define MPS_YMAX 1.5
#define MPS_SIGMA 24.5

// ---- flow (ins) macros ----
#define insInitialConditions2D(nu, t, x, y, u, v, p)                          \
{                                                                             \
  *(u) = 0.0;                                                                 \
  *(v) = 0.0;                                                                 \
  *(p) = MPS_RHO_LIQUID*(MPS_GX*((x)-MPS_XMAX) + MPS_GY*((y)-MPS_YMAX))       \
       + ((((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC)*((y)-MPS_YC) < MPS_R*MPS_R)   \
          ? MPS_SIGMA/MPS_R : 0.0);                                           \
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

#define mpsExactLevelSet2D(t, x, y, phi)                                   \
{                                                                          \
  *(phi) = sqrt(((x)-MPS_XC)*((x)-MPS_XC) + ((y)-MPS_YC)*((y)-MPS_YC)) - MPS_R; \
}

#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                                          \
  *(phiB) = phiM;                                                          \
}
