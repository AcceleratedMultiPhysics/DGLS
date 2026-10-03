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

// Hydrostatic balance of two fluids at rest: interface y = 0 (phi = y), the fluid
// inside (phi < 0, below) set by DENSITY INSIDE, gravity from GRAVITY Y.
// Walls on all sides. Exact solution: u = 0, grad p = rho g.

// Initial pressure: hydrostatic with the sharp densities, p = -rho g y (continuous at
// y = 0). These constants must match DENSITY INSIDE/OUTSIDE and GRAVITY Y in the setup;
// set MPS_HYDRO_P0 to 0 to start from p = 0 (impulsive start).
#ifndef MPS_HYDRO_P0
#define MPS_HYDRO_P0 1
#endif
#define MPS_HYDRO_RHO_IN  1000.0
#define MPS_HYDRO_RHO_OUT 1.0
#define MPS_HYDRO_G      (-1.0)

// ---- flow (ins) macros ----
#define insInitialConditions2D(nu, t, x, y, u, v, p)                          \
{                                                                             \
  *(u) = 0.0;                                                                 \
  *(v) = 0.0;                                                                 \
  *(p) = MPS_HYDRO_P0*((y) < 0.0 ? MPS_HYDRO_RHO_IN : MPS_HYDRO_RHO_OUT)      \
         *MPS_HYDRO_G*(y);                                                    \
}

/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insVelocityDirichletConditions2D(bc, nu, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
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

#define mpsExactLevelSet2D(t, x, y, phi) \
{                                        \
  *(phi) = (y);                          \
}

#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                                          \
  *(phiB) = phiM;                                                          \
}
