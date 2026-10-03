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

// Reinitialization test of Karakus, Chalmers & Warburton (CAMWA 2022, eq. 22):
// perturbed signed distance to the unit circle on [-2,2]^2, no flow.
// The exact signed distance has a kink at the origin.

#define mpsPrescribedVelocity2D(t, x, y, u, v) \
{                                              \
  *(u) = 0.0;                                  \
  *(v) = 0.0;                                  \
}

#define mpsInitialLevelSet2D(t, x, y, phi)                                  \
{                                                                           \
  *(phi) = (((x)-1.0)*((x)-1.0) + ((y)-1.0)*((y)-1.0) + 0.1)                \
          *(sqrt((x)*(x) + (y)*(y)) - 1.0);                                 \
}

#define mpsExactLevelSet2D(t, x, y, phi)                                    \
{                                                                           \
  *(phi) = sqrt((x)*(x) + (y)*(y)) - 1.0;                                   \
}

#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB)  \
{                                                                           \
  *(phiB) = phiM;                                                           \
}
