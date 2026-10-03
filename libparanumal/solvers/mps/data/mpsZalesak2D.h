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

// Zalesak's slotted disk: domain [-50,50]^2 (the classic [0,100]^2 shifted),
// disk centered at (0,25) with radius 15, slot width 5 reaching y = 35,
// one revolution in t = 628. Signed distance is negative inside the disk.

#define MPS_OMEGA (2.0*M_PI/628.0)
#define MPS_XC 0.0
#define MPS_YC 25.0
#define MPS_R  15.0
#define MPS_SLOT_HW 2.5
#define MPS_SLOT_TOP 35.0

// signed distance to an axis-aligned box with center (cx,cy) and half widths (hx,hy)
dfloat mpsBoxDistance(const dfloat x, const dfloat y,
                      const dfloat cx, const dfloat cy,
                      const dfloat hx, const dfloat hy){
  const dfloat qx = fabs(x-cx) - hx;
  const dfloat qy = fabs(y-cy) - hy;
  const dfloat ox = (qx > 0.0) ? qx : 0.0;
  const dfloat oy = (qy > 0.0) ? qy : 0.0;
  const dfloat outside = sqrt(ox*ox + oy*oy);
  const dfloat mq = (qx > qy) ? qx : qy;
  const dfloat inside = (mq < 0.0) ? mq : 0.0;
  return outside + inside;
}

dfloat mpsZalesakPhi0(const dfloat x, const dfloat y){
  const dfloat dDisk = sqrt((x-MPS_XC)*(x-MPS_XC) + (y-MPS_YC)*(y-MPS_YC)) - MPS_R;
  // slot: x in [-hw,hw], y from below the disk up to MPS_SLOT_TOP
  const dfloat ybot = MPS_YC - MPS_R - 5.0;
  const dfloat dSlot = mpsBoxDistance(x, y, MPS_XC, 0.5*(ybot+MPS_SLOT_TOP),
                                      MPS_SLOT_HW, 0.5*(MPS_SLOT_TOP-ybot));
  // disk minus slot
  return (dDisk > -dSlot) ? dDisk : -dSlot;
}

#define mpsPrescribedVelocity2D(t, x, y, u, v) \
{                                         \
  *(u) = -MPS_OMEGA*(y);                  \
  *(v) =  MPS_OMEGA*(x);                  \
}

#define mpsExactLevelSet2D(t, x, y, phi)          \
{                                                 \
  const dfloat ct = cos(MPS_OMEGA*(t));           \
  const dfloat st = sin(MPS_OMEGA*(t));           \
  const dfloat xr =  ct*(x) + st*(y);             \
  const dfloat yr = -st*(x) + ct*(y);             \
  *(phi) = mpsZalesakPhi0(xr, yr);                \
}

#define mpsLevelSetDirichletConditions2D(bc, t, x, y, nx, ny, phiM, phiB) \
{                                                 \
  mpsExactLevelSet2D(t, x, y, phiB);              \
}
