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

// Extend the polynomial of element src into element e: project src onto tensor Legendre
// modes of degree <= p and evaluate them at the nodes of e, mapped into the reference
// coordinates of src. Assumes affine (parallelogram) quads with the standard vertex order.
// The result keeps the sign of e's current (clipped) value and is clamped to
// lo <= |phi| <= clip: an element entering layer J lies at least ~(J-1)h from the interface,
// so the extension must not create a zero crossing there.
static void ExtendPolynomial(const mesh_t& mesh, const memory<dfloat>& invV,
                             memory<dfloat>& phi, const dlong src, const dlong e,
                             const int p, const dfloat lo, const dfloat clip){

  const int Nq = mesh.Nq, Np = mesh.Np, Nv = mesh.Nverts;

  // modal coefficients c[n][m] (m along r, n along s)
  memory<dfloat> t(Np), c(Np);
  for (int j=0;j<Nq;++j)
    for (int m=0;m<Nq;++m) {
      dfloat r = 0.0;
      for (int i=0;i<Nq;++i) r += invV[m*Nq+i]*phi[src*Np + j*Nq + i];
      t[j*Nq+m] = r;
    }
  for (int n=0;n<Nq;++n)
    for (int m=0;m<Nq;++m) {
      dfloat r = 0.0;
      for (int j=0;j<Nq;++j) r += invV[n*Nq+j]*t[j*Nq+m];
      c[n*Nq+m] = r;
    }

  // affine map of src: x = x0 + (r+1)/2 (x1-x0) + (s+1)/2 (x3-x0)
  const dfloat x0 = mesh.EX[src*Nv+0], y0 = mesh.EY[src*Nv+0];
  const dfloat ax = 0.5*(mesh.EX[src*Nv+1]-x0), ay = 0.5*(mesh.EY[src*Nv+1]-y0);
  const dfloat bx = 0.5*(mesh.EX[src*Nv+3]-x0), by = 0.5*(mesh.EY[src*Nv+3]-y0);
  const dfloat det = ax*by - ay*bx;

  dfloat mean = 0.0;
  for (int k=0;k<Np;++k) mean += phi[e*Np+k];
  const dfloat sgn = (mean >= 0.0) ? 1.0 : -1.0;

  memory<dfloat> Pr(Nq), Ps(Nq);
  for (int k=0;k<Np;++k) {
    const dfloat dx = mesh.x[e*Np+k]-x0, dy = mesh.y[e*Np+k]-y0;
    const dfloat r = ( by*dx - bx*dy)/det - 1.0;
    const dfloat s = (-ay*dx + ax*dy)/det - 1.0;
    for (int m=0;m<=p;++m) {
      mesh_t::OrthonormalBasis1D(r, m, Pr[m]);
      mesh_t::OrthonormalBasis1D(s, m, Ps[m]);
    }
    dfloat v = 0.0;
    for (int n=0;n<=p;++n)
      for (int m=0;m<=p;++m)
        v += c[n*Nq+m]*Pr[m]*Ps[n];
    phi[e*Np+k] = sgn*std::max(lo, std::min(clip, sgn*v));
  }
}

// Build the band of element layers around the cut elements (L0) by breadth first
// search over vertex neighbors (elements sharing a node), so J layers reach about J*h
// from the interface in every direction. Elements that leave the band are clipped to +/-bandClip.
// Elements that enter the band are initialized with the polynomial extension (degree
// extensionDegree) of an adjacent local element that was already in the band, clamped
// to +/-bandClip. Only local elements serve as sources (host coordinates of halo
// elements are not available); an element without a source keeps its clipped value.
void mps_t::BandBuild(const bool initial){

  const dlong Nelements = mesh.Nelements;
  const dlong Ntotal = Nelements + mesh.totalHaloPairs;
  const int Np = mesh.Np;
  const int Nfaces = mesh.Nfaces;

  if (!bandEnabled) {
    NbandElements = Nelements;
    return;
  }

  // full element data including halo
  o_phi.copyTo(phi);
  mesh.halo.Exchange(phi, Np);

  memory<int> oldLayer(Ntotal);
  for (dlong e=0;e<Ntotal;++e) oldLayer[e] = elementLayer[e];

  // L0: elements with a sign change over their nodes
  memory<int> layer(Ntotal);
  for (dlong e=0;e<Ntotal;++e) layer[e] = -1;
  for (dlong e=0;e<Nelements;++e) {
    dfloat pmin = phi[e*Np], pmax = phi[e*Np];
    for (int n=1;n<Np;++n) {
      pmin = std::min(pmin, phi[e*Np+n]);
      pmax = std::max(pmax, phi[e*Np+n]);
    }
    if (pmin <= 0.0 && pmax >= 0.0) layer[e] = 0;
  }

  // L1..LJ: flag the nodes of layer k-1, max over shared nodes, collect touched elements
  memory<dfloat> flag(Nelements*Np);
  for (int k=1;k<=bandLayers;++k) {
    for (dlong e=0;e<Nelements;++e)
      for (int n=0;n<Np;++n)
        flag[e*Np+n] = (layer[e] == k-1) ? 1.0 : 0.0;

    nodeOgs.GatherScatter(flag, 1, ogs::Max, ogs::Sym);

    for (dlong e=0;e<Nelements;++e) {
      if (layer[e] >= 0) continue;
      for (int n=0;n<Np;++n)
        if (flag[e*Np+n] > 0.0) { layer[e] = k; break; }
    }
  }
  mesh.halo.Exchange(layer, 1);

  bool changed = false;

  // activate new elements, in order of increasing layer so sources exist
  if (!initial) {
    memory<int> active(Ntotal);
    for (dlong e=0;e<Ntotal;++e) active[e] = (oldLayer[e] >= 0) ? 1 : 0;

    for (int k=0;k<=bandLayers;++k) {
      for (dlong e=0;e<Nelements;++e) {
        if (layer[e] != k || active[e]) continue;

        // adjacent source already in the band, with the lowest layer
        hlong src = -1;
        int srcLayer = bandLayers+1;
        for (int f=0;f<Nfaces;++f) {
          const hlong eP = mesh.EToE[e*Nfaces+f];
          if (eP < 0 || eP >= Nelements || !active[eP]) continue;
          const int lP = (layer[eP] >= 0) ? layer[eP] : bandLayers+1;
          if (lP < srcLayer) { src = eP; srcLayer = lP; }
        }

        if (src >= 0)
          ExtendPolynomial(mesh, invV1D, phi, src, e, extensionDegree,
                           std::min((k-1)*hmin, bandClip), bandClip);
        active[e] = 1;
        changed = true;
      }
    }
  }

  // clip elements outside the band
  for (dlong e=0;e<Nelements;++e) {
    if (layer[e] >= 0) continue;
    if (!initial && oldLayer[e] < 0) continue; // already clipped
    dfloat mean = 0.0;
    for (int n=0;n<Np;++n) mean += phi[e*Np+n];
    const dfloat sgn = (mean >= 0.0) ? 1.0 : -1.0;
    for (int n=0;n<Np;++n) phi[e*Np+n] = sgn*bandClip;
    changed = true;
  }

  if (changed) o_phi.copyFrom(phi, Nelements*Np);

  for (dlong e=0;e<Ntotal;++e) elementLayer[e] = layer[e];
  o_elementLayer.copyFrom(elementLayer);

  NbandElements = 0;
  for (dlong e=0;e<Nelements;++e)
    if (layer[e] >= 0) bandElements[NbandElements++] = e;
  if (NbandElements)
    o_bandElements.copyFrom(bandElements, NbandElements);
}

// rebuild once a cut element lies beyond rebuildLayer
void mps_t::BandCheck(){

  if (!bandEnabled) return;

  bandCutKernel(NbandElements, o_bandElements, o_phi, o_elementCut);
  o_elementCut.copyTo(elementCut);

  int rebuild = 0;
  for (dlong l=0;l<NbandElements;++l) {
    const dlong e = bandElements[l];
    if (elementCut[e] && elementLayer[e] > rebuildLayer) rebuild = 1;
    // a sign change in the outermost layer: rebuild around it (it becomes L0) and report
    // where it happened, instead of aborting (thin films and splashes can do this)
    if (elementCut[e] && elementLayer[e] == bandLayers) {
      NedgeCuts++;
      printf("mps warning: zero level in the outermost band layer at (%g, %g), band rebuilt\n",
             mesh.x[e*mesh.Np], mesh.y[e*mesh.Np]);
      if (getenv("MPS_DEBUG_EDGE") && NedgeCuts==1) { // TEMP DEBUG
        o_phi.copyTo(phi);
        dfloat pmin=1e30, pmax=-1e30; int nz=0;
        for (int n=0;n<mesh.Np;++n) { const dfloat p=phi[e*mesh.Np+n]; pmin=std::min(pmin,p); pmax=std::max(pmax,p); if (p==0.0) nz++; }
        printf("DEBUG edge cut: e %d phi [%g, %g], exact zeros %d, reinits %d, mass shift %g\n", (int)e, pmin, pmax, nz, Nreinits, lastMassShift);
        for (int f=0;f<mesh.Nfaces;++f) { hlong eP = mesh.EToE[e*mesh.Nfaces+f]; if (eP<0||eP>=mesh.Nelements) continue;
          dfloat a=1e30,b=-1e30; for (int n=0;n<mesh.Np;++n){a=std::min(a,phi[eP*mesh.Np+n]); b=std::max(b,phi[eP*mesh.Np+n]);}
          printf("   nbr %lld layer %d phi [%g, %g]\n", (long long)eP, elementLayer[eP], a, b); }
        fflush(stdout);
        memory<dfloat> lay(mesh.Nelements); for (dlong ee=0;ee<mesh.Nelements;++ee) lay[ee]=elementLayer[ee];
        o_elementMark.copyTo(elementMark);
        PlotFields(phi, U, elementMark, lay, std::string(getenv("MPS_DEBUG_EDGE")));
      }
    }
  }
  comm.Allreduce(rebuild, Comm::Max);

  if (rebuild) {
    BandBuild(false);
    Nrebuilds++;
  }
}
