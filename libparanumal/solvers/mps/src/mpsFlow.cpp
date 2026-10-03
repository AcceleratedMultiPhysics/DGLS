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
#include <sstream>
#include <iomanip>

// EXT/BDF coefficients (as in libs/timeStepper extbdf3), row = order-1
static const dfloat extA[3][3] = {{1., 0., 0.}, {2., -1., 0.}, {3., -3., 1.}};
static const dfloat bdfB[3][4] = {{1., 1., 0., 0.}, {1.5, 2., -0.5, 0.}, {11./6., 3., -1.5, 1./3.}};

static std::string toString(const dfloat v){
  std::ostringstream s; s << std::setprecision(17) << v; return s.str();
}

// Set up the ins operators and solvers with the constant viscosity nu0 of the split
// viscous term, plus the two-phase fields and kernels.
void mps_t::FlowSetup(properties_t& kernelInfo){

  settings.getSetting("DENSITY INSIDE", rhoIn);
  settings.getSetting("DENSITY OUTSIDE", rhoOut);
  settings.getSetting("VISCOSITY INSIDE", muIn);
  settings.getSetting("VISCOSITY OUTSIDE", muOut);
  settings.getSetting("SURFACE TENSION", sigma);
  settings.getSetting("GRAVITY X", gravX);
  settings.getSetting("GRAVITY Y", gravY);
  settings.getSetting("FLOW ORDER", flowOrder);
  settings.getSetting("PRESSURE EXTRAPOLATION ORDER", pressureOrder);
  settings.getSetting("TIME STEP", flowDt);
  explicitViscous = settings.compareSetting("EXPLICIT VISCOUS REMAINDER","TRUE") ? 1 : 0;
  settings.getSetting("FLOW FILTER STRENGTH", filterStrength);
  pressureIncremental = settings.compareSetting("PRESSURE SPLITTING","INCREMENTAL") ? 1 : 0;
  pressureVariable = settings.compareSetting("PRESSURE SOLVER","VARIABLE") ? 1 : 0;
  if (pressureVariable) pressureIncremental = 1; // pressure increment on top of p~
  settings.getSetting("PRESSURE ITERATIONS", pressureIterations);
  settings.getSetting("PRESSURE ITERATION TOLERANCE", pressureIterationTol);
  settings.getSetting("PRESSURE UPDATE RELAXATION", pressureRelaxation);
  settings.getSetting("ADVECTION SUBCYCLES", advectionSubcycles);
  if (pressureIterations < 1) pressureIterations = 1;
  NpressureIters = 0;
  pressureIterationChange = 0.0;
  // the incremental form with a second order extrapolated pressure is unstable
  if (pressureIncremental) pressureOrder = 1;

  dfloat thickness, factor;
  settings.getSetting("PROPERTY THICKNESS", thickness);
  settings.getSetting("VISCOUS SPLIT FACTOR", factor);
  propertyEps = thickness*hmin/mesh.N;
  harmonicDensity = settings.compareSetting("DENSITY INTERPOLATION","HARMONIC") ? 1 : 0;
  densityScaling = settings.compareSetting("SURFACE FORCE DENSITY SCALING","TRUE") ? 1 : 0;
  curvatureLSQ = settings.compareSetting("CURVATURE METHOD","LSQ") ? 1 : 0;
  curvatureWallMirror = settings.compareSetting("CURVATURE WALL MIRROR","TRUE") ? 1 : 0;
  if (curvatureWallMirror && (settings.compareSetting("CURVATURE WALL MIRROR GRADIENT SCALING","TRUE") ||
                              settings.compareSetting("CURVATURE WALL MIRROR GRADIENT SCALING","ELEMENT"))) curvatureWallMirror = 2;
  if (curvatureWallMirror && settings.compareSetting("CURVATURE WALL MIRROR GRADIENT SCALING","CONTACT")) curvatureWallMirror = 3;
  settings.getSetting("CONTACT ANGLE", contactAngle);
  cosContactAngle = cos(contactAngle*M_PI/180.0);
  settings.getSetting("CURVATURE WALL MIRROR WEIGHT", curvatureMirrorWeight);
  contactAngleModel = settings.compareSetting("CONTACT ANGLE MODEL","LEVELSET") ? 1 : 0;
  wallReinit = settings.compareSetting("WALL REINITIALIZATION","TRUE") ? 1 : 0;
  reinitGhostWall = wallReinit;   // the contact angle ghost condition is then also used by the global reinitialization
  settings.getSetting("WALL REINITIALIZATION INTERVAL", wallReinitInterval);
  if (wallReinitInterval < 1) wallReinitInterval = 1;
  settings.getSetting("WALL REINITIALIZATION RELAXATION", wallReinitRelax);
  settings.getSetting("CURVATURE FIT NEIGHBOUR WEIGHT", curvatureNbrWeight);
  settings.getSetting("CURVATURE LIMIT", kappaMax);
  if (kappaMax < 0.0) kappaMax = 0.5*mesh.N/hmin;
  surfaceForceSubcell = settings.compareSetting("SURFACE FORCE SUBCELL","TRUE") ? 1 : 0;
  settings.getSetting("SURFACE FORCE SUBCELL THRESHOLD", surfaceForceSubcellThr);
  surfaceForceSubcellThr *= kappaMax;
  viscousImplicit = settings.compareSetting("VISCOUS TREATMENT","IMPLICIT") ? 1 : 0;
  settings.getSetting("IMPLICIT VISCOUS TOLERANCE", viscousImplicitTol);
  settings.getSetting("IMPLICIT VISCOUS MAX ITERATIONS", viscousImplicitMaxIt);

  // Dodd & Ferrante: rho0 = min(rho); implicit viscosity nu0 >= max(mu/rho) for stability
  rho0 = std::min(rhoIn, rhoOut);
  nu0 = factor*std::max(muIn/rhoIn, muOut/rhoOut);

  if      (settings.compareSetting("LEVEL SET MOTION","STATIC")) levelSetMotion = 1;
  else if (settings.compareSetting("LEVEL SET MOTION","FLOW"))   levelSetMotion = 2;
  else                                                           levelSetMotion = 0;

  flowSettings.changeSetting("TIME INTEGRATOR", "EXTBDF3");
  flowSettings.changeSetting("PRESSURE INCREMENT", pressureIncremental ? "TRUE" : "FALSE");
  flowSettings.changeSetting("VISCOSITY", toString(nu0));
  flowSettings.changeSetting("OUTPUT TO FILE", "FALSE");
  if (mesh.rank==0) flowSettings.report();

  flow.Setup(platform, mesh, flowSettings);

  // variable coefficient pressure operator on the ins pressure space
  if (pressureVariable) {
    LIBP_ABORT("mps: PRESSURE SOLVER VARIABLE needs PRESSURE DISCRETIZATION CONTINUOUS",
               !flow.pDisc_c0);
    static_cast<elliptic_t&>(varPressure) = flow.pSolver;
    std::string fileName = std::string(DMPS "/okl/") + "mpsVarPoisson" + mesh.elementSuffix() + ".okl";
    varPressure.varAxKernel = platform.buildKernel(fileName, "mpsVarPartialAx" + mesh.elementSuffix(), kernelInfo);
    massWeightKernel = platform.buildKernel(fileName, "mpsMassWeight", kernelInfo);
    inverseKernel    = platform.buildKernel(fileName, "mpsInverse", kernelInfo);
    varProjectKernel = platform.buildKernel(fileName, "mpsVarProject", kernelInfo);
  }

  const dlong Ntot = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
  o_rho   = platform.malloc<dfloat>(Ntot);
  o_mu    = platform.malloc<dfloat>(Ntot);
  o_Hf    = platform.malloc<dfloat>(Ntot);
  o_Gf    = platform.malloc<dfloat>(Ntot);
  o_kappa = platform.malloc<dfloat>(Ntot);
  o_wallPhiOld = platform.malloc<dfloat>(Ntot);
  {
    // nodal distance to the nearest wall face (bc 1, 4, 5), for the taper of the local wall reinitialization
    std::vector<dfloat> seg;
    for (dlong e=0;e<mesh.Nelements;++e)
      for (int f=0;f<mesh.Nfaces;++f) {
        const int bc = mesh.EToB[e*mesh.Nfaces+f];
        if (!(bc==1 || bc==4 || bc==5)) continue;
        const int n0 = mesh.faceNodes[f*mesh.Nfp], n1 = mesh.faceNodes[f*mesh.Nfp+mesh.Nfp-1];
        seg.push_back(mesh.x[e*mesh.Np+n0]); seg.push_back(mesh.y[e*mesh.Np+n0]);
        seg.push_back(mesh.x[e*mesh.Np+n1]); seg.push_back(mesh.y[e*mesh.Np+n1]);
      }
    memory<dfloat> wd(Ntot, (dfloat)1e30);
    for (dlong n=0;n<mesh.Nelements*mesh.Np;++n) {
      const dfloat px = mesh.x[n], py = mesh.y[n];
      dfloat dmin = 1e30;
      for (size_t k=0;k<seg.size();k+=4) {
        const dfloat ax = seg[k], ay = seg[k+1], bx = seg[k+2], by = seg[k+3];
        const dfloat ex = bx-ax, ey = by-ay, l2 = ex*ex + ey*ey;
        dfloat t = (l2 > 0.0) ? ((px-ax)*ex + (py-ay)*ey)/l2 : 0.0;
        t = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);
        const dfloat dx = px-ax-t*ex, dy = py-ay-t*ey;
        dmin = std::min(dmin, (dfloat)sqrt(dx*dx + dy*dy));
      }
      wd[n] = dmin;
    }
    o_wallDist = platform.malloc<dfloat>(wd);
  }
  o_Uhist = platform.malloc<dfloat>(3*Ntot*NVfields);
  o_Nhist = platform.malloc<dfloat>(3*Ntot*NVfields);
  o_pPrev  = platform.malloc<dfloat>(Ntot);
  o_pTilde = platform.malloc<dfloat>(Ntot);
  o_Pnew   = platform.malloc<dfloat>(Ntot);
  o_GP      = platform.malloc<dfloat>(Ntot*NVfields);
  o_Visc    = platform.malloc<dfloat>(Ntot*NVfields);
  o_Fst     = platform.malloc<dfloat>(Ntot*NVfields);
  o_stFlag  = platform.malloc<int>(mesh.Nelements);
  if (viscousImplicit) {
    LIBP_ABORT("mps: VISCOUS TREATMENT IMPLICIT needs EXPLICIT VISCOUS REMAINDER TRUE", !explicitViscous);
    const int m = viscousImplicitMaxIt;
    o_gmV = platform.malloc<dfloat>((m+1)*Ntot*NVfields);
    o_gmZ = platform.malloc<dfloat>(m*Ntot*NVfields);
    o_gmW = platform.malloc<dfloat>(Ntot*NVfields);
    // preconditioner solves: own PCG objects with a zero initial guess, so that the Krylov vectors do not
    // enter the extrapolation history of the main velocity solve
    const dlong Nl = flow.uSolver.Ndofs, Nh = flow.uSolver.Nhalo;
    precULinearSolver.Setup<LinearSolver::pcg<dfloat>>(Nl, Nh, platform, flow.vSettings, mesh.comm);
    precVLinearSolver.Setup<LinearSolver::pcg<dfloat>>(Nl, Nh, platform, flow.vSettings, mesh.comm);
    precULinearSolver.SetupInitialGuess<InitialGuess::Zero<dfloat>>(Nl, platform, flow.vSettings, mesh.comm);
    precVLinearSolver.SetupInitialGuess<InitialGuess::Zero<dfloat>>(Nl, platform, flow.vSettings, mesh.comm);
  }
  o_gradPhi = platform.malloc<dfloat>(Ntot*NVfields);
  o_Utilde  = platform.malloc<dfloat>(Ntot*NVfields);
  o_Ustar   = platform.malloc<dfloat>(Ntot*NVfields);
  o_flowRhs = platform.malloc<dfloat>(Ntot*NVfields);
  o_rhsP    = platform.malloc<dfloat>(Ntot);
  o_dP      = platform.malloc<dfloat>(Ntot);
  if (pressureVariable) {
    o_irho = platform.malloc<dfloat>(Ntot);
    o_dPL  = platform.malloc<dfloat>(Ntot);
    o_GdP  = platform.malloc<dfloat>(Ntot*NVfields);
    varPressure.o_irho = o_irho;
  }
  o_S8      = platform.malloc<dfloat>(Ntot*8);
  o_gradMu  = platform.malloc<dfloat>(Ntot*NVfields);
  s8TraceHalo = mesh.HaloTraceSetup(8);

  std::string oklFilePrefix = DMPS "/okl/";
  std::string suffix = mesh.elementSuffix();
  std::string fileName = oklFilePrefix + "mpsFlow.okl";
  propertiesKernel   = platform.buildKernel(fileName, "mpsProperties", kernelInfo);
  flowRhsKernel      = platform.buildKernel(fileName, "mpsFlowRhs", kernelInfo);
  combine3Kernel     = platform.buildKernel(fileName, "mpsCombine3", kernelInfo);
  normalizeKernel    = platform.buildKernel(fileName, "mpsNormalize", kernelInfo);
  surfaceForceKernel = platform.buildKernel(fileName, "mpsSurfaceForce", kernelInfo);
  surfaceForceSubcellKernel = platform.buildKernel(fileName, "mpsSurfaceForceSubcell" + suffix, kernelInfo);
  filterKernel       = platform.buildKernel(fileName, "mpsFilter" + suffix, kernelInfo);

  // 1D modal filter F = V diag(sigma) V^{-1}, top mode(s) above kc = N-1 attenuated by alpha
  {
    const int Nq = mesh.Nq;
    memory<dfloat> V; mesh.Vandermonde1D(mesh.N, mesh.gllz, V);
    memory<dfloat> Vi(Nq*Nq); for (int n=0;n<Nq*Nq;++n) Vi[n] = V[n];
    linAlg_t::matrixInverse(Nq, Vi);
    const int kc = std::max(0, mesh.N-1);
    memory<dfloat> F(Nq*Nq);
    for (int i=0;i<Nq;++i)
      for (int j=0;j<Nq;++j) {
        dfloat s = 0.0;
        for (int k=0;k<Nq;++k) {
          dfloat sig = 1.0;
          if (k > kc) { const dfloat r = (dfloat)(k-kc)/(mesh.N-kc); sig = 1.0 - filterStrength*r*r; }
          s += V[i*Nq+k]*sig*Vi[k*Nq+j];
        }
        F[i*Nq+j] = s;
      }
    o_filter1D = platform.malloc<dfloat>(F);
  }
  fileName = oklFilePrefix + "mpsViscous" + suffix + ".okl";
  viscousGradKernel = platform.buildKernel(fileName, "mpsViscousGrad" + suffix, kernelInfo);
  viscousDivKernel  = platform.buildKernel(fileName, "mpsViscousDiv" + suffix, kernelInfo);
  curvatureKernel   = platform.buildKernel(fileName, "mpsCurvature" + suffix, kernelInfo);
  curvatureLSQKernel = platform.buildKernel(fileName, "mpsCurvatureLSQ" + suffix, kernelInfo);
  wallBlendKernel = platform.buildKernel(fileName, "mpsWallBlend", kernelInfo);
  contactAngleKernel = platform.buildKernel(std::string(DMPS "/okl/") + "mpsContactAngle" + suffix + ".okl",
                                            "mpsContactAngle" + suffix, kernelInfo);
  o_gllz = platform.malloc<dfloat>(mesh.gllz);

  // subcycling: the ins subcycle advection kernels (advecting velocity Ue, advected field U)
  if (advectionSubcycles > 0) {
    LIBP_ABORT("mps: ADVECTION SUBCYCLES is implemented for Quad2D only",
               !(mesh.elementType==Mesh::QUADRILATERALS && mesh.dim==2));
    properties_t subInfo = kernelInfo;
    const int maxNodes = std::max(mesh.Np, mesh.Nfp*mesh.Nfaces);
    subInfo["defines/" "p_maxNodes"] = maxNodes;
    subInfo["defines/" "p_NblockV"] = std::max(1, 256/mesh.Np);
    subInfo["defines/" "p_NblockS"] = std::max(1, 256/maxNodes);
    const std::string subFile = std::string(DMPS "/../ins/okl/") + "insSubcycleAdvection" + suffix + ".okl";
    subAdvVolumeKernel  = platform.buildKernel(subFile, "insSubcycleAdvectionVolume" + suffix, subInfo);
    subAdvSurfaceKernel = platform.buildKernel(subFile, "insSubcycleAdvectionSurface" + suffix, subInfo);
    const dlong NUs = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
    o_Usub   = platform.malloc<dfloat>(3*NUs);
    o_Ue     = platform.malloc<dfloat>(NUs);
    o_subRes = platform.malloc<dfloat>(NUs);
    o_subRhs = platform.malloc<dfloat>(NUs);
  }

  flowStep = 0;
  flowShift = 0;
  lsFlowT0 = 0.0; lsFlowDt = 0.0; NlsSubsteps = 0;
  NiterU = NiterV = NiterP = 0;
}

// rho, mu, H from phi; curvature kappa = div(grad phi/|grad phi|) and the surface force
// sigma kappa grad G with the same DG gradient as the pressure (balanced force), G the
// density weighted Heaviside (the force sits on the heavy side and f/rho does not blow
// up in the light phase at large density ratios)
void mps_t::UpdateProperties(const dfloat T){

  const dlong N = mesh.Nelements*mesh.Np;
  propertiesKernel(N, densityScaling, harmonicDensity, propertyEps, rhoIn, rhoOut, muIn, muOut, o_phi, o_rho, o_mu, o_Hf, o_Gf);

  if (pressureVariable) inverseKernel(N, o_rho, o_irho);

  // viscosity gradient for the (grad u)^T grad mu term
  flow.Gradient(1.0, o_mu, 0.0, o_gradMu, T);

  const dlong Ntot = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
  if (sigma != 0.0) {
    // normals and curvature; the divergence uses the interior trace on boundary faces
    // (a 90 degree contact angle) instead of the ins wall velocity condition
    if (curvatureLSQ) {
      platform.linAlg().set(Ntot, (dfloat)0.0, o_kappa);
      curvatureLSQKernel(NbandElements, o_bandElements, mesh.Nelements, o_EToE, o_elementLayer,
                         mesh.o_x, mesh.o_y, mesh.o_EToB, mesh.o_vmapM, mesh.o_sgeo, curvatureWallMirror, cosContactAngle, curvatureMirrorWeight,
                         hmin, kappaMax, curvatureNbrWeight, o_phi, o_kappa);
    } else {
      flow.Gradient(1.0, o_phi, 0.0, o_gradPhi, T);
      normalizeKernel(mesh.Nelements, o_gradPhi);
      flow.vTraceHalo.Exchange(o_gradPhi, 1);
      curvatureKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_sgeo, mesh.o_D,
                      mesh.o_vmapP, mesh.o_EToB, kappaMax, o_gradPhi, o_kappa);
    }

    flow.Gradient(1.0, o_Gf, 0.0, o_Fst, T);
    // thin structures: gradient of G on the subcell nodes instead of the DG gradient
    if (surfaceForceSubcell)
      surfaceForceSubcellKernel(mesh.Nelements, surfaceForceSubcellThr, mesh.o_vgeo, o_gllz, o_kappa, o_Gf, o_Fst, o_stFlag);
    surfaceForceKernel(mesh.Nelements, sigma, o_kappa, o_Fst);
  } else {
    platform.linAlg().set(Ntot*NVfields, (dfloat)0.0, o_Fst);
  }
}

void mps_t::ViscousRemainder(deviceMemory<dfloat>& o_Uin, deviceMemory<dfloat>& o_Vout, const dfloat T){

  flow.vTraceHalo.Exchange(o_Uin, 1);
  viscousGradKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_sgeo, mesh.o_D,
                    mesh.o_vmapM, mesh.o_vmapP, mesh.o_EToB, T,
                    mesh.o_x, mesh.o_y, nu0, o_mu, o_Uin, o_S8);
  s8TraceHalo.Exchange(o_S8, 1);
  viscousDivKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_sgeo, mesh.o_D,
                   mesh.o_vmapP, mesh.o_EToB, nu0, o_rho, o_gradMu, o_S8, o_Vout);
}

// One step of the JCP 2019 splitting (eq. 17) with variable density, Dodd & Ferrante
// pressure force grad p/rho0 + (1/rho - 1/rho0) grad p~:
//  (a) explicit: rhs = sum b_i U/dt + sum a_i F(U) + V(u~) + g + Fst/rho - (1/rho - s/rho0) grad p~
//  (b) gamma U* - nu0 Lap U* = rhs                     (ins velocity Helmholtz, nu0)
//  INCREMENTAL (s = 0, default; p~ = p^n, as the first order pressure increment of
//  JCP 2019 -- with p~ = 2p^n - p^{n-1} this form is unstable):
//  (c) -Lap dp = -(gamma rho0) div U*,  p = p~ + dp
//  (d) U = U* - (1/(gamma rho0)) grad dp
//  Pressure iterations (INCREMENTAL): (a)-(d) are repeated with p~ = the latest pressure
//  until ||p - p~|| is small, which removes the split error (1/rho0 - 1/rho) grad(p - p~)
//  (each pass contracts it by about 1 - rho0/rho_max, so it helps most at moderate
//  density ratios). Advection, the viscous remainder and the surface force are kept.
//  DODD-FERRANTE (s = 1, testing only): in this velocity-first order U* already meets
//  the wall condition, so the full pressure gets dp/dn = 0 at walls instead of the
//  physical rho g.n (Saini et al. solve the pressure before the viscous step and use
//  the consistent Neumann condition, their eq. 74); an O(1) wall error results.
//  (c) -Lap p = -(gamma rho0) div U*
//  (d) U = U* - (1/(gamma rho0)) grad p
void mps_t::FlowStep(const dfloat time, const dfloat _dt){

  const dlong Ntot = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np;
  const dlong NU = Ntot*NVfields;
  const int order = std::min(flowOrder, flowStep+1);
  const dfloat *a = extA[order-1];
  const dfloat *b = bdfB[order-1];
  const dfloat T = time + _dt;

  // history slices: U^{n-i} and F(U^{n-i})
  deviceMemory<dfloat> U0 = o_Uhist + ((flowShift+0)%3)*NU;
  deviceMemory<dfloat> U1 = o_Uhist + ((flowShift+1)%3)*NU;
  deviceMemory<dfloat> U2 = o_Uhist + ((flowShift+2)%3)*NU;
  deviceMemory<dfloat> N0 = o_Nhist + ((flowShift+0)%3)*NU;
  deviceMemory<dfloat> N1 = o_Nhist + ((flowShift+1)%3)*NU;
  deviceMemory<dfloat> N2 = o_Nhist + ((flowShift+2)%3)*NU;

  // advection: explicit EXT (N^{n-i} histories) or OIFS subcycling (u~_i in the BDF sum, no N term)
  if (advectionSubcycles > 0)
    SubcycleAdvection(time, _dt, order);
  else
    flow.Advection(1.0, U0, 0.0, N0, time);
  const dlong NUs = NU;
  deviceMemory<dfloat> B0 = (advectionSubcycles > 0) ? o_Usub + 0*NUs : U0;
  deviceMemory<dfloat> B1 = (advectionSubcycles > 0) ? o_Usub + 1*NUs : U1;
  deviceMemory<dfloat> B2 = (advectionSubcycles > 0) ? o_Usub + 2*NUs : U2;
  const dfloat ea0 = (advectionSubcycles > 0) ? 0.0 : a[0];
  const dfloat ea1 = (advectionSubcycles > 0) ? 0.0 : a[1];
  const dfloat ea2 = (advectionSubcycles > 0) ? 0.0 : a[2];

  // extrapolated velocity and the explicit viscous remainder at t^{n+1}
  combine3Kernel(NU, a[0], U0, a[1], U1, a[2], U2, o_Utilde);
  if (explicitViscous)
    ViscousRemainder(o_Utilde, o_Visc, T);
  else
    platform.linAlg().set(NU, (dfloat)0.0, o_Visc);

  // extrapolated pressure
  if (pressureOrder==2 && flowStep>0)
    combine3Kernel(Ntot, 2.0, flow.o_p, -1.0, o_pPrev, 0.0, o_pPrev, o_pTilde);
  else
    o_pTilde.copyFrom(flow.o_p);

  const dfloat sRho0 = pressureIncremental ? 0.0 : 1.0/rho0;
  const dfloat gamma = b[0]/_dt;
  const dfloat c = 1.0/(gamma*rho0);
  const int maxPasses = (pressureIncremental && !pressureVariable) ? pressureIterations : 1;

  for (int pass=0;pass<maxPasses;++pass) {

    // p~ from the previous pass
    if (pass>0) o_pTilde.copyFrom(o_Pnew);
    flow.Gradient(1.0, o_pTilde, 0.0, o_GP, T);

    flowRhsKernel(mesh.Nelements, _dt, b[1], b[2], b[3], ea0, ea1, ea2,
                  B0, B1, B2, N0, N1, N2,
                  o_Visc, o_GP, o_Fst, o_rho, sRho0, gravX, gravY, o_flowRhs);

    // velocity Helmholtz with the constant nu0 of the ins solver
    if (pass==0) o_Ustar.copyFrom(o_Utilde);
    flow.VelocitySolve(o_Ustar, o_flowRhs, gamma, T);
    NiterU = flow.NiterU; NiterV = flow.NiterV;
    // implicit variable viscosity: correct U* from the explicit V(u~) to V(U*)
    if (viscousImplicit) ImplicitViscousSolve(o_Ustar, gamma, T);

    // constant coefficient pressure (increment) solve: -Lap p = -(gamma rho0) div U*
    flow.Divergence(-1.0, o_Ustar, 0.0, o_rhsP, T);
    o_Pnew.copyFrom(o_pTilde);
    if (pressureVariable) {
      // -div((1/rho) grad dp) = -gamma div U*: mass weighted local rhs, gather, solve
      const dlong Nloc = mesh.Nelements*mesh.Np;
      deviceMemory<dfloat> o_rhsL = platform.reserve<dfloat>(Nloc);
      massWeightKernel(Nloc, gamma, mesh.o_wJ, o_rhsP, o_rhsL);
      deviceMemory<dfloat> o_Grhs = platform.reserve<dfloat>(varPressure.Ndofs+varPressure.Nhalo);
      deviceMemory<dfloat> o_Gdp  = platform.reserve<dfloat>(varPressure.Ndofs+varPressure.Nhalo);
      varPressure.ogsMasked.Gather(o_Grhs, o_rhsL, 1, ogs::Add, ogs::Trans);
      platform.linAlg().set(varPressure.Ndofs+varPressure.Nhalo, (dfloat)0.0, o_Gdp);
      NiterP = varPressure.Solve(flow.pLinearSolver, o_Gdp, o_Grhs, flow.presTOL, 5000, 0);
      varPressure.ogsMasked.Scatter(o_dPL, o_Gdp, 1, ogs::NoTrans);
      o_Gdp.free(); o_Grhs.free(); o_rhsL.free();
      // p = p~ + dp, U = U* - (1/gamma)(1/rho) grad dp
      platform.linAlg().axpy(Nloc, (dfloat)1.0, o_dPL, (dfloat)1.0, o_Pnew);
      flow.Gradient(1.0, o_dPL, 0.0, o_GdP, T);
      varProjectKernel(mesh.Nelements, (dfloat)(1.0/gamma), o_irho, o_GdP, o_Ustar);
    } else if (pressureIncremental)
      flow.PressureIncrementSolve(o_Pnew, o_rhsP, c, T, _dt);
    else
      flow.PressureSolve(o_Pnew, o_rhsP, c, T);
    if (!pressureVariable) NiterP = flow.NiterP;
    NpressureIters = pass+1;

    if (maxPasses > 1) {
      const dlong Nloc = mesh.Nelements*mesh.Np;
      combine3Kernel(Nloc, 1.0, o_Pnew, -1.0, o_pTilde, 0.0, o_pTilde, o_dP);
      const dfloat d2 = platform.linAlg().innerProd(Nloc, o_dP, o_dP, mesh.comm);
      const dfloat p2 = platform.linAlg().innerProd(Nloc, o_Pnew, o_Pnew, mesh.comm);
      pressureIterationChange = (p2 > 0.0) ? sqrt(d2/p2) : sqrt(d2);
      if (pressureIterationChange < pressureIterationTol) break;
    }
  }

  // projection (done above for the variable coefficient pressure)
  if (!pressureVariable) {
    flow.Gradient(-c, o_Pnew, 1.0, o_Ustar, T);
    if (pressureIncremental)
      flow.Gradient( c, o_pTilde, 1.0, o_Ustar, T);
  }

  // under-relaxed pressure update p = p~ + theta (p - p~); the velocity keeps the full projection
  if (pressureIncremental && pressureRelaxation != 1.0)
    combine3Kernel(Ntot, pressureRelaxation, o_Pnew, 1.0 - pressureRelaxation, o_pTilde, 0.0, o_pTilde, o_Pnew);

  // divergence and continuity penalty of the projected velocity (ins option, off by default)
  flow.PenaltyStep(o_Ustar, _dt);

  // modal filter of the new velocity (spectral-vanishing-viscosity-like stabilization)
  if (filterStrength > 0.0) filterKernel(mesh.Nelements, o_filter1D, o_Ustar);

  // rotate histories: U^{n+1} goes to the slot of U^{n-2}
  U2.copyFrom(o_Ustar, NU);
  flowShift = (flowShift+2)%3;
  flow.o_u.copyFrom(o_Ustar, NU);
  o_pPrev.copyFrom(flow.o_p);
  flow.o_p.copyFrom(o_Pnew);

  flowStep++;
}

void mps_t::RunFlow(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  // level set and band
  levelSetInitialKernel(mesh.Nelements, startTime, mesh.o_x, mesh.o_y, mesh.o_z, o_phi);
  BandBuild(true);
  if (reinitEnabled && reinitInitial) Reinitialize();
  if (massCorrection) { dfloat Lint; MassIntegrals(0.0, massTarget, Lint); }
  MarkElements(o_phi);

  // flow initial condition (ins data macros), all history levels equal
  flow.initialConditionKernel(mesh.Nelements, startTime, mesh.o_x, mesh.o_y, mesh.o_z,
                              nu0, flow.o_u, flow.o_p);
  const dlong NU = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
  for (int i=0;i<3;++i) (o_Uhist + i*NU).copyFrom(flow.o_u, NU);
  o_pPrev.copyFrom(flow.o_p);

  UpdateProperties(startTime);

  // time step: given, or advective CFL and the capillary limit (Brackbill et al.)
  dfloat cfl = 1.0;
  settings.getSetting("CFL NUMBER", cfl);
  if (flowDt <= 0.0) {
    // ins wave speed estimate already carries the 1/h scaling
    dfloat vmax = std::max(flow.MaxWaveSpeed(flow.o_u, startTime), (dfloat)(1.0/hmin));
    flowDt = cfl/(vmax*(mesh.N+1.)*(mesh.N+1.));
    if (advectionSubcycles > 0) flowDt *= advectionSubcycles;   // advection is subcycled
    if (sigma > 0.0) {
      const dfloat hN = hmin/mesh.N;
      flowDt = std::min(flowDt, (dfloat)sqrt((rhoIn+rhoOut)*hN*hN*hN/(4.0*M_PI*sigma)));
    }
  }
  dt = flowDt;

  if (mesh.rank==0)
    printf("two-phase flow: rho in/out %g/%g, mu in/out %g/%g, sigma %g, g (%g,%g), rho0 %g, nu0 %g, dt %g, %s split, p~ order %d, EXT/BDF order %d\n",
           rhoIn, rhoOut, muIn, muOut, sigma, gravX, gravY, rho0, nu0, flowDt,
           pressureIncremental ? "incremental" : "Dodd-Ferrante", pressureOrder, flowOrder);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat time = startTime;
  int tstep = 0;
  Report(time, tstep);
  dfloat outputTime = time + outputInterval;
  const dfloat tol = 1e-10*flowDt;

  while (time < finalTime - tol) {

    // level set first, then properties at t^{n+1}, then the flow (Saini et al., Alg. 1)
    if (levelSetMotion == 2) {
      LevelSetAdvance(time, flowDt, tstep+1);
      // local reinitialization of the wall elements with the contact angle (ghost condition through the LSQ mirror)
      if (wallReinit && NbandElements && (tstep+1) % wallReinitInterval == 0) {
        WallRedistance();
        if (massCorrection) MassCorrection();
        BandCheck();
      }
      // level set contact angle correction in wall elements near the contact line
      if (contactAngleModel==1 && NbandElements)
        contactAngleKernel(NbandElements, o_bandElements, o_elementLayer, mesh.o_EToB, mesh.o_vgeo,
                           mesh.o_sgeo, mesh.o_D, o_gllz, cosContactAngle,
                           (dfloat)sin(contactAngle*M_PI/180.0), o_phi);
      UpdateProperties(time + flowDt);
    } else if (levelSetMotion == 0) {
      LIBP_ABORT("mps: use LEVEL SET MOTION STATIC or FLOW with TWO PHASE FLOW", true);
    }

    FlowStep(time, flowDt);
    time += flowDt;
    tstep++;

    // flag slow pressure solves; stop (with a final snapshot) once the solver hits its limit
    if (NiterP >= 1000 && mesh.rank==0)
      printf("   %.5f (%6d): pressure solve took %d iterations (U %d, V %d)\n",
             time, tstep, NiterP, NiterU, NiterV);
    if (NiterP >= 5000) {
      Report(time, tstep);
      LIBP_FORCE_ABORT("mps: pressure solver did not converge in 5000 iterations");
    }

    if (time >= outputTime - tol) {
      Report(time, tstep);
      outputTime += outputInterval;
    }
  }
}

// Level set over one flow step: LSERK4 substeps on the band with the flow velocity
// extrapolated from U^n, U^{n-1}; then the band check, adaptive/interval reinitialization
// and the mass correction, as in the level set only run.
void mps_t::LevelSetAdvance(const dfloat t0, const dfloat _dt, const int tstep){

  lsFlowT0 = t0;
  lsFlowDt = _dt;

  // level set CFL from the current flow velocity (ins estimate carries 1/h)
  dfloat cfl = 1.0;
  settings.getSetting("CFL NUMBER", cfl);
  dfloat vmax = flow.MaxWaveSpeed(flow.o_u, t0);
  dfloat dtLS = (vmax > 0.0) ? cfl/(vmax*(mesh.N+1.)*(mesh.N+1.)) : _dt;
  NlsSubsteps = std::max(1, (int)ceil(_dt/dtLS - 1e-12));
  if (NlsSubsteps > 10 && mesh.rank==0)
    printf("   %.5f (%6d): %d level set substeps (wave speed %g)\n", t0, tstep, NlsSubsteps, vmax);
  if (NlsSubsteps > 200) {
    Report(t0, tstep);
    LIBP_FORCE_ABORT("mps: level set substeps exceed 200, the flow velocity has blown up");
  }
  const dfloat h = _dt/NlsSubsteps;

  for (int s=0;s<NlsSubsteps;++s) {
    Step(t0 + s*h, h);
    BandCheck();
  }

  if (reinitEnabled) {
    bool doReinit = (reinitInterval > 0 && tstep % reinitInterval == 0);
    if (!doReinit && reinitTolerance > 0.0 && tstep % reinitCheckInterval == 0)
      doReinit = (GradientDeviation() > reinitTolerance);
    if (doReinit) {
      Reinitialize();
      if (massCorrection) MassCorrection();
      BandCheck();
    }
  }
}

// A q = -div((1/rho) grad q) on the gathered C0 space (cf. elliptic_t::Operator)
void mpsVarPressure_t::Operator(deviceMemory<double>& o_q, deviceMemory<double>& o_Aq){

  deviceMemory<double> o_AqL = platform.reserve<double>(mesh.Np*mesh.Nelements);

  gHalo.Exchange(o_q, 1);

  if (mesh.NlocalGatherElements)
    varAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList, o_GlobalToLocal,
                mesh.o_ggeo, mesh.o_D, o_irho, o_q, o_AqL);
  if (mesh.NglobalGatherElements)
    varAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList, o_GlobalToLocal,
                mesh.o_ggeo, mesh.o_D, o_irho, o_q, o_AqL);

  ogsMasked.Gather(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);
}

// OIFS subcycling (Maday, Patera & Ronquist; as in the ins SSBDF scheme): for every BDF history level i,
// integrate u~_t = -div(u~ (x) v) from t^{n+1-i} to t^{n+1} with i*ADVECTION SUBCYCLES LSERK4 substeps,
// starting from u^{n+1-i}; v(t) is the linear extrapolation of u^n, u^{n-1}. The results replace
// u^{n+1-i} in the BDF sum, and the explicit advection term is dropped.
void mps_t::SubcycleAdvection(const dfloat time, const dfloat _dt, const int order){

  const dlong NU = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
  deviceMemory<dfloat> Uh[3];
  for (int i=0;i<3;++i) Uh[i] = o_Uhist + ((flowShift+i)%3)*NU;
  const dfloat dts = _dt/advectionSubcycles;

  for (int i=1;i<=order;++i) {
    deviceMemory<dfloat> o_S = o_Usub + (i-1)*NU;
    o_S.copyFrom(Uh[i-1], NU);
    platform.linAlg().set(NU, (dfloat)0.0, o_subRes);
    dfloat t = time - (i-1)*_dt;
    const int Nsteps = i*advectionSubcycles;
    for (int k=0;k<Nsteps;++k) {
      for (int rk=0;rk<5;++rk) {
        const dfloat ts = t + rkc[rk]*dts;
        // advecting velocity at ts
        const dfloat th = (ts - time)/_dt;
        const dfloat c1 = (flowStep > 0) ? -th : 0.0;
        combine3Kernel(NU, 1.0 - c1, Uh[0], c1, Uh[1], 0.0, Uh[1], o_Ue);
        flow.vTraceHalo.Exchange(o_Ue, 1);
        flow.vTraceHalo.Exchange(o_S, 1);
        subAdvVolumeKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_Ue, o_S, o_subRhs);
        subAdvSurfaceKernel(mesh.Nelements, mesh.o_sgeo, mesh.o_LIFT, mesh.o_vmapM, mesh.o_vmapP,
                            mesh.o_EToB, ts, mesh.o_x, mesh.o_y, mesh.o_z, nu0, o_Ue, o_S, o_subRhs);
        combine3Kernel(NU, rka[rk], o_subRes, dts, o_subRhs, 0.0, o_subRhs, o_subRes);
        combine3Kernel(NU, 1.0, o_S, rkb[rk], o_subRes, 0.0, o_subRes, o_S);
      }
      t += dts;
    }
  }
}

// local reinitialization near walls: the flow of time reinitialization with the contact angle ghost condition on
// the walls, blended into phi near the walls only, phi <- phi + beta w (phi_reinit - phi), w = 1 within h of a wall,
// cosine taper to 0 at 2h
void mps_t::WallRedistance(){
  const dlong Nlocal = mesh.Nelements*mesh.Np;
  o_wallPhiOld.copyFrom(o_phi);
  Reinitialize();
  Nreinits--;
  wallBlendKernel(Nlocal, wallReinitRelax, hmin, o_wallDist, o_wallPhiOld, o_phi);
  NwallReinits++;
}
