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

#ifndef MPS_HPP
#define MPS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "ins.hpp"

#define DMPS LIBP_DIR"/solvers/mps/"

using namespace libp;

class mpsSettings_t: public settings_t {
public:
  mpsSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     insSettings_t& flowSettings,
                     const std::string filename);
};

// Variable coefficient pressure operator -div((1/rho) grad p) on the C0 pressure space of
// the ins pressure solver: a copy of that elliptic_t whose Operator is replaced, so the
// linear solver, the (constant coefficient) preconditioner, masking and null space
// handling are reused unchanged.
class mpsVarPressure_t: public elliptic_t {
public:
  kernel_t varAxKernel;
  deviceMemory<dfloat> o_irho;   // local nodal 1/rho

  void Operator(deviceMemory<double>& o_q, deviceMemory<double>& o_Aq) override;
};

// Two-phase incompressible solver with a narrow-band DG level set.
// Level set transport by a prescribed velocity field (flow coupling comes later).
class mps_t: public solver_t {
public:
  mesh_t mesh;

  // low storage RK4 (Carpenter-Kennedy) restricted to the band element list
  memory<dfloat> rka, rkb, rkc;
  deviceMemory<dfloat> o_rhsPhi, o_resPhi;
  dfloat dt;

  int NVfields;            // number of velocity components
  dfloat hmin;             // minimum element characteristic length

  ogs::halo_t phiTraceHalo;
  ogs::halo_t vTraceHalo;

  // level set field
  memory<dfloat> phi;
  deviceMemory<dfloat> o_phi;

  // velocity, element-major: U[e*Np*NVfields + fld*Np + n]
  memory<dfloat> U;
  deviceMemory<dfloat> o_U;

  // subcell finite volume marking (one flag per element, 1 = subcell FV)
  int subcellMode;         // 0 = none, 1 = all elements, 2 = modal detector
  dfloat detectorOn;       // mark if decay exponent s < detectorOn
  dfloat detectorOff;      // keep a marked element while s < detectorOff
  memory<int> elementMark;
  deviceMemory<int> o_elementMark;
  memory<dfloat> decayRate;
  deviceMemory<dfloat> o_decayRate;
  memory<dfloat> invV1D;          // 1D inverse Vandermonde (orthonormal Legendre on GLL)
  deviceMemory<dfloat> o_invV1D;
  deviceMemory<dfloat> o_gllw;    // GLL weights = subcell widths in reference space
  int subcellOrder;               // 1 = piecewise constant, 2 = MUSCL (minmod)

  // reconstructed subcell edge states, element-major with 4 fields:
  //  0: left edge in r (face 3), 1: right edge in r (face 1),
  //  2: left edge in s (face 0), 3: right edge in s (face 2)
  // unmarked elements store their nodal values (the DG trace)
  deviceMemory<dfloat> o_phiEdge;
  ogs::halo_t edgeTraceHalo;
  deviceMemory<int> o_EToF;

  // narrow band of element layers around the interface
  int bandEnabled;
  int bandLayers;          // J: layers L0 (cut) .. LJ
  int rebuildLayer;        // rebuild when a cut element lies beyond this layer
  dfloat bandClip;         // phi = +/- bandClip outside the band
  int bandTaper;           // taper the velocity near the band edge (needs reinitialization)
  dfloat taperStart;       // velocity taper c(|phi|) = 1 below, 0 above taperEnd
  dfloat taperEnd;
  memory<int> elementLayer;       // with halo, -1 outside the band
  deviceMemory<int> o_elementLayer;
  memory<int> elementCut;
  deviceMemory<int> o_elementCut;
  dlong NbandElements;
  memory<dlong> bandElements;
  deviceMemory<dlong> o_bandElements;
  int Nrebuilds;
  int NedgeCuts;           // sign changes found in the outermost layer
  ogs::ogs_t nodeOgs;             // gather-scatter over shared nodes (vertex adjacency)
  int extensionDegree;            // degree of the polynomial extension into new band elements

  // reinitialization (flow of time Eikonal, CAMWA 2022) on the band
  int reinitEnabled;
  int reinitInterval;      // steps between reinitializations (0: never during the run)
  int reinitInitial;       // reinitialize the initial condition
  dfloat reinitDistance;   // pseudo time (= distance) covered
  dfloat reinitCFL;
  int reinitMarkMode;      // 0 none, 1 all, 2 detector on u, v every pseudo step
  dfloat reinitDetectorOn, reinitDetectorOff;
  int reinitPreserve;      // cut elements: 0 reinitialized, 1 kept, 2 normalized
  dfloat reinitTolerance;  // adaptive trigger: reinitialize when the mean ||grad phi|-1| in L0,L1 exceeds this
  int reinitCheckInterval; // steps between adaptive checks
  int massCorrection;      // restore the enclosed area after each reinitialization
  dfloat massTarget;       // enclosed (smoothed Heaviside) area to restore
  dfloat lastMassShift;
  kernel_t gradientDeviationKernel, massIntegralsKernel, massShiftKernel;
  deviceMemory<dfloat> o_gradDev, o_massVL;
  kernel_t reinitPreserveKernel;
  static constexpr int Nhist = 4;
  ogs::halo_t reinitTraceHalo;
  deviceMemory<dfloat> o_Q, o_Qrhs, o_Qres, o_Qhist, o_phi0, o_arrival;
  deviceMemory<int> o_found, o_reinitMark, o_reinitMarkV;
  deviceMemory<dfloat> o_reinitDecay;
  int Nreinits;
  kernel_t reinitInitKernel, reinitRhsKernel, reinitUpdateKernel, reinitSaveHistoryKernel;
  kernel_t reinitArrivalKernel, reinitFinalizeKernel, reinitClearFoundKernel;
  kernel_t reinitElementMinKernel, reinitBoundKernel;
  deviceMemory<dfloat> o_Qmin, o_Qprev;
  deviceMemory<hlong> o_EToE;
  kernel_t levelSetInitialKernel;

  // two-phase flow: libParanumal ins operators and solvers (JCP 2019 splitting, eq. 17)
  // with the variable density adapted as in Dodd & Ferrante / Saini et al.
  int twoPhase;
  int levelSetMotion;      // 0 prescribed, 1 static, 2 flow velocity
  insSettings_t flowSettings;
  ins_t flow;
  dfloat rhoIn, rhoOut, muIn, muOut, sigma, gravX, gravY;
  dfloat rho0, nu0;        // constant density and viscosity of the split operators
  dfloat propertyEps;      // property transition half width
  dfloat kappaMax;         // curvature limit
  int flowOrder, pressureOrder;
  int explicitViscous;
  dfloat filterStrength;          // modal filter of the velocity after each step (0: off)
  deviceMemory<dfloat> o_filter1D;
  kernel_t filterKernel;
  int pressureIncremental;
  int pressureVariable;           // 1: variable coefficient pressure Poisson (no split)
  mpsVarPressure_t varPressure;
  deviceMemory<dfloat> o_irho, o_dPL, o_GdP;
  kernel_t massWeightKernel, inverseKernel, varProjectKernel;
  int pressureIterations;         // max split passes per step
  dfloat pressureIterationTol;
  dfloat pressureRelaxation;   // p^{n+1} = p~ + theta dp
  int NpressureIters;             // passes taken in the last step
  dfloat pressureIterationChange; // last relative pressure change
  deviceMemory<dfloat> o_dP;
  dfloat flowDt;
  int flowStep, flowShift; // steps taken, history ring index of U^n
  deviceMemory<dfloat> o_rho, o_mu, o_Hf, o_Gf, o_kappa;
  int densityScaling;
  int harmonicDensity;     // 1/rho (instead of rho) linear in H
  int curvatureWallMirror; // LSQ fit with mirrored nodes at walls (contact angle)
  dfloat curvatureMirrorWeight;   // weight of the mirrored points in the LSQ fit
  dfloat contactAngle, cosContactAngle; // static contact angle [degrees], measured in the phi < 0 fluid
  deviceMemory<dfloat> o_Uhist, o_Nhist;   // 3 levels each, NVfields
  deviceMemory<dfloat> o_pPrev, o_pTilde, o_Pnew;
  deviceMemory<dfloat> o_GP, o_Visc, o_Fst, o_gradPhi, o_Utilde, o_Ustar, o_flowRhs, o_rhsP;
  deviceMemory<dfloat> o_S8;               // mu grad u and grad u (8 fields)
  deviceMemory<dfloat> o_gradMu;
  ogs::halo_t s8TraceHalo;
  kernel_t propertiesKernel, flowRhsKernel, combine3Kernel, normalizeKernel, surfaceForceKernel;
  kernel_t viscousGradKernel, viscousDivKernel, curvatureKernel, curvatureLSQKernel;
  kernel_t contactAngleKernel;
  kernel_t wallBlendKernel;
  kernel_t surfaceForceSubcellKernel;
  int surfaceForceSubcell=0;
  dfloat surfaceForceSubcellThr=0;
  deviceMemory<int> o_stFlag;
  kernel_t wallContactPointsKernel;
  deviceMemory<dfloat> o_contactPts;   // per element: contact point on a wall face (ghost condition of the reinitialization)
  int wallReinit, wallReinitInterval, NwallReinits=0;
  dfloat wallReinitRelax;
  int reinitGhostWall=0;              // contact angle ghost condition in the reinitialization
  deviceMemory<dfloat> o_wallDist;     // nodal distance to the nearest wall

  deviceMemory<dfloat> o_wallPhiOld;  // phi before the local wall reinitialization
  // OIFS subcycling of the momentum advection (ADVECTION SUBCYCLES > 0)
  int advectionSubcycles;
  kernel_t subAdvVolumeKernel, subAdvSurfaceKernel;
  deviceMemory<dfloat> o_Usub, o_Ue, o_subRes, o_subRhs;
  void SubcycleAdvection(const dfloat time, const dfloat _dt, const int order);
  int contactAngleModel;           // 0 CURVATURE (mirror in the LSQ fit only), 1 LEVELSET (also correct phi at walls)
  deviceMemory<dfloat> o_gllz;
  int curvatureLSQ;        // 1: least squares patch fit curvature
  dfloat curvatureNbrWeight;
  int NiterU, NiterV, NiterP;
  dfloat lsFlowT0, lsFlowDt;       // flow step the level set is advanced over
  int NlsSubsteps;

  // advance the level set from t0 to t0+dt with the flow velocity, then band, reinit, mass
  void LevelSetAdvance(const dfloat t0, const dfloat dt, const int tstep);

  kernel_t levelSetVolumeKernel;
  kernel_t levelSetUpdateKernel;
  kernel_t bandCutKernel;
  kernel_t bandTaperKernel;
  kernel_t levelSetDetectorKernel;
  kernel_t levelSetReconstructKernel;
  kernel_t levelSetSurfaceKernel;
  kernel_t levelSetExactKernel;
  kernel_t velocityFieldKernel;

  mps_t() = default;
  mps_t(platform_t &_platform, mesh_t &_mesh,
        mpsSettings_t& _settings, insSettings_t& _flowSettings):
    flowSettings(_flowSettings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& platform, mesh_t& mesh,
             mpsSettings_t& settings);

  // two-phase flow setup (called from Setup when TWO PHASE FLOW is on)
  void FlowSetup(properties_t& kernelInfo);

  // density, viscosity, curvature and surface tension force from phi at time T
  void UpdateProperties(const dfloat T);

  // explicit viscous remainder (1/rho) div(mu (grad u + grad u^T)) - nu0 Lap u
  void ViscousRemainder(deviceMemory<dfloat>& o_Uin, deviceMemory<dfloat>& o_Vout, const dfloat T);

  // implicit variable viscosity: FGMRES on (gamma - nu0 Lap - V) U = f, preconditioned by the nu0 Helmholtz solve
  int viscousImplicit=0, viscousImplicitMaxIt=20, NiterVisc=0;
  dfloat viscousImplicitTol=1e-6, viscousImplicitRes=0;
  linearSolver_t<dfloat> precULinearSolver, precVLinearSolver;
  deviceMemory<dfloat> o_gmV, o_gmZ, o_gmW;
  void PrecVelocitySolve(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_out, const dfloat gamma, const dfloat T);
  void ImplicitViscousSolve(deviceMemory<dfloat>& o_Uc, const dfloat gamma, const dfloat T);

  // one flow step from time to time+dt
  void FlowStep(const dfloat time, const dfloat dt);

  // flow run loop
  void RunFlow();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Phi, memory<dfloat> V, memory<int> Mark,
                  memory<dfloat> Decay, const std::string fileName,
                  const std::vector<std::pair<std::string, memory<dfloat>>>& extra = {});

  // level set rhs: dphi/dt = -u.grad(phi)
  void rhsf(deviceMemory<dfloat>& o_phi, deviceMemory<dfloat>& o_rhs, const dfloat time);

  // one LSERK4 step over the band
  void Step(const dfloat time, const dfloat dt);

  // (re)build the band from the current level set; initial = clip everything outside
  void BandBuild(const bool initial);

  // rebuild the band if the interface has moved too far into it
  void BandCheck();

  // replace phi by the signed distance to its zero level set, on the band
  void Reinitialize();
  void WallRedistance();

  // mean ||grad phi| - 1| over layers L0 and L1 (global, wJ weighted)
  dfloat GradientDeviation();

  // enclosed area and interface length (smoothed) after shifting band phi by s
  void MassIntegrals(const dfloat s, dfloat& V, dfloat& Lint);

  // shift phi on the band so that the enclosed area equals massTarget
  void MassCorrection();

  // update subcell marks from the current level set
  void MarkElements(deviceMemory<dfloat>& o_Phi);

  // number of marked elements in the band (global)
  hlong NumberOfMarkedElements();

  // fill o_U with the prescribed velocity field at time T
  void SetVelocity(const dfloat T);

  // max velocity magnitude over all nodes
  dfloat MaxVelocity();

  // L2 error, interface (L1 Heaviside) error, and enclosed area
  void LevelSetErrors(const dfloat T, dfloat& errL2, dfloat& errInf,
                      dfloat& errL1H, dfloat& area, dfloat& areaExact);
};

#endif
