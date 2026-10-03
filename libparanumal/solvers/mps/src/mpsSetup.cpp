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

void mps_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  mpsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = mesh.comm;
  settings = _settings;

  LIBP_ABORT("mps: only quadrilateral meshes are supported for now",
             mesh.elementType!=Mesh::QUADRILATERALS);

  NVfields = mesh.dim;

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max"});

  /*setup trace halo exchanges */
  phiTraceHalo = mesh.HaloTraceSetup(1); //one field
  vTraceHalo   = mesh.HaloTraceSetup(NVfields);

  // low storage RK4 coefficients (Carpenter & Kennedy), as in libs/timeStepper
  const dfloat _rka[5] = {0.0,
       -567301805773.0/1357537059087.0 ,
       -2404267990393.0/2016746695238.0 ,
       -3550918686646.0/2091501179385.0  ,
       -1275806237668.0/842570457699.0};
  const dfloat _rkb[5] = { 1432997174477.0/9575080441755.0 ,
        5161836677717.0/13612068292357.0 ,
        1720146321549.0/2090206949498.0  ,
        3134564353537.0/4481467310338.0  ,
        2277821191437.0/14882151754819.0};
  const dfloat _rkc[6] = {0.0  ,
       1432997174477.0/9575080441755.0 ,
       2526269341429.0/6820363962896.0 ,
       2006345519317.0/3224310063776.0 ,
       2802321613138.0/2924317926251.0 ,
       1.0};
  rka.malloc(5); rka.copyFrom(_rka);
  rkb.malloc(5); rkb.copyFrom(_rkb);
  rkc.malloc(6); rkc.copyFrom(_rkc);

  o_rhsPhi = platform.malloc<dfloat>(Nlocal+Nhalo);
  o_resPhi = platform.malloc<dfloat>(Nlocal+Nhalo);

  // level set and velocity storage (with halo)
  phi.malloc(Nlocal+Nhalo);
  o_phi = platform.malloc<dfloat>(phi);

  U.malloc((Nlocal+Nhalo)*NVfields);
  o_U = platform.malloc<dfloat>(U);

  hmin = mesh.MinCharacteristicLength();

  // narrow band
  bandEnabled = settings.compareSetting("NARROW BAND","TRUE") ? 1 : 0;
  settings.getSetting("BAND LAYERS", bandLayers);
  settings.getSetting("BAND REBUILD LAYER", rebuildLayer);
  bandTaper = settings.compareSetting("BAND TAPER","TRUE") ? 1 : 0;
  settings.getSetting("BAND EXTENSION DEGREE", extensionDegree);
  if (extensionDegree < 0 || extensionDegree > mesh.N) extensionDegree = mesh.N;
  settings.getSetting("BAND TAPER START", taperStart);
  settings.getSetting("BAND TAPER END", taperEnd);
  taperStart *= hmin;
  taperEnd   *= hmin;
  // clip value outside the band, above anything the band holds: a layer J element can
  // lie up to ~(J+1) sqrt(2) h from the interface along a diagonal. (Clipping inside the
  // band as well, to a plateau below the band depth, was tried: the kink it creates
  // costs more than the band edge jump.)
  settings.getSetting("BAND CLIP DISTANCE", bandClip);
  if (bandClip < 0.0) bandClip = (bandLayers+2)*sqrt(2.0);
  bandClip *= hmin;
  Nrebuilds = 0;
  NedgeCuts = 0;

  LIBP_ABORT("mps: BAND LAYERS must exceed BAND REBUILD LAYER + 1",
             bandEnabled && bandLayers < rebuildLayer+2);

  elementLayer.malloc(mesh.Nelements+mesh.totalHaloPairs);
  for (dlong e=0;e<mesh.Nelements+mesh.totalHaloPairs;++e) elementLayer[e] = 0;
  o_elementLayer = platform.malloc<int>(elementLayer);

  elementCut.malloc(mesh.Nelements);
  o_elementCut = platform.malloc<int>(elementCut);

  // gather-scatter over continuous node numbering: elements sharing any node are neighbors
  {
    memory<hlong> ids(Nlocal);
    for (dlong n=0;n<Nlocal;++n) ids[n] = mesh.globalIds[n];
    nodeOgs.Setup(Nlocal, ids, comm, ogs::Signed, ogs::Auto, false, false, platform);
  }

  // reinitialization
  reinitEnabled = settings.compareSetting("REINITIALIZATION","TRUE") ? 1 : 0;
  settings.getSetting("REINIT INTERVAL", reinitInterval);
  reinitInitial = settings.compareSetting("REINIT INITIAL","TRUE") ? 1 : 0;
  settings.getSetting("REINIT DISTANCE", reinitDistance);
  if (reinitDistance < 0.0) reinitDistance = bandLayers+1;
  reinitDistance *= hmin;
  if (bandEnabled) reinitDistance = std::min(reinitDistance, bandClip);
  settings.getSetting("REINIT CFL NUMBER", reinitCFL);
  if      (settings.compareSetting("REINIT SUBCELL MARKING","ALL"))      reinitMarkMode = 1;
  else if (settings.compareSetting("REINIT SUBCELL MARKING","DETECTOR")) reinitMarkMode = 2;
  else                                                                    reinitMarkMode = 0;
  settings.getSetting("REINIT DETECTOR THRESHOLD", reinitDetectorOn);
  settings.getSetting("REINIT DETECTOR RELEASE THRESHOLD", reinitDetectorOff);
  if      (settings.compareSetting("REINIT CUT ELEMENTS","KEEP"))      reinitPreserve = 1;
  else if (settings.compareSetting("REINIT CUT ELEMENTS","NORMALIZE")) reinitPreserve = 2;
  else                                                                  reinitPreserve = 0;
  settings.getSetting("REINIT GRADIENT TOLERANCE", reinitTolerance);
  settings.getSetting("REINIT CHECK INTERVAL", reinitCheckInterval);
  if (reinitCheckInterval < 1) reinitCheckInterval = 1;
  massCorrection = settings.compareSetting("MASS CORRECTION","TRUE") ? 1 : 0;
  massTarget = 0.0;
  lastMassShift = 0.0;
  o_gradDev = platform.malloc<dfloat>(2*mesh.Nelements);
  o_massVL  = platform.malloc<dfloat>(2*mesh.Nelements);
  Nreinits = 0;

  reinitTraceHalo = mesh.HaloTraceSetup(2);
  o_Q     = platform.malloc<dfloat>((Nlocal+Nhalo)*2);
  o_Qrhs  = platform.malloc<dfloat>((Nlocal+Nhalo)*2);
  o_Qres  = platform.malloc<dfloat>((Nlocal+Nhalo)*2);
  o_Qhist = platform.malloc<dfloat>((Nlocal+Nhalo)*2*Nhist);
  o_phi0  = platform.malloc<dfloat>(Nlocal+Nhalo);
  o_arrival = platform.malloc<dfloat>(Nlocal+Nhalo);
  o_found = platform.malloc<int>(Nlocal+Nhalo);
  o_Qmin  = platform.malloc<dfloat>((mesh.Nelements+mesh.totalHaloPairs)*2);
  o_Qprev = platform.malloc<dfloat>((Nlocal+Nhalo)*2);
  o_EToE  = platform.malloc<hlong>(mesh.EToE);
  {
    memory<int> m(mesh.Nelements);
    for (dlong e=0;e<mesh.Nelements;++e) m[e] = (reinitMarkMode==1) ? 1 : 0;
    o_reinitMark  = platform.malloc<int>(m);
    o_reinitMarkV = platform.malloc<int>(m);
  }
  o_reinitDecay = platform.malloc<dfloat>(mesh.Nelements);

  // start with all elements active; BandBuild narrows this down
  NbandElements = mesh.Nelements;
  bandElements.malloc(mesh.Nelements);
  for (dlong e=0;e<mesh.Nelements;++e) bandElements[e] = e;
  o_bandElements = platform.malloc<dlong>(bandElements);

  // subcell marking
  if      (settings.compareSetting("SUBCELL MARKING","ALL"))      subcellMode = 1;
  else if (settings.compareSetting("SUBCELL MARKING","DETECTOR")) subcellMode = 2;
  else                                                             subcellMode = 0;
  settings.getSetting("DETECTOR THRESHOLD", detectorOn);
  settings.getSetting("DETECTOR RELEASE THRESHOLD", detectorOff);

  elementMark.malloc(mesh.Nelements);
  for (dlong e=0;e<mesh.Nelements;++e) elementMark[e] = (subcellMode==1) ? 1 : 0;
  o_elementMark = platform.malloc<int>(elementMark);

  decayRate.malloc(mesh.Nelements);
  for (dlong e=0;e<mesh.Nelements;++e) decayRate[e] = 0.0;
  o_decayRate = platform.malloc<dfloat>(decayRate);

  // modal coefficients c = invV*phi along each direction
  memory<dfloat> V1D;
  mesh.Vandermonde1D(mesh.N, mesh.gllz, V1D);
  linAlg_t::matrixInverse(mesh.Nq, V1D);
  invV1D = V1D;
  o_invV1D = platform.malloc<dfloat>(V1D);

  o_gllw = platform.malloc<dfloat>(mesh.gllw);

  subcellOrder = settings.compareSetting("SUBCELL RECONSTRUCTION","MUSCL") ? 2 : 1;

  edgeTraceHalo = mesh.HaloTraceSetup(4);
  o_phiEdge = platform.malloc<dfloat>((Nlocal+Nhalo)*4);
  o_EToF = platform.malloc<int>(mesh.EToF);

  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_NVfields"]= NVfields;
  kernelInfo["defines/" "p_blockSize"]= 256;
  kernelInfo["defines/" "p_subcellOrder"]= subcellOrder;
  kernelInfo["defines/" "p_Nhist"]= Nhist;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DMPS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName   = oklFilePrefix + "mpsLevelSetVolume" + suffix + oklFileSuffix;
  kernelName = "mpsLevelSetVolume" + suffix;
  levelSetVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsLevelSetUpdate" + oklFileSuffix;
  kernelName = "mpsLevelSetUpdate";
  levelSetUpdateKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsReinit" + oklFileSuffix;
  reinitInitKernel        = platform.buildKernel(fileName, "mpsReinitInit", kernelInfo);
  reinitUpdateKernel      = platform.buildKernel(fileName, "mpsReinitUpdate", kernelInfo);
  reinitSaveHistoryKernel = platform.buildKernel(fileName, "mpsReinitSaveHistory", kernelInfo);
  reinitArrivalKernel     = platform.buildKernel(fileName, "mpsReinitArrival", kernelInfo);
  reinitFinalizeKernel    = platform.buildKernel(fileName, "mpsReinitFinalize", kernelInfo);
  reinitClearFoundKernel  = platform.buildKernel(fileName, "mpsReinitClearFound", kernelInfo);
  reinitElementMinKernel  = platform.buildKernel(fileName, "mpsReinitElementMin", kernelInfo);
  reinitBoundKernel       = platform.buildKernel(fileName, "mpsReinitBound", kernelInfo);
  wallContactPointsKernel = platform.buildKernel(fileName, "mpsWallContactPoints", kernelInfo);

  reinitPreserveKernel = platform.buildKernel(fileName, "mpsReinitPreserveCut" + suffix, kernelInfo);

  fileName   = oklFilePrefix + "mpsGradientDeviation" + suffix + oklFileSuffix;
  gradientDeviationKernel = platform.buildKernel(fileName, "mpsGradientDeviation" + suffix, kernelInfo);

  fileName   = oklFilePrefix + "mpsMass" + oklFileSuffix;
  massIntegralsKernel = platform.buildKernel(fileName, "mpsMassIntegrals", kernelInfo);
  massShiftKernel     = platform.buildKernel(fileName, "mpsMassShift", kernelInfo);

  fileName   = oklFilePrefix + "mpsReinitRhs" + suffix + oklFileSuffix;
  kernelName = "mpsReinitRhs" + suffix;
  reinitRhsKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsBand" + oklFileSuffix;
  bandCutKernel   = platform.buildKernel(fileName, "mpsBandCut", kernelInfo);
  bandTaperKernel = platform.buildKernel(fileName, "mpsBandTaper", kernelInfo);

  fileName   = oklFilePrefix + "mpsLevelSetDetector" + suffix + oklFileSuffix;
  kernelName = "mpsLevelSetDetector" + suffix;
  levelSetDetectorKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsLevelSetReconstruct" + suffix + oklFileSuffix;
  kernelName = "mpsLevelSetReconstruct" + suffix;
  levelSetReconstructKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsLevelSetSurface" + suffix + oklFileSuffix;
  kernelName = "mpsLevelSetSurface" + suffix;
  levelSetSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  std::string dimSuffix = (mesh.dim==2) ? "2D" : "3D";

  fileName   = oklFilePrefix + "mpsLevelSetExact" + dimSuffix + oklFileSuffix;
  kernelName = "mpsLevelSetExact" + dimSuffix;
  levelSetExactKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsLevelSetInitial" + dimSuffix + oklFileSuffix;
  kernelName = "mpsLevelSetInitial" + dimSuffix;
  levelSetInitialKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "mpsVelocityField" + dimSuffix + oklFileSuffix;
  kernelName = "mpsVelocityField" + dimSuffix;
  velocityFieldKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  // two-phase flow
  twoPhase = settings.compareSetting("TWO PHASE FLOW","TRUE") ? 1 : 0;
  if (twoPhase) FlowSetup(kernelInfo);
}
