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

//settings for mps solver
mpsSettings_t::mpsSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/mpsZalesak2D.h",
             "Velocity field, boundary and initial conditions header");

  newSetting("TIME INTEGRATOR",
             "LSERK4",
             "Time integration method for the level set",
             {"LSERK4"});

  newSetting("CFL NUMBER",
             "1.0",
             "Multiplier for timestep stability bound");

  newSetting("START TIME",
             "0",
             "Start time for time integration");

  newSetting("FINAL TIME",
             "10",
             "End time for time integration");

  newSetting("SUBCELL MARKING",
             "DETECTOR",
             "Elements solved with subcell finite volumes",
             {"NONE", "ALL", "DETECTOR"});

  newSetting("SUBCELL RECONSTRUCTION",
             "MUSCL",
             "Subcell finite volume reconstruction",
             {"FIRST", "MUSCL"});

  newSetting("DETECTOR THRESHOLD",
             "2.0",
             "Mark element when modal decay exponent drops below this value");

  newSetting("DETECTOR RELEASE THRESHOLD",
             "2.5",
             "Unmark element only when modal decay exponent exceeds this value");

  newSetting("NARROW BAND",
             "TRUE",
             "Solve the level set only in a band of element layers around the interface",
             {"TRUE", "FALSE"});

  newSetting("BAND LAYERS",
             "6",
             "Number of element layers around the cut elements");

  newSetting("BAND REBUILD LAYER",
             "1",
             "Rebuild the band when a cut element lies beyond this layer");

  newSetting("BAND EXTENSION DEGREE",
             "1",
             "Degree of the polynomial extension into elements entering the band (-1: N; high degrees amplify noise)");

  newSetting("BAND CLIP DISTANCE",
             "-1",
             "phi outside the band is clipped to +/- this distance (in units of h; -1: (BAND LAYERS+2) sqrt(2))");

  newSetting("BAND TAPER",
             "FALSE",
             "Taper the velocity to zero near the band edge (Peng et al.); only with reinitialization",
             {"TRUE", "FALSE"});

  newSetting("BAND TAPER START",
             "2.0",
             "Velocity taper starts at |phi| = value*h");

  newSetting("BAND TAPER END",
             "3.0",
             "Velocity taper reaches zero at |phi| = value*h");

  newSetting("REINITIALIZATION",
             "TRUE",
             "Reinitialize the level set with the flow of time Eikonal equation",
             {"TRUE", "FALSE"});

  newSetting("REINIT INTERVAL",
             "0",
             "Time steps between reinitializations (0: none during the run)");

  newSetting("REINIT INITIAL",
             "FALSE",
             "Reinitialize the initial condition",
             {"TRUE", "FALSE"});

  newSetting("REINIT DISTANCE",
             "-1",
             "Distance (in units of h) covered by reinitialization (-1: BAND LAYERS+1)");

  newSetting("REINIT CFL NUMBER",
             "1.0",
             "Pseudo time step multiplier");

  newSetting("REINIT SUBCELL MARKING",
             "DETECTOR",
             "Elements solved with subcell differences during reinitialization",
             {"NONE", "ALL", "DETECTOR"});

  newSetting("REINIT CUT ELEMENTS",
             "REINITIALIZE",
             "Treatment of interface (cut) elements: reinitialize, keep phi, or phi/|grad phi|",
             {"REINITIALIZE", "KEEP", "NORMALIZE"});

  newSetting("REINIT GRADIENT TOLERANCE",
             "0.2",
             "Reinitialize when the mean ||grad phi|-1| in the two innermost layers exceeds this (0: off)");

  newSetting("REINIT CHECK INTERVAL",
             "10",
             "Time steps between checks of the gradient tolerance");

  newSetting("MASS CORRECTION",
             "TRUE",
             "Shift phi after each reinitialization to restore the initial enclosed area",
             {"TRUE", "FALSE"});

  newSetting("REINIT DETECTOR THRESHOLD",
             "3.0",
             "Mark element during reinitialization when the modal decay exponent drops below this value");

  newSetting("REINIT DETECTOR RELEASE THRESHOLD",
             "3.5",
             "Unmark element during reinitialization only when the decay exponent exceeds this value");

  newSetting("TWO PHASE FLOW",
             "FALSE",
             "Solve the variable density incompressible flow",
             {"TRUE", "FALSE"});

  newSetting("LEVEL SET MOTION",
             "PRESCRIBED",
             "Velocity that moves the level set",
             {"PRESCRIBED", "STATIC", "FLOW"});

  newSetting("DENSITY INSIDE",  "1.0", "Density where phi < 0");
  newSetting("DENSITY OUTSIDE", "1.0", "Density where phi > 0");
  newSetting("VISCOSITY INSIDE",  "0.01", "Dynamic viscosity where phi < 0");
  newSetting("VISCOSITY OUTSIDE", "0.01", "Dynamic viscosity where phi > 0");
  newSetting("SURFACE TENSION", "0.0", "Surface tension coefficient");
  newSetting("GRAVITY X", "0.0", "Body force per unit mass, x");
  newSetting("GRAVITY Y", "0.0", "Body force per unit mass, y");

  newSetting("PROPERTY THICKNESS",
             "3.0",
             "Half width of the property transition in units of h/N");

  newSetting("CURVATURE WALL MIRROR",
             "TRUE",
             "LSQ curvature: add the mirror images of wall element nodes to the fit (90 degree contact angle)",
             {"TRUE", "FALSE"});

  newSetting("CURVATURE WALL MIRROR WEIGHT",
             "4.0",
             "Weight of the mirrored wall points in the LSQ curvature fit (the element's own points have weight 1)");

  newSetting("CURVATURE WALL MIRROR GRADIENT SCALING",
             "FALSE",
             "Scale the contact angle shift of the mirrored values with the local |grad phi|, d phi/dn = -|grad phi| cos(theta): ELEMENT (or TRUE) = plane fit on the element nodes, CONTACT = |d phi/ds|/sin(theta) at the zero of phi on the wall face",
             {"FALSE", "TRUE", "ELEMENT", "CONTACT"});

  newSetting("SURFACE FORCE SUBCELL",
             "FALSE",
             "Thin structures: in elements whose largest |kappa| exceeds THRESHOLD*CURVATURE LIMIT, compute the gradient of the smoothed Heaviside of the surface force by compact differences on the GLL subcell nodes instead of the DG gradient",
             {"TRUE", "FALSE"});

  newSetting("VISCOUS TREATMENT",
             "EXPLICIT",
             "EXPLICIT: the variable viscous remainder V(u~) is extrapolated; IMPLICIT: (gamma - nu0 Lap - V) U = f is solved by FGMRES preconditioned with the nu0 velocity Helmholtz solve (homogeneous velocity boundary data)",
             {"EXPLICIT", "IMPLICIT"});

  newSetting("IMPLICIT VISCOUS TOLERANCE",
             "1e-6",
             "Relative residual reduction of the implicit viscous FGMRES");

  newSetting("IMPLICIT VISCOUS MAX ITERATIONS",
             "20",
             "Krylov dimension (no restart) of the implicit viscous FGMRES");

  newSetting("SURFACE FORCE SUBCELL THRESHOLD",
             "0.25",
             "Fraction of the curvature limit above which an element uses the subcell surface force");

  newSetting("WALL REINITIALIZATION",
             "FALSE",
             "Local reinitialization near walls: the flow of time reinitialization with the contact angle ghost condition d phi/dn = -cos(theta) on wall faces near the interface, blended into phi within 2h of the walls (the ghost condition is then also used by the global reinitialization)",
             {"TRUE", "FALSE"});

  newSetting("WALL REINITIALIZATION RELAXATION",
             "1.0",
             "Relaxation of the local wall reinitialization: phi <- phi + beta w (phi_reinit - phi), w = 1 within h of a wall, 0 beyond 2h");

  newSetting("WALL REINITIALIZATION INTERVAL",
             "50",
             "Flow steps between local wall reinitializations");

  newSetting("CONTACT ANGLE MODEL",
             "CURVATURE",
             "CURVATURE: static contact angle only through the wall mirror of the LSQ curvature; LEVELSET: in addition, rotate the level set gradient at walls to the contact angle in L0/L1 wall elements after every step",
             {"CURVATURE", "LEVELSET"});

  newSetting("CONTACT ANGLE",
             "90.0",
             "Static contact angle in degrees, measured in the phi < 0 fluid, imposed through the wall mirror of the LSQ curvature");

  newSetting("DENSITY INTERPOLATION",
             "ARITHMETIC",
             "ARITHMETIC: rho linear in H; HARMONIC: 1/rho linear in H (smooth 1/rho across the transition)",
             {"ARITHMETIC", "HARMONIC"});

  newSetting("SURFACE FORCE DENSITY SCALING",
             "TRUE",
             "Weight the surface force with rho/rho_avg (Brackbill et al. 1992)",
             {"TRUE", "FALSE"});

  newSetting("CURVATURE METHOD",
             "LSQ",
             "LSQ: least squares cubic patch fit of phi (element and face neighbours), kappa at the closest interface point; DG: div of the DG normals",
             {"LSQ", "DG"});

  newSetting("CURVATURE FIT NEIGHBOUR WEIGHT",
             "0.1",
             "Least squares weight of the face neighbours' nodes (the element's own nodes have weight 1)");

  newSetting("CURVATURE LIMIT",
             "-1",
             "Clip |kappa| to this value (-1: N/(2h), the curvature the mesh resolves; 0: no limit)");

  newSetting("VISCOUS SPLIT FACTOR",
             "1.0",
             "Implicit viscosity nu0 = factor * max(mu/rho)");

  newSetting("EXPLICIT VISCOUS REMAINDER",
             "TRUE",
             "Include (1/rho) div(mu (grad u + grad u^T)) - nu0 Lap u explicitly",
             {"TRUE", "FALSE"});

  newSetting("FLOW FILTER STRENGTH",
             "0",
             "Attenuation of the highest velocity mode per step (modal filter; 0: off)");

  newSetting("TIME STEP",
             "0",
             "Flow time step (0: from the CFL number and the capillary limit)");

  newSetting("FLOW ORDER",
             "2",
             "Order of the EXT/BDF flow time stepping",
             {"1", "2", "3"});

  newSetting("PRESSURE SOLVER",
             "SPLIT",
             "SPLIT: constant coefficient pressure with the density split; VARIABLE: -div((1/rho) grad dp) = -gamma div U* solved with the ins pressure preconditioner",
             {"SPLIT", "VARIABLE"});

  // (option names must not contain each other: settings compare by substring)
  newSetting("PRESSURE SPLITTING",
             "INCREMENTAL",
             "INCREMENTAL: pressure increment with p~ = p^n (JCP 2019 eq. 17 with sigma^1 p); DODD-FERRANTE: full pressure solve -- inconsistent at walls in the velocity-first order (dp/dn = 0 instead of rho g.n), for testing only",
             {"INCREMENTAL", "DODD-FERRANTE"});

  newSetting("PRESSURE ITERATIONS",
             "1",
             "Maximum passes of the velocity/pressure split per step, each using the latest pressure as p~");

  newSetting("ADVECTION SUBCYCLES",
             "0",
             "OIFS subcycling of the momentum advection: number of LSERK4 substeps per flow step (0: explicit EXT advection)");

  newSetting("PRESSURE UPDATE RELAXATION",
             "1.0",
             "theta in p^{n+1} = p~ + theta dp (the projection always uses the full dp); theta < 1 under-relaxes the pressure update, as the ins EXT/BDF step does with theta = 1/b0");

  newSetting("PRESSURE ITERATION TOLERANCE",
             "1e-5",
             "Stop the pressure iterations when ||p_new - p~||/||p_new|| drops below this");

  newSetting("PRESSURE EXTRAPOLATION ORDER",
             "2",
             "Order of the explicit pressure in the density splitting",
             {"1", "2"});

  newSetting("INTERFACE THICKNESS",
             "1.0",
             "Half-width of smoothed Heaviside in units of h/N (used for diagnostics)");

  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "mps");
}

void mpsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "MPS Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("TIME INTEGRATOR");
    reportSetting("CFL NUMBER");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("SUBCELL MARKING");
    reportSetting("SUBCELL RECONSTRUCTION");
    reportSetting("DETECTOR THRESHOLD");
    reportSetting("DETECTOR RELEASE THRESHOLD");
    reportSetting("NARROW BAND");
    reportSetting("BAND LAYERS");
    reportSetting("BAND REBUILD LAYER");
    reportSetting("BAND EXTENSION DEGREE");
    reportSetting("BAND CLIP DISTANCE");
    reportSetting("BAND TAPER");
    reportSetting("BAND TAPER START");
    reportSetting("BAND TAPER END");
    reportSetting("REINITIALIZATION");
    reportSetting("REINIT INTERVAL");
    reportSetting("REINIT INITIAL");
    reportSetting("REINIT DISTANCE");
    reportSetting("REINIT CFL NUMBER");
    reportSetting("REINIT SUBCELL MARKING");
    reportSetting("REINIT CUT ELEMENTS");
    reportSetting("REINIT GRADIENT TOLERANCE");
    reportSetting("REINIT CHECK INTERVAL");
    reportSetting("MASS CORRECTION");
    reportSetting("REINIT DETECTOR THRESHOLD");
    reportSetting("REINIT DETECTOR RELEASE THRESHOLD");
    reportSetting("TWO PHASE FLOW");
    if (compareSetting("TWO PHASE FLOW","TRUE")) {
      reportSetting("LEVEL SET MOTION");
      reportSetting("DENSITY INSIDE");
      reportSetting("DENSITY OUTSIDE");
      reportSetting("VISCOSITY INSIDE");
      reportSetting("VISCOSITY OUTSIDE");
      reportSetting("SURFACE TENSION");
      reportSetting("GRAVITY X");
      reportSetting("GRAVITY Y");
      reportSetting("PROPERTY THICKNESS");
      reportSetting("DENSITY INTERPOLATION");
      reportSetting("CURVATURE WALL MIRROR");
      reportSetting("CONTACT ANGLE");
      reportSetting("CONTACT ANGLE MODEL");
      reportSetting("CURVATURE WALL MIRROR WEIGHT");
      reportSetting("CURVATURE WALL MIRROR GRADIENT SCALING");
      reportSetting("SURFACE FORCE SUBCELL");
      reportSetting("SURFACE FORCE SUBCELL THRESHOLD");
      reportSetting("VISCOUS TREATMENT");
      reportSetting("IMPLICIT VISCOUS TOLERANCE");
      reportSetting("IMPLICIT VISCOUS MAX ITERATIONS");
      reportSetting("WALL REINITIALIZATION");
      reportSetting("WALL REINITIALIZATION INTERVAL");
      reportSetting("WALL REINITIALIZATION RELAXATION");
      reportSetting("SURFACE FORCE DENSITY SCALING");
      reportSetting("CURVATURE METHOD");
      reportSetting("CURVATURE FIT NEIGHBOUR WEIGHT");
      reportSetting("CURVATURE LIMIT");
      reportSetting("VISCOUS SPLIT FACTOR");
      reportSetting("TIME STEP");
      reportSetting("FLOW ORDER");
      reportSetting("PRESSURE SOLVER");
      reportSetting("PRESSURE SPLITTING");
      reportSetting("PRESSURE EXTRAPOLATION ORDER");
      reportSetting("PRESSURE ITERATIONS");
      reportSetting("PRESSURE ITERATION TOLERANCE");
      reportSetting("PRESSURE UPDATE RELAXATION");
      reportSetting("ADVECTION SUBCYCLES");
    }
    reportSetting("INTERFACE THICKNESS");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

// Setting routing: DATA FILE goes to both mps and the flow (one header holds the
// level set and the ins macros); names mps defines go to mps; "FLOW <name>" goes to the
// ins settings as <name>; other ins names (VELOCITY ..., PRESSURE ...) go to ins.
void mpsSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                  meshSettings_t& meshSettings,
                                  insSettings_t& flowSettings,
                                  const std::string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

  for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();
    const std::string val = set.getVal<std::string>();
    if (platformSettings.hasSetting(name))
      platformSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (name == "DATA FILE") {
      changeSetting(name, val);
      flowSettings.changeSetting(name, val);
    } else if (name.rfind("FLOW ", 0) == 0 && flowSettings.hasSetting(name.substr(5)))
      flowSettings.changeSetting(name.substr(5), val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else if (flowSettings.hasSetting(name))
      flowSettings.changeSetting(name, val);
    else  {
      LIBP_FORCE_ABORT("Unknown setting: [" << name << "] requested");
    }
  }
}
