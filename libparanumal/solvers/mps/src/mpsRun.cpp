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

void mps_t::Run(){

  if (twoPhase) { RunFlow(); return; }

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  levelSetInitialKernel(mesh.Nelements,
                      startTime,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      o_phi);

  // band around the initial interface, clip outside
  BandBuild(true);

  if (reinitEnabled && reinitInitial) Reinitialize();

  // area to preserve
  if (massCorrection) {
    dfloat Lint;
    MassIntegrals(0.0, massTarget, Lint);
  }

  // initial marks
  MarkElements(o_phi);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step from the initial (untapered) velocity field
  SetVelocity(startTime);
  dfloat vmax = MaxVelocity();

  if (vmax <= 0.0) vmax = 1.0; // no flow (e.g. reinitialization tests)
  dt = cfl*hmin/(vmax*(mesh.N+1.)*(mesh.N+1.));

  if (mesh.rank==0)
    printf("hmin = %g, vmax = %g, dt = %g\n", hmin, vmax, dt);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat time = startTime;
  int tstep = 0;
  Report(time, tstep);
  dfloat outputTime = time + outputInterval;

  const dfloat tol = 1e-10*dt;
  while (time < finalTime - tol) {
    dfloat stepdt = std::min(dt, finalTime - time);
    stepdt = std::min(stepdt, outputTime - time);

    Step(time, stepdt);
    time += stepdt;
    tstep++;

    if (reinitEnabled) {
      bool doReinit = (reinitInterval > 0 && tstep % reinitInterval == 0);
      if (!doReinit && reinitTolerance > 0.0 && tstep % reinitCheckInterval == 0)
        doReinit = (GradientDeviation() > reinitTolerance);
      if (doReinit) {
        Reinitialize();
        if (massCorrection) MassCorrection();
      }
    }

    BandCheck();

    if (time >= outputTime - tol) {
      Report(time, tstep);
      outputTime += outputInterval;
    }
  }

  // final errors against the exact (or initial, for periodic tests) solution
  dfloat errL2, errInf, errL1H, area, areaExact;
  LevelSetErrors(finalTime, errL2, errInf, errL1H, area, areaExact);

  if(mesh.rank==0) {
    printf("\nFinal: N = %d, h = %g\n", mesh.N, hmin);
    printf("  L2 error        = %.6e\n", errL2);
    printf("  Linf error      = %.6e\n", errInf);
    printf("  L1 Heaviside    = %.6e\n", errL1H);
    printf("  area / exact    = %.10f\n", area/areaExact);
    printf("  reinits         = %d\n", Nreinits);
  }
}
