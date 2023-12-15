/*----------------------------------------------------------------------*/
/*! \file

\brief Solve fluid problems on ALE mesh

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_fsi_fluid_ale.H"

#include "baci_inpar_validparameters.H"
#include "baci_io_pstream.H"
#include "baci_lib_globalproblem.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAleAlgorithm::FluidAleAlgorithm(const Epetra_Comm& comm)
    : FluidMovingBoundaryBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(), "FSICoupling"),
      comm_(comm)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID() == 0) DRT::INPUT::PrintDefaultParameters(IO::cout, fsidyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::ReadRestart(int step)
{
  time_ = MBFluidField()->ReadRestart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  if (Comm().MyPID() == 0)
    std::cout << "\n"
              << "TIME:  " << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_ << "     STEP = " << std::setw(4) << step_
              << "/" << std::setw(4) << nstep_ << "\n\n";

  MBFluidField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Solve() { MBFluidField()->NonlinearSolve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Update() { MBFluidField()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Output() { MBFluidField()->Output(); }

BACI_NAMESPACE_CLOSE
