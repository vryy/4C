#ifdef CCADISCRET

#include "fsi_fluid_xfem.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidXFEMAlgorithm::FluidXFEMAlgorithm(Epetra_Comm& comm)
  : FluidMovingBoundaryBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams(),"FSICoupling"),
    comm_(comm)
{
  const Teuchos::ParameterList& fluiddyn   = DRT::Problem::Instance()->FluidDynamicParams();

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fluiddyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fluiddyn.get<double>("TIMESTEP");
  nstep_ = fluiddyn.get<int>("NUMSTEP");
  maxtime_ = fluiddyn.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidXFEMAlgorithm::~FluidXFEMAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::Timeloop()
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
void FSI::FluidXFEMAlgorithm::ReadRestart(int step)
{
  time_ = MBFluidField().ReadRestart(step);
  step_ = step;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  if (Comm().MyPID()==0)
    std::cout << "\n"
              << "TIME:  "    << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_
              << "     STEP = " YELLOW_LIGHT << setw(4) << step_ << END_COLOR "/" << setw(4) << nstep_
              << "\n\n";

  MBFluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::Solve()
{
  MBFluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::Update()
{
  MBFluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::Output()
{
  MBFluidField().Output();
}

#endif
