
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "fsi_fluid_ale.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAleAdapter::FluidAleAdapter(const Teuchos::ParameterList& prbdyn)
  : fluid_(prbdyn,true),
    ale_()
{
  FSI::Coupling& coupfa = FluidAleCoupling();

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *fluidnodemap,
                       *alenodemap);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  // the ale matrix is build just once
  AleField().BuildSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAdapter::PrepareTimeStep()
{
  FluidField().PrepareTimeStep();
  AleField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAdapter::Update()
{
  FluidField().Update();
  AleField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAdapter::Output()
{
  FluidField().Output();
  AleField().Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAdapter::NonlinearSolve()
{
  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FSI::FluidAleAdapter::CreateFieldTest()
{
  return FluidField().CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GeneralFluidBaseAlgorithm::GeneralFluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  // here we could do some decision what kind of generalized fluid to build
  fluid_ = Teuchos::rcp(new FluidAleAdapter(prbdyn));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAleAlgorithm::FluidAleAlgorithm(Epetra_Comm& comm)
  : GeneralFluidBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams()),
    comm_(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fsidyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAleAlgorithm::~FluidAleAlgorithm()
{
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
void FSI::FluidAleAlgorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  if (Comm().MyPID()==0)
    std::cout << "\n"
              << "TIME:  "    << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_
              << "     STEP = " YELLOW_LIGHT << setw(4) << step_ << END_COLOR "/" << setw(4) << nstep_
              << "\n\n";

  FluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Solve()
{
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Update()
{
  FluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Output()
{
  FluidField().Output();
}

#endif
