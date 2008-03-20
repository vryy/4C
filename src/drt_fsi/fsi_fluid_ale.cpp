
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
  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                   FluidField().Interface(),
                                  *AleField().Discretization(),
                                   AleField().Interface(),
                                  "FSICoupling");

  //FSI::Coupling& coupfa = FluidAleFieldCoupling();

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa_.SetupCoupling(*FluidField().Discretization(),
                        *AleField().Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  FluidField().SetMeshMap(coupfa_.MasterDofMap());

  // the ale matrix is build just once
  AleField().BuildSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> FSI::FluidAleAdapter::Discretization()
{
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const LINALG::MapExtractor& FSI::FluidAleAdapter::Interface() const
{
  return FluidField().Interface();
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

  FluidField().LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidAleAdapter::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);
  return FluidField().Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAdapter::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                          Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));
    FluidField().ApplyInterfaceVelocities(ivel);
  }

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                                  double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  fluiddisp->Scale(1./dt);

  FluidField().ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1./dt);

  return FluidField().RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::ExtractInterfaceForces()
{
  return FluidField().ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::IntegrateInterfaceShape()
{
  return FluidField().IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FSI::FluidAleAdapter::CreateFieldTest()
{
  return FluidField().CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAleAdapter::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
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
