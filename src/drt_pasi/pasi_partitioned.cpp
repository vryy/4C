/*---------------------------------------------------------------------------*/
/*! \file
\brief partitioned algorithm for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned.H"

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle_algorithm/particle_algorithm.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PASI::PartitionedAlgo::PartitionedAlgo(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : AlgorithmBase(comm, params), isinit_(false), issetup_(false)
{
  // empty constructor
}

void PASI::PartitionedAlgo::Init()
{
  // reset setup flag
  SetIsSetup(false);

  // init structure field
  InitStructureField();

  // init particle algorithm
  InitParticleAlgorithm();

  // set communication object at the interface
  interface_ = structurefield_->Interface();

  // construct interface states
  intfdispnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
  intfvelnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
  intfaccnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);

  // set init flag
  SetIsInit(true);
}

void PASI::PartitionedAlgo::Setup()
{
  // check correct initialization
  CheckIsInit();

  // setup particle algorithm
  particlealgorithm_->Setup();

  // set setup flag
  SetIsSetup(true);
}

void PASI::PartitionedAlgo::ReadRestart(int restartstep)
{
  // read restart information for structure field
  structurefield_->ReadRestart(restartstep);

  // read restart information for particle algorithm
  particlealgorithm_->ReadRestart(restartstep);

  // set time and step after restart
  SetTimeStep(structurefield_->TimeOld(), restartstep);

  // extract interface states
  ExtractInterfaceStates();

  // set interface states
  SetInterfaceStates(intfdispnp_, intfvelnp_, intfaccnp_);
}

void PASI::PartitionedAlgo::TestResults(const Epetra_Comm& comm)
{
  // get instance of global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // add structure field specific result test object
  problem->AddFieldTest(structurefield_->CreateFieldTest());

  // create particle field specific result test objects
  std::vector<std::shared_ptr<DRT::ResultTest>> allresulttests =
      particlealgorithm_->CreateResultTests();

  // add particle field specific result test objects
  for (auto& resulttest : allresulttests)
    if (resulttest) problem->AddFieldTest(Teuchos::rcp(resulttest));

  // perform all tests
  problem->TestAll(comm);
}

void PASI::PartitionedAlgo::PrepareTimeStep(bool printheader)
{
  // increment time and step
  IncrementTimeAndStep();

  // print header
  if (printheader) PrintHeader();

  // prepare time step of structure field and particle algorithm
  structurefield_->PrepareTimeStep();
  particlealgorithm_->PrepareTimeStep(false);
}

void PASI::PartitionedAlgo::StructStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::StructStep");

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    printf("-------------------- STRUCTURE SOLVER --------------------\n");

  // integrate structural time step
  structurefield_->Solve();
}

void PASI::PartitionedAlgo::ParticleStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::ParticleStep");

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    printf("-------------------- PARTICLE SOLVER ---------------------\n");

  // integrate time step
  particlealgorithm_->IntegrateTimeStep();
}

void PASI::PartitionedAlgo::PostEvaluateTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::PostEvaluateTimeStep");

  // post evaluate time step
  particlealgorithm_->PostEvaluateTimeStep();
}

void PASI::PartitionedAlgo::ExtractInterfaceStates()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::ExtractInterfaceStates");

  intfdispnp_ = interface_->ExtractPASICondVector(structurefield_->Dispnp());
  intfvelnp_ = interface_->ExtractPASICondVector(structurefield_->Velnp());
  intfaccnp_ = interface_->ExtractPASICondVector(structurefield_->Accnp());
}

void PASI::PartitionedAlgo::SetInterfaceStates(Teuchos::RCP<const Epetra_Vector> intfdispnp,
    Teuchos::RCP<const Epetra_Vector> intfvelnp, Teuchos::RCP<const Epetra_Vector> intfaccnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::SetInterfaceStates");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->GetParticleWallHandlerInterface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef DEBUG
  if (walldatastate->GetDispRow() == Teuchos::null or walldatastate->GetDispCol() == Teuchos::null)
    dserror("wall displacements not initialized!");
  if (walldatastate->GetVelCol() == Teuchos::null) dserror("wall velocities not initialized!");
  if (walldatastate->GetAccCol() == Teuchos::null) dserror("wall accelerations not initialized!");
#endif

  // export displacement, velocity and acceleration states
  LINALG::Export(*intfdispnp, *walldatastate->GetMutableDispCol());
  LINALG::Export(*intfvelnp, *walldatastate->GetMutableVelCol());
  LINALG::Export(*intfaccnp, *walldatastate->GetMutableAccCol());

  // export column to row displacements (no communication)
  LINALG::Export(*walldatastate->GetDispCol(), *walldatastate->GetMutableDispRow());

  // print norm of interface displacement to the screen
  if (PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    double normintfdisp(0.0);
    intfdispnp->Norm2(&normintfdisp);

    if (Comm().MyPID() == 0) printf("--> Norm of interface displacement: %10.5E\n", normintfdisp);
  }
}

void PASI::PartitionedAlgo::StructOutput()
{
  // calculate stresses, strains, energies
  structurefield_->PrepareOutput();

  // update all single field solvers
  structurefield_->Update();

  // write output to files
  structurefield_->Output();
}

void PASI::PartitionedAlgo::ParticleOutput()
{
  // write output to files
  particlealgorithm_->Output();
}

void PASI::PartitionedAlgo::InitStructureField()
{
  // get instance of global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter list
  const Teuchos::ParameterList& params = problem->StructuralDynamicParams();

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");

  // build structure
  if (params.get<std::string>("INT_STRATEGY") == "Standard")
  {
    // create and init structure base algorithm
    struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(params);
    struct_adapterbase_ptr_->Init(params, const_cast<Teuchos::ParameterList&>(params), structdis);
  }
  else if (params.get<std::string>("INT_STRATEGY") == "Old")
    dserror(
        "Old time integration not supported in particle structure interaction!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");
  else
    dserror(
        "Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");

  // build and register structure model evaluator
  BuildStructureModelEvaluator();
}

void PASI::PartitionedAlgo::InitParticleAlgorithm()
{
  // get instance of global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter list
  const Teuchos::ParameterList& params = problem->ParticleParams();

  // reference to vector of initial particles
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles = problem->Particles();

  // create and init particle algorithm
  particlealgorithm_ = Teuchos::rcp(new PARTICLEALGORITHM::ParticleAlgorithm(Comm(), params));
  particlealgorithm_->Init(initialparticles);
}

void PASI::PartitionedAlgo::BuildStructureModelEvaluator()
{
  // if adapter base has not already been set up outside.
  if (not struct_adapterbase_ptr_->IsSetup())
  {
    // build and register pasi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> pasi_model_ptr =
        Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedPASI());

    struct_adapterbase_ptr_->RegisterModelEvaluator("Partitioned Coupling Model", pasi_model_ptr);

    // call Setup() on structure base algorithm (wrapper is created inside)
    struct_adapterbase_ptr_->Setup();

    // get wrapper and cast it to specific type
    structurefield_ = Teuchos::rcp_dynamic_cast<::ADAPTER::PASIStructureWrapper>(
        struct_adapterbase_ptr_->StructureField());

    if (structurefield_ == Teuchos::null)
      dserror("No valid pointer to ADAPTER::PASIStructureWrapper set!");

    // set pointer to model evaluator in PASIStructureWrapper
    structurefield_->SetModelEvaluatorPtr(
        Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedPASI>(pasi_model_ptr));
  }
}
