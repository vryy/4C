/*---------------------------------------------------------------------------*/
/*! \file
\brief partitioned algorithm for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_pasi_partitioned.hpp"

#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_pasiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_particle_algorithm.hpp"
#include "4C_particle_wall_datastate.hpp"
#include "4C_particle_wall_interface.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PaSI::PartitionedAlgo::PartitionedAlgo(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : AlgorithmBase(comm, params), isinit_(false), issetup_(false)
{
  // empty constructor
}

void PaSI::PartitionedAlgo::Init()
{
  // reset setup flag
  set_is_setup(false);

  // init structure field
  init_structure_field();

  // init particle algorithm
  init_particle_algorithm();

  // set communication object at the interface
  interface_ = structurefield_->Interface();

  // construct interface states
  intfdispnp_ = Core::LinAlg::CreateVector(*interface_->PASICondMap(), true);
  intfvelnp_ = Core::LinAlg::CreateVector(*interface_->PASICondMap(), true);
  intfaccnp_ = Core::LinAlg::CreateVector(*interface_->PASICondMap(), true);

  // set init flag
  set_is_init(true);
}

void PaSI::PartitionedAlgo::setup()
{
  // check correct initialization
  check_is_init();

  // setup particle algorithm
  particlealgorithm_->setup();

  // write initial output
  structurefield_->Output();

  // set setup flag
  set_is_setup(true);
}

void PaSI::PartitionedAlgo::read_restart(int restartstep)
{
  // read restart information for structure field
  structurefield_->read_restart(restartstep);

  // read restart information for particle algorithm
  particlealgorithm_->read_restart(restartstep);

  // set time and step after restart
  SetTimeStep(structurefield_->TimeOld(), restartstep);

  // extract interface states
  extract_interface_states();

  // set interface states
  set_interface_states(intfdispnp_, intfvelnp_, intfaccnp_);
}

void PaSI::PartitionedAlgo::TestResults(const Epetra_Comm& comm)
{
  // get instance of global problem
  Global::Problem* problem = Global::Problem::Instance();

  // add structure field specific result test object
  problem->AddFieldTest(structurefield_->CreateFieldTest());

  // create particle field specific result test objects
  std::vector<std::shared_ptr<Core::UTILS::ResultTest>> allresulttests =
      particlealgorithm_->CreateResultTests();

  // add particle field specific result test objects
  for (auto& resulttest : allresulttests)
    if (resulttest) problem->AddFieldTest(Teuchos::rcp(resulttest));

  // perform all tests
  problem->TestAll(comm);
}

void PaSI::PartitionedAlgo::prepare_time_step(bool printheader)
{
  // increment time and step
  increment_time_and_step();

  // print header
  if (printheader) print_header();

  // prepare time step of structure field and particle algorithm
  structurefield_->prepare_time_step();
  particlealgorithm_->prepare_time_step(false);
}

void PaSI::PartitionedAlgo::pre_evaluate_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::pre_evaluate_time_step");

  // pre evaluate time step
  particlealgorithm_->pre_evaluate_time_step();
}

void PaSI::PartitionedAlgo::struct_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::StructStep");

  if ((Comm().MyPID() == 0) and print_screen_evry() and (Step() % print_screen_evry() == 0))
    printf("-------------------- STRUCTURE SOLVER --------------------\n");

  // integrate structural time step
  structurefield_->Solve();
}

void PaSI::PartitionedAlgo::particle_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::particle_step");

  if ((Comm().MyPID() == 0) and print_screen_evry() and (Step() % print_screen_evry() == 0))
    printf("-------------------- PARTICLE SOLVER ---------------------\n");

  // integrate time step
  particlealgorithm_->IntegrateTimeStep();
}

void PaSI::PartitionedAlgo::post_evaluate_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::post_evaluate_time_step");

  // post evaluate time step
  particlealgorithm_->post_evaluate_time_step();
}

void PaSI::PartitionedAlgo::extract_interface_states()
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::extract_interface_states");

  intfdispnp_ = interface_->ExtractPASICondVector(structurefield_->Dispnp());
  intfvelnp_ = interface_->ExtractPASICondVector(structurefield_->Velnp());
  intfaccnp_ = interface_->ExtractPASICondVector(structurefield_->Accnp());
}

void PaSI::PartitionedAlgo::set_interface_states(Teuchos::RCP<const Epetra_Vector> intfdispnp,
    Teuchos::RCP<const Epetra_Vector> intfvelnp, Teuchos::RCP<const Epetra_Vector> intfaccnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("PaSI::PartitionedAlgo::set_interface_states");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->get_particle_wall_handler_interface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (walldatastate->GetDispRow() == Teuchos::null or walldatastate->GetDispCol() == Teuchos::null)
    FOUR_C_THROW("wall displacements not initialized!");
  if (walldatastate->GetVelCol() == Teuchos::null) FOUR_C_THROW("wall velocities not initialized!");
  if (walldatastate->GetAccCol() == Teuchos::null)
    FOUR_C_THROW("wall accelerations not initialized!");
#endif

  // export displacement, velocity and acceleration states
  Core::LinAlg::Export(*intfdispnp, *walldatastate->GetDispCol());
  Core::LinAlg::Export(*intfvelnp, *walldatastate->GetVelCol());
  Core::LinAlg::Export(*intfaccnp, *walldatastate->GetAccCol());

  // export column to row displacements (no communication)
  Core::LinAlg::Export(*walldatastate->GetDispCol(), *walldatastate->GetDispRow());

  // print norm of interface displacement to the screen
  if (print_screen_evry() and (Step() % print_screen_evry() == 0))
  {
    double normintfdisp(0.0);
    intfdispnp->Norm2(&normintfdisp);

    if (Comm().MyPID() == 0) printf("--> Norm of interface displacement: %10.5E\n", normintfdisp);
  }
}

void PaSI::PartitionedAlgo::struct_output()
{
  // calculate stresses, strains, energies
  constexpr bool force_prepare = false;
  structurefield_->prepare_output(force_prepare);

  // update all single field solvers
  structurefield_->Update();

  // write output
  structurefield_->Output();
}

void PaSI::PartitionedAlgo::particle_output()
{
  // write output
  particlealgorithm_->WriteOutput();

  // write restart information
  particlealgorithm_->write_restart();
}

void PaSI::PartitionedAlgo::init_structure_field()
{
  // get instance of global problem
  Global::Problem* problem = Global::Problem::Instance();

  // get parameter list
  const Teuchos::ParameterList& params = problem->structural_dynamic_params();

  // access the structural discretization
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");

  // build structure
  if (params.get<std::string>("INT_STRATEGY") == "Standard")
  {
    // create and init structure base algorithm
    struct_adapterbase_ptr_ = Adapter::build_structure_algorithm(params);
    struct_adapterbase_ptr_->Init(params, const_cast<Teuchos::ParameterList&>(params), structdis);
  }
  else if (params.get<std::string>("INT_STRATEGY") == "Old")
    FOUR_C_THROW(
        "Old time integration not supported in particle structure interaction!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");
  else
    FOUR_C_THROW(
        "Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");

  // build and register structure model evaluator
  build_structure_model_evaluator();
}

void PaSI::PartitionedAlgo::init_particle_algorithm()
{
  // get instance of global problem
  Global::Problem* problem = Global::Problem::Instance();

  // get parameter list
  const Teuchos::ParameterList& params = problem->ParticleParams();

  // reference to vector of initial particles
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles = problem->Particles();

  // create and init particle algorithm
  particlealgorithm_ = Teuchos::rcp(new PARTICLEALGORITHM::ParticleAlgorithm(Comm(), params));
  particlealgorithm_->Init(initialparticles);
}

void PaSI::PartitionedAlgo::build_structure_model_evaluator()
{
  // if adapter base has not already been set up outside.
  if (not struct_adapterbase_ptr_->is_setup())
  {
    // build and register pasi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> pasi_model_ptr =
        Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedPASI());

    struct_adapterbase_ptr_->register_model_evaluator("Partitioned Coupling Model", pasi_model_ptr);

    // call setup() on structure base algorithm (wrapper is created inside)
    struct_adapterbase_ptr_->setup();

    // get wrapper and cast it to specific type
    structurefield_ = Teuchos::rcp_dynamic_cast<Adapter::PASIStructureWrapper>(
        struct_adapterbase_ptr_->structure_field());

    if (structurefield_ == Teuchos::null)
      FOUR_C_THROW("No valid pointer to Adapter::PASIStructureWrapper set!");

    // set pointer to model evaluator in PASIStructureWrapper
    structurefield_->set_model_evaluator_ptr(
        Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedPASI>(pasi_model_ptr));
  }
}

FOUR_C_NAMESPACE_CLOSE
