/*!----------------------------------------------------------------------
\file pasi_partitioned.cpp

\brief partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_partitioned.H"
#include "pasi_str_model_evaluator_partitioned.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

/*----------------------------------------------------------------------*
 | constructor                                           sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
PASI::PartitionedAlgo::PartitionedAlgo(
    const Epetra_Comm&              comm,        //! communicator
    const Teuchos::ParameterList&   pasi_params  //! input parameters for particle structure interaction
    ) :
    // instantiate base class
    AlgorithmBase(comm,pasi_params),

    structure_(Teuchos::null),
    particles_(Teuchos::null),
    struct_adapterbase_ptr_(Teuchos::null),
    issetup_(false),
    isinit_(false)
{

  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.

} // PASI::PartitionedAlgo::PartitionedAlgo()

/*----------------------------------------------------------------------*
 | Init this class                                       sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::Init(
    const Epetra_Comm& comm  //! communicator
    )
{
  // reset the setup flag
  SetIsSetup(false);

  // get the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter lists
  const Teuchos::ParameterList& pasi_params = problem->PASIDynamicParams();
  const Teuchos::ParameterList& struct_params = problem->StructuralDynamicParams();

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");

  // build structure
  if (struct_params.get<std::string>("INT_STRATEGY")=="Standard")
  {
    struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(struct_params);

    // initialize structure base algorithm
    struct_adapterbase_ptr_->Init(
        pasi_params,
        const_cast<Teuchos::ParameterList&>(struct_params),
        structdis);
  }
  else if (struct_params.get<std::string>("INT_STRATEGY")=="Old")
    dserror("Old time integration not supported in particle structure interaction!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");
  else
    dserror("Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");

  // build particle algorithm with pasi_params for time settings
  particles_ = Teuchos::rcp(new PARTICLE::Algorithm(comm,pasi_params));

  // init particle simulation
  particles_->Init(false);

  // set isinit_ flag true
  SetIsInit(true);

  return;
} // PASI::PartitionedAlgo::Init()

/*----------------------------------------------------------------------*
 | Setup this class                                      sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // if adapter base has not already been set up outside.
  if(not struct_adapterbase_ptr_->IsSetup())
  {
    // build and register pasi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> pasi_model_ptr =
        Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedPASI());

    struct_adapterbase_ptr_->RegisterModelEvaluator("Partitioned Coupling Model",pasi_model_ptr);

    // call Setup() on structure base algorithm (wrapper is created inside)
    struct_adapterbase_ptr_->Setup();

    // get wrapper and cast it to specific type
    structure_ =
        Teuchos::rcp_dynamic_cast< ::ADAPTER::PASIStructureWrapper>(struct_adapterbase_ptr_->StructureField());

    if(structure_ == Teuchos::null)
      dserror("No valid pointer to ADAPTER::PASIStructureWrapper set!");

    // set pointer to model evaluator in PASIStructureWrapper
    structure_->SetModelEvaluatorPtr(Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedPASI>(pasi_model_ptr));
  }

  // set flag issetup true
  SetIsSetup(true);

  return;
} // PASI::PartitionedAlgo::Setup()

/*----------------------------------------------------------------------*
 | read restart data                                     sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ReadRestart(
    int step  //! time step for restart
    )
{
  // read structure restart variables
  structure_->ReadRestart(step);

  // read particle restart variables
  particles_->ReadRestart(step);

  // set time and time step
  SetTimeStep(structure_->TimeOld(),step);

  return;
} // PASI::PartitionedAlgo::ReadRestart()

/*----------------------------------------------------------------------*
 | prepare time step                                     sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();

  if(printheader)
    PrintHeader();

  particles_->PrepareTimeStep(false);

  structure_->PrepareTimeStep();

  return;
} // PASI::PartitionedAlgo::PrepareTimeStep()

/*----------------------------------------------------------------------*
 | single fields are tested                              sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(particles_->AdapterParticle()->CreateFieldTest());
  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->TestAll(comm);

  return;
} // PASI::PartitionedAlgo::TestResults()
