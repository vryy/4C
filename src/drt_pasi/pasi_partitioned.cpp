/*---------------------------------------------------------------------------*/
/*!
\brief partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

#include "../linalg/linalg_mapextractor.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
PASI::PartitionedAlgo::PartitionedAlgo(
    const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params)
    : AlgorithmBase(comm, pasi_params),
      structure_(Teuchos::null),
      particles_(Teuchos::null),
      struct_adapterbase_ptr_(Teuchos::null),
      isinit_(false),
      issetup_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init pasi algorithm                                        sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::Init()
{
  // reset setup flag
  SetIsSetup(false);

  // get the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter lists
  const Teuchos::ParameterList& struct_params = problem->StructuralDynamicParams();
  const Teuchos::ParameterList& particle_params = problem->ParticleParamsOld();

  // safety check
  if ((bool)DRT::INPUT::IntegralValue<int>(particle_params, "MOVING_WALLS") == false)
  {
    const_cast<Teuchos::ParameterList&>(particle_params).set<std::string>("MOVING_WALLS", "Yes");
    std::cout << "WARNING: MOVING_WALLS flag overwritten to 'Yes' for Particle Structure "
                 "Interaction simulations!!!"
              << std::endl;
  }

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");

  // build structure
  if (struct_params.get<std::string>("INT_STRATEGY") == "Standard")
  {
    struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(struct_params);

    // initialize structure base algorithm
    struct_adapterbase_ptr_->Init(
        struct_params, const_cast<Teuchos::ParameterList&>(struct_params), structdis);
  }
  else if (struct_params.get<std::string>("INT_STRATEGY") == "Old")
    dserror(
        "Old time integration not supported in particle structure interaction!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");
  else
    dserror(
        "Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!");

  // build and register structure model evaluator
  BuildStructureModelEvaluator();

  // build particle algorithm with pasi_params for time settings
  particles_ = Teuchos::rcp(new PARTICLE::Algorithm(Comm(), particle_params));

  // init particle simulation
  particles_->Init(false);

  // set init flag
  SetIsInit(true);
}

/*---------------------------------------------------------------------------*
 | setup pasi algorithm                                       sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::Setup()
{
  // check correct initialization
  CheckIsInit();

  // nothing to do

  // set setup flag
  SetIsSetup(true);
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ReadRestart(int restartstep)
{
  // read structure restart variables
  structure_->ReadRestart(restartstep);

  // read particle restart variables
  particles_->ReadRestart(restartstep);

  // set time and step after restart
  SetTimeStep(structure_->TimeOld(), restartstep);
}

/*---------------------------------------------------------------------------*
 | perform result tests                                       sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(particles_->AdapterParticle()->CreateFieldTest());
  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->TestAll(comm);
}

/*---------------------------------------------------------------------------*
 | prepare time step                                          sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::PrepareTimeStep(bool printheader)
{
  // increment time and step
  IncrementTimeAndStep();

  // print header
  if (printheader) PrintHeader();

  structure_->PrepareTimeStep();
  particles_->PrepareTimeStep(false);
}

/*---------------------------------------------------------------------------*
 | structural time step                                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::StructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << " STRUCTURE SOLVER" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::StructStep");

  // integrate structural time step
  structure_->Solve();
}

/*---------------------------------------------------------------------------*
 | particle time step                                         sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ParticleStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << " PARTICLE SOLVER" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::ParticleStep");

  // apply external forces or accelerations and to erase outdated vectors
  particles_->AdapterParticle()->UpdateExtActions();

  // integrate particle time step
  particles_->AdapterParticle()->IntegrateStep();
}

/*---------------------------------------------------------------------------*
 | set structural states                                      sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::SetStructStates()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::SetStructStates");

  // extract wall states from structural states
  Teuchos::RCP<const Epetra_Vector> walldispn =
      particles_->GetWallExtractor()->ExtractCondVector(structure_->Dispn());
  Teuchos::RCP<const Epetra_Vector> walldispnp =
      particles_->GetWallExtractor()->ExtractCondVector(structure_->Dispnp());
  Teuchos::RCP<const Epetra_Vector> wallvelnp =
      particles_->GetWallExtractor()->ExtractCondVector(structure_->Velnp());

  double normbdrydispnp;
  walldispnp->Norm2(&normbdrydispnp);

  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << " Norm of boundary displacements:  " << std::setprecision(7) << normbdrydispnp
              << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }

  // hand wall states to particle field
  particles_->SetWallStates(walldispn, walldispnp, wallvelnp);

  // set wall states to particle wall discretization
  particles_->SetUpWallDiscret();
}

/*---------------------------------------------------------------------------*
 | output of structure field                                  sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::StructOutput()
{
  // calculate stresses, strains, energies
  structure_->PrepareOutput();

  // update all single field solvers
  structure_->Update();

  // write output to files
  structure_->Output();
}

/*---------------------------------------------------------------------------*
 | output of particle field                                   sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ParticleOutput()
{
  // calculate all output quantities
  particles_->AdapterParticle()->PrepareOutput();

  // update all single field solvers
  particles_->AdapterParticle()->Update();

  // write output to files
  particles_->Output();
}

/*---------------------------------------------------------------------------*
 | build and register structure model evaluator               sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
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
    structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::PASIStructureWrapper>(
        struct_adapterbase_ptr_->StructureField());

    if (structure_ == Teuchos::null)
      dserror("No valid pointer to ADAPTER::PASIStructureWrapper set!");

    // set pointer to model evaluator in PASIStructureWrapper
    structure_->SetModelEvaluatorPtr(
        Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedPASI>(pasi_model_ptr));
  }
}
