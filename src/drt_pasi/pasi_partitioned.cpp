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

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

#include <Teuchos_TimeMonitor.hpp>

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
  // note that time parameters for subproblems are set in PASI::UTILS::ChangeTimeParameter
  const Teuchos::ParameterList& struct_params = problem->StructuralDynamicParams();
  const Teuchos::ParameterList& particle_params = problem->ParticleParams();

  // safety check: moving_walls_ flag in particle algorithm is always set to 'Yes' for pasi simulations
  if ((bool)DRT::INPUT::IntegralValue<int>(particle_params,"MOVING_WALLS") == false)
  {
    const_cast<Teuchos::ParameterList&>(particle_params).set<std::string>("MOVING_WALLS","Yes");
    std::cout << "WARNING: MOVING_WALLS flag overwritten to 'Yes' for Particle Structure Interaction simulations!!!" << std::endl;
  }

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");

  // build structure
  if (struct_params.get<std::string>("INT_STRATEGY")=="Standard")
  {
    struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(struct_params);

    // initialize structure base algorithm
    struct_adapterbase_ptr_->Init(
        struct_params,
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
  particles_ = Teuchos::rcp(new PARTICLE::Algorithm(comm,particle_params));

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
 | structural time step                                  sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::StructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "                           STRUCTURE SOLVER                           " << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::StructStep");

  // Newton-Raphson iteration
  structure_->Solve();

  return;
} // PASI::PartitionedAlgo::StructStep()

/*----------------------------------------------------------------------*
 | particle time step                                    sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ParticleStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "                           PARTICLE SOLVER                            " << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::ParticleStep");

  // apply external forces or accelerations and to erase outdated vectors
  particles_->AdapterParticle()->UpdateExtActions();

  // integrate particle time step
  particles_->AdapterParticle()->IntegrateStep();

  // adaptions for Normal_DEM_thermo
  particles_->NormDemThermoAdapt();

  return;
} // PASI::PartitionedAlgo::ParticleStep()

/*----------------------------------------------------------------------*
 | set structural displacements and velocities           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::SetStructDispVel(
    Teuchos::RCP<const Epetra_Vector> dispnp,  //! structural displacements \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> velnp    //! structural velocities at \f$t_{n+1}\f$
    )
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PartitionedAlgo::SetStructDispVel");

  double normbdrydispnp;
  double normbdryvelnp;
  dispnp->Norm2(&normbdrydispnp);
  velnp->Norm2(&normbdryvelnp);

  if(Comm().MyPID() == 0)
  {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << " Norm of boundary displacements:  " << std::setprecision(7) << normbdrydispnp <<std::endl;
    std::cout << " Norm of boundary velocities:     " << std::setprecision(7) << normbdryvelnp <<std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
  }

  // hand structural states to particle field
  particles_->SetStructStates(structure_->Dispn(),dispnp,velnp);

  // set structural displacements to particle wall
  particles_->SetUpWallDiscret();

  return;
} // PASI::PartitionedAlgo::UpdateParticleWall()

/*----------------------------------------------------------------------*
 | structural output                                     sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::StructOutput()
{
  // calculate stresses, strains, energies
  structure_->PrepareOutput();

  // update all single field solvers
  structure_->Update();

  // write output to files
  structure_->Output();

  return;
} // PASI::PartitionedAlgo::StructOutput()

/*----------------------------------------------------------------------*
 | particle output                                       sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PartitionedAlgo::ParticleOutput()
{
  // calculate all output quantities
  particles_->AdapterParticle()->PrepareOutput();

  // update particle connectivity
  particles_->UpdateConnectivity();

  // update all single field solvers
  particles_->AdapterParticle()->Update();

  // write output to files
  particles_->Output();

  return;
} // PASI::PartitionedAlgo::ParticleOutput()

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
