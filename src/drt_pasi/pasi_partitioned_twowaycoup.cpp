/*---------------------------------------------------------------------------*/
/*! \file
\brief two way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned_twowaycoup.H"

#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle_algorithm/particle_algorithm.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup::PASI_PartTwoWayCoup(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PartitionedAlgo(comm, params),
      itmax_(params.get<int>("ITEMAX")),
      ittol_(params.get<double>("CONVTOL")),
      ignoreconvcheck_(DRT::INPUT::IntegralValue<bool>(params, "IGNORE_CONV_CHECK")),
      writerestartevery_(params.get<int>("RESTARTEVRY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init pasi algorithm                                        sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Init()
{
  // call base class init
  PASI::PartitionedAlgo::Init();

  // construct state and increment vectors
  forcenp_ = LINALG::CreateVector(*structurefield_->DofRowMap(), true);
  dispincnp_ = LINALG::CreateVector(*structurefield_->DofRowMap(), true);
  forceincnp_ = LINALG::CreateVector(*structurefield_->DofRowMap(), true);
}

/*---------------------------------------------------------------------------*
 | setup pasi algorithm                                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Setup()
{
  // call base class setup
  PASI::PartitionedAlgo::Setup();

  // safety check
  {
    // get interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
        particlealgorithm_->GetParticleWallHandlerInterface();

    // get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
        particlewallinterface->GetWallDataState();

    if (walldatastate->GetDispRow() == Teuchos::null or
        walldatastate->GetDispCol() == Teuchos::null)
      dserror("wall displacements not initialized!");
    if (walldatastate->GetVelCol() == Teuchos::null) dserror("wall velocities not initialized!");
    if (walldatastate->GetAccCol() == Teuchos::null) dserror("wall accelerations not initialized!");
    if (walldatastate->GetForceCol() == Teuchos::null) dserror("wall forces not initialized!");
  }
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PartitionedAlgo::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structurefield_->Discretization(), restartstep);
  if (restartstep != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  // get forcenp_ from restart
  reader.ReadVector(forcenp_, "forcenp_");
}

/*---------------------------------------------------------------------------*
 | partitioned two way coupled timeloop                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  while (NotFinished())
  {
    // counter and print header
    PrepareTimeStep();

    // iteration loop between coupled fields
    Outerloop();

    // output of fields
    Output();
  }
}

/*---------------------------------------------------------------------------*
 | iteration loop between coupled fields                      sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+-------------------------------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS                                                                |\n");
    printf("+-------------------------------------------------------------------------------------------------------+\n");
    // clang-format on
  }

  // save particle states
  SaveParticleStates();

  while (stopnonliniter == false)
  {
    // increment number of iteration
    ++itnum;

    // update the current states in every iteration
    IterUpdateStates(structurefield_->Dispnp(), forcenp_);

    // set wall forces
    SetWallForces(forcenp_);

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // reset particle states
    ResetParticleStates();

    // clear wall forces
    ClearWallForces();

    // particle time step
    ParticleStep();

    // get wall forces
    GetWallForces();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);
  }
}

/*---------------------------------------------------------------------------*
 | output of fields                                           sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Output()
{
  // output of structure field
  StructOutput();

  // write interface force in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
    structurefield_->Discretization()->Writer()->WriteVector("forcenp_", forcenp_);

  // output of particle field
  ParticleOutput();
}

/*---------------------------------------------------------------------------*
 | update the current states in every iteration               sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::IterUpdateStates(
    Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<const Epetra_Vector> forcenp)
{
  // store last solutions (current states)
  // will be compared in ConvergenceCheck to the solutions
  // obtained from the next Structure and Particle steps
  dispincnp_->Update(1.0, *dispnp, 0.0);
  forceincnp_->Update(1.0, *forcenp, 0.0);
}

/*---------------------------------------------------------------------------*
 | set wall forces                                            sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::SetWallForces(Teuchos::RCP<const Epetra_Vector> forcenp)
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::SetWallForces");

  double normwallforce(0.0);
  forcenp->Norm2(&normwallforce);

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << " Norm of wall forces: " << std::setprecision(7) << normwallforce << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }

  // apply wall force on structure discretization
  structurefield_->ApplyInterfaceForce(
      structurefield_->Interface()->ExtractPASICondVector(forcenp));
}

/*---------------------------------------------------------------------------*
 | reset particle states                                      sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::ResetParticleStates()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::ResetParticleStates");

  // get interface to particle engine
  std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface =
      particlealgorithm_->GetParticleEngineInterface();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface->GetParticleContainerBundle();

  // iterate over particle types
  for (const auto& type : particlecontainerbundle->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(type, PARTICLEENGINE::Owned);

    // get particle states stored in container
    const std::set<PARTICLEENGINE::StateEnum>& particlestates = container->GetStoredStates();

    // reset particle position, velocity and acceleration states of all particles
    container->UpdateState(0.0, PARTICLEENGINE::Position, 1.0, PARTICLEENGINE::LastIterPosition);
    container->UpdateState(0.0, PARTICLEENGINE::Velocity, 1.0, PARTICLEENGINE::LastIterVelocity);
    container->UpdateState(
        0.0, PARTICLEENGINE::Acceleration, 1.0, PARTICLEENGINE::LastIterAcceleration);

    // reset angular velocity state of all particles
    if (particlestates.count(PARTICLEENGINE::AngularVelocity))
      container->UpdateState(
          0.0, PARTICLEENGINE::AngularVelocity, 1.0, PARTICLEENGINE::LastIterAngularVelocity);

    // reset angular acceleration state of all particles
    if (particlestates.count(PARTICLEENGINE::AngularAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::AngularAcceleration, 1.0,
          PARTICLEENGINE::LastIterAngularAcceleration);

    // reset modified acceleration state of all particles
    if (particlestates.count(PARTICLEENGINE::ModifiedAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedAcceleration, 1.0,
          PARTICLEENGINE::LastIterModifiedAcceleration);
  }
}

/*---------------------------------------------------------------------------*
 | clear wall forces                                          sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::ClearWallForces()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::ClearWallForces");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->GetParticleWallHandlerInterface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef DEBUG
  if (walldatastate->GetForceCol() == Teuchos::null) dserror("wall forces not initialized!");
#endif

  // clear wall forces
  walldatastate->GetMutableForceCol()->PutScalar(0.0);
}

/*---------------------------------------------------------------------------*
 | get wall forces                                            sfuchs 04/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::GetWallForces()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::GetWallForces");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->GetParticleWallHandlerInterface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef DEBUG
  if (walldatastate->GetForceCol() == Teuchos::null) dserror("wall forces not initialized!");
#endif

  // clear wall forces
  forcenp_->PutScalar(0.0);

  // assemble wall forces
  Epetra_Export exporter(walldatastate->GetForceCol()->Map(), forcenp_->Map());
  int err = forcenp_->Export(*walldatastate->GetForceCol(), exporter, Add);
  if (err) dserror("export of wall forces failed with err=%d", err);
}

/*---------------------------------------------------------------------------*
 | convergence check for structure and particles fields       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
bool PASI::PASI_PartTwoWayCoup::ConvergenceCheck(int itnum)
{
  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  // variables to save different L2-Norms
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);
  double forceincnorm_L2(0.0);
  double forcenorm_L2(0.0);

  // build the current displacement increment
  dispincnp_->Update(1.0, *(structurefield_->Dispnp()), -1.0);

  // build the L2-norm of the displacement increment and the displacement
  dispincnp_->Norm2(&dispincnorm_L2);
  structurefield_->Dispnp()->Norm2(&dispnorm_L2);

  // build the current force increment
  forceincnp_->Update(1.0, *forcenp_, -1.0);

  // build the L2-norm of the force increment and the force
  forceincnp_->Norm2(&forceincnorm_L2);
  forcenp_->Norm2(&forcenorm_L2);

  // care for the case that there is (almost) zero scalar
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (forcenorm_L2 < 1e-6) forcenorm_L2 = 1.0;

  // absolute and relative displacement increment
  double abs_disp_inc = dispincnorm_L2;
  double rel_disp_inc = dispincnorm_L2 / dispnorm_L2;

  // absolute and relative force increment
  double abs_force_inc = forceincnorm_L2;
  double rel_force_inc = forceincnorm_L2 / forcenorm_L2;

  // print the incremental based convergence check to the screen
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+------------+--------------------+----------------+----------------+-----------------+-----------------+\n");
    printf("|  step/max  |  tol       [norm]  |  abs-disp-inc  |  rel-disp-inc  |  abs-force-inc  |  rel-force-inc  |\n");
    printf("|   %3d/%3d  | %10.3E [L_2 ]  |    %10.3E  |    %10.3E  |     %10.3E  |     %10.3E  |\n", itnum, itmax_, ittol_, abs_disp_inc, rel_disp_inc, abs_force_inc, rel_force_inc);
    printf("+------------+--------------------+----------------+----------------+-----------------+-----------------+\n");
    // clang-format on
  }

  // converged
  if (abs_disp_inc <= ittol_ and rel_disp_inc <= ittol_ and abs_force_inc <= ittol_ and
      rel_force_inc <= ittol_)
  {
    stopnonliniter = true;

    if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    {
      // clang-format off
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                                             |\n", itnum, itmax_);
      printf("+-------------------------------------------------------------------------------------------------------+\n");
      // clang-format on
    }
  }

  // stop if itemax is reached without convergence
  if ((itnum == itmax_) and (abs_disp_inc > ittol_ or rel_disp_inc > ittol_ or
                                abs_force_inc > ittol_ or rel_force_inc > ittol_))
  {
    stopnonliniter = true;

    // ignore convergence check and proceed simulation
    if (ignoreconvcheck_)
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      {
        // clang-format off
        printf("|  ATTENTION: Outer Iteration loop not converged in itemax = %3d steps!                                 |\n", itmax_);
        printf("+-------------------------------------------------------------------------------------------------------+\n");
        // clang-format on
      }
    }
    // abort the simulation
    else
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      {
        // clang-format off
        printf("|  STOP: Outer Iteration loop not converged in itemax = %3d steps                                       |\n", itmax_);
        printf("+-------------------------------------------------------------------------------------------------------+\n");
        // clang-format on
      }
      dserror("The partitioned PASI solver did not converge in ITEMAX steps!");
    }
  }

  return stopnonliniter;
}

/*---------------------------------------------------------------------------*
 | save particle states                                       sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::SaveParticleStates()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::SaveParticleStates");

  // get interface to particle engine
  std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface =
      particlealgorithm_->GetParticleEngineInterface();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface->GetParticleContainerBundle();

  // iterate over particle types
  for (const auto& type : particlecontainerbundle->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(type, PARTICLEENGINE::Owned);

    // get particle states stored in container
    const std::set<PARTICLEENGINE::StateEnum>& particlestates = container->GetStoredStates();

    // reset particle position, velocity and acceleration states of all particles
    container->UpdateState(0.0, PARTICLEENGINE::LastIterPosition, 1.0, PARTICLEENGINE::Position);
    container->UpdateState(0.0, PARTICLEENGINE::LastIterVelocity, 1.0, PARTICLEENGINE::Velocity);
    container->UpdateState(
        0.0, PARTICLEENGINE::LastIterAcceleration, 1.0, PARTICLEENGINE::Acceleration);

    // reset angular velocity state of all particles
    if (particlestates.count(PARTICLEENGINE::AngularVelocity))
      container->UpdateState(
          0.0, PARTICLEENGINE::LastIterAngularVelocity, 1.0, PARTICLEENGINE::AngularVelocity);

    // reset angular acceleration state of all particles
    if (particlestates.count(PARTICLEENGINE::AngularAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::LastIterAngularAcceleration, 1.0,
          PARTICLEENGINE::AngularAcceleration);

    // reset modified acceleration state of all particles
    if (particlestates.count(PARTICLEENGINE::ModifiedAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::LastIterModifiedAcceleration, 1.0,
          PARTICLEENGINE::ModifiedAcceleration);
  }
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup_ForceRelax::PASI_PartTwoWayCoup_ForceRelax(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PASI_PartTwoWayCoup(comm, params), omega_(params.get<double>("STARTOMEGA"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | iteration loop between coupled fields with relaxed forces  sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+-------------------------------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS WITH RELAXED FORCES                                            |\n");
    printf("+-------------------------------------------------------------------------------------------------------+\n");
    // clang-format on
  }

  // init the relaxed input
  Teuchos::RCP<Epetra_Vector> forcenp = LINALG::CreateVector(*(structurefield_->DofRowMap()), true);

  forcenp->Update(1.0, *forcenp_, 0.0);

  // save particle states
  SaveParticleStates();

  while (stopnonliniter == false)
  {
    // increment number of iteration
    ++itnum;

    // update the states to the last solutions obtained
    IterUpdateStates(structurefield_->Dispnp(), forcenp);

    // set wall forces
    SetWallForces(forcenp);

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // reset particle states
    ResetParticleStates();

    // clear wall forces
    ClearWallForces();

    // particle time step
    ParticleStep();

    // get wall forces
    GetWallForces();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);

    // calculate relaxation parameter
    CalcOmega(omega_, itnum);

    // force relaxation
    forcenp->Update(omega_, *forceincnp_, 1.0);
  }
}

/*---------------------------------------------------------------------------*
 | calculate relaxation parameter                             sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::CalcOmega(double& omega, const int itnum)
{
  // output constant relaxation parameter
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    std::cout << "Fixed relaxation parameter: " << omega << std::endl;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::PASI_PartTwoWayCoup_ForceRelaxAitken(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PASI_PartTwoWayCoup_ForceRelax(comm, params),
      maxomega_(params.get<double>("MAXOMEGA")),
      minomega_(params.get<double>("MINOMEGA"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init pasi algorithm                                        sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::Init()
{
  // call base class init
  PASI::PASI_PartTwoWayCoup_ForceRelax::Init();

  // construct old increment vector
  forceincnpold_ = LINALG::CreateVector(*structurefield_->DofRowMap(), true);
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PASI_PartTwoWayCoup_ForceRelax::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structurefield_->Discretization(), restartstep);
  if (restartstep != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  // get omega_ from restart
  omega_ = reader.ReadDouble("omega_");
}

/*---------------------------------------------------------------------------*
 | output of fields                                           sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::Output()
{
  // output of structure field
  StructOutput();

  // write interface force and relaxation parameter in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
  {
    structurefield_->Discretization()->Writer()->WriteVector("forcenp_", forcenp_);
    structurefield_->Discretization()->Writer()->WriteDouble("omega_", omega_);
  }

  // output of particle field
  ParticleOutput();
}

/*---------------------------------------------------------------------------*
 | calculate relaxation parameter                             sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::CalcOmega(double& omega, const int itnum)
{
  // Aitken relaxation following PhD thesis U. Kuettler, equation (3.5.29)

  Teuchos::RCP<Epetra_Vector> forceincnpdiff =
      LINALG::CreateVector(*structurefield_->DofRowMap(), true);
  forceincnpdiff->Update(1.0, *forceincnp_, (-1.0), *forceincnpold_, 0.0);

  double forceincnpdiffnorm(0.0);
  forceincnpdiff->Norm2(&forceincnpdiffnorm);

  if (forceincnpdiffnorm <= 1e-06)
  {
    if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      std::cout << "Warning: The norm of force increment is to small to use it for Aitken "
                   "relaxation. Reuse previous Aitken relaxation parameter instead!"
                << std::endl;
  }

  // in first iteration reuse Aitken relaxation parameter from previous step
  if (itnum != 1 and forceincnpdiffnorm > 1e-06)
  {
    double forceincsdot(0.0);
    forceincnpdiff->Dot(*forceincnp_, &forceincsdot);

    // update Aitken relaxation parameter
    omega = omega * (1.0 - (forceincsdot) / (forceincnpdiffnorm * forceincnpdiffnorm));

    // allowed range for Aitken relaxation parameter
    if (omega < minomega_)
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
        std::cout << "Warning: The calculation of the relaxation parameter via Aitken did lead to "
                     "a value smaller than MINOMEGA!"
                  << std::endl;
      omega = minomega_;
    }
    if (omega > maxomega_)
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
        std::cout << "Warning: The calculation of the relaxation parameter via Aitken did lead to "
                     "a value bigger than MAXOMEGA!"
                  << std::endl;
      omega = maxomega_;
    }
  }

  // update force increment
  forceincnpold_->Update(1.0, *forceincnp_, 0.0);

  // output Aitken relaxation parameter
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    std::cout << "Aitken relaxation parameter: " << omega << std::endl;
}
