/*---------------------------------------------------------------------------*/
/*! \file
\brief two way coupled partitioned algorithm for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned_twowaycoup.H"

#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle_algorithm/particle_algorithm.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup::PASI_PartTwoWayCoup(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PartitionedAlgo(comm, params),
      itmax_(params.get<int>("ITEMAX")),
      convtolrelativedisp_(params.get<double>("CONVTOLRELATIVEDISP")),
      convtolscaleddisp_(params.get<double>("CONVTOLSCALEDDISP")),
      convtolrelativeforce_(params.get<double>("CONVTOLRELATIVEFORCE")),
      convtolscaledforce_(params.get<double>("CONVTOLSCALEDFORCE")),
      ignoreconvcheck_(DRT::INPUT::IntegralValue<bool>(params, "IGNORE_CONV_CHECK")),
      writerestartevery_(params.get<int>("RESTARTEVRY"))
{
  // empty constructor
}

void PASI::PASI_PartTwoWayCoup::Init()
{
  // call base class init
  PASI::PartitionedAlgo::Init();

  // construct interface force
  intfforcenp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);

  // construct interface increment states
  intfdispincnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
  intfforceincnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);

  // safety check
  if (convtolrelativedisp_ < 0.0 and convtolscaleddisp_ < 0.0 and convtolrelativeforce_ < 0.0 and
      convtolscaledforce_ < 0.0)
    dserror("no convergence tolerance for partitioned iterations set!");
}

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

void PASI::PASI_PartTwoWayCoup::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PartitionedAlgo::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structurefield_->Discretization(), restartstep);
  if (restartstep != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  // get interface force from restart
  reader.ReadVector(intfforcenp_, "intfforcenp");
}

void PASI::PASI_PartTwoWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  while (NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // iteration loop between coupled fields
    Outerloop();

    // post evaluate time step
    PostEvaluateTimeStep();

    // output of fields
    Output();
  }
}

void PASI::PASI_PartTwoWayCoup::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS                                       |\n");
    printf("+------------------------------------------------------------------------------+\n");
    // clang-format on
  }

  // save particle states
  SaveParticleStates();

  while (stopnonliniter == false)
  {
    // increment number of iteration
    ++itnum;

    // reset increment states
    ResetIncrementStates(intfdispnp_, intfforcenp_);

    // reset particle states
    ResetParticleStates();

    // clear interface forces
    ClearInterfaceForces();

    // particle time step
    ParticleStep();

    // get interface forces
    GetInterfaceForces();

    // set interface forces
    SetInterfaceForces(intfforcenp_);

    // structural time step
    StructStep();

    // extract interface states
    ExtractInterfaceStates();

    // build increment states
    BuildIncrementStates();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);

    // set interface states
    SetInterfaceStates(intfdispnp_, intfvelnp_, intfaccnp_);
  }
}

void PASI::PASI_PartTwoWayCoup::Output()
{
  // output of structure field
  StructOutput();

  // write interface force in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
    structurefield_->Discretization()->Writer()->WriteVector("intfforcenp", intfforcenp_);

  // output of particle field
  ParticleOutput();
}

void PASI::PASI_PartTwoWayCoup::ResetIncrementStates(
    Teuchos::RCP<const Epetra_Vector> intfdispnp, Teuchos::RCP<const Epetra_Vector> intfforcenp)
{
  intfdispincnp_->Update(1.0, *intfdispnp, 0.0);
  intfforceincnp_->Update(1.0, *intfforcenp, 0.0);
}

void PASI::PASI_PartTwoWayCoup::BuildIncrementStates()
{
  intfdispincnp_->Update(1.0, *intfdispnp_, -1.0);
  intfforceincnp_->Update(1.0, *intfforcenp_, -1.0);
}

void PASI::PASI_PartTwoWayCoup::SetInterfaceForces(Teuchos::RCP<const Epetra_Vector> intfforcenp)
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::SetInterfaceForces");

  // apply interface force on structure discretization
  structurefield_->ApplyInterfaceForce(intfforcenp);

  // print norm of interface force to the screen
  if (PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    double normintfforce(0.0);
    intfforcenp->Norm2(&normintfforce);

    if (Comm().MyPID() == 0) printf("--> Norm of interface force: %10.5E\n", normintfforce);
  }
}

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

    // reset position, velocity and acceleration states of all particles
    container->UpdateState(0.0, PARTICLEENGINE::Position, 1.0, PARTICLEENGINE::LastIterPosition);
    container->UpdateState(0.0, PARTICLEENGINE::Velocity, 1.0, PARTICLEENGINE::LastIterVelocity);
    container->UpdateState(
        0.0, PARTICLEENGINE::Acceleration, 1.0, PARTICLEENGINE::LastIterAcceleration);

    // reset angular velocity state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::AngularVelocity))
      container->UpdateState(
          0.0, PARTICLEENGINE::AngularVelocity, 1.0, PARTICLEENGINE::LastIterAngularVelocity);

    // reset angular acceleration state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::AngularAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::AngularAcceleration, 1.0,
          PARTICLEENGINE::LastIterAngularAcceleration);

    // reset modified acceleration state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedAcceleration, 1.0,
          PARTICLEENGINE::LastIterModifiedAcceleration);

    // reset density state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::DensityDot))
      container->UpdateState(0.0, PARTICLEENGINE::Density, 1.0, PARTICLEENGINE::LastIterDensity);

    // reset temperature state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::TemperatureDot))
      container->UpdateState(
          0.0, PARTICLEENGINE::Temperature, 1.0, PARTICLEENGINE::LastIterTemperature);
  }
}

void PASI::PASI_PartTwoWayCoup::ClearInterfaceForces()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::ClearInterfaceForces");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->GetParticleWallHandlerInterface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef DEBUG
  if (walldatastate->GetForceCol() == Teuchos::null) dserror("wall forces not initialized!");
#endif

  // clear interface forces
  walldatastate->GetMutableForceCol()->PutScalar(0.0);
}

void PASI::PASI_PartTwoWayCoup::GetInterfaceForces()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::GetInterfaceForces");

  // get interface to particle wall handler
  std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
      particlealgorithm_->GetParticleWallHandlerInterface();

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface->GetWallDataState();

#ifdef DEBUG
  if (walldatastate->GetForceCol() == Teuchos::null) dserror("wall forces not initialized!");
#endif

  // clear interface forces
  intfforcenp_->PutScalar(0.0);

  // assemble interface forces
  Epetra_Export exporter(walldatastate->GetForceCol()->Map(), intfforcenp_->Map());
  int err = intfforcenp_->Export(*walldatastate->GetForceCol(), exporter, Add);
  if (err) dserror("export of interface forces failed with err=%d", err);
}

bool PASI::PASI_PartTwoWayCoup::ConvergenceCheck(int itnum)
{
  bool stopnonliniter = false;

  // variables to save different L2-Norms
  double intfdispincnorm_L2(0.0);
  double intfdispnorm_L2(0.0);
  double intfforceincnorm_L2(0.0);
  double intfforcenorm_L2(0.0);

  // build L2-norm of interface displacement increment and interface displacement
  intfdispincnp_->Norm2(&intfdispincnorm_L2);
  intfdispnp_->Norm2(&intfdispnorm_L2);

  // build L2-norm of interface force increment and interface force
  intfforceincnp_->Norm2(&intfforceincnorm_L2);
  intfforcenp_->Norm2(&intfforcenorm_L2);

  // care for the case that there is (almost) zero scalar
  if (intfdispnorm_L2 < 1e-6) intfdispnorm_L2 = 1.0;
  if (intfforcenorm_L2 < 1e-6) intfforcenorm_L2 = 1.0;

  // scaled and relative interface displacement increment
  double scaled_disp_inc = intfdispincnorm_L2 / (Dt() * sqrt(intfdispincnp_->GlobalLength()));
  double relative_disp_inc = intfdispincnorm_L2 / intfdispnorm_L2;

  // scaled and relative interface force increment
  double scaled_force_inc = intfforceincnorm_L2 / (Dt() * sqrt(intfforceincnp_->GlobalLength()));
  double relative_force_inc = intfforceincnorm_L2 / intfforcenorm_L2;

  // print the incremental based convergence check to the screen
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+----------+-----------------+--------------+------------------+---------------+\n");
    printf("| step/max | scaled-disp-inc | rel-disp-inc | scaled-force-inc | rel-force-inc |\n");
    printf("|  %3d/%3d |      %10.3E |   %10.3E |       %10.3E |    %10.3E |\n", itnum, itmax_, scaled_disp_inc, relative_disp_inc, scaled_force_inc, relative_force_inc);
    printf("+----------+-----------------+--------------+------------------+---------------+\n");
    // clang-format on
  }

  bool isconverged = true;

  // check convergence of scaled interface displacement increment
  if (convtolscaleddisp_ > 0.0) isconverged &= scaled_disp_inc <= convtolscaleddisp_;

  // check convergence of relative interface displacement increment
  if (convtolrelativedisp_ > 0.0) isconverged &= relative_disp_inc <= convtolrelativedisp_;

  // check convergence of scaled interface force increment
  if (convtolscaledforce_ > 0.0) isconverged &= scaled_force_inc <= convtolscaledforce_;

  // check convergence of relative interface force increment
  if (convtolrelativeforce_ > 0.0) isconverged &= relative_force_inc <= convtolrelativeforce_;

  // converged
  if (isconverged)
  {
    stopnonliniter = true;

    if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    {
      // clang-format off
      printf("|  Outer iteration loop converged after iteration %3d/%3d !                    |\n", itnum, itmax_);
      printf("+------------------------------------------------------------------------------+\n");
      // clang-format on
    }
  }

  // stop if maximum iteration number is reached without convergence
  if ((itnum == itmax_) and (not isconverged))
  {
    stopnonliniter = true;

    // ignore convergence check and proceed simulation
    if (ignoreconvcheck_)
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      {
        // clang-format off
        printf("|  ATTENTION: Outer iteration loop not converged in itemax = %3d steps!        |\n", itmax_);
        printf("+------------------------------------------------------------------------------+\n");
        // clang-format on
      }
    }
    // abort the simulation
    else
    {
      if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      {
        // clang-format off
        printf("|  STOP: Outer iteration loop not converged in itemax = %3d steps              |\n", itmax_);
        printf("+------------------------------------------------------------------------------+\n");
        // clang-format on
      }
      dserror("The partitioned PASI solver did not converge in ITEMAX steps!");
    }
  }

  return stopnonliniter;
}

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

    // save position, velocity and acceleration states of all particles
    container->UpdateState(0.0, PARTICLEENGINE::LastIterPosition, 1.0, PARTICLEENGINE::Position);
    container->UpdateState(0.0, PARTICLEENGINE::LastIterVelocity, 1.0, PARTICLEENGINE::Velocity);
    container->UpdateState(
        0.0, PARTICLEENGINE::LastIterAcceleration, 1.0, PARTICLEENGINE::Acceleration);

    // save angular velocity state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::AngularVelocity))
      container->UpdateState(
          0.0, PARTICLEENGINE::LastIterAngularVelocity, 1.0, PARTICLEENGINE::AngularVelocity);

    // save angular acceleration state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::AngularAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::LastIterAngularAcceleration, 1.0,
          PARTICLEENGINE::AngularAcceleration);

    // save modified acceleration state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedAcceleration))
      container->UpdateState(0.0, PARTICLEENGINE::LastIterModifiedAcceleration, 1.0,
          PARTICLEENGINE::ModifiedAcceleration);

    // save density state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::DensityDot))
      container->UpdateState(0.0, PARTICLEENGINE::LastIterDensity, 1.0, PARTICLEENGINE::Density);

    // save temperature state of all particles
    if (container->HaveStoredState(PARTICLEENGINE::TemperatureDot))
      container->UpdateState(
          0.0, PARTICLEENGINE::LastIterTemperature, 1.0, PARTICLEENGINE::Temperature);
  }
}

PASI::PASI_PartTwoWayCoup_DispRelax::PASI_PartTwoWayCoup_DispRelax(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PASI_PartTwoWayCoup(comm, params), omega_(params.get<double>("STARTOMEGA"))
{
  // empty constructor
}

void PASI::PASI_PartTwoWayCoup_DispRelax::Init()
{
  // call base class init
  PASI::PASI_PartTwoWayCoup::Init();

  // construct relaxed interface states
  relaxintfdispnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
  relaxintfvelnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
  relaxintfaccnp_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
}

void PASI::PASI_PartTwoWayCoup_DispRelax::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
  {
    // clang-format off
    printf("+------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS WITH RELAXED DISPLACEMENTS            |\n");
    printf("+------------------------------------------------------------------------------+\n");
    // clang-format on
  }

  // init relaxation of interface states
  InitRelaxationInterfaceStates();

  // set interface states
  SetInterfaceStates(relaxintfdispnp_, relaxintfvelnp_, relaxintfaccnp_);

  // save particle states
  SaveParticleStates();

  while (stopnonliniter == false)
  {
    // increment number of iteration
    ++itnum;

    // reset increment states
    ResetIncrementStates(relaxintfdispnp_, intfforcenp_);

    // reset particle states
    ResetParticleStates();

    // clear interface forces
    ClearInterfaceForces();

    // particle time step
    ParticleStep();

    // get interface forces
    GetInterfaceForces();

    // set interface forces
    SetInterfaceForces(intfforcenp_);

    // structural time step
    StructStep();

    // extract interface states
    ExtractInterfaceStates();

    // build increment states
    BuildIncrementStates();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);

    // calculate relaxation parameter
    CalcOmega(omega_, itnum);

    // perform relaxation of interface states
    PerformRelaxationInterfaceStates();

    // set interface states
    SetInterfaceStates(relaxintfdispnp_, relaxintfvelnp_, relaxintfaccnp_);
  }
}

void PASI::PASI_PartTwoWayCoup_DispRelax::CalcOmega(double& omega, const int itnum)
{
  // output constant relaxation parameter
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    std::cout << "Fixed relaxation parameter: " << omega << std::endl;
}

void PASI::PASI_PartTwoWayCoup_DispRelax::InitRelaxationInterfaceStates()
{
  relaxintfdispnp_->Update(1.0, *intfdispnp_, 0.0);
  relaxintfvelnp_->Update(1.0, *intfvelnp_, 0.0);
  relaxintfaccnp_->Update(1.0, *intfaccnp_, 0.0);
}

void PASI::PASI_PartTwoWayCoup_DispRelax::PerformRelaxationInterfaceStates()
{
  relaxintfdispnp_->Update(omega_, *intfdispincnp_, 1.0);

  relaxintfvelnp_->Update(1.0, *intfdispnp_, 0.0);
  relaxintfvelnp_->Update(1.0 / Dt(), *relaxintfdispnp_, -1.0 / Dt());

  relaxintfaccnp_->Update(1.0, *intfvelnp_, 0.0);
  relaxintfaccnp_->Update(1.0 / Dt(), *relaxintfvelnp_, -1.0 / Dt());
}

PASI::PASI_PartTwoWayCoup_DispRelaxAitken::PASI_PartTwoWayCoup_DispRelaxAitken(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PASI_PartTwoWayCoup_DispRelax(comm, params),
      maxomega_(params.get<double>("MAXOMEGA")),
      minomega_(params.get<double>("MINOMEGA"))
{
  // empty constructor
}

void PASI::PASI_PartTwoWayCoup_DispRelaxAitken::Init()
{
  // call base class init
  PASI::PASI_PartTwoWayCoup_DispRelax::Init();

  // construct old interface increment state
  intfdispincnpold_ = LINALG::CreateVector(*interface_->PASICondMap(), true);
}

void PASI::PASI_PartTwoWayCoup_DispRelaxAitken::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PASI_PartTwoWayCoup_DispRelax::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structurefield_->Discretization(), restartstep);
  if (restartstep != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  // get relaxation parameter from restart
  omega_ = reader.ReadDouble("omega");
}

void PASI::PASI_PartTwoWayCoup_DispRelaxAitken::Output()
{
  // output of structure field
  StructOutput();

  // write interface force and relaxation parameter in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
  {
    structurefield_->Discretization()->Writer()->WriteVector("intfforcenp", intfforcenp_);
    structurefield_->Discretization()->Writer()->WriteDouble("omega", omega_);
  }

  // output of particle field
  ParticleOutput();
}

void PASI::PASI_PartTwoWayCoup_DispRelaxAitken::CalcOmega(double& omega, const int itnum)
{
  Teuchos::RCP<Epetra_Vector> intfdispincnpdiff =
      LINALG::CreateVector(*interface_->PASICondMap(), true);
  intfdispincnpdiff->Update(1.0, *intfdispincnp_, (-1.0), *intfdispincnpold_, 0.0);

  double dispincnpdiffnorm(0.0);
  intfdispincnpdiff->Norm2(&dispincnpdiffnorm);

  if (dispincnpdiffnorm <= 1e-06)
  {
    if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
      std::cout << "Warning: The norm of displacement increment is to small to use it for Aitken "
                   "relaxation. Reuse previous Aitken relaxation parameter instead!"
                << std::endl;
  }

  // in first iteration reuse Aitken relaxation parameter from previous step
  if (itnum != 1 and dispincnpdiffnorm > 1e-06)
  {
    double dispincsdot(0.0);
    intfdispincnpdiff->Dot(*intfdispincnp_, &dispincsdot);

    // update Aitken relaxation parameter
    omega = omega * (1.0 - (dispincsdot) / (dispincnpdiffnorm * dispincnpdiffnorm));

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

  // output Aitken relaxation parameter
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0))
    std::cout << "Aitken relaxation parameter: " << omega << std::endl;

  // store current interface displacement increment for next iteration
  intfdispincnpold_->Update(1.0, *intfdispincnp_, 0.0);
}
