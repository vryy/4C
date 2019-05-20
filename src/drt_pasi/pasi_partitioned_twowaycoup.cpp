/*---------------------------------------------------------------------------*/
/*!
\brief two way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned_twowaycoup.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FEVector.h>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup::PASI_PartTwoWayCoup(
    const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params)
    : PartitionedAlgo(comm, pasi_params),
      forcenp_(Teuchos::null),
      dispincnp_(Teuchos::null),
      forceincnp_(Teuchos::null),
      itmax_(-1),
      ittol_(-1.0),
      ignoreconvcheck_(false),
      writerestartevery_(-1)
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
  forcenp_ = LINALG::CreateVector(*structure_->DofRowMap(), true);
  dispincnp_ = LINALG::CreateVector(*structure_->DofRowMap(), true);
  forceincnp_ = LINALG::CreateVector(*structure_->DofRowMap(), true);

  // get the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter lists
  const Teuchos::ParameterList& pasi_params = problem->PASIDynamicParams();
  const Teuchos::ParameterList& pasi_params_part =
      problem->PASIDynamicParams().sublist("PARTITIONED");
  const Teuchos::ParameterList& particle_params = problem->ParticleParamsOld();

  // get the parameters for the ConvergenceCheck
  itmax_ = pasi_params.get<int>("ITEMAX");
  ittol_ = pasi_params_part.get<double>("CONVTOL");
  ignoreconvcheck_ = (bool)DRT::INPUT::IntegralValue<int>(pasi_params_part, "IGNORE_CONV_CHECK");

  // write restart every n steps
  writerestartevery_ = pasi_params.get<int>("RESTARTEVRY");

  // safety check: two way coupled pasi not implemented for all Particle Interaction Types
  if (particles_->ParticleInteractionType() != INPAR::PARTICLEOLD::Normal_DEM and
      particles_->ParticleInteractionType() != INPAR::PARTICLEOLD::NormalAndTang_DEM)
  {
    dserror(
        "Two way coupled partitioned PASI not implemented yet for ParticleInteractionType: 'None' "
        "and 'Normal_MD'!");
  }

  // safety check: two way coupled pasi currently just implemented for CentrDiff time integration
  // scheme in particle field
  if (DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::DynamicType>(particle_params, "DYNAMICTYP") !=
      INPAR::PARTICLEOLD::dyna_centrdiff)
  {
    dserror("Two way coupled partitioned PASI just implemented for DYNAMICTYP: 'CentrDiff'");
  }
}

/*---------------------------------------------------------------------------*
 | setup pasi algorithm                                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Setup()
{
  // call base class setup
  PASI::PartitionedAlgo::Setup();

  // construct vector of particle forces on structural field
  Teuchos::RCP<Epetra_FEVector> f_structure =
      Teuchos::rcp(new Epetra_FEVector(*structure_->DofRowMap(), true));

  // set vector of particle forces on structural field in particle field
  particles_->AdapterParticle()->SetFstructure(f_structure);
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PartitionedAlgo::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structure_->Discretization(), restartstep);
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
    // redistribute load in parallel
    particles_->DynamicLoadBalancing();

    // counter and print header
    PrepareTimeStep();

    // iteration loop between coupled fields
    Outerloop();

    // update and output
    UpdateOutput();
  }
}

/*---------------------------------------------------------------------------*
 | iteration loop between coupled fields                      sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    printf("+--------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS                                         |\n");
    printf("+--------------------------------------------------------------------------------+\n");
  }

  while (stopnonliniter == false)
  {
    // increment number of iteration
    itnum++;

    // update the current states in every iteration
    IterUpdateStates(structure_->Dispnp(), forcenp_);

    // set particle forces
    SetParticleForces(forcenp_);

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // reset particle states
    ResetParticleStates();

    // particle time step
    ParticleStep();

    // get particle forces
    GetParticleForces();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);
  }
}

/*---------------------------------------------------------------------------*
 | update and output                                          sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::UpdateOutput()
{
  // output of structure field
  StructOutput();

  // write interface force in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
    structure_->Discretization()->Writer()->WriteVector("forcenp_", forcenp_);

  // output of particle field
  ParticleOutput();
}

/*---------------------------------------------------------------------------*
 | reset particle states                                      sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::ResetParticleStates()
{
  // reset displacements, velocities and accelerations
  particles_->AdapterParticle()->WriteAccessDispnp()->Update(
      1.0, *(particles_->AdapterParticle()->Dispn()), 0.0);
  particles_->AdapterParticle()->WriteAccessVelnp()->Update(
      1.0, *(particles_->AdapterParticle()->Veln()), 0.0);
  particles_->AdapterParticle()->WriteAccessAccnp()->Update(
      1.0, *(particles_->AdapterParticle()->Accn()), 0.0);

  // reset angular velocities and accelerations
  if (particles_->ParticleInteractionType() != INPAR::PARTICLEOLD::None)
  {
    particles_->AdapterParticle()->WriteAccessAngVelnp()->Update(
        1.0, *(particles_->AdapterParticle()->AngVeln()), 0.0);
    particles_->AdapterParticle()->WriteAccessAngAccnp()->Update(
        1.0, *(particles_->AdapterParticle()->AngAccn()), 0.0);
  }

  // clear vector of particle forces on structural discretization
  particles_->AdapterParticle()->WriteAccessFstructure()->PutScalar(0.0);
}

/*---------------------------------------------------------------------------*
 | get particle forces                                        sfuchs 04/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::GetParticleForces()
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::GetParticleForces");

  // get vector of particle forces on structural boundary
  Teuchos::RCP<Epetra_FEVector> f_structure =
      particles_->AdapterParticle()->WriteAccessFstructure();

  if (f_structure == Teuchos::null)
    dserror("GetParticleForces() returned no vector of particle forces on structural boundaries!");

  // call global assemble for vector of particle forces on structural boundary
  const int err = f_structure->GlobalAssemble(Add, false);
  if (err < 0) dserror("global assemble into structforces failed");

  // save vector of particle forces on structural boundary as Epetra_Vector
  forcenp_ = Teuchos::rcp(new Epetra_Vector(Copy, *f_structure, 0));
}

/*---------------------------------------------------------------------------*
 | set particle forces                                        sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup::SetParticleForces(Teuchos::RCP<const Epetra_Vector> forcenp)
{
  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartTwoWayCoup::SetParticleForces");

  double normbdryforce;
  forcenp->Norm2(&normbdryforce);

  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << " Norm of boundary forces:         " << std::setprecision(7) << normbdryforce
              << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }

  // apply particle force on structure discretization
  structure_->ApplyInterfaceForce(structure_->Interface()->ExtractPASICondVector(forcenp));
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
  dispincnp_->Update(1.0, *(structure_->Dispnp()), -1.0);

  // build the L2-norm of the displacement increment and the displacement
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

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
  if (Comm().MyPID() == 0)
  {
    printf(
        "+------------+--------------------+----------------+----------------+-----------------+---"
        "--------------+\n");
    printf(
        "|  step/max  |  tol       [norm]  |  abs-disp-inc  |  rel-disp-inc  |  abs-force-inc  |  "
        "rel-force-inc  |\n");
    printf("|   %3d/%3d  | %10.3E [L_2 ]  |    %10.3E  |    %10.3E  |     %10.3E  |     %10.3E  |",
        itnum, itmax_, ittol_, abs_disp_inc, rel_disp_inc, abs_force_inc, rel_force_inc);
    printf("\n");
    printf(
        "+------------+--------------------+----------------+----------------+-----------------+---"
        "--------------+\n");
  }

  // converged
  if (abs_disp_inc <= ittol_ and rel_disp_inc <= ittol_ and abs_force_inc <= ittol_ and
      rel_force_inc <= ittol_)
  {
    stopnonliniter = true;

    if (Comm().MyPID() == 0)
    {
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                             "
          "                |\n",
          itnum, itmax_);
      printf(
          "+---------------------------------------------------------------------------------------"
          "----------------+\n");
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
      if (Comm().MyPID() == 0)
      {
        printf(
            "|  ATTENTION: Outer Iteration loop not converged in itemax steps!                     "
            "                  |\n");
        printf(
            "+-------------------------------------------------------------------------------------"
            "------------------+\n");
      }
    }
    // abort the simulation
    else
    {
      if (Comm().MyPID() == 0)
      {
        printf(
            "|  STOP: Outer Iteration loop not converged in itemax steps!                          "
            "                  |\n");
        printf(
            "+-------------------------------------------------------------------------------------"
            "------------------+\n");
      }
      dserror("The partitioned PASI solver did not converge in ITEMAX steps!");
    }
  }

  return stopnonliniter;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup_ForceRelax::PASI_PartTwoWayCoup_ForceRelax(
    const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params)
    : PASI_PartTwoWayCoup(comm, pasi_params), omega_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init pasi algorithm                                        sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::Init()
{
  // call base class init
  PASI::PASI_PartTwoWayCoup::Init();

  // get parameter lists
  const Teuchos::ParameterList& pasi_params_part =
      DRT::Problem::Instance()->PASIDynamicParams().sublist("PARTITIONED");

  omega_ = pasi_params_part.get<double>("STARTOMEGA");
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::ReadRestart(int restartstep)
{
  // call base class read restart
  PASI::PartitionedAlgo::ReadRestart(restartstep);

  IO::DiscretizationReader reader(structure_->Discretization(), restartstep);
  if (restartstep != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  // get forcenp_ from restart
  reader.ReadVector(forcenp_, "forcenp_");

  // get omega_ from restart
  omega_ = reader.ReadDouble("omega_");
}

/*---------------------------------------------------------------------------*
 | iteration loop between coupled fields with relaxed forces  sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::Outerloop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    printf("+--------------------------------------------------------------------------------+\n");
    printf("|  ITERATION LOOP BETWEEN COUPLED FIELDS WITH RELAXED FORCES                     |\n");
    printf("+--------------------------------------------------------------------------------+\n");
  }

  // init the relaxed input
  Teuchos::RCP<Epetra_Vector> forcenp = LINALG::CreateVector(*(structure_->DofRowMap()), true);

  forcenp->Update(1.0, *forcenp_, 0.0);

  while (stopnonliniter == false)
  {
    // increment number of iteration
    itnum++;

    // update the states to the last solutions obtained
    IterUpdateStates(structure_->Dispnp(), forcenp);

    // set particle forces
    SetParticleForces(forcenp);

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // reset particle states
    ResetParticleStates();

    // particle time step
    ParticleStep();

    // get particle forces
    GetParticleForces();

    // convergence check for structure and particles fields
    stopnonliniter = ConvergenceCheck(itnum);

    // calculate relaxation parameter
    CalcOmega(omega_, itnum);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    forcenp->Update(omega_, *forceincnp_, 1.0);
  }
}

/*---------------------------------------------------------------------------*
 | update and output                                          sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::UpdateOutput()
{
  // output of structure field
  StructOutput();

  // write interface force and relaxation parameter in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
  {
    structure_->Discretization()->Writer()->WriteVector("forcenp_", forcenp_);
    structure_->Discretization()->Writer()->WriteDouble("omega_", omega_);
  }

  // output of particle field
  ParticleOutput();
}

/*---------------------------------------------------------------------------*
 | calculate relaxation parameter                             sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelax::CalcOmega(double& omega, const int itnum)
{
  // output constant relaxation parameter omega
  if (Comm().MyPID() == 0) std::cout << "Fixed relaxation parameter omega: " << omega << std::endl;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::PASI_PartTwoWayCoup_ForceRelaxAitken(
    const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params)
    : PASI_PartTwoWayCoup_ForceRelax(comm, pasi_params),
      forceincnpold_(Teuchos::null),
      maxomega_(0.0),
      minomega_(0.0)
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
  forceincnpold_ = LINALG::CreateVector(*structure_->DofRowMap(), true);

  // get parameter lists
  const Teuchos::ParameterList& pasi_params_part =
      DRT::Problem::Instance()->PASIDynamicParams().sublist("PARTITIONED");

  maxomega_ = pasi_params_part.get<double>("MAXOMEGA");
  minomega_ = pasi_params_part.get<double>("MINOMEGA");
}

/*---------------------------------------------------------------------------*
 | calculate relaxation parameter                             sfuchs 03/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartTwoWayCoup_ForceRelaxAitken::CalcOmega(double& omega, const int itnum)
{
  // forceincnpdiff =  r^{i+1}_{n+1} - r^i_{n+1}
  Teuchos::RCP<Epetra_Vector> forceincnpdiff = LINALG::CreateVector(*structure_->DofRowMap(), true);
  forceincnpdiff->Update(1.0, *forceincnp_, (-1.0), *forceincnpold_, 0.0);

  double forceincnpdiffnorm = 0.0;
  forceincnpdiff->Norm2(&forceincnpdiffnorm);

  if (forceincnpdiffnorm <= 1e-06 and Comm().MyPID() == 0)
    std::cout << "Warning: The scalar increment is to small in order to use it for Aitken "
                 "relaxation. Using the previous omega instead!"
              << std::endl;

  // calculate dot product
  double forceincsdot = 0.0;  // delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  forceincnpdiff->Dot(*forceincnp_, &forceincsdot);

  // in first iteration reuse omega from previous step
  if (itnum != 1 and forceincnpdiffnorm > 1e-06)
  {
    // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) /
    // |r^{i+1} - r^{i}|^2 results in
    omega = omega * (1.0 - (forceincsdot) /
                               (forceincnpdiffnorm *
                                   forceincnpdiffnorm));  // compare e.g. PhD thesis U. Kuettler

    // we force omega to be in the range defined in the input file
    if (omega < minomega_)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value smaller than MINOMEGA!"
                  << std::endl;
      omega = minomega_;
    }
    if (omega > maxomega_)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value bigger than MAXOMEGA!"
                  << std::endl;
      omega = maxomega_;
    }
  }

  // output Aitken relaxation parameter omega
  if (Comm().MyPID() == 0) std::cout << "Aitken relaxation parameter omega: " << omega << std::endl;

  // update history vector old increment r^i_{n+1}
  forceincnpold_->Update(1.0, *forceincnp_, 0.0);
}
