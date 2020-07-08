/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of Two Phase algorithms utilizing standard fluid and levelset

\level 3

            089/28915236

*/
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_fluid/fluid_timint_two_phase.H"
#include "../drt_fluid/fluid_timint_two_phase_genalpha.H"
#include "../drt_fluid/fluid_timint_two_phase_ost.H"
#include "../drt_scatra/scatra_timint_ost.H"

#include "two_phase_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TWOPHASEFLOW::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams)
    : ScaTraFluidCouplingAlgorithm(comm, prbdyn, false, "scatra", solverparams),
      dt_(0.0),
      maxtime_(0.0),
      stepmax_(0),
      itmax_(0),
      ittol_(1.0),
      write_center_of_mass_(false),
      turbinflow_(false),
      numinflowsteps_(0),
      velnpi_(Teuchos::null),
      phinpi_(Teuchos::null),
      prbdyn_(prbdyn),
      smoothedgradphitype_(INPAR::TWOPHASE::smooth_grad_phi_l2_projection),
      scalesmoothedgradients_(false)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TWOPHASEFLOW::Algorithm::~Algorithm() { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  // call Init() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(prbdyn, scatradyn, solverparams, disname, isale);

  // time-step length, maximum time and maximum number of steps
  dt_ = prbdyn_.get<double>("TIMESTEP");
  maxtime_ = prbdyn_.get<double>("MAXTIME");
  stepmax_ = prbdyn_.get<int>("NUMSTEP");

  // Output specific criterions
  write_center_of_mass_ = DRT::INPUT::IntegralValue<bool>(prbdyn_, "WRITE_CENTER_OF_MASS");

  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_ = prbdyn_.get<double>("CONVTOL");
  itmax_ = prbdyn_.get<int>("ITEMAX");

  // flag for special flow and start of sampling period from fluid parameter list
  const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();

  // flag for turbulent inflow
  turbinflow_ =
      DRT::INPUT::IntegralValue<int>(fluiddyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW");
  // number of inflow steps
  numinflowsteps_ = fluiddyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP");

  smoothedgradphitype_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SmoothGradPhi>(
      prbdyn.sublist("SURFACE TENSION"), "SMOOTHGRADPHI");
  scalesmoothedgradients_ = DRT::INPUT::IntegralValue<bool>(
      prbdyn.sublist("SURFACE TENSION"), "SCALE_SMOOTHED_GRADIENTS");

  // Instantiate vectors contatining outer loop increment data
  fsvelincnorm_.reserve(itmax_);
  fspressincnorm_.reserve(itmax_);
  fsphiincnorm_.reserve(itmax_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::Setup()
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  // Fluid-Scatra Iteration vectors are initialized
  velnpi_ = Teuchos::rcp(new Epetra_Vector(FluidField()->Velnp()->Map()), true);
  velnpi_->Update(1.0, *FluidField()->Velnp(), 0.0);
  phinpi_ = Teuchos::rcp(new Epetra_Vector(ScaTraField()->Phinp()->Map()), true);
  phinpi_->Update(1.0, *ScaTraField()->Phinp(), 0.0);

  // Values of velocity field are transferred to ScaTra field.
  // This function overwrites the previous initialization in the
  // constructor of ScaTraFluidCouplingAlgorithm().
  // This is necessary to correctly initialize a particle algorithm.
  // old particle framework removed -> todo: requires clean up
  SetFluidValuesInScaTra(true);

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::TimeLoop()
{
  OutputInitialField();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop for particular type of algorithm
    OuterLoop();

    // update for next time step
    TimeUpdate();

    // write output to screen and files
    Output();

  }  // time loop

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary TPF problem                                    winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID() == 0)
  {
    printf(
        "------Stationary-Two-Phase-Flow------  time step "
        "----------------------------------------\n");
  }

  // check if ScaTraField()->initialvelset == true
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // Give Scatra Values to fluid.
  SetScaTraValuesInFluid();

  // solve nonlinear Navier-Stokes equations
  FluidField()->Solve();

  // solve level set equation
  if (Comm().MyPID() == 0)
    std::cout << "/!\\ warning === Level-set field not solved for stationary problems" << std::endl;

  // write output to screen and files
  Output();

  return;
}


/*---------------------------------------------------------------------------------------*
| Prepares values and variables needed in the outer iteration                            |
*----------------------------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::PrepareOuterIteration()
{
  // Update phi for outer loop convergence check
  phinpi_->Update(1.0, *ScaTraField()->Phinp(), 0.0);
  velnpi_->Update(1.0, *FluidField()->Velnp(), 0.0);

  // Clear the vectors containing the data for the partitioned increments
  fsvelincnorm_.clear();
  fspressincnorm_.clear();
  fsphiincnorm_.clear();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField()->PrepareTimeStep();

  // set scatra values required in fluid
  SetScaTraValuesInFluid();

  // prepare fluid time step, among other things, predict velocity field
  FluidField()->PrepareTimeStep();

  // synchronicity check between fluid algorithm and level set algorithms
  if (FluidField()->Time() != Time())
    dserror("Time in Fluid time integration differs from time in two phase flow algorithm");
  if (ScaTraField()->Time() != Time())
    dserror("Time in ScaTra time integration differs from time in two phase flow algorithm");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::OuterLoop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n", Time(), maxtime_, dt_,
        ScaTraField()->MethodTitle().c_str(), Step(), stepmax_);
  }

  // initially solve scalar transport equation
  // (values for intermediate time steps were calculated at the end of PerpareTimeStep)
  DoScaTraField();

  // Prepare variables for convergence check.
  PrepareOuterIteration();

  while (stopnonliniter == false)
  {
    itnum++;

    DoFluidField();

    DoScaTraField();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::DoFluidField()
{
  if (Comm().MyPID() == 0)
    std::cout << "\n****************************************\n              FLUID "
                 "SOLVER\n****************************************\n";


  // Set relevant ScaTra values in Fluid field.
  SetScaTraValuesInFluid();

  // Print time step information
  Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidField())->PrintTimeStepInfo();

  // Solve the Fluid field.
  FluidField()->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::DoScaTraField()
{
  if (Comm().MyPID() == 0)
    std::cout << "\n****************************************\n        SCALAR TRANSPORT "
                 "SOLVER\n****************************************\n";


  // Set relevant Fluid values in ScaTra field.
  SetFluidValuesInScaTra(false);

  // Solve the ScaTra field.
  ScaTraField()->Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::SetScaTraValuesInFluid()
{
  // set level set in fluid field
  // derivatives and discretization based on time-integration scheme

  // check assumption
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_one_step_theta and
      ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("OST for level-set field assumed!");

  switch (FluidField()->TimIntScheme())
  {
    case INPAR::FLUID::timeint_stationary:
    {
      Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidField())
          ->SetScalarFields(ScaTraField()->Phinp(),
              Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())
                  ->GetNodalCurvature(
                      ScaTraField()->Phinp(), GetSmoothedLevelSetGradient(ScaTraField()->Phinp())),
              GetSmoothedLevelSetGradient(ScaTraField()->Phinp()), ScaTraField()->Residual(),
              ScaTraField()->Discretization(), -1);
    }
    break;
    case INPAR::FLUID::timeint_afgenalpha:
    {
      if (Comm().MyPID() == 0)
      {
        std::cout << "============================================================================="
                     "==================================="
                  << std::endl;
        std::cout << "||LEVELSET USES ONE_STEP_THETA CURRENTLY. LINEAR INTERPOLATION IS DONE FOR "
                     "VALUES AT TIME ALPHA_F and ALPHA_M.||"
                  << std::endl;
        std::cout << "============================================================================="
                     "==================================="
                  << std::endl;
      }

      // TimIntParam() provides 1 - alphaf
      const double alphaf = 1.0 - FluidField()->TimIntParam();
      const double alpham =
          Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhaseGenAlpha>(FluidField())->AlphaM();
      // get phi at intermediate time level n+alpha_F
      const Teuchos::RCP<const Epetra_Vector> phiaf =
          Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())
              ->Phinptheta(alphaf);

      // As The Levelset algorithm does not have support for gen-alpha timeintegration ->
      // OST values from the ScaTra field will have to be used in the fluid field.
      Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidField())
          ->SetIterScalarFields(phiaf,
              Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())
                  ->Phinptheta(alpham),
              Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())
                  ->Phidtnptheta(alpham),
              Teuchos::null,  // ScaTraField()->FsPhi()
              Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())
                  ->GetNodalCurvature(phiaf, GetSmoothedLevelSetGradient(phiaf)),
              GetSmoothedLevelSetGradient(phiaf), ScaTraField()->Discretization());


      //    // As The Levelset algorithm does not have support for gen-alpha timeintegration ->
      //    // OST values from the ScaTra field will have to be used in the fluid field.
      //    // Here the gen-alpha values have to be alpha_M=alpha_F=1
      //    Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidField())->SetIterScalarFields(ScaTraField()->Phinp(),
      //                                       ScaTraField()->Phinp(),
      //                                       ScaTraField()->Phidtnp(),
      //                                       ScaTraField()->FsPhi(), //Teuchos::null
      //                                       Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(),
      //                                       ScaTraField()->ReconstructGradientAtNodes(),
      //                                       ScaTraField()->Discretization());
    }
    break;
    case INPAR::FLUID::timeint_one_step_theta:
    {
      Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidField())
          ->SetIterScalarFields(ScaTraField()->Phinp(), ScaTraField()->Phin(),
              ScaTraField()->Phidtnp(), ScaTraField()->FsPhi(),
              Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())
                  ->GetNodalCurvature(
                      ScaTraField()->Phinp(), GetSmoothedLevelSetGradient(ScaTraField()->Phinp())),
              GetSmoothedLevelSetGradient(ScaTraField()->Phinp()), ScaTraField()->Discretization());

      Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhaseOst>(FluidField())
          ->SetIterScalarFieldsn(Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())
                                     ->GetNodalCurvature(ScaTraField()->Phin(),
                                         GetSmoothedLevelSetGradient(ScaTraField()->Phin())),
              GetSmoothedLevelSetGradient(ScaTraField()->Phin()), ScaTraField()->Discretization());
    }
    break;
    default:
      dserror("Time integration scheme not supported");
      break;
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool TWOPHASEFLOW::Algorithm::ConvergenceCheck(int itnum)
{
  // define flags for Fluid and ScaTra convergence check
  bool fluidstopnonliniter = false;
  bool scatrastopnonliniter = false;

  bool notconverged = false;

  if (itmax_ <= 0) dserror("Set iterations to something reasonable!!!");

  double fsvelincnorm = 1.0;
  double fspressincnorm = 1.0;
  double fsphiincnorm = 1.0;
  // Get increment for outer loop of Fluid and ScaTra
  GetOuterLoopIncFluid(fsvelincnorm, fspressincnorm, itnum);
  GetOuterLoopIncScaTra(fsphiincnorm, itnum);

  fsvelincnorm_[itnum - 1] = fsvelincnorm;
  fspressincnorm_[itnum - 1] = fspressincnorm;
  fsphiincnorm_[itnum - 1] = fsphiincnorm;

  if (Comm().MyPID() == 0)
  {
    printf(
        "\n|+------ TWO PHASE FLOW CONVERGENCE CHECK:  time step %2d, outer iteration %2d ------+|",
        Step(), itnum);
    printf(
        "\n|- iter/itermax -|----tol-[Norm]---|-- fluid-inc --|-- press inc --|-- levset inc --|");
  }

  for (int k_itnum = 0; k_itnum < itnum; k_itnum++)
  {
    if (k_itnum == 0)
    {
      if (Comm().MyPID() == 0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |       -       |       -       |   %10.3E   |",
            (k_itnum + 1), itmax_, ittol_, fsphiincnorm_[k_itnum]);
      }  // end if processor 0 for output
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |  %10.3E   |  %10.3E   |   %10.3E   |",
            (k_itnum + 1), itmax_, ittol_, fsvelincnorm_[k_itnum], fspressincnorm_[k_itnum],
            fsphiincnorm_[k_itnum]);
      }  // end if processor 0 for output
    }
  }
  if (Comm().MyPID() == 0)
    printf(
        "\n|+---------------------------------------------------------------------------------+|"
        "\n");

  if ((fsvelincnorm <= ittol_) and (fspressincnorm <= ittol_) and itnum > 1)
    fluidstopnonliniter = true;

  if ((fsphiincnorm <= ittol_)) scatrastopnonliniter = true;


  // If tolerance or number of maximum iterations are reached
  if ((fluidstopnonliniter and scatrastopnonliniter) or (itnum >= itmax_))
  {
    notconverged = true;
  }

  if (Comm().MyPID() == 0)
  {
    if ((itnum == stepmax_) and (notconverged == true))
    {
      printf("|+---------------- not converged ----------------------+|");
      printf("\n|+-----------------------------------------------------+|\n");
    }
  }

  return notconverged;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::TimeUpdate()
{
  // update scalar
  ScaTraField()->Update();

  // update fluid
  FluidField()->Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::Output()
{
  // set scalar and thermodynamic pressure at n+1 and SCATRA true residual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  FluidField()->SetScalarFields(
      ScaTraField()->Phinp(), 0.0, ScaTraField()->TrueResidual(), ScaTraField()->Discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField()->StatisticsAndOutput();

  ScaTraField()->Output();

  // TODO: This is not the best way to include additional scatra output
  //       Moving this out entirely to the level-set algorithm and combining it with
  //       the respective function for sharp interfaces, which should move from the old
  //       combustion module to the new level-set algorithm, seem to be a smarter solution!
  if (write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/* -------------------------------------------------------------------------------*
 | Restart a two phase problem                                              winter|
 | remark: Levelset is solved before fluid                                        |
 |         switching the order would affect the restart                           |
 * -------------------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::Restart(int step, const bool restartscatrainput)
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "##################################################" << std::endl;
    std::cout << "#                                                #" << std::endl;
    std::cout << "#     Restart of T(wo) P(hase) F(low) problem    #" << std::endl;
    if (restartscatrainput)
    {
      std::cout << "#                                                #" << std::endl;
      std::cout << "#   -Restart with scalar field from input file   #" << std::endl;
    }
    std::cout << "#                                                #" << std::endl;
    std::cout << "##################################################" << std::endl;
  }

  // read level-set field from input file instead of restart file
  Teuchos::RCP<Epetra_Vector> oldphidtn = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphin = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphinp = Teuchos::null;
  if (restartscatrainput)
  {
    oldphin = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phin())));
    oldphinp = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())));
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      dserror("Check restart from scatra input for gen alpha!");
      // do we have to keep it?
      oldphidtn = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phidtn())));
    }
  }

  // for particle level-set method, we cannot call the level-set restart, since there are no
  // particles yet (restart from flow generation via fluid)
  // old particle framework removed -> todo: requires clean up
  if (!restartscatrainput)
    ScaTraField()->ReadRestart(step);
  else
  {
    // call restart function of base time integrator only
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
      Teuchos::rcp_dynamic_cast<SCATRA::TimIntOneStepTheta>(ScaTraField())
          ->SCATRA::TimIntOneStepTheta::ReadRestart(step);
    else  // particles are not yet supported by other time integration schemes
      ScaTraField()->ReadRestart(step);
  }

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> scatradis = ScaTraField()->Discretization();

  // restart of fluid field
  FluidField()->ReadRestart(step);

  // read level-set field from input file instead of restart file
  if (restartscatrainput)
  {
    if (Comm().MyPID() == 0)
      std::cout << "---  overwriting scalar field with field from input file... " << std::endl;
    // now overwrite restart phis w/ the old phis
    ScaTraField()->Phinp()->Update(1.0, *(oldphinp), 0.0);
    ScaTraField()->Phin()->Update(1.0, *(oldphin), 0.0);
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      ScaTraField()->Phidtn()->Update(1.0, *(oldphidtn), 0.0);
      ScaTraField()->ComputeIntermediateValues();
    }
    else if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
    {
      ScaTraField()->Phidtnp()->PutScalar(0.0);
      // calls CalcInitialTimeDerivative()
      // ApplyDirichletBC() called within this function doesn't pose any problem since we
      // don't have any DBCs for level-set problems
      ScaTraField()->PrepareFirstTimeStep();
      // CalcInitialTimeDerivative() copies phidtnp to phidtn
    }
    else
      dserror("Unknown time-integration scheme!");
    if (Comm().MyPID() == 0) std::cout << "done" << std::endl;
  }

  // set time in scalar transport time integration scheme
  // this has also to be done for restartscatrainput, since a particle field may now be included
  // old particle framework removed -> todo: requires clean up
  if (restartscatrainput) ScaTraField()->SetTimeStep(FluidField()->Time(), step);

  SetTimeStep(FluidField()->Time(), step);

  // TODO: Give TPF the capability to use generating flow  as seen in COMBUST module.
  // assign the fluid velocity field to the levelset function as convective velocity field
  // bool true allows for setting old convective velocity required for particle coupling
  //  if (not gen_flow_)
  //    SetVelocityLevelSet(true);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::ReadInflowRestart(int restart)
{
  // in case a inflow generation in the inflow section has been performed,
  // there are not any scatra results available and the initial field is used
  // caution: if AVM3Preparation is called ,e.g., for multifractal subgrid-scale
  //          modeling the physical parameters (dens, visc, diff) are required
  //          to obtain non-zero values which otherwise cause troubles when dividing by them
  // set initial scalar field
  FluidField()->SetScalarFields(
      ScaTraField()->Phinp(), 0.0, Teuchos::null, ScaTraField()->Discretization());
  FluidField()->ReadRestart(restart);
  // as ReadRestart is only called for the FluidField
  // time and step have not been set in the superior class and the ScaTraField
  SetTimeStep(FluidField()->Time(), FluidField()->Step());
  ScaTraField()->SetTimeStep(FluidField()->Time(), FluidField()->Step());
  return;
}


void TWOPHASEFLOW::Algorithm::SetFluidValuesInScaTra(bool init)
{
  // get convel at the correct time
  const Teuchos::RCP<const Epetra_Vector> convel =
      (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha) ? (FluidField()->Velaf())
                                                                        : (FluidField()->Velnp());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())
        ->SetVelocityField(
            convel, Teuchos::null, Teuchos::null, FluidField()->FsVel(), false, init);
  else  // temporary solution, since level-set algorithm does not yet support gen-alpha
  {
    if (Comm().MyPID() == 0)
      std::cout << "CORRECT THESE LINES WHEN LEVEL-SET ALOGITHM SUPPORTS GEN-ALPHA" << std::endl;

    ScaTraField()->SetVelocityField(convel, Teuchos::null, Teuchos::null, FluidField()->FsVel(), 1);
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output of initial field                                                winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::OutputInitialField()
{
  if (Step() == 0)
  {
    // output fluid initial state
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary) FluidField()->Output();

    // output Levelset function initial state
    if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary) ScaTraField()->Output();
  }

  if (write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform result test                                     winter 06/14 |
 *----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::TestResults()
{
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
  {
    // DRT::Problem::Instance()->TestAll() is called in level-set field
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
  }
  else
  {
    DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(Comm());
  }

  return;
}

void TWOPHASEFLOW::Algorithm::GetOuterLoopIncFluid(
    double& fsvelincnorm, double& fspressincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> velnpip = FluidField()->Velnp();  // Contains Fluid and Pressure

  Teuchos::RCP<const Epetra_Vector> onlyvel = FluidField()->ExtractVelocityPart(velnpip);
  Teuchos::RCP<const Epetra_Vector> onlypress = FluidField()->ExtractPressurePart(velnpip);

  Teuchos::RCP<Epetra_Vector> onlyveli = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()), true);
  Teuchos::RCP<Epetra_Vector> onlypressi = Teuchos::rcp(new Epetra_Vector(onlypress->Map()), true);

  if (itnum > 1)
  {
    onlyveli->Update(1.0, *(FluidField()->ExtractVelocityPart(velnpi_)), 0.0);
    onlypressi->Update(1.0, *(FluidField()->ExtractPressurePart(velnpi_)), 0.0);
    //    onlyveli   = FluidField()->ExtractVelocityPart(velnpi_);
    //    onlypressi = FluidField()->ExtractPressurePart(velnpi_);
  }

  double velnormL2 = 1.0;
  double pressnormL2 = 1.0;

  onlyvel->Norm2(&velnormL2);
  onlypress->Norm2(&pressnormL2);

  if (velnormL2 < 1e-5) velnormL2 = 1.0;
  if (pressnormL2 < 1e-5) pressnormL2 = 1.0;

  double fsvelnormL2 = 1.0;
  double fspressnormL2 = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()), true);
  incvel->Update(1.0, *onlyvel, -1.0, *onlyveli, 0.0);
  incvel->Norm2(&fsvelnormL2);

  Teuchos::RCP<Epetra_Vector> incpress = Teuchos::rcp(new Epetra_Vector(onlypress->Map()), true);
  incpress->Update(1.0, *onlypress, -1.0, *onlypressi, 0.0);
  incpress->Norm2(&fspressnormL2);
  //-----------------------------------------------------

  fsvelincnorm = fsvelnormL2 / velnormL2;
  fspressincnorm = fspressnormL2 / pressnormL2;

#if DEBUG
  //-------------------------
  std::cout << "fsvelnormL2: " << fsvelnormL2 << std::endl;
  std::cout << "velnormL2: " << velnormL2 << std::endl << std::endl;

  std::cout << "fspressnormL2: " << fspressnormL2 << std::endl;
  std::cout << "pressnormL2: " << pressnormL2 << std::endl << std::endl;
//-------------------------
#endif

  velnpi_->Update(1.0, *velnpip, 0.0);
}

void TWOPHASEFLOW::Algorithm::GetOuterLoopIncScaTra(double& fsphiincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> phinpip = ScaTraField()->Phinp();

  double phinormL2 = 1.0;

  phinpip->Norm2(&phinormL2);
  if (phinormL2 < 1e-5) phinormL2 = 1.0;
  double fsphinormL2 = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incphi = Teuchos::rcp(new Epetra_Vector(phinpip->Map(), true));
  incphi->Update(1.0, *phinpip, -1.0, *phinpi_, 0.0);
  incphi->Norm2(&fsphinormL2);
  //-----------------------------------------------------

  fsphiincnorm = fsphinormL2 / phinormL2;

#if DEBUG
  //-------------------------
  std::cout << "fsphinormL2: " << fsphinormL2 << std::endl;
  std::cout << "phinormL2: " << phinormL2 << std::endl << std::endl;
//-------------------------
#endif


  phinpi_->Update(1.0, *phinpip, 0.0);
}

/*----------------------------------------------------------------------*
 | Set relevant values from ScaTra field in the Fluid field.            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> TWOPHASEFLOW::Algorithm::GetSmoothedLevelSetGradient(
    const Teuchos::RCP<const Epetra_Vector>& phinp)
{
  bool normalize_gradients = scalesmoothedgradients_;

  switch (smoothedgradphitype_)
  {
    case INPAR::TWOPHASE::smooth_grad_phi_l2_projection:
    {
      return ScaTraField()->ReconstructGradientAtNodesL2Projection(phinp, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_superconvergent_patch_recovery_3D:
    {
      return ScaTraField()->ReconstructGradientAtNodesPatchRecon(phinp, 3, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_superconvergent_patch_recovery_2Dz:
    {
      return ScaTraField()->ReconstructGradientAtNodesPatchRecon(phinp, 2, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_meanvalue:
    {
      return ScaTraField()->ReconstructGradientAtNodesMeanAverage(phinp, normalize_gradients);
      break;
    }
    default:
      dserror("The chosen smoothing is not as of yet supported!");
      break;
  }

  return Teuchos::null;
}
