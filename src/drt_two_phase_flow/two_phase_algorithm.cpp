/*----------------------------------------------------------------------*/
/*!
\file two_phase_algorithm.cpp

\brief Basis of Two Phase algorithms utilizing standard fluid and levelset

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>
*/
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_fluid/fluid_timint_two_phase.H"
#include "../drt_fluid/fluid_timint_two_phase_genalpha.H"
#include "../drt_scatra/scatra_timint_ost.H"

#include "two_phase_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TWOPHASEFLOW::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,"scatra",solverparams)
{
  // time-step length, maximum time and maximum number of steps
  dt_      = prbdyn.get<double>("TIMESTEP");
  maxtime_ = prbdyn.get<double>("MAXTIME");
  stepmax_ = prbdyn.get<int>("NUMSTEP");

  //Output specific criterions
  write_center_of_mass_ = DRT::INPUT::IntegralValue<bool>(prbdyn,"WRITE_CENTER_OF_MASS");

  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_ = prbdyn.get<double>("CONVTOL");
  itmax_ = prbdyn.get<int>("ITEMAX");

  // flag for special flow and start of sampling period from fluid parameter list
  const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
  //const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

  //Values of velocity field are transfered to ScaTra field. This function overwrites the previous initialization in the constructor
  //of ScaTraFluidCouplingAlgorithm(). This is necessary to correctly initialize a particle algorithm.
  SetFluidValuesInScaTra(true);

  // flag for turbulent inflow
  turbinflow_ = DRT::INPUT::IntegralValue<int>(fluiddyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW");
  // number of inflow steps
  numinflowsteps_ = fluiddyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TWOPHASEFLOW::Algorithm::~Algorithm()
{
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

  } // time loop

return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary TPF problem                                    winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("------Stationary-Two-Phase-Flow------  time step ----------------------------------------\n");
  }

  // check if ScaTraField()->initialvelset == true
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  //Give Scatra Values to fluid.
  SetScaTraValuesInFluid();

  // solve nonlinear Navier-Stokes equations
  FluidField().Solve();

  // solve level set equation
  if (Comm().MyPID()==0)
    std::cout << "/!\\ warning === Level-set field not solved for stationary problems" << std::endl;

  // write output to screen and files
  Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField()->PrepareTimeStep();

  // prepare fluid time step, among other things, predict velocity field
  FluidField().PrepareTimeStep();

  // synchronicity check between fluid algorithm and level set algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid time integration differs from time in two phase flow algorithm");
  if (ScaTraField()->Time() != Time())
    dserror("Time in ScaTra time integration differs from time in two phase flow algorithm");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           Time(),maxtime_,dt_,ScaTraField()->MethodTitle().c_str(),Step(),stepmax_);
  }

  // set fluid values required in scatra
  SetFluidValuesInScaTra(false);

  // initially solve scalar transport equation
  // (values for intermediate time steps were calculated at the end of PerpareTimeStep)
  if (Comm().MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
  ScaTraField()->Solve();

  while (stopnonliniter==false)
  {
    itnum++;

    // set scatra values required in fluid
    SetScaTraValuesInFluid();

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) std::cout<<"\n****************************************\n              FLUID SOLVER\n****************************************\n";
    FluidField().Solve();

    // set fluid values required in scatra
    SetFluidValuesInScaTra(false);

    // solve scalar transport equation
    if (Comm().MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
    ScaTraField()->Solve();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::SetScaTraValuesInFluid()
{
  // set level set in fluid field
  // derivatives and discretization based on time-integration scheme

  switch(FluidField().TimIntScheme())
  {
  case INPAR::FLUID::timeint_stationary:
  {
    Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidFieldrcp())->SetScalarFields(ScaTraField()->Phinp(),
                                       Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(),
                                       ScaTraField()->GetGradientAtNodes(),
                                       ScaTraField()->Residual(),
                                       ScaTraField()->Discretization(),
                                       -1);
  }
  break;
  case INPAR::FLUID::timeint_afgenalpha:
  {

    if (Comm().MyPID()==0)
    {
      std::cout << "================================================================================================================" << std::endl;
      std::cout << "||LEVELSET USES ONE_STEP_THETA CURRENTLY. LINEAR INTERPOLATION IS DONE FOR VALUES AT TIME ALPHA_F and ALPHA_M.||" << std::endl;
      std::cout << "================================================================================================================" << std::endl;
    }

    //TimIntParam() provides 1 - alphaf
    const double alphaf = 1.0 - FluidField().TimIntParam();
    const double alpham = Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhaseGenAlpha>(FluidFieldrcp())->AlphaM();
    //Necessary function for the structure of GetNodalCurvature() and GetGradientAtNodes()
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())->SetAlphaF(alphaf); //Set to -1.0 for Gradients and NodalCurv calculations at time (n+1)


    // As The Levelset algorithm does not have support for gen-alpha timeintegration ->
    // OST values from the ScaTra field will have to be used in the fluid field.
    Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidFieldrcp())->SetIterScalarFields(Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())->PhiafOst(alphaf),
                             Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())->PhiamOst(alpham),
                                               Teuchos::rcp_dynamic_cast<SCATRA::LevelSetTimIntOneStepTheta>(ScaTraField())->PhidtamOst(alpham),
                             Teuchos::null, //ScaTraField()->FsPhi()
                             Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(),
                             ScaTraField()->GetGradientAtNodes(),
                             ScaTraField()->Discretization());


//    // As The Levelset algorithm does not have support for gen-alpha timeintegration ->
//    // OST values from the ScaTra field will have to be used in the fluid field.
//    // Here the gen-alpha values have to be alpha_M=alpha_F=1
//    Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidFieldrcp())->SetIterScalarFields(ScaTraField()->Phinp(),
//                                       ScaTraField()->Phinp(),
//                                       ScaTraField()->Phidtnp(),
//                                       ScaTraField()->FsPhi(), //Teuchos::null
//                                       Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(),
//                                       ScaTraField()->GetGradientAtNodes(),
//                                       ScaTraField()->Discretization());

  }
  break;
  case INPAR::FLUID::timeint_one_step_theta:
  {
    Teuchos::rcp_dynamic_cast<FLD::TimIntTwoPhase>(FluidFieldrcp())->SetIterScalarFields(ScaTraField()->Phinp(),
                                   ScaTraField()->Phin(),
                                   ScaTraField()->Phidtnp(),
                                   ScaTraField()->FsPhi(), //Teuchos::null
                                   Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(),
                                   ScaTraField()->GetGradientAtNodes(),
                                   ScaTraField()->Discretization());
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
  // define flags for fluid and scatra convergence check
  bool fluidstopnonliniter  = false;
  bool scatrastopnonliniter = false;

  // fluid convergence check
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n  CONVERGENCE CHECK FOR ITERATION STEP\n****************************************\n";
    std::cout<<"\n****************************************\n              FLUID CHECK\n****************************************\n";
  }
  fluidstopnonliniter  = FluidField().ConvergenceCheck(itnum,itmax_,ittol_);

  //Levelset/ScaTra convergence check.
  if (Comm().MyPID()==0) std::cout<<"\n****************************************\n         SCALAR TRANSPORT CHECK\n****************************************\n";
  scatrastopnonliniter = ScaTraField()->ConvergenceCheck(itnum,itmax_,ittol_);

  if (fluidstopnonliniter == true and scatrastopnonliniter == true) return true;
  else                                                              return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::TimeUpdate()
{
  // update scalar
  ScaTraField()->Update();

  // update fluid
  FluidField().Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TWOPHASEFLOW::Algorithm::Output()
{
  // set scalar and thermodynamic pressure at n+1 and SCATRA true residual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  FluidField().SetScalarFields(ScaTraField()->Phinp(),
                                 0.0,
                                 ScaTraField()->TrueResidual(),
                                 ScaTraField()->Discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();

  ScaTraField()->Output();

  if(write_center_of_mass_)
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
  if (Comm().MyPID()==0)
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
  Teuchos::RCP<Epetra_Vector> oldphidtn  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphin    = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphinp   = Teuchos::null;
  if (restartscatrainput)
  {
    oldphin  = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phin())));
    oldphinp = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())));
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      dserror("Check restart from scatra input for gen alpha!");
      // do we have to keep it?
      oldphidtn= Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phidtn())));
    }
  }

  // for particle level-set method, we cannot call the level-set restart, since there are no particles yet
  // (restart from flow generation via fluid)
  if (!restartscatrainput)
    ScaTraField()->ReadRestart(step);
  else
  {
    // call restart function of base time integrator only
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
      Teuchos::rcp_dynamic_cast<SCATRA::TimIntOneStepTheta>(ScaTraField())->SCATRA::TimIntOneStepTheta::ReadRestart(step);
    else // particles are not yet supported by other time integration schemes
      ScaTraField()->ReadRestart(step);
  }

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> scatradis = ScaTraField()->Discretization();

  // restart of fluid field
  FluidField().ReadRestart(step);

  // read level-set field from input file instead of restart file
  if (restartscatrainput)
  {
    if (Comm().MyPID()==0)
      std::cout << "---  overwriting scalar field with field from input file... " << std::endl;
    // now overwrite restart phis w/ the old phis
    ScaTraField()->Phinp()->Update(1.0, *(oldphinp), 0.0);
    ScaTraField()->Phin() ->Update(1.0, *(oldphin),  0.0);
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      ScaTraField()->Phidtn() ->Update(1.0, *(oldphidtn),  0.0);
      ScaTraField()->ComputeIntermediateValues();
    }
    else if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
    {
      ScaTraField()->Phidtn()->PutScalar(0.0);
      // calls CalcInitialPhidt()
      // ApplyDirichletBC() called within this function doesn't pose any problem since we
      // don't have any DBCs for level-set problems
      ScaTraField()->PrepareFirstTimeStep();
      // CalcInitialPhidt() copies phidtn to phidtnp
    }
    else
      dserror("Unknown time-integration scheme!");
    if (Comm().MyPID()==0)
      std::cout << "done" << std::endl;

    // Gmsh post processing files, can be added for debugging.
    // See combust_algorithm.cpp (COMBUST::Algorithm::Restart(*))
  }

  // set time in scalar transport time integration scheme
  // this has also to be done for restartscatrainput, since a particle field may now be included
  if(restartscatrainput)
    ScaTraField()->SetTimeStep(FluidField().Time(),step);

  SetTimeStep(FluidField().Time(),step);

  //TODO: Give TPF the capability to use generating flow  as seen in COMBUST module.
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
  FluidField().SetScalarFields(ScaTraField()->Phinp(),
                                 0.0,
                                 Teuchos::null,
                                 ScaTraField()->Discretization());
  FluidField().ReadRestart(restart);
  // as ReadRestart is only called for the FluidField
  // time and step have not been set in the superior class and the ScaTraField
  SetTimeStep(FluidField().Time(),FluidField().Step());
  ScaTraField()->SetTimeStep(FluidField().Time(),FluidField().Step());
  return;
}


void TWOPHASEFLOW::Algorithm::SetFluidValuesInScaTra(bool init)
{
  // get convel at the correct time
  const Teuchos::RCP<const Epetra_Vector> convel = (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
                                                       ?(FluidField().Velaf())
                                                       :(FluidField().Velnp());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->SetVelocityField(convel,
        Teuchos::null,
        Teuchos::null,
        FluidField().FsVel(),
        FluidField().Discretization()->GetDofSetProxy(),
        //FluidField().DofSet(),
        FluidField().Discretization(),
        init);
  else // temporary solution, since level-set algorithm does not yet support gen-alpha
  {
    if (Comm().MyPID()==0)
      std::cout << "CORRECT THESE LINES WHEN LEVEL-SET ALOGITHM SUPPORTS GEN-ALPHA" << std::endl;

    ScaTraField()->SetVelocityField(convel,
        Teuchos::null,
        Teuchos::null,
        FluidField().FsVel(),
        FluidField().Discretization()->GetDofSetProxy(),
//        FluidField().DofSet(),
        FluidField().Discretization());
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
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
      FluidField().Output();

    // output Levelset function initial state
    if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
      ScaTraField()->Output();
  }

  if(write_center_of_mass_)
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
  DRT::Problem::Instance()->AddFieldTest(FluidField().CreateFieldTest());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
  {
    // DRT::Problem::Instance()->TestAll() is called in level-set field after adding particles
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
  }
  else
  {
    DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(Comm());
  }

  return;
}
