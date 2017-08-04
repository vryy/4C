/*----------------------------------------------------------------------*/
/*!
\file scatra_algorithm.cpp

\brief Transport of passive scalars in Navier-Stokes velocity field

\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/

#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "scatra_timint_implicit.H"
#include "scatra_algorithm.H"
#include "../drt_adapter/adapter_coupling_volmortar.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraAlgorithm::ScaTraAlgorithm(
  const Epetra_Comm& comm,                     ///< communicator
  const Teuchos::ParameterList& scatradyn,     ///< scatra parameter list
  const Teuchos::ParameterList& fdyn,          ///< fluid parameter list
  const std::string scatra_disname,                   ///< scatra discretization name
  const Teuchos::ParameterList& solverparams   ///< solver parameter list
)
  :  ScaTraFluidCouplingAlgorithm(comm,scatradyn,false,scatra_disname,solverparams),
     natconv_(DRT::INPUT::IntegralValue<int>(scatradyn,"NATURAL_CONVECTION")),
     natconvitmax_(scatradyn.sublist("NONLINEAR").get<int>("ITEMAX_OUTER")),
     natconvittol_(scatradyn.sublist("NONLINEAR").get<double>("CONVTOL_OUTER")),
     velincnp_ (Teuchos::null),
     phiincnp_ (Teuchos::null),
     samstart_(fdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_START")),
     samstop_(fdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_STOP"))
{
  // no stuff to add here at the moment
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraAlgorithm::~ScaTraAlgorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Setup()
{
  // call setup in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  // create vectors
  velincnp_ = Teuchos::rcp(new Epetra_Vector(*(FluidField()->ExtractVelocityPart(FluidField()->Velnp()))));
  phiincnp_ = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())));

  if(velincnp_==Teuchos::null)
    dserror("velincnp_ == Teuchos::null");
  if(phiincnp_==Teuchos::null)
      dserror("phiincnp_ == Teuchos::null");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Init(
    const Teuchos::ParameterList&   prbdyn,         ///< parameter list for global problem
    const Teuchos::ParameterList&   scatradyn,      ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList&   solverparams,   ///< parameter list for scalar transport solver
    const std::string&              disname,        ///< name of scalar transport discretization
    const bool                      isale           ///< ALE flag
)
{
  // call init in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(
      prbdyn,
      scatradyn,
      solverparams,
      disname,
      isale);

  return;
}


/*----------------------------------------------------------------------*
 | Provide required time loop                                fang 07/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::TimeLoop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // provide information about initial field
  PrepareTimeLoop();

  if (natconv_ == false)
    TimeLoopOneWay();   // one-way coupling
  else
    TimeLoopTwoWay();   // two-way coupling (natural convection)

  return;
}


/*----------------------------------------------------------------------*
 | Time loop for one-way coupled problems                    fang 07/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::TimeLoopOneWay()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // write velocities since they may be needed in PrepareFirstTimeStep() of the scatra field
  SetVelocityField();
  // write initial output to file
  if (Step()==0)
    Output();

  // time loop (no-subcycling at the moment)
  while(NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // solve nonlinear Navier-Stokes system
    DoFluidStep();

    // solve transport equations
    DoTransportStep();

    // update all field solvers
    Update();

    // compute error for problems with analytical solution
    ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

    // write output to screen and files
    Output();
  } // while(NotFinished())

  return;
}


/*--------------------------------------------------------------------------*
 | Time loop for two-way coupled problems (natural convection)   fang 07/14 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::TimeLoopTwoWay()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  PrepareTimeLoopTwoWay();

  // time loop
  while(NotFinished())
  {
    IncrementTimeAndStep();

    // prepare next time step
    PrepareTimeStepConvection();

    // Outer iteration loop
    OuterIterationConvection();

    // update all single field solvers and density
    UpdateConvection();

    // write output to screen and files
    Output();
  } // while(NotFinished())

  return;
}


/*----------------------------------------------------------------------*
 | Initial calculations for two-way coupling time loop       fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::PrepareTimeLoopTwoWay()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // safety check
  switch((FluidField()->TimIntScheme()))
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_bdf2:
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_stationary:
    break;
  default:
  {
    dserror("Selected time integration scheme is not available!");
    break;
  }
  }

  // compute initial mean concentrations and load densification coefficients
  ScaTraField()->SetupNatConv();

  // compute initial density
  ScaTraField()->ComputeDensity();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::PrepareTimeStep()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  IncrementTimeAndStep();

  FluidField()->PrepareTimeStep();

  // prepare time step
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField()->PrepareTimeStep();

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n******************\n   TIME STEP     \n******************\n";
    std::cout<<"\nStep:   " << Step() << " / " << NStep() << "\n";
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::PrepareTimeStepConvection()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // OST/BDF2 time integration schemes are implemented
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST and BDF2 (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to 0.0

  switch((FluidField()->TimIntScheme()))
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_bdf2:
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_stationary:
  {
    FluidField()->SetIterScalarFields(
        ScaTraField()->Densafnp(),
        ScaTraField()->Densafnp(), // not needed, provided as dummy vector
        Teuchos::null,
        ScaTraField()->Discretization());
    break;
  }
  default:
  {
    dserror("Selected time integration scheme is not available!");
    break;
  }
  }

  FluidField()->PrepareTimeStep();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step()==1)
    ScaTraField()->SetVelocityField(
        FluidField()->Velnp(),
        FluidField()->Hist(),
        Teuchos::null,
        Teuchos::null,
        1);

  // prepare time step (+ initialize one-step-theta scheme correctly with
  // velocity given above)
  ScaTraField()->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*
 | Print scatra solver type to screen                        fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::PrintScaTraSolver()
{
  if (Comm().MyPID()==0)
    std::cout<<"\n****************************\n      TRANSPORT SOLVER\n****************************\n";

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::DoFluidStep()
{
  // solve nonlinear Navier-Stokes system
  if (Comm().MyPID()==0)
    std::cout<<"\n************************\n      FLUID SOLVER\n************************\n";

  // currently only required for forced homogeneous isotropic turbulence with
  // passive scalar transport; does nothing otherwise
  FluidField()->CalcIntermediateSolution();

  FluidField()->Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::DoTransportStep()
{
  PrintScaTraSolver();

  // transfer velocities to scalar transport field solver
  SetVelocityField();

  // solve the transport equation(s)
  ScaTraField()->Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::SetVelocityField()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // transfer velocities to scalar transport field solver
  // NOTE: so far, the convective velocity is chosen to equal the fluid velocity
  //       since it is not yet clear how the grid velocity should be interpolated
  //       properly -> hence, ScaTraAlgorithm does not support moving
  //       meshes yet

  //this is ugly, but FsVel() may give a Null pointer which we canNOT give to the volmortart framework
  //TODO (thon): make this somehow prettier..
  Teuchos::RCP<const Epetra_Vector> fsvel = FluidField()->FsVel();
  if (fsvel!=Teuchos::null)
    fsvel=FluidToScatra(fsvel);

  switch(FluidField()->TimIntScheme())
  {
  case INPAR::FLUID::timeint_npgenalpha:
  case INPAR::FLUID::timeint_afgenalpha:
  {
    ScaTraField()->SetVelocityField(
        FluidToScatra( FluidField()->Velaf() ),
        FluidToScatra( FluidField()->Accam() ),
        FluidToScatra( FluidField()->Velaf() ),
        fsvel,
        1);
    break;
  }
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_bdf2:
  case INPAR::FLUID::timeint_stationary:
  {
    ScaTraField()->SetVelocityField(
        FluidToScatra( FluidField()->Velnp() ),
        FluidToScatra( FluidField()->Hist() ),
        FluidToScatra( FluidField()->Velnp() ),
        fsvel,
        1);
    break;
  }
  default:
  {
    dserror("Time integration scheme not supported");
    break;
  }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::OuterIterationConvection()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  int  natconvitnum = 0;
  bool stopnonliniter = false;

  // Outer Iteration loop starts
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n";
    std::cout<<"**************************************************************\n";
    std::cout<<"      OUTER ITERATION LOOP ("<<ScaTraField()->MethodTitle()<<")\n";
    printf("      Time Step %3d/%3d \n",Step(), ScaTraField()->NStep());
    std::cout<<"**************************************************************\n";
  }

#ifdef OUTPUT
  // Output after each Outer Iteration step
  const int numdim = 3;
  //create output file name
  std::stringstream temp;
  temp << DRT::Problem::Instance()->OutputControlFile()->FileName()<<"_nonliniter_step"<<Step();
  std::string outname = temp.str();
  std::string probtype = DRT::Problem::Instance()->ProblemName();

  Teuchos::RCP<IO::OutputControl> myoutputcontrol = Teuchos::rcp(new IO::OutputControl(ScaTraField().Discretization()->Comm(),probtype,"Polynomial","myinput",outname,numdim,0,1000));
  // create discretization writer with my own control settings
  Teuchos::RCP<IO::DiscretizationWriter> myoutput = ScaTraField().Discretization()->Writer();
  myoutput->SetOutput(myoutputcontrol);
  // write mesh at step 0
  myoutput->WriteMesh(0,0.0);
#endif

  while (stopnonliniter==false)
  {
    natconvitnum ++;

    phiincnp_->Update(1.0,*ScaTraField()->Phinp(),0.0);
    velincnp_->Update(1.0,*FluidField()->ExtractVelocityPart(FluidField()->Velnp()),0.0);

    // solve nonlinear Navier-Stokes system with body forces
    DoFluidStep();

    // solve scalar transport equation
    DoTransportStep();

    // compute new density field and pass it to the fluid discretization
    ScaTraField()->ComputeDensity();
    FluidField()->SetScalarFields(
        ScaTraField()->Densafnp(),
        0.0,
        Teuchos::null,
        ScaTraField()->Discretization());

    // convergence check based on incremental values
    stopnonliniter = ConvergenceCheck(natconvitnum,natconvitmax_,natconvittol_);

    // Test: Output of the flux across the boundary into the output text file
    // after each outer Iteration step
    // print mean concentration
#ifdef OUTPUT
    if (stopnonliniter==false)
    {
    printf("\n");
    printf("Flux: Outer Iterations step: %3d \n", natconvitnum);
    ScaTraField().OutputProblemSpecific();
    ScaTraField().OutputTotalAndMeanScalars();
    printf("\n");
    }

//#ifdef OUTPUT
   // iteration number (only after that data output is possible)
   myoutput->NewStep(natconvitnum,natconvitnum);
   myoutput->WriteVector("phinp", ScaTraField().Phinp());
   myoutput->WriteVector("convec_velocity", ScaTraField().ConVel());
   // myoutput->WriteVector("velnp", FluidField()->Velnp());
#endif

  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Update(const int num)
{
  FluidField()->Update();
  ScaTraField()->Update(num);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::UpdateConvection()
{
  // update scatra and fluid fields
  ScaTraField()->Update();
  FluidField()->Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.

  if ((Step()>=samstart_) and (Step()<=samstop_))
  {
    // if statistics for one-way coupled problems is performed, provide
    // the field for the first scalar!
    FluidField()->SetScalarFields(
        ScaTraField()->Phinp(),
        0.0,
        Teuchos::null,
        ScaTraField()->Discretization(),
        0 // do statistics for FIRST dof at every node!!
    );
  }

  FluidField()->StatisticsAndOutput();
  ScaTraField()->Output();

  // we have to call the output of averaged fields for scatra separately
  if (  FluidField()->TurbulenceStatisticManager() != Teuchos::null)
    FluidField()->TurbulenceStatisticManager()
        ->DoOutputForScaTra(*ScaTraField()->DiscWriter(),ScaTraField()->Step());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScaTraAlgorithm::ConvergenceCheck(
  int            natconvitnum,
  const int      natconvitmax,
  const double   natconvittol
  )
{
  // convergence check based on phi and velocity increment

  //   | phi increment |_2
  //  -------------------------- < Tolerance
  //     | phi_n+1 |_2

  bool stopnonliniter = false;
  Teuchos::RCP<LINALG::MapExtractor> phisplitter = ScaTraField()->Splitter();
  // Variables to save different L2 - Norms

  double velincnorm_L2(0.0);
  double velnorm_L2(0.0);
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);

  // Calculate velocity increment and velocity L2 - Norm
  // velincnp_ = 1.0 * convelnp_ - 1.0 * conveln_

  velincnp_->Update(1.0,*FluidField()->ExtractVelocityPart(FluidField()->Velnp()),-1.0);
  velincnp_->Norm2(&velincnorm_L2); // Estimation of the L2 - norm save values to both variables (velincnorm_L2 and velnorm_L2)
  FluidField()->ExtractVelocityPart(FluidField()->Velnp())->Norm2(&velnorm_L2);

  // Calculate phi increment and phi L2 - Norm
  // tempincnp_ includes the concentration and the potential increment
  // tempincnp_ = 1.0 * phinp_ - 1.0 * phin_

  phiincnp_->Update(1.0,*ScaTraField()->Phinp(),-1.0);
  phiincnp_->Norm2(&phiincnorm_L2);
  ScaTraField()->Phinp()->Norm2(&phinorm_L2);

  // care for the case that there is (almost) zero temperature or velocity
  // (usually not required for temperature)
  if (velnorm_L2 < 1e-6) velnorm_L2 = 1.0;
  if (phinorm_L2 < 1e-6) phinorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
    if (natconvitnum != 1)
    {
      if (Comm().MyPID() == 0)
      {
      std::cout<<"\n";
      std::cout<<"*****************************************************************************\n";
      std::cout<<"                          OUTER ITERATION STEP\n";
      std::cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E  | %10.3E  |",
          natconvitnum, natconvitmax, natconvittol, phiincnorm_L2/phinorm_L2, velincnorm_L2/velnorm_L2);
      printf("\n");
      printf("+------------+-------------------+-------------+-------------+\n");
      }

      // Converged or Not
      if ((phiincnorm_L2/phinorm_L2 <= natconvittol)  &&
          (velincnorm_L2/velnorm_L2 <= natconvittol))
        //if ((incconnorm_L2/connorm_L2 <= natconvittol))
      {
        stopnonliniter=true;
        if (Comm().MyPID() == 0)
        {
        printf("| Outer Iteration loop converged after iteration %3d/%3d                    |\n", natconvitnum,natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
        }
      }
      else
      {
        if (Comm().MyPID() == 0)
        {
        printf("| Outer Iteration loop is not converged after iteration %3d/%3d             |\n", natconvitnum,natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
        }
      }

    }
    else
    {
      // first outer iteration loop: fluid solver has not got the new density yet
      // => minimum two outer iteration loops
      stopnonliniter=false;
      if (Comm().MyPID() == 0)
      {
      std::cout<<"\n";
      std::cout<<"*****************************************************************************\n";
      std::cout<<"                          OUTER ITERATION STEP\n";
      std::cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |       -     |      -      |",
          natconvitnum, natconvitmax, natconvittol);
      printf("\n");
      printf("+------------+-------------------+-------------+-------------+\n");
    }
  }

    // warn if natconvitemax is reached without convergence, but proceed to next timestep
    // natconvitemax = 1 is also possible for segregated coupling approaches (not fully implicit)
    if (natconvitnum == natconvitmax)
    {
      if (((phiincnorm_L2/phinorm_L2 > natconvittol) ||
          (velincnorm_L2/velnorm_L2 > natconvittol)) or (natconvitmax==1))
      {
        stopnonliniter=true;
        if ((Comm().MyPID() == 0))
        {
          printf("|     >>>>>> not converged in itemax steps!     |\n");
          printf("+-----------------------------------------------+\n");
          printf("\n");
          printf("\n");
        }
      }
    }

  return stopnonliniter;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::ReadInflowRestart(int restart)
{
  // in case a inflow generation in the inflow section has been performed,
  // there are not any scatra results available and the initial field is used
  FluidField()->ReadRestart(restart);
  // as ReadRestart is only called for the FluidField
  // time and step have not been set in the superior class and the ScaTraField
  SetTimeStep(FluidField()->Time(),FluidField()->Step());
  ScaTraField()->SetTimeStep(FluidField()->Time(),FluidField()->Step());
  return;
}
