/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.cpp

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "loma_algorithm.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::Algorithm(
    Epetra_Comm&                  comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false)
{
  // flag for monolithic solver
  monolithic_ = (DRT::INPUT::IntegralValue<int>(prbdyn,"MONOLITHIC"));

  // time-step length, maximum time and maximum number of steps
  dt_      = prbdyn.get<double>("TIMESTEP");
  maxtime_ = prbdyn.get<double>("MAXTIME");
  stepmax_ = prbdyn.get<int>("NUMSTEP");

  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_    = prbdyn.get<double>("CONVTOL");
  itmaxpre_ = prbdyn.get<int>("ITEMAX");

  // flag for constant thermodynamic pressure
  consthermpress_ = prbdyn.get<string>("CONSTHERMPRESS");

  // flag for special flow and start of sampling period from fluid parameter list
  special_flow_ = prbdyn.get<string>("CANONICAL_FLOW");
  samstart_     = prbdyn.get<int>("SAMPLING_START");

  // check scatra solver type, which should be incremental, for the time being
  if (ScaTraField().Incremental() == false)
    dserror("Incremental ScaTra formulation required for low-Mach-number flow");

  // preparatives for monolithic solver
  if (monolithic_)
  {
    // check whether (fluid) linearization scheme is a fixed-point-like scheme,
    // which is the only one enabled for monolithic solver, for the time being
    const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
    INPAR::FLUID::LinearisationAction linearization = DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fluiddyn, "NONLINITER");
    if (linearization != INPAR::FLUID::fixed_point_like)
      dserror("Only a fixed-point-like iteration scheme is enabled for monolithic low-Mach-number solver, for the time being!");

    // generate proxy of fluid dof set to be used by scatra field
    Teuchos::RCP<DRT::DofSet> fluiddofset = FluidField().Discretization()->GetDofSetProxy();
    // generate proxy of scatra dof set to be used by fluid field
    Teuchos::RCP<DRT::DofSet> scatradofset = ScaTraField().Discretization()->GetDofSetProxy();

    // check number of dof sets in respective fields
    if (FluidField().Discretization()->AddDofSet(scatradofset)!=1)
      dserror("Incorrect number of dof sets in fluid field!");
    if (ScaTraField().Discretization()->AddDofSet(fluiddofset)!=1)
      dserror("Incorrect number of dof sets in scatra field!");

    // create combined map for loma problem
    std::vector<Teuchos::RCP<const Epetra_Map> > dofrowmaps;

    // insert actual (zeroth) map of the discretization: first fluid, then scatra
    dofrowmaps.push_back(FluidField().DofRowMap(0));
    dofrowmaps.push_back(ScaTraField().DofRowMap(0));

    // check existence of elements
    if (dofrowmaps[0]->NumGlobalElements()==0) dserror("No fluid elements!");
    if (dofrowmaps[1]->NumGlobalElements()==0) dserror("No scatra elements!");

    Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(dofrowmaps);

    // full loma block dofrowmap
    lomablockdofrowmap_.Setup(*fullmap,dofrowmaps);

    // create loma solver
    const Teuchos::ParameterList& lomasolverparams = DRT::Problem::Instance()->CoupledFluidAndScalarTransportSolverParams();

    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(lomasolverparams,"SOLVER");
    if (solvertype != INPAR::SOLVER::aztec_msr)
      dserror("Aztec solver expected!");

    const int azprectype = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(lomasolverparams,"AZPREC");
    if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
      dserror("Block Gauss-Seidel preconditioner expected!");

    // use loma solver object
    lomasolver_ = rcp(new LINALG::Solver(lomasolverparams,FluidField().Discretization()->Comm(),DRT::Problem::Instance()->ErrorFile()->Handle()));

    lomasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->FluidSolverParams());
    lomasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->ScalarTransportFluidSolverParams());

    FluidField().Discretization()->ComputeNullSpaceIfNecessary(lomasolver_->Params().sublist("Inverse1"));
    ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(lomasolver_->Params().sublist("Inverse2"));

    // create loma block matrix
    lomasystemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(lomablockdofrowmap_,lomablockdofrowmap_,135,false,true));

    // create loma rhs vector
    lomarhs_ = rcp(new Epetra_Vector(*lomablockdofrowmap_.FullMap(),true));

    // create loma increment vector
    lomaincrement_ = rcp(new Epetra_Vector(*lomablockdofrowmap_.FullMap(),true));

    // create vector of zeros for enforcing zero Dirichlet boundary conditions
    zeros_ = rcp(new Epetra_Vector(*lomablockdofrowmap_.FullMap(),true));

    // create combined Dirichlet boundary condition map
    const Teuchos::RCP<const Epetra_Map > fdbcmap = FluidField().GetDBCMapExtractor()->CondMap();
    const Teuchos::RCP<const Epetra_Map > sdbcmap = ScaTraField().DirichMaps()->CondMap();
    lomadbcmap_ = LINALG::MergeMap(fdbcmap,sdbcmap,false);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeLoop()
{
  // do initial calculations
  InitialCalculations();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop for particular type of algorithm
    if (monolithic_) MonoLoop();
    else             OuterLoop();

    // update for next time step
    TimeUpdate();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::InitialCalculations()
{
  // set initial velocity field for evaluation of initial scalar time
  // derivative in SCATRA
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 Teuchos::null,
                                 Teuchos::null,
                                 Teuchos::null,
                                 FluidField().Discretization());

  // set initial value of thermodynamic pressure in SCATRA
  ScaTraField().SetInitialThermPressure();

  // energy conservation: compute initial time derivative of therm. pressure
  // mass conservation: compute initial mass (initial time deriv. assumed zero)
  if (consthermpress_=="No_energy")
    ScaTraField().ComputeInitialThermPressureDeriv();
  else if (consthermpress_=="No_mass")
    ScaTraField().ComputeInitialMass();

  // set initial scalar field and thermodynamic pressure for evaluation of
  // Neumann boundary conditions in FLUID at beginning of first time step
  FluidField().SetTimeLomaFields(ScaTraField().Phinp(),
                                 ScaTraField().ThermPressNp(),
                                 null,
                                 ScaTraField().Discretization());

  // write initial fields
  //Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField().PrepareTimeStep();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy") ScaTraField().PredictThermPressure();

  // prepare fluid time step, among other things, predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           ScaTraField().Time(),maxtime_,dt_,ScaTraField().MethodTitle().c_str(),ScaTraField().Step(),stepmax_);
  }

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  // set fluid values required in scatra
  SetFluidValuesInScaTra();

  // initially solve scalar transport equation
  // (values for intermediate time steps were calculated at the end of PerpareTimeStep)
  if (Comm().MyPID()==0) cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
  ScaTraField().Solve();

  while (stopnonliniter==false)
  {
    itnum++;

    // in case of non-constant thermodynamic pressure: compute
    // (either based on energy conservation or based on mass conservation)
    if (consthermpress_=="No_energy")
      ScaTraField().ComputeThermPressure();
    else if (consthermpress_=="No_mass")
      ScaTraField().ComputeThermPressureFromMassCons();

    // set scatra values required in fluid
    SetScaTraValuesInFluid();

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n****************************************\n              FLUID SOLVER\n****************************************\n";
    FluidField().MultiCorrector();

    // set fluid values required in scatra
    SetFluidValuesInScaTra();

    // solve scalar transport equation
    if (Comm().MyPID()==0) cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
    ScaTraField().Solve();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::MonoLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    cout<<"\n****************************************\n       MONOLITHIC ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           ScaTraField().Time(),maxtime_,dt_,ScaTraField().MethodTitle().c_str(),ScaTraField().Step(),stepmax_);
  }

  // maximum number of iterations tolerance for monolithic iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  while (stopnonliniter==false)
  {
    itnum++;

    // set fluid values required in scatra
    SetFluidValuesInScaTra();

    // in case of non-constant thermodynamic pressure: compute
    // (either based on energy conservation or based on mass conservation)
    if (consthermpress_=="No_energy")
      ScaTraField().ComputeThermPressure();
    else if (consthermpress_=="No_mass")
      ScaTraField().ComputeThermPressureFromMassCons();

    // set scatra values required in fluid
    SetScaTraValuesInFluid();

    // preparatives for scalar transport and fluid solver
    ScaTraField().PrepareLinearSolve();
    FluidField().PrepareSolve();

    // set up matrix and right-hand-side for monolithic low-Mach-number system
    SetupMonoLomaMatrix();
    SetupMonoLomaRHS();

    // solve monolithic low-Mach-number system
    MonoLomaSystemSolve();

    // update for next iteration step
    IterUpdate();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SetFluidValuesInScaTra()
{
  // set respective field vectors for velocity/pressure, acceleration
  // and discretization based on time-integration scheme
  if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    ScaTraField().SetVelocityField(FluidField().Velaf(),
                                   FluidField().Accam(),
                                   Teuchos::null,
                                   Teuchos::null,
                                   FluidField().Discretization());
  else
    ScaTraField().SetVelocityField(FluidField().Velnp(),
                                   FluidField().Hist(),
                                   Teuchos::null,
                                   Teuchos::null,
                                   FluidField().Discretization());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SetScaTraValuesInFluid()
{
  // set scalar and thermodynamic pressure values as well as time
  // derivatives and discretization based on time-integration scheme
  if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    FluidField().SetIterLomaFields(ScaTraField().Phiaf(),
                                   ScaTraField().Phiam(),
                                   ScaTraField().Phidtam(),
                                   ScaTraField().ThermPressAf(),
                                   ScaTraField().ThermPressAm(),
                                   ScaTraField().ThermPressDtAf(),
                                   ScaTraField().ThermPressDtAm(),
                                   ScaTraField().Discretization());
  else
    FluidField().SetIterLomaFields(ScaTraField().Phinp(),
                                   ScaTraField().Phin(),
                                   ScaTraField().Phidtnp(),
                                   ScaTraField().ThermPressNp(),
                                   ScaTraField().ThermPressN(),
                                   ScaTraField().ThermPressDtNp(),
                                   ScaTraField().ThermPressDtNp(),
                                   ScaTraField().Discretization());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SetupMonoLomaMatrix()
{
  // set loma block matrix to zero
  lomasystemmatrix_->Zero();

  //----------------------------------------------------------------------
  // 1st diagonal block (upper left): fluid weighting - fluid solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ff = FluidField().SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_ff->UnComplete();

  // assign matrix block
  lomasystemmatrix_->Assign(0,0,View,*mat_ff);

  //----------------------------------------------------------------------
  // 2nd diagonal block (lower right): scatra weighting - scatra solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ss = ScaTraField().SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_ss->UnComplete();

  // assign matrix block
  lomasystemmatrix_->Assign(1,1,View,*mat_ss);

  // complete loma block matrix
  lomasystemmatrix_->Complete();

  //----------------------------------------------------------------------
  // 1st off-diagonal block (upper right): fluid weighting - scatra solution
  //----------------------------------------------------------------------
  // create matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_fs = Teuchos::null;
  mat_fs = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField().Discretization()->DofRowMap(0)),27,true,true));

  // evaluate loma off-diagonal matrix block in fluid
  EvaluateLomaODBlockMatFluid(mat_fs);

  // uncomplete matrix block (appears to be required in certain cases)
  mat_fs->UnComplete();

  // assign matrix block
  lomasystemmatrix_->Assign(0,1,View,*mat_fs);

  //----------------------------------------------------------------------
  // 2nd off-diagonal block (lower left): scatra weighting - fluid solution
  //----------------------------------------------------------------------
  // create matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_sf = Teuchos::null;
  mat_sf = Teuchos::rcp(new LINALG::SparseMatrix(*(ScaTraField().Discretization()->DofRowMap(0)),108,true,true));

  // evaluate loma off-diagonal matrix block in scatra
  // (for present fixed-point-like iteration: no entries)
  //EvaluateLomaODBlockMatScaTra(mat_sf);

  // uncomplete matrix block (appears to be required in certain cases)
  mat_sf->UnComplete();

  // assign matrix block
  lomasystemmatrix_->Assign(1,0,View,*mat_sf);

  // complete loma block matrix
  lomasystemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::EvaluateLomaODBlockMatFluid(
  Teuchos::RCP<LINALG::SparseMatrix> mat_fs
  )
{
  // create parameters for fluid discretization
  ParameterList fparams;

  // set action type
  fparams.set("action","calc_loma_mono_odblock");

  // set general vector values needed by elements
  FluidField().Discretization()->ClearState();
  FluidField().Discretization()->SetState(0,"hist" ,FluidField().Hist());
  FluidField().Discretization()->SetState(0,"accam",FluidField().Accam());
  FluidField().Discretization()->SetState(0,"scaaf",FluidField().Scaaf());
  FluidField().Discretization()->SetState(0,"scaam",FluidField().Scaam());

  // set time-integration-scheme-specific element parameters and vector values
  if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
  {
    // set thermodynamic pressures
    fparams.set("thermpress at n+alpha_F/n+1",ScaTraField().ThermPressAf());
    fparams.set("thermpress at n+alpha_M/n",ScaTraField().ThermPressAm());
    fparams.set("thermpressderiv at n+alpha_F/n+1",ScaTraField().ThermPressDtAf());
    fparams.set("thermpressderiv at n+alpha_M/n+1",ScaTraField().ThermPressDtAm());

    // set velocity vector
    FluidField().Discretization()->SetState(0,"velaf",FluidField().Velaf());
  }
  else
  {
    // set thermodynamic pressures
    fparams.set("thermpress at n+alpha_F/n+1",ScaTraField().ThermPressNp());
    fparams.set("thermpress at n+alpha_M/n",ScaTraField().ThermPressN());
    fparams.set("thermpressderiv at n+alpha_F/n+1",ScaTraField().ThermPressDtNp());
    fparams.set("thermpressderiv at n+alpha_M/n+1",ScaTraField().ThermPressDtNp());

    // set velocity vector
    FluidField().Discretization()->SetState(0,"velaf",FluidField().Velnp());
 }

  // build specific assemble strategy for this off-diagonal matrix block,
  // which is assembled in fluid solver
  // fluid dof set = 0, scatra dof set = 1
  DRT::AssembleStrategy fluidstrategy(
                          0,  // rows: fluid dof set
                          1,  // columns: scatra dof set
                          mat_fs,
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null
                          );

  // evaluate off-diagonal matrix block entries for fluid element
  FluidField().Discretization()->Evaluate(fparams,fluidstrategy);
  FluidField().Discretization()->ClearState();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SetupMonoLomaRHS()
{
  // define fluid and scatra residual vectors
  Teuchos::RCP<const Epetra_Vector> fluidres  = FluidField().RHS();
  Teuchos::RCP<const Epetra_Vector> scatrares = ScaTraField().Residual();

  // insert fluid and scatra residual vectors into loma residual vector
  lomablockdofrowmap_.InsertVector(*fluidres, 0,*lomarhs_);
  lomablockdofrowmap_.InsertVector(*scatrares,1,*lomarhs_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::MonoLomaSystemSolve()
{
  // set incremental solution vector to zero
  lomaincrement_->PutScalar(0.0);

  // apply Dirichlet boundary conditions to system
  LINALG::ApplyDirichlettoSystem(lomasystemmatrix_,lomaincrement_,lomarhs_,Teuchos::null,zeros_,*lomadbcmap_);

  // solve monolithic low-Mach-number system
  lomasolver_->Solve(lomasystemmatrix_->EpetraOperator(),lomaincrement_,lomarhs_,true,true);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::IterUpdate()
{
  // define incremental fluid and scatra solution vectors
  Teuchos::RCP<const Epetra_Vector> incfluid;
  Teuchos::RCP<const Epetra_Vector> incscatra;

  // extract incremental fluid and scatra solution vectors
  // from incremental low-Mach-number solution vector
  incfluid  = lomablockdofrowmap_.ExtractVector(lomaincrement_,0);
  incscatra = lomablockdofrowmap_.ExtractVector(lomaincrement_,1);

  // add incremental fluid and scatra solution vectors to
  // respective solution vectors from last iteration step
  FluidField().IterUpdate(incfluid);
  ScaTraField().UpdateIter(incscatra);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool LOMA::Algorithm::ConvergenceCheck(int itnum)
{
  // define flags for fluid and scatra convergence check
  bool fluidstopnonliniter  = false;
  bool scatrastopnonliniter = false;

  // fluid convergence check
  if (Comm().MyPID()==0)
  {
    cout<<"\n****************************************\n  CONVERGENCE CHECK FOR ITERATION STEP\n****************************************\n";
    cout<<"\n****************************************\n              FLUID CHECK\n****************************************\n";
  }
  fluidstopnonliniter  = FluidField().ConvergenceCheck(itnum,itmax_,ittol_);

  // scatra convergence check
  if (Comm().MyPID()==0) cout<<"\n****************************************\n         SCALAR TRANSPORT CHECK\n****************************************\n";
  scatrastopnonliniter = ScaTraField().ConvergenceCheck(itnum,itmax_,ittol_);

  if (fluidstopnonliniter == true and scatrastopnonliniter == true) return true;
  else                                                              return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeUpdate()
{
  // update scalar
  ScaTraField().Update();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    ScaTraField().UpdateThermPressure();

  // update fluid
  FluidField().Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Output()
{
  // set scalar and thermodynamic pressure at n+1 and SCATRA trueresidual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  FluidField().SetTimeLomaFields(ScaTraField().Phinp(),
                                 ScaTraField().ThermPressNp(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().Discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();

  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
