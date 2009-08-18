/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_implicit.cpp
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_timint_implicit.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_lib/drt_function.H"
#include "../drt_fluid/fluid_utils.H" // for splitter
#include "scatra_utils.H" // for splitstrategy

// for the condition writer output
/*
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
*/
/*
// for output of intermediate states in Newton loop
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
*/


/*----------------------------------------------------------------------*
 |  Constructor (public)                                        vg 05/07|
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::ScaTraTimIntImpl(
    RCP<DRT::Discretization>      actdis,
    RCP<LINALG::Solver>           solver,
    RCP<ParameterList>            params,
    RCP<IO::DiscretizationWriter> output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  myrank_ (discret_->Comm().MyPID()),
  time_   (0.0),
  step_   (0),
  prbtype_  (params_->get<string>("problem type")),
  solvtype_ (params_->get<string>("solver type")),
  reaction_ (params_->get<string>("reaction")),
  isale_    (params_->get<bool>("isale")),
  stepmax_  (params_->get<int>("max number timesteps")),
  maxtime_  (params_->get<double>("total time")),
  timealgo_ (params_->get<INPAR::SCATRA::TimeIntegrationScheme>("time int algo")),
  upres_    (params_->get<int>("write solution every")),
  uprestart_(params_->get<int>("write restart every")),
  writeflux_(params_->get<string>("write flux")),
  outmean_  (params_->get<bool>("write mean values")),
  dta_      (params_->get<double>("time step size")),
  dtp_      (params_->get<double>("time step size")),
  cdvel_    (params_->get<int>("velocity field")),
  convform_ (params_->get<string>("form of convective term")),
  neumannin_(params_->get<string>("Neumann inflow")),
  fssgd_    (params_->get<string>("fs subgrid diffusivity")),
  frt_      (96485.3399/(8.314472 * params_->get<double>("TEMPERATURE",298.15))),
  tpn_      (1.0),
  errfile_  (params_->get<FILE*>("err file")),
  initialvelset_(false),
  lastfluxoutputstep_(-1)
{
  // -------------------------------------------------------------------
  // determine whether linear incremental or nonlinear solver
  // -------------------------------------------------------------------
  if (solvtype_ == "nonlinear")
  {
    incremental_ = true;
    nonlinear_   = true;
  }
  else if (solvtype_ == "linear_incremental")
  {
    incremental_ = true;
    nonlinear_   = false;
  }
  else
  {
    incremental_ = false;
    nonlinear_   = false;
  }

  // -------------------------------------------------------------------
  // determine whether Neumann inflow terms need to be accounted for
  // -------------------------------------------------------------------
  neumanninflow_ = false;
  if (neumannin_ == "yes") neumanninflow_ = true;

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  pbc_ = rcp(new PeriodicBoundaryConditions (discret_));
  pbc_->UpdateDofsForPeriodicBoundaryConditions();

  pbcmapmastertoslave_ = pbc_->ReturnAllCoupledNodesOnThisProc();

  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  numscal_ = discret_->NumDof(discret_->lRowNode(0));
  if (prbtype_ == "elch")
  {
    // number of concentrations transported is numdof-1
    numscal_ -= 1;
    // set up the concentration-el.potential splitter
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_,splitter_);
  }
  else if (prbtype_ == "loma" and numscal_ > 1)
  {
    // set up a species-temperature splitter
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_-1,splitter_);
  }

  if (params_->get<int>("BLOCKPRECOND",0) )
  {
    // we need a block sparse matrix here
    if (prbtype_ != "elch")
      dserror("Block-Preconditioning is only for ELCH problems");
    // initial guess for non-zeros per row: 27 neighboring nodes for hex8 times (numscal_+1) dofs

    Teuchos::RCP<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(splitter_,splitter_,27*(numscal_+1),false,true));
    blocksysmat->SetNumScal(numscal_);

    // the fluid assembly strategy still works, but it does not make use
    // of the ELCH-specific sparsity pattern
    //Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
    //Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(splitter_,splitter_,27*(numscal_+1),false,true));
    //blocksysmat->SetNumdim(numscal_);

    sysmat_ = blocksysmat;
  }
  else
  {
    // initialize standard (stabilized) system matrix (and save its graph!)
    // in standard case, but do not save the graph if fine-scale subgrid
    // diffusivity is used
    if (fssgd_ == "No") sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
    else sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));
  }

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1 and n
  phinp_ = LINALG::CreateVector(*dofrowmap,true);
  phin_  = LINALG::CreateVector(*dofrowmap,true);

  // density at time n+1
  densnp_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  hist_ = LINALG::CreateVector(*dofrowmap,true);

  // convective velocity (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  convel_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // subgrid-scale velocity (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  sgvel_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  if (isale_)
  {
    // displacement field for moving mesh applications using ALE
    // (get noderowmap of discretization for creating this multivector)
    dispnp_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));
  }

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  // -------------------------------------------------------------------
  // create vectors associated to solution process
  // -------------------------------------------------------------------
  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // the residual vector --- more or less the rhs
  residual_ = LINALG::CreateVector(*dofrowmap,true);

  // residual vector containing the normal boundary fluxes
  trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

  // incremental solution vector
  increment_ = LINALG::CreateVector(*dofrowmap,true);

  // subgrid-diffusivity(-scaling) vector
  // (used either for AVM3 approach or temperature equation
  //  with all-scale subgrid-diffusivity model)
  subgrdiff_ = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // set parameters associated to potential statistical flux evaluations
  // -------------------------------------------------------------------
  // get fluid turbulence sublist
  ParameterList * turbparams =&(params_->sublist("TURBULENCE PARAMETERS"));

  // parameters for statistical evaluation of normal fluxes
  samstart_  = turbparams->get<int>("SAMPLING_START");
  samstop_   = turbparams->get<int>("SAMPLING_STOP" );
  dumperiod_ = turbparams->get<int>("DUMPING_PERIOD");

  // initialize vector for statistics (assume a maximum of 10 conditions)
  sumnormfluxintegral_ = Teuchos::rcp(new Epetra_SerialDenseVector(10));

  // -------------------------------------------------------------------
  // necessary only for AVM3 approach:
  // initialize subgrid-diffusivity matrix + respective ouptput
  // -------------------------------------------------------------------
  if (fssgd_ != "No")
  {
    sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

    // Output
    if (myrank_ == 0)
    {
      cout << "Fine-scale subgrid-diffusivity approach based on AVM3: ";
      cout << fssgd_;
      cout << &endl << &endl;
    }
  }

  // -------------------------------------------------------------------
  // create specific vectors for low-Mach-number case
  // (i.e., scalar variable is temperature)
  // get also turbulence model and parameters
  // -------------------------------------------------------------------
  turbmodel_ = false;
  if (prbtype_ == "loma")
  {
    // temperature and velocity increment at time n+1
    tempincnp_ = LINALG::CreateVector(*dofrowmap,true);
    //velincnp_  = rcp(new Epetra_MultiVector(*noderowmap,3,true));

    // potential turbulence model
    if (turbparams->get<string>("PHYSICAL_MODEL") != "no_model")
      turbmodel_ = true;

    // warning if classical (all-scale) turbulence model and fine-scale
    // subgrid-viscosity approach are intended to be used simultaneously
    if (turbmodel_ and fssgd_ != "No")
      dserror("No combination of classical (all-scale) turbulence model and fine-scale subgrid-diffusivity approach currently possible!");

    // turbulent Prandtl number
    tpn_ = turbparams->get<double>("C_TURBPRANDTL",1.0);

    // Output
    if (turbmodel_ and myrank_ == 0)
    {
      cout << "All-scale subgrid-diffusivity model: ";
      cout << turbparams->get<string>("PHYSICAL_MODEL");
      cout << &endl << &endl;
    }
  }

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(params_->get<int>("scalar initial field"), params_->get<int>("scalar initial field func number"));

  // -------------------------------------------------------------------
  // set initial density to prescribed value (default = 1.0):
  // - used throughout simulation for non-temperature case
  // - used as good initial guess for stationary temperature case
  // -------------------------------------------------------------------
  const double initdens = params_->get<double>("initial density",1.0);
  densnp_->PutScalar(initdens);

  // initializes variables for natural convection (ELCH) if necessary
  SetupElchNatConv();

  // screen output (has to come after SetInitialField)
  if (prbtype_ == "elch")
  {
    double sigma = ComputeConductivity(); // every processor has to do this call
    if (myrank_==0)
    {
      cout<<"\nSetup of splitter: numscal = "<<numscal_<<endl;
      cout<<"Temperature value T (Kelvin)     = "<<params_->get<double>("TEMPERATURE",298.0)<<endl;
      cout<<"Constant F/RT                    = "<<frt_<<endl;
      cout<<"Conductivity of electrolyte      = "<<sigma<<endl<<endl;
    }
  }

  // sysmat might be singular (some modes are defined only up to a constant)
  // in this case, we need basis vectors for the nullspace/kernel
  vector<DRT::Condition*> KSPCond;
  discret_->GetCondition("KrylovSpaceProjection",KSPCond);
  int nummodes = KSPCond.size();

  if (nummodes > 0)
  {
    project_ = true;
    w_       = rcp(new Epetra_MultiVector(*dofrowmap,nummodes,true));
    c_       = rcp(new Epetra_MultiVector(*dofrowmap,nummodes,true));
    if (myrank_ == 0)
      cout<<"\nSetup of KrylovSpaceProjection:\n"<<
      "    => number of kernel basis vectors: "<<nummodes<<endl<<endl;
  }
  else
  {
    project_ = false;
    w_       = Teuchos::null;
    c_       = Teuchos::null;
  }

  // -------------------------------------------------------------------
  // preliminaries for reactive flow using progress-variable approach
  // -------------------------------------------------------------------
  if (reaction_ == "Arrhenius_pv")
  {
    ParameterList eleparams;
    eleparams.set("action","get_density_values");
    eleparams.set("isale",isale_);
    discret_->Evaluate(eleparams,null,null,null,null,null);
    unbdens_ = eleparams.get("unburnt density", 1.161);
    const double burdens = eleparams.get("burnt density", 0.29);
    if (burdens > unbdens_)
      dserror("density in burnt phase larger than in unburnt phase!");
    densdiff_ = burdens - unbdens_;
  }

  return;

} // ScaTraTimIntImpl::ScaTraTimIntImpl


/*----------------------------------------------------------------------*
| returns matching string for each time integration scheme   gjb 08/08 |
*----------------------------------------------------------------------*/
std::string SCATRA::ScaTraTimIntImpl::MapTimIntEnumToString
(
   const enum INPAR::SCATRA::TimeIntegrationScheme term
)
{
  // length of return string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPAR::SCATRA::timeint_one_step_theta :
    return "One-Step-Theta";
    break;
  case INPAR::SCATRA::timeint_bdf2 :
    return "    BDF2      ";
    break;
  case INPAR::SCATRA::timeint_stationary :
    return "  Stationary  ";
    break;
  case INPAR::SCATRA::timeint_gen_alpha :
    return "  Gen. Alpha  ";
    break;
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }
}


/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::TimeLoop()
{
  // write out inital state
  // Output();

  // compute error for problems with analytical solution (initial field!)
  EvaluateErrorComparedToAnalyticalSol();

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  while ((step_<stepmax_) and ((time_+ EPS12) < maxtime_))
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    PrepareTimeStep();

    // -------------------------------------------------------------------
    // compute values at intermediate time steps (only for gen.-alpha)
    // -------------------------------------------------------------------
    ComputeIntermediateValues();

    // -------------------------------------------------------------------
    //                  solve nonlinear / linear equation
    // -------------------------------------------------------------------
    if (nonlinear_) NonlinearSolve();
    else            Solve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();

    // -------------------------------------------------------------------
    // evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

  } // while

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // ScaTraTimIntImpl::TimeLoop


/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                    vg 08/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareTimeStep()
{
  // time measurement: prepare time step
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + prepare time step");

  // -------------------------------------------------------------------
  //                       initialization
  // -------------------------------------------------------------------
  if (step_ == 0)
  {
    // if initial velocity field has not been set here, the initial time derivative of phi will be
    // calculated wrongly for some time integration schemes
    if (initialvelset_)
      PrepareFirstTimeStep();
    else
      dserror("Initial velocity field has not been set");
  }

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  IncrementTimeAndStep();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  ApplyDirichletBC(time_,phinp_,Teuchos::null);
  ApplyNeumannBC(time_,phinp_,neumann_loads_);

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if (step_==1 and fssgd_ != "No") AVM3Preparation();

  return;

} // ScaTraTimIntImpl::PrepareTimeStep


/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}           gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> phinp,
  Teuchos::RCP<Epetra_Vector> phidt
)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:      + apply dirich cond.");

  // needed parameters
  ParameterList p;
  p.set("total time",time);  // actual time t_{n+1}

  // predicted Dirichlet values
  // \c  phinp then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->EvaluateDirichlet(p,phinp,phidt,Teuchos::null,Teuchos::null,dbcmaps_);
  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions at t_{n+1}             gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyNeumannBC
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> phinp,
  Teuchos::RCP<Epetra_Vector> neumann_loads
)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // set time for evaluation of Neumann boundary conditions as parameter
  // depending on time-integration scheme
  ParameterList p;
  SetTimeForNeumannEvaluation(p);

  discret_->ClearState();
  // evaluate Neumann conditions at actual time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateNeumann(p,*neumann_loads);
  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 | export multivector to column map & add it to parameter list gjb 06/09|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddMultiVectorToParameterList
(Teuchos::ParameterList& p,
    const std::string name,
    Teuchos::RCP<Epetra_MultiVector> vec
)
{
  //provide data in node-based multi-vector for usage on element level
  // -> export to column map is necessary for parallel evaluation
  //SetState cannot be used since this multi-vector is nodebased and not dofbased!
  const Epetra_Map* nodecolmap = discret_->NodeColMap();
  int numcol = vec->NumVectors();
  RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,numcol));
  LINALG::Export(*vec,*tmp);
  p.set(name,tmp);

  return;
}


/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                       gjb 09/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::NonlinearSolve()
{
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:   + nonlin. iteration/lin. solve");

  // out to screen
  PrintTimeStepInfo();
  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- con-res ---|-- pot-res ---|-- con-inc ---|-- pot-inc ---|\n");
  }

  // ---------------------------------------------- nonlinear iteration
  //stop nonlinear iteration when both increment-norms are below this bound
  const double  ittol = params_->sublist("NONLINEAR").get<double>("CONVTOL");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = (getIntegralValue<int>(params_->sublist("NONLINEAR"),"ADAPTCONV") == 1);
  const double adaptolbetter = params_->sublist("NONLINEAR").get<double>("ADAPTCONV_BETTER");
  double       actresidual(0.0);

  int   itnum = 0;
  int   itemax = params_->sublist("NONLINEAR").get<int>("ITEMAX");
  bool  stopnonliniter = false;

  // perform explicit predictor step (-> better starting point for nonlinear solver)
  const bool explpredictor = (getIntegralValue<int>(params_->sublist("NONLINEAR"),"EXPLPREDICT") == 1);
  if (explpredictor)
    ExplicitPredictor();

/*
  const int numdim = 3;
  //create output file name
  stringstream temp;
  temp<< DRT::Problem::Instance()->OutputControlFile()->FileName()<<".nonliniter_step"<<step_;
  string outname = temp.str();
  string probtype = DRT::Problem::Instance()->ProblemType(); // = "elch"

  RCP<IO::OutputControl> myoutputcontrol = rcp(new IO::OutputControl(discret_->Comm(),probtype,"Polynomial","myinput",outname,numdim,0,1000));
  // create discretization writer with my own control settings
  RCP<IO::DiscretizationWriter> myoutput =
    rcp(new IO::DiscretizationWriter(discret_,myoutputcontrol));
  // write mesh at step 0
  myoutput->WriteMesh(0,0.0);
*/

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and rhs and assemble
    // -------------------------------------------------------------------
    AssembleMatAndRHS();

    // -------------------------------------------------------------------
    // potential residual scaling and potential addition of Neumann terms
    // -------------------------------------------------------------------
    ScalingAndNeumann();

    // add contributions due to electrode kinetics conditions
    EvaluateElectrodeKinetics(sysmat_,residual_);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the Dirichlet positions
    // are not used anyway.
    // We could avoid this though, if the dofrowmap would not include
    // the Dirichlet values as well. But it is expensive to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    // abort nonlinear iteration if desired
    if (AbortNonlinIter(itnum,itemax,ittol,actresidual))
       break;

    //--------- Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    increment_->PutScalar(0.0);

    // Apply dirichlet boundary conditions to system matrix
    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }

    //------------------------------------------------solve
    {
      // get cpu time
      const double tcpusolve=ds_cputime();

      // time measurement: call linear solver
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + call linear solver");

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        solver_->AdaptTolerance(ittol,actresidual,adaptolbetter);
      }

/*
      // matrix printing options (DEBUGGING!)
      RCP<LINALG::SparseMatrix> A = SystemMatrix();
      if (A != Teuchos::null)
      {
        // print to file in matlab format
        const std::string fname = "sparsematrix.mtl";
        LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
        // print to screen
        (A->EpetraMatrix())->Print(cout);
        // print sparsity pattern to file
        LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
      }
      else
      {
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = BlockSystemMatrix();
        const std::string fname = "sparsematrix.mtl";
        LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
      }
      */

      PrepareKrylovSpaceProjection();
      solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,itnum==1,w_,c_,project_);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpusolve;
    }

    //------------------------------------------------ update solution vector
    phinp_->Update(1.0,*increment_,1.0);

    // iteration number (only after that data output is possible)
  /*
    myoutput->NewStep(itnum,itnum);
    myoutput->WriteVector("phinp", phinp_);
   */

  } // nonlinear iteration
  return;
}


/*----------------------------------------------------------------------*
 | check if to stop the nonlinear iteration                    gjb 09/08|
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::AbortNonlinIter(
    const int itnum,
    const int itemax,
    const double ittol,
    double& actresidual)
{
  //----------------------------------------------------- compute norms
  double incconnorm_L2(0.0);
  double incpotnorm_L2(0.0);

  double connorm_L2(0.0);
  double potnorm_L2(0.0);

  double conresnorm(0.0);
  double potresnorm(0.0);

  if (prbtype_ == "elch")
  {
    Teuchos::RCP<Epetra_Vector> onlycon = splitter_.ExtractOtherVector(residual_);
    onlycon->Norm2(&conresnorm);

    splitter_.ExtractOtherVector(increment_,onlycon);
    onlycon->Norm2(&incconnorm_L2);

    splitter_.ExtractOtherVector(phinp_,onlycon);
    onlycon->Norm2(&connorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypot = splitter_.ExtractCondVector(residual_);
    onlypot->Norm2(&potresnorm);

    splitter_.ExtractCondVector(increment_,onlypot);
    onlypot->Norm2(&incpotnorm_L2);

    splitter_.ExtractCondVector(phinp_,onlypot);
    onlypot->Norm2(&potnorm_L2);
  }
  else
  {
    residual_ ->Norm2(&conresnorm);
    increment_->Norm2(&incconnorm_L2);
    phinp_    ->Norm2(&connorm_L2);
  }

  // care for the case that nothing really happens in the concentration
  // or potential field
  if (connorm_L2 < 1e-5)
  {
    connorm_L2 = 1.0;
  }
  if (potnorm_L2 < 1e-5)
  {
    potnorm_L2 = 1.0;
  }

  //-------------------------------------------------- output to screen
  /* special case of very first iteration step:
      - solution increment is not yet available
      - do not do a solver call when the initial residuals are < EPS15*/
  if (itnum == 1)
  {
    if (myrank_ == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
          itnum,itemax,ittol,conresnorm,potresnorm);
      printf(" (      --     ,te=%10.3E",dtele_);
      printf(")\n");
    }
    // abort iteration for ELCH, when there's nothing to do
    if ((conresnorm < EPS15) && (potresnorm < EPS15))
    {
      // print 'finish line'
      if (myrank_ == 0)
      {
        printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
      }
      return true;
    }
  }
  /* ordinary case later iteration steps:
      - solution increment can be printed
      - convergence check should be done*/
  else
  {
    // print the screen info
    if (myrank_ == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
          itnum,itemax,ittol,conresnorm,potresnorm,
          incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
      printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
    }

    // this is the convergence check
    // We always require at least one solve. We test the L_2-norm of the
    // current residual. Norm of residual is just printed for information
    if (conresnorm <= ittol and potresnorm <= ittol and
        incconnorm_L2/connorm_L2 <= ittol and incpotnorm_L2/potnorm_L2 <= ittol)
    {
      if (myrank_ == 0)
      {
        // print 'finish line'
        printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
        // write info to error file
        FILE* errfile = params_->get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"elch solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
              itnum,itemax,ittol,conresnorm,potresnorm,
              incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
        }
      }
      // yes, we stop the iteration
      return true;
    }
    // if not yet converged go on...
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if ((itnum == itemax))
  {
    if (myrank_ == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");

      FILE* errfile = params_->get<FILE*>("err file",NULL);
      if (errfile!=NULL)
      {
        fprintf(errfile,"elch divergent solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
            itnum,itemax,ittol,conresnorm,potresnorm,
            incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
      }
    }
    // yes, we stop the iteration
    return true;
  }

  // return the maximum residual value -> used for adaptivity of linear solver tolarance
  actresidual = max(conresnorm,potresnorm);
  actresidual = max(actresidual,incconnorm_L2/connorm_L2);
  actresidual = max(actresidual,incpotnorm_L2/potnorm_L2);

  return false;
}


/*----------------------------------------------------------------------*
 | contains the linear solver                                  vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Solve()
{
  // -------------------------------------------------------------------
  //                        output to screen
  // -------------------------------------------------------------------
  PrintTimeStepInfo();

  // -------------------------------------------------------------------
  // call elements to calculate system matrix and rhs and assemble
  // -------------------------------------------------------------------
  AssembleMatAndRHS();

  // -------------------------------------------------------------------
  // potential residual scaling and potential addition of Neumann terms
  // -------------------------------------------------------------------
  ScalingAndNeumann();

  // -------------------------------------------------------------------
  // Apply Dirichlet boundary conditions to system matrix and solve
  // system in incremental or non-incremental case
  // -------------------------------------------------------------------
  if (incremental_)
  {
    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the Dirichlet positions
    // are not used anyway.
    // We could avoid this though, if the dofrowmap would not include
    // the Dirichlet values as well. But it is expensive to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    //--------- Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    increment_->PutScalar(0.0);

    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }

    //-------solve
    {
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

      // get cpu time
      const double tcpusolve=ds_cputime();

      solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,true);

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpusolve;
    }

    //------------------------------------------------ update solution vector
    phinp_->Update(1.0,*increment_,1.0);

    //--------------------------------------------- compute norm of increment
    double incnorm_L2(0.0);
    double scalnorm_L2(0.0);
    increment_->Norm2(&incnorm_L2);
    phinp_    ->Norm2(&scalnorm_L2);

    if (myrank_ == 0)
    {
      printf("+-------------+-------------+\n");
      printf("|  increment  | %10.3E  |\n",incnorm_L2/scalnorm_L2);
      printf("+-------------+-------------+\n");
    }
  }
  else
  {
    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,phinp_,*(dbcmaps_->CondMap()));
    }

    //-------solve
    {
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

      // get cpu time
      const double tcpusolve=ds_cputime();

      solver_->Solve(sysmat_->EpetraOperator(),phinp_,residual_,true,true);

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpusolve;
    }
  }

  return;
} // ScaTraTimIntImpl::Solve


/*----------------------------------------------------------------------*
 | contains the assembly process for matrix and rhs            vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AssembleMatAndRHS()
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

  // get cpu time
  const double tcpuele = ds_cputime();

  // zero out matrix entries
  sysmat_->Zero();

  // reset the residual vector
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set("action","calc_condif_systemmat_and_residual");

  // other parameters that might be needed by the elements
  eleparams.set("time-step length",dta_);
  eleparams.set("problem type",prbtype_);
  eleparams.set("incremental solver",incremental_);
  eleparams.set("reaction",reaction_);
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid diffusivity",fssgd_);
  eleparams.set("turbulence model",turbmodel_);
  eleparams.set("frt",frt_);// ELCH specific factor F/RT

  // provide velocity field (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"subgrid-scale velocity field",sgvel_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist",hist_);
  if (turbmodel_) discret_->SetState("subgrid diffusivity",subgrdiff_);

  // AVM3 separation
  if (incremental_ and fssgd_ != "No")
  {
    discret_->SetState("subgrid diffusivity",subgrdiff_);
    AVM3Separation();
  }

  // add element parameters and density state according to time-int. scheme
  AddSpecificTimeIntegrationParameters(eleparams);

  // call loop over elements with subgrid-diffusivity(-scaling) vector
  discret_->Evaluate(eleparams,sysmat_,null,residual_,subgrdiff_,null);
  discret_->ClearState();

  // AVM3 scaling
  if (not incremental_ and fssgd_ != "No") AVM3Scaling(eleparams);

  // finalize the complete matrix
  sysmat_->Complete();

  // end time measurement for element
  dtele_=ds_cputime()-tcpuele;

  return;
} // ScaTraTimIntImpl::AssembleMatAndRHS


/*----------------------------------------------------------------------*
 | contains the residual scaling and addition of Neumann terms vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ScalingAndNeumann()
{
  // scaling to get true residual vector for all time integration schemes
  // in incremental case: boundary flux values can be computed from trueresidual
  if (incremental_) trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  AddNeumannToResidual();

  // add potential Neumann inflow
  if (neumanninflow_) ComputeNeumannInflow(sysmat_,residual_);

  return;
} // ScaTraTimIntImpl::ScalingAndNeumann


/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                          gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Output()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // solution output and potentially restart data and/or flux data
  if ((step_%upres_==0 )or (step_%uprestart_==0))
  {
    // step number and time (only after that data output is possible)
    output_->NewStep(step_,time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_==upres_) output_->WriteElementData();

    // write state vectors
    OutputState();

    // add restart data
    if (step_%uprestart_==0) OutputRestart();

    // write flux vector field
    if (writeflux_!="No") OutputFlux();

    // write mean values of scalars and density
    if (outmean_)
    {
      OutputMeanTempAndDens();
      OutputElectrodeInfo();
    }
  }
  else
  {
    // calculation of statistics for normal fluxes (no output to file)
    if (step_>=samstart_ and step_<=samstop_ and writeflux_!="No") CalcFlux();
  }

  return;
} // ScaTraTimIntImpl::Output


/*----------------------------------------------------------------------*
 |  write current state to BINIO                             gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputState()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  // density (only required in low-Mach-number case)
  if (prbtype_ == "loma") output_->WriteVector("densnp", densnp_);

  // velocity
  output_->WriteVector("convec_velocity", convel_,IO::DiscretizationWriter::nodevector);

  // displacement field
  if (isale_) output_->WriteVector("dispnp", dispnp_);

  return;
}


/*----------------------------------------------------------------------*
 | update the velocity field                                  gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField()
{
  if (cdvel_ == 0) // zero
    convel_->PutScalar(0); // just to be sure!
  else if ((cdvel_ == 1))  // function
  {
    const int numdim = 3; // the velocity field is always 3D
    const int velfuncno = params_->get<int>("velocity function number");
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value=DRT::UTILS::FunctionManager::Instance().Funct(velfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);
        convel_->ReplaceMyValue (lnodeid, index, value);
      }
    }
  }
  else
    dserror("Wrong SetVelocity() action for velocity field type %d!",cdvel_);

  // initial velocity field has now been set
  if (step_ == 0) initialvelset_ = true;

  return;

} // ScaTraImplicitTimeInt::SetVelocityField


/*----------------------------------------------------------------------*
 | set convective velocity field, (subgrid velocity field and subgrid   |
 | viscosity field)                                           gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(
Teuchos::RCP<const Epetra_Vector> fluidvel,
Teuchos::RCP<const Epetra_Vector> fluidsgvelvisc,
Teuchos::RCP<DRT::Discretization> fluiddis)
{
  if (cdvel_ != 2)
    dserror("Wrong SetVelocityField() called for velocity field type %d!",cdvel_);

  TEUCHOS_FUNC_TIME_MONITOR("SCATRA: set convective velocity field");

  // boolean indicating whether subgrid velocity-diffusivity vector does exist
  bool subgridswitch = fluidsgvelvisc!=Teuchos::null;

  // variable to store subgrid diffusivity temporarily
  double subgriddiff;

  // get dofrowmap of fluid discretization
  const Epetra_Map* fluiddofrowmap = fluiddis->DofRowMap();

  // loop over all local nodes of scatra discretization
  for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // Here we rely on the fact that the scatra discretization
    // is a clone of the fluid mesh. => a scatra node has the same
    // local (and global) ID as its corresponding fluid node!

    // get the processor's local fluid node with the same lnodeid
    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this fluid node
    vector<int> fluidnodedofs = fluiddis->Dof(fluidlnode);
    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;

    // now we transfer velocity dofs only
    for(int index=0;index < numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[index];
      const int flid = fluiddofrowmap->LID(fgid);

      // get value of corresponding velocity component
      double velocity = (*fluidvel)[flid];
      // insert velocity value into node-based vector
      convel_->ReplaceMyValue(lnodeid, index, velocity);

      if (subgridswitch)
      {
        // get subgrid-scale velocity component
        double sgvelocity = (*fluidsgvelvisc)[flid];
        // insert subgrid-scale velocity value into node-based vector
        sgvel_->ReplaceMyValue(lnodeid, index, sgvelocity);
      }
    }

    if (subgridswitch)
    {
      // now we transfer "pressure" dofs
      // remark: the "pressure" dof is used to store the nodal values of
      // the subgrid viscosity

      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[numdim];
      const int flid = fluiddofrowmap->LID(fgid);

      // get subgrid viscosity for this processor's local fluid dof and
      // divide by turbulent Prandtl number to get subgrid diffusivity
      subgriddiff  = (*fluidsgvelvisc)[flid]/tpn_;
      // insert subgrid diffusivity value into node-based vector
      subgrdiff_->ReplaceMyValues(1,&subgriddiff,&lnodeid);
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index=numdim; index < 3; ++index)
    {
      convel_->ReplaceMyValue(lnodeid, index, 0.0);
      if (subgridswitch) sgvel_->ReplaceMyValue(lnodeid, index, 0.0);
    }
  }

  // initial velocity field has now been set
  if (step_ == 0) initialvelset_ = true;

  return;

} // ScaTraTimIntImpl::SetVelocityField


/*----------------------------------------------------------------------*
 | apply moving mesh data                                     gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> dispnp,
    Teuchos::RCP<DRT::Discretization> fluiddis
)
{
  if (isale_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA: apply mesh movement");

    if (dispnp == Teuchos::null) dserror("Got null pointer for displacements");

    // get dofrowmap of fluid discretization
    const Epetra_Map* fluiddofrowmap = fluiddis->DofRowMap();

    // loop over all local nodes of scatra discretization
    for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      // Here we rely on the fact that the scatra discretization
      // is a clone of the fluid mesh. => a scatra node has the same
      // local (and global) ID as its corresponding fluid node!

      // get the processor's local fluid node with the same lnodeid
      DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
      // get the degrees of freedom associated with this fluid node
      vector<int> fluidnodedofs = fluiddis->Dof(fluidlnode);
      // determine number of space dimensions (numdof - pressure dof)
      const int numdim = ((int) fluidnodedofs.size()) -1;

      // now we transfer velocity dofs only
      for(int index=0;index < numdim; ++index)
      {
        // global and processor's local fluid dof ID
        const int fgid = fluidnodedofs[index];
        const int flid = fluiddofrowmap->LID(fgid);

        // get value of corresponding velocity component
        double disp = (*dispnp)[flid];
        // insert velocity value into node-based vector
        dispnp_->ReplaceMyValue(lnodeid, index, disp);
      }

      // for security reasons in 1D or 2D problems:
      // set zeros for all unused velocity components
      for (int index=numdim; index < 3; ++index)
      {
        dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
      }

    } // for lnodid
  } // if (isale_)

  return;

} // ScaTraTimIntImpl::SetDisplacementField


/*----------------------------------------------------------------------*
 |  set initial field for phi                                 gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetInitialField(int init, int startfuncno)
{
  if (init == 0) // zero_field
  {
    phin_-> PutScalar(0.0);
    phinp_-> PutScalar(0.0);
  }
  else if (init == 1 || init == 3)  // (disturbed_)field_by_function
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function
        double initialval = DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(k,lnode->X(),0.0,NULL);
        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }

    // add random perturbation for initial field of turbulent flows
    if(init==3)
    {
      int err = 0;

      // random noise is relative to difference of max-min values of initial profile
      double perc = params_->sublist("TURBULENCE PARAMETERS").get<double>("CHAN_AMPL_INIT_DIST",0.1);

      // out to screen
      if (myrank_==0)
      {
        cout << "Disturbed initial scalar profile:   max. " << perc*100 << "% random perturbation\n";
        cout << "\n\n";
      }

      // get overall max and min values and range between min and max
      double maxphi(0.0);
      double minphi(0.0);
      err = phinp_->MaxValue(&maxphi);
      if (err > 0) dserror("Error during evaluation of maximum value.");
      err = phinp_->MinValue(&minphi);
      if (err > 0) dserror("Error during evaluation of minimum value.");
      double range = abs(maxphi - minphi);

      // disturb initial field for all degrees of freedom
      for (int k=0; k < phinp_->MyLength(); ++k)
      {
        double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);
        double noise = perc * range * randomnumber;
        err += phinp_->SumIntoMyValues(1,&noise,&k);
        err += phin_ ->SumIntoMyValues(1,&noise,&k);
        if (err!=0) dserror("Error while disturbing initial field.");
      }
    }
  }
  else if (init==2) // field_by_condition
  {
    // access the initial field condition
    vector<DRT::Condition*> cond;
    discret_->GetCondition("InitialField", cond);

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    for (unsigned i=0; i<cond.size(); ++i)
    {
      cout<<"Applied InitialField Condition "<<i<<endl;

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        // get the processor local node
        DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

        vector<DRT::Condition*> mycond;
        lnode->GetCondition("InitialField",mycond);

        if (mycond.size()>0)
        {
          // the set of degrees of freedom associated with the node
          vector<int> nodedofset = discret_->Dof(lnode);

          int numdofs = nodedofset.size();
          for (int k=0;k< numdofs;++k)
          {
            // set 1.0 as initial value if node belongs to condition
            double phi0 = 1.0;
            // set initial value
            const int dofgid = nodedofset[k];
            int doflid = dofrowmap->LID(dofgid);
            phin_->ReplaceMyValues(1,&phi0,&doflid);
            // initialize also the solution vector. These values are a pretty good guess for the
            // solution after the first time step (much better than starting with a zero vector)
            phinp_->ReplaceMyValues(1,&phi0,&doflid);
          }
        }
      }
    }
  }
  // discontinuous 0-1 field for progress variable in 1-D
  else if (init == 4)
  {
    // check whether it is indeed a progress-variable approach
    if (reaction_ != "Arrhenius_pv") dserror("This initial field is supposed to be used only in the context of a reactive flow with progress-variable approach!");

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get coordinate
      const double x = lnode->X()[0];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        double initialval = 0.0;
        if (x > -EPS10) initialval = 1.0;

        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }
  }
  // analytical reactive profile in x2-direction due to Ferziger and Echekki (1993)
  // for two-dimensional flame-vortex interaction problem (x2=0-200)
  else if (init == 5)
  {
    // check whether it is indeed a progress-variable approach
    if (reaction_ != "Arrhenius_pv") dserror("This initial field is supposed to be used only in the context of a reactive flow with progress-variable approach!");

    // get flame parameter beta and diffusive flame thickness
    ParameterList eleparams;
    eleparams.set("action","get_flame_parameters");
    eleparams.set("isale",isale_);
    discret_->Evaluate(eleparams,null,null,null,null,null);
    const double beta  = eleparams.get("flame parameter beta", 1.0);
    const double delta = eleparams.get("diffusive flame thickness", 0.0);
    if (delta < EPS15) dserror("Diffusive flame thickness is zero or negative!");

    // define flame location in x2-direction
    const double floc = 100.0;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x2-coordinate
      const double x2 = lnode->X()[1];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        double initialval = 0.0;
        if (x2 < floc-EPS10) initialval = (1.0-(1.0/beta))*exp((x2-floc)/delta);
        else                 initialval = 1.0-(exp((1.0-beta)*(x2-floc)/delta)/beta);

        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }
  }
  else
    dserror("Unknown option for initial field: %d", init);

  return;
} // ScaTraTimIntImpl::SetInitialField


/*----------------------------------------------------------------------*
 | prepare Krylov space projection                            gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareKrylovSpaceProjection()
{
  if (project_)
  {
    if (isale_ or (step_<=1)) // fixed grid: compute w_,c_ only once!
    {
      vector<DRT::Condition*> KSPcond;
      discret_->GetCondition("KrylovSpaceProjection",KSPcond);
      int nummodes = KSPcond.size();

      for (int imode = 0; imode < nummodes; ++imode)
      {
        // zero w and c completely
        if (imode == 0)
        {
          w_->PutScalar(0.0);
          c_->PutScalar(0.0);
        }

        // in this case, we want to project out some zero pressure modes
        const string* definition = KSPcond[imode]->Get<string>("weight vector definition");

        if(*definition == "pointvalues")
        {
          dserror("option pointvalues not implemented");
        }
        else if(*definition == "integration")
        {
          // get rigid body modes
          const vector<double>* mode = KSPcond[imode]->Get<vector<double> >("mode");

          ParameterList mode_params;

          // set parameters for elements
          mode_params.set("action","integrate_shape_functions");

          int numdof = 0;
          Epetra_IntSerialDenseVector dofids(6);
          for(int rr=0;rr<6;rr++)
          {
            if(abs((*mode)[rr])>1e-14)
            {
              numdof++;
              dofids(rr)=rr;
            }
            else
              dofids(rr)=-1;
          }

        //  if (numdof > 1) dserror("only basis vectors with 1 dof active are allowed");
          mode_params.set("dofids",dofids);

          mode_params.set("isale",isale_);
          if (isale_)
            AddMultiVectorToParameterList(mode_params,"dispnp",dispnp_);

          /* evaluate KrylovSpaceProjection condition in order to get
    // integrated nodal basis functions w_
    // Note that in the case of definition integration based,
    // the average pressure will vanish in an integral sense
    //
    //                    /              /                      /
    //   /    \          |              |  /          \        |  /    \
    //  | w_*p | = p_i * | N_i(x) dx =  | | N_i(x)*p_i | dx =  | | p(x) | dx = 0
    //   \    /          |              |  \          /        |  \    /
    //                   /              /                      /
           */

          // get an RCP of the current column Epetra_Vector of the MultiVector
          Teuchos::RCP<Epetra_Vector> wi = rcp((*w_)(imode),false);

          // compute integral of shape functions
          discret_->EvaluateCondition
              (mode_params       ,
              Teuchos::null      ,
              Teuchos::null      ,
              wi                 ,
              Teuchos::null      ,
              Teuchos::null      ,
              "KrylovSpaceProjection");

          // set the current kernel basis vector
          for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
          {
            DRT::Node* node = discret_->lRowNode(inode);
            vector<int> gdof = discret_->Dof(node);
            int numdof = gdof.size();
            if (numdof > 6) dserror("only up to 6 dof per node supported");
            for(int rr=0;rr<numdof;++rr)
            {
              const double val = (*mode)[rr];
              int err = c_->ReplaceGlobalValue(gdof[rr],imode,val);
              if (err != 0) dserror("error while inserting value into c_");
            }
          }
        }
        else
        {
          dserror("unknown definition of weight vector w for restriction of Krylov space");
        }

      } // loop over nummodes
    }
  }

  return;

} // ScaTraTimIntImpl::PrepareKrylovSpaceProjection


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                   gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl()
{
  return;
}


#endif /* CCADISCRET       */
