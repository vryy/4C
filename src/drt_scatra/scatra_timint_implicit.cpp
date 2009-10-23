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
#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"
//REINHARD
#include "../drt_geometry/element_volume.H"
//end REINHARD

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
    RCP<ParameterList>            extraparams,
    RCP<IO::DiscretizationWriter> output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  extraparams_(extraparams),
  output_ (output),
  myrank_ (discret_->Comm().MyPID()),
  time_   (0.0),
  step_   (0),
  prbtype_  (extraparams_->get<string>("problem type")),
  solvtype_ (Teuchos::getIntegralValue<INPAR::SCATRA::SolverType>(*params,"SOLVERTYPE")),
  isale_    (extraparams_->get<bool>("isale")),
  scatratype_  (Teuchos::getIntegralValue<INPAR::SCATRA::ScaTraType>(*params,"SCATRATYPE")),
  reinitaction_(INPAR::SCATRA::reinitaction_none),
  masscalc_    (INPAR::SCATRA::masscalc_none),
  reinitswitch_(false),
  stepmax_  (params_->get<int>("NUMSTEP")),
  maxtime_  (params_->get<double>("MAXTIME")),
  timealgo_ (Teuchos::getIntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(*params,"TIMEINTEGR")),
  upres_    (params_->get<int>("UPRES")),
  uprestart_(params_->get<int>("RESTARTEVRY")),
  writeflux_(Teuchos::getIntegralValue<INPAR::SCATRA::FluxType>(*params,"WRITEFLUX")),
  outmean_  (Teuchos::getIntegralValue<int>(*params,"OUTMEAN")),
  dta_      (params_->get<double>("TIMESTEP")),
  dtp_      (params_->get<double>("TIMESTEP")),
  cdvel_    (Teuchos::getIntegralValue<INPAR::SCATRA::VelocityField>(*params,"VELOCITYFIELD")),
  convform_ (Teuchos::getIntegralValue<INPAR::SCATRA::ConvForm>(*params,"CONVFORM")),
  neumanninflow_(Teuchos::getIntegralValue<int>(*params,"NEUMANNINFLOW")),
  fssgd_    (Teuchos::getIntegralValue<INPAR::SCATRA::FSSUGRDIFF>(*params,"FSSUGRDIFF")),
  frt_      (0.0),
  tpn_      (1.0),
  errfile_  (extraparams_->get<FILE*>("err file")),
  initialvelset_(false),
  lastfluxoutputstep_(-1)
{
  // -------------------------------------------------------------------
  // determine whether linear incremental or nonlinear solver
  // -------------------------------------------------------------------
  switch(solvtype_)
  {
  case INPAR::SCATRA::solvertype_nonlinear:
  {
    incremental_ = true;
    nonlinear_   = true;
  }
  break;
  case INPAR::SCATRA::solvertype_linear_incremental:
  {
    incremental_ = true;
    nonlinear_   = false;
  }
  break;
  case INPAR::SCATRA::solvertype_linear_full:
  {
    incremental_ = false;
    nonlinear_   = false;
  }
  break;
  default:
    dserror("Received illegal scatra solvertype enum.");
  }

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
    // set up a species-temperature splitter (if more than one scalar)
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_-1,splitter_);
  }

  if (scatratype_ == INPAR::SCATRA::scatratype_levelset)
  {
    reinitaction_ = Teuchos::getIntegralValue<INPAR::SCATRA::ReinitializationAction>(params_->sublist("LEVELSET"),"REINITIALIZATION");
    masscalc_     = Teuchos::getIntegralValue<INPAR::SCATRA::MassCalculation>(params_->sublist("LEVELSET"),"MASSCALCULATION");
  }

  if (Teuchos::getIntegralValue<int>(*params_,"BLOCKPRECOND"))
  {
    // we need a block sparse matrix here
    if (prbtype_ != "elch")
      dserror("Block-Preconditioning is only for ELCH problems");
    // initial guess for non-zeros per row: 27 neighboring nodes for hex8 times (numscal_+1) dofs
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    Teuchos::RCP<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(splitter_,splitter_,27*(numscal_+1),false,true));
    blocksysmat->SetNumScal(numscal_);

    sysmat_ = blocksysmat;
  }
  else
  {
    // initialize standard (stabilized) system matrix (and save its graph!)
    // in standard case, but do not save the graph if fine-scale subgrid
    // diffusivity is used in non-incremental case
    if (fssgd_ != INPAR::SCATRA::fssugrdiff_no and not incremental_)
         sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));
    else sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
  }

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1 and n
  phinp_ = LINALG::CreateVector(*dofrowmap,true);
  phin_  = LINALG::CreateVector(*dofrowmap,true);

  // temporal solution derivative at time n+1
  phidtnp_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  hist_ = LINALG::CreateVector(*dofrowmap,true);

  // convective velocity (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  convel_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // acceleration and pressure required for computation of subgrid-scale
  // velocity (always four components per node)
  accpre_ = rcp(new Epetra_MultiVector(*noderowmap,4,true));

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
  // ensure that the Transport string was removed from conditions
  // -------------------------------------------------------------------
  {
    DRT::Condition* cond = discret_->GetCondition("TransportDirichlet");
    if (cond) dserror("Found a Transport Dirichlet condition. Remove Transport string!");
    cond = discret_->GetCondition("TransportNeumann");
    if (cond) dserror("Found a Transport Neumann condition. Remove Transport string!");
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
  ParameterList * turbparams =&(extraparams_->sublist("TURBULENCE PARAMETERS"));

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
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
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
    // scalar and velocity increment at time n+1
    phiincnp_ = LINALG::CreateVector(*dofrowmap,true);
    //velincnp_  = rcp(new Epetra_MultiVector(*noderowmap,3,true));

    // potential turbulence model
    if (turbparams->get<string>("PHYSICAL_MODEL") != "no_model")
      turbmodel_ = true;

    // warning No. 1: if classical (all-scale) turbulence model other than
    // constant-coefficient Smagorinsky is intended to be used
    if (turbparams->get<string>("PHYSICAL_MODEL") == "Smagorinsky_with_van_Driest_damping" or
        turbparams->get<string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
      dserror("No classical (all-scale) turbulence model other than constant-coefficient Smagorinsky model currently possible!");

    // warning No. 2: if classical (all-scale) turbulence model and fine-scale
    // subgrid-viscosity approach are intended to be used simultaneously
    if (turbmodel_ and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
      dserror("No combination of classical (all-scale) turbulence model and fine-scale subgrid-diffusivity approach currently possible!");

    // Smagorinsky constant
    Cs_ = turbparams->get<double>("C_SMAGORINSKY",0.0);

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
  SetInitialField(Teuchos::getIntegralValue<INPAR::SCATRA::InitialField>(*params_,"INITIALFIELD"),
      params_->get<int>("INITFUNCNO"));

  // initializes variables for natural convection (ELCH) if necessary
  SetupElchNatConv();

  // screen output (has to come after SetInitialField)
  if (prbtype_ == "elch")
  {
    frt_ = 96485.3399/(8.314472 * extraparams_->get<double>("TEMPERATURE"));

    double sigma = ComputeConductivity(); // every processor has to do this call
    if (myrank_==0)
    {
      cout<<"\nSetup of splitter: numscal = "<<numscal_<<endl;
      cout<<"Temperature value T (Kelvin)     = "<<extraparams_->get<double>("TEMPERATURE")<<endl;
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

//REINHARD
    // -------------------------------------------------------------------
    //                reinitialize the level set function
    // -------------------------------------------------------------------
    if ((scatratype_ == INPAR::SCATRA::scatratype_levelset) and (reinitaction_ != INPAR::SCATRA::reinitaction_none))
      Reinitialize();
//end REINHARD

    // -------------------------------------------------------------------
    //                update time derivative after solution
    // -------------------------------------------------------------------
    UpdateTimeDerivative();

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

//REINHARD
    // -------------------------------------------------------------------
    // compute mass loss for level set function and print it to display
    // -------------------------------------------------------------------
    if ((scatratype_ == INPAR::SCATRA::scatratype_levelset) and (masscalc_ != INPAR::SCATRA::masscalc_none))
      CalculateMassLoss();
//end REINHARD

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
    if (initialvelset_) PrepareFirstTimeStep();
    else dserror("Initial velocity field has not been set");
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
  ApplyDirichletBC(time_,phinp_,phidtnp_);
  ApplyNeumannBC(time_,phinp_,neumann_loads_);

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if (step_==1 and fssgd_ != INPAR::SCATRA::fssugrdiff_no) AVM3Preparation();

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
  const bool   isadapttol    = (getIntegralValue<int>(params_->sublist("NONLINEAR"),"ADAPTCONV"));
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
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid diffusivity",fssgd_);
  eleparams.set("turbulence model",turbmodel_);
  eleparams.set("Smagorinsky constant",Cs_);
  eleparams.set("turbulent Prandtl number",tpn_);
  eleparams.set("frt",frt_);// ELCH specific factor F/RT

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set type of scalar transport problem
  eleparams.set("scatratype",scatratype_);

  // set switch for reinitialization
  eleparams.set("reinitswitch",reinitswitch_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  // set vector values needed by elements
  discret_->ClearState();
  if (turbmodel_) discret_->SetState("subgrid diffusivity",subgrdiff_);

  // AVM3 separation
  if (incremental_ and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
  {
    discret_->SetState("subgrid diffusivity",subgrdiff_);
    AVM3Separation();
  }

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(eleparams);

  // call loop over elements with subgrid-diffusivity(-scaling) vector
  discret_->Evaluate(eleparams,sysmat_,null,residual_,subgrdiff_,null);
  discret_->ClearState();

  // AVM3 scaling
  if (not incremental_ and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    AVM3Scaling(eleparams);

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
    if (writeflux_!=INPAR::SCATRA::flux_no) OutputFlux();

    // write mean values of scalar(s)
    if (outmean_)
    {
      OutputMeanScalars();
      OutputElectrodeInfo();
    }
  }
  else
  {
    // calculation of statistics for normal fluxes (no output to file)
    if (step_>=samstart_ and step_<=samstop_ and writeflux_!=INPAR::SCATRA::flux_no) CalcFlux();
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
  if (cdvel_ == INPAR::SCATRA::velocity_zero) // zero
    convel_->PutScalar(0); // just to be sure!
  else if ((cdvel_ == INPAR::SCATRA::velocity_function))  // function
  {
    const int numdim = 3; // the velocity field is always 3D
    const int velfuncno = params_->get<int>("VELFUNCNO");
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value=DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);
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
 | set convective velocity field (+ pressure and acceleration field     |
 | if required)                                               gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(
Teuchos::RCP<const Epetra_Vector> fluidvel,
Teuchos::RCP<const Epetra_Vector> fluidacc,
Teuchos::RCP<const DRT::DofSet> dofset,
Teuchos::RCP<DRT::Discretization> fluiddis)
{
  if (cdvel_ != INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Wrong SetVelocityField() called for velocity field type %d!",cdvel_);

  TEUCHOS_FUNC_TIME_MONITOR("SCATRA: set convective velocity field");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid mesh
  // at the beginning, the fluid nodal distribution can have changed meanwhile
  // (e.g., caused by periodic boundary conditions applied only on the fluid side)!
  // We have to be sure, that everything is still matching.
  if (not fluiddis->NodeRowMap()->SameAs(*(discret_->NodeRowMap())))
    dserror("Fluid and Scatra noderowmaps are NOT identical. Emergency!");
#endif

  // boolean indicating whether acceleration vector exists
  // -> if yes, subgrid-scale velocity may need to be computed on element level
  bool sgvelswitch = fluidacc != Teuchos::null;

  // loop over all local nodes of scatra discretization
  for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // Here we rely on the fact that the scatra discretization
    // is a clone of the fluid mesh. => a scatra node has the same
    // local (and global) ID as its corresponding fluid node!

    // get the processor's local fluid node with the same lnodeid
    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);

    // care for the slave nodes of rotationally symm. periodic boundary conditions
    double rotangle(0.0);
    bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(fluidlnode,rotangle);

    // get the degrees of freedom associated with this fluid node
    vector<int> fluidnodedofs;
    if (dofset == Teuchos::null)
      fluidnodedofs = fluiddis->Dof(fluidlnode);
    else // ask a different dofset e.g. for a XFEM fluid
      fluidnodedofs = (*dofset).Dof(fluidlnode);

    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;

    // now we transfer velocity dofs only
    for(int index=0;index < numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[index];
//      const int flid = fluiddofrowmap->LID(fgid);
      const int flid = fluidvel->Map().LID(fgid);
      if (flid < 0) dserror("lid not found in map for given gid");

      // get value of corresponding velocity component
      double velocity = (*fluidvel)[flid];
      if (havetorotate)
      {
        // this is the desired component of the rotated vector field
        velocity = FLD::GetComponentOfRotatedVectorField(index,fluidvel,flid,rotangle);
      }

      // insert velocity value into node-based vector
      convel_->ReplaceMyValue(lnodeid, index, velocity);

      if (sgvelswitch)
      {
        // get value of corresponding acceleration component
        double acceleration = (*fluidacc)[flid];

        if (havetorotate)
        {
          // this is the desired component of the rotated vector field
          acceleration = FLD::GetComponentOfRotatedVectorField(index,fluidacc,flid,rotangle);
        }

        // insert acceleration value into node-based vector
        accpre_->ReplaceMyValue(lnodeid, index, acceleration);
      }
    }

    // now we transfer pressure dofs only
    if (sgvelswitch)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[numdim];
      // const int flid = fluiddofrowmap->LID(fgid);
      const int flid = fluidvel->Map().LID(fgid);
      if (flid < 0) dserror("lid not found in map for given gid");

      // get value of corresponding pressure component
      double pressure = (*fluidvel)[flid];
      // insert pressure value into node-based vector
      accpre_->ReplaceMyValue(lnodeid, numdim, pressure);
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index=numdim; index < 3; ++index)
    {
      convel_->ReplaceMyValue(lnodeid, index, 0.0);
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
void SCATRA::ScaTraTimIntImpl::SetInitialField(
    const INPAR::SCATRA::InitialField init,
    const int startfuncno)
{
  switch(init)
  {
  case INPAR::SCATRA::initfield_zero_field:
  {
    phin_-> PutScalar(0.0);
    phinp_-> PutScalar(0.0);
    break;
  }
  case INPAR::SCATRA::initfield_field_by_function:
  case INPAR::SCATRA::initfield_disturbed_field_by_function:
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
        double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(k,lnode->X(),0.0,NULL);
        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }

    // add random perturbation for initial field of turbulent flows
    if(init==INPAR::SCATRA::initfield_disturbed_field_by_function)
    {
      int err = 0;

      // random noise is relative to difference of max-min values of initial profile
      double perc = extraparams_->sublist("TURBULENCE PARAMETERS").get<double>("CHAN_AMPL_INIT_DIST",0.1);

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
    break;
  }
  case INPAR::SCATRA::initfield_field_by_condition:
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
    break;
  }
  // discontinuous 0-1 field for progress variable in 1-D
  case INPAR::SCATRA::initfield_DISCONTPV_1D:
  {
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
    break;
  }
  // analytical reactive profile in x2-direction due to Ferziger and Echekki (1993)
  // for two-dimensional flame-vortex interaction problem (x2=0-200)
  case INPAR::SCATRA::initfield_FVI_FERECHPRO:
  {
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
    break;
  }
  // initial mixture-fraction profile for Rayleigh-Taylor instability
  case INPAR::SCATRA::initfield_RAYTAYMIXFRAC:
  {
    // define interface thickness, sinusoidal disturbance wave amplitude and pi
    const double delta = 0.002;
    const double alpha = 0.001;
    const double pi = 3.141592654;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x1- and x2-coordinate
      const double x1 = lnode->X()[0];
      const double x2 = lnode->X()[1];

      // interface disturbance
      //double x2_int = 0.05*cos(pi*(x1+0.5));
      //double x2_int = 0.05*cos(2.0*pi*x1);
      double x2_int = 0.0;
      x2_int -= cos(4*pi*x1);
      x2_int -= cos(14*pi*x1);
      x2_int -= cos(23*pi*x1);
      x2_int -= cos(28*pi*x1);
      x2_int -= cos(33*pi*x1);
      x2_int -= cos(42*pi*x1);
      x2_int -= cos(51*pi*x1);
      x2_int -= cos(59*pi*x1);
      x2_int *= alpha;

      const double value = (x2_int-x2)/(2.0*delta);

      // values required for tanh-distribution
      const double vp = exp(value);
      const double vm = exp(-value);

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        // compute tanh-distribution
        double initialval = 0.0;
        initialval = 0.5*(1.0+(vp-vm)/(vp+vm));

        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }
    break;
  }
  default:
    dserror("Unknown option for initial field: %d", init);
  } // switch(init)

  return;
} // ScaTraTimIntImpl::SetInitialField


//REINHARD
/*----------------------------------------------------------------------*
 | do reinitialization by function
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Reinitialize()
{
	if (reinitaction_ == INPAR::SCATRA::reinitaction_directdistance)
	{
		//Reinitialization as distance calculation to some interfacial points
		//written for hex8 elements
		cout << "Reinitialize as distance calculation" << endl;
		const Epetra_Map* elemrowmap = discret_->ElementRowMap();
		const int nsd = 3;

		//first: a set of known interfacial points has to be computed
		std::vector<std::vector<double> > interpointscollection;
		//loop over all elements
		for (int j=0;j<discret_->NumMyColElements();++j)
		{
			//tell, whether this element is intersected
			//nearly identical to void COMBUST::RefinementCell::IdentifyIntersectionStatus()
			DRT::Element* actele = discret_->gElement(elemrowmap->GID(j));
			bool intersected = false;
			int counter = 0;
			int numnode = actele->NumNode();
		vector<int> lm(numnode);
		// get vector of node GIDs for this element
		const int* nodeids = actele->NodeIds();
		for (unsigned inode=0; inode < lm.size(); inode++)
			lm[inode] = nodeids[inode];
		// create vector "ephinp" holding scalar phi values for this element
		vector<double> ephinp(numnode);
		//do extraction at last
		DRT::UTILS::ExtractMyValues(*phinp_, ephinp, lm);

			if (ephinp[0] == 0.0)
				intersected = true;
			else if (ephinp[0] < 0.0)
			{
				while ((counter<(numnode-1)) and (intersected == false))
				{
					counter++;
					if (ephinp[counter]>=0.0)
						intersected = true;
				}
			}
			else if (ephinp[0]>0.0)
			{
				while((counter<(numnode-1)) and  (intersected == false))
				{
					counter++;
					if (ephinp[counter]<=0.0)
						intersected = true;
				}
			}

			//the following is only done, if the current element is intersected...
			//similar to void COMBUST::FlameFront::FindIntersectionPoints(...)

			if (intersected == true)
			{
				//get global coordinates of nodes of current elements
				LINALG::SerialDenseMatrix xyze(nsd,numnode);
				xyze = GEO::InitialPositionArray(actele);
				if(actele->Shape()!=DRT::Element::hex8)
					dserror("The reinitialization is only supported for hex8 elements");
				std::vector<std::vector<int> > lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex8);
				//loop over all lines
				for(std::size_t ilines=0; ilines<lines.size(); ilines++)
				{
					int node1 = lines[ilines][0];
					int node2 = lines[ilines][1];
					double phival1 = ephinp[node1];
					double phival2 = ephinp[node2];
					vector<double> coordinates(nsd);

					if (phival1*phival2 < 0)	//real intersection point
					{
						for (int dim=0; dim<nsd; dim++)
						{
							if(xyze(dim,node1)==xyze(dim,node2))
								coordinates[dim]=xyze(dim,node1);
							else 	//calculate intersectionpoint
								coordinates[dim]=xyze(dim,node1)-phival1/(phival2-phival1)*(xyze(dim,node2)-xyze(dim,node1));
						}
						//now vector coordinates contains global coordinates of intersection point
						//interpointscollection[interpointscollection.size()]=coordinates;
						//compare with existing points and write in interpointscollection

						//flag is true, if coordinates shall be written in interpointscollection
						bool flag = true;
						for (unsigned i=0; i<interpointscollection.size(); i++)
						{
							if (interpointscollection[i] == coordinates)
							{
								flag = false;
								break;
							}
						}
						if (flag == true)
							interpointscollection.push_back(coordinates);
					}
					else if ((phival1 == 0.0) or (phival2 == 0.0))
					{
						if (phival1 == 0.0)
						{
							for (int dim=0; dim<nsd; dim++)
								coordinates[dim]=xyze(dim,node1);
							bool flag = true;
							for (unsigned i=0; i<interpointscollection.size(); i++)
							{
								if (interpointscollection[i] == coordinates)
								{
									flag = false;
									break;
								}
							}
							if (flag == true)
								interpointscollection.push_back(coordinates);
						}
						if (phival2 == 0.0)
						{
							for (int dim=0; dim<nsd; dim++)
								coordinates[dim]=xyze(dim,node2);
							bool flag = true;
							for (unsigned i=0; i<interpointscollection.size(); i++)
							{
								if (interpointscollection[i] == coordinates)
								{
									flag = false;
									break;
								}
							}
							if (flag == true)
								interpointscollection.push_back(coordinates);

						}
					}
					//else do nothing
				}//end loop over all lines
			}//end of calculation of intersection points if intersected == true
		}//end loop over all elements for calculation of intersection points

		//now, intersection points are in std::vector<std::vector<double> > interpointscollection
		//and reinitialization for all points shall be done

		//loop over all nodes
		for (int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
		{
			//get current node
			DRT::Node* actnode = discret_->lRowNode(lnodeid);
			const double* nodecoords = actnode->X();

			//new_phi will become new signed distance function
			vector<double> new_phi(1);
			//distance to current intersection point
			double curdistance = 0.0;
			for (int k=0; k<nsd; k++)
				curdistance = curdistance + DSQR(nodecoords[k]-interpointscollection[0][k]);
			curdistance = sqrt(curdistance);
			new_phi[0] = curdistance;

			for (unsigned j=1; j<interpointscollection.size(); j++)
			{
				curdistance = 0;
				for (int k=0; k<nsd; k++)
					curdistance = curdistance + DSQR(nodecoords[k]-interpointscollection[j][k]);
				curdistance = sqrt(curdistance);
				if (curdistance<new_phi[0])
					new_phi[0] = curdistance;
			}
			//the right sign still has to be adjusted...
			vector<double> phi_old(1);
			vector<int> lm_node(1);
			lm_node[0] = lnodeid;
			DRT::UTILS::ExtractMyValues(*phinp_,phi_old,lm_node);

			if (phi_old[0]<0.0)
				new_phi[0] = -new_phi[0];

			phinp_->ReplaceMyValues(1,&new_phi[0],&lnodeid);

		}//end loop over all nodes

	}//end of reinitialization with distance calculation to some interfacial points

	else if (reinitaction_ == INPAR::SCATRA::reinitaction_sussman)
	{
		//Reinitialization according to Sussman 1994 or the Hamilton-Jacobi formulation in Marchandise
		cout << "Reinitialize according to Sussman" << endl;
		//compare suggestion of Sussman: timestepsize = 0.9*characteristic length
		double timestepsize_reinit = 0.009;

		RCP<Epetra_Vector> phi_old = Teuchos::rcp(new Epetra_Vector(*phin_));
		RCP<Epetra_Vector> phidtn_old = Teuchos::rcp(new Epetra_Vector(*phidtn_));

		//set a new time step size and save the old time step size
		// this is up to now directly adapted for my problem of Zalesaks disk
		double timestepsize_old = params_->get<double> ("time step size");
		params_->set<double> ("time step size" , timestepsize_reinit);
		double dta_save = dta_;
		double dtp_save = dtp_;
		dta_ = timestepsize_reinit;
		dtp_ = timestepsize_reinit;


		//maximal number the time loop is done
		int nmax = 3;
		int loopindex=0;
		//flag to decide on end of reinitialization
		bool reinitflag = false;

		//is switch finally for Sysmat and BodyForce, as slightly different actions have to be done
		//in this reinitialization compared to the normal level set solving
		reinitswitch_ = true;

		//time loop for reinitialization
		//for (int loopindex=0; loopindex<nmax; loopindex++)
		while ((loopindex<nmax) and (reinitflag==false))
		{
			Update();

			//set phidtn_ to zero
			phidtn_->PutScalar(0.0);

			SetReinitVelocityField();
			SetOldPartOfRighthandside();

			if(nonlinear_) NonlinearSolve();
			else Solve();

			dtp_=dta_;

			//now, phinp_ and phin_ are compared to decide whether another reinit. step is necessary
			//sum of abs values of changes during this reinitialization step
			double change = 0.0;
			vector<double> ephinp(discret_->NumMyRowNodes());
			vector<double> ephin(discret_->NumMyRowNodes());
			vector<int> lm(discret_->NumMyRowNodes());
			for (int lnodeid=0; lnodeid<discret_->NumMyRowNodes();lnodeid++)
				lm[lnodeid]=lnodeid;
			DRT::UTILS::ExtractMyValues(*phinp_,ephinp,lm);
			DRT::UTILS::ExtractMyValues(*phin_,ephin,lm);
			for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
				change = change + abs(ephinp[lnodeid]-ephin[lnodeid]);

			if (change<1.0)
				reinitflag = true;

			loopindex = loopindex + 1;
		} //end time loop for reinitialization

		//set back phin_ and phidtn_
		*phin_ = *phi_old;
		*phidtn_ = *phidtn_old;

		//set back the param "Type of reinitialization" in sublist LEVELSET
		reinitswitch_ = false;

		//end reinitialization by setting back the velocity field
		SetVelocityField();
		//set back the time step
		params_->set<double> ("time step size" , timestepsize_old);
		dta_ = dta_save;
		dtp_ = dtp_save;


	}	//end of Reinitialization according to Sussman 1994

	else if (reinitaction_ == INPAR::SCATRA::reinitaction_interfaceprojection)
	{
		//Reinitialization according to Parolini/Burman
		//here, extrapolation is still done for demonstrating its structure
		//Reinitialization by function according to:
		//A local projection reinitialization procedure for the level set equation on unstructured grids
		//by Nicola Parolini and Erik Burman
		cout << "Reinitialize as InterfaceProjection according to Parolini" << endl;
		const Epetra_Map* dofrowmap = discret_->DofRowMap();
		// phi at times n+1
		RCP<Epetra_Vector> phinpnew_ = LINALG::CreateVector(*dofrowmap,true);
		phinpnew_->PutScalar(0.0);

		//calculation of expol
		Epetra_SerialDenseMatrix expol(8,8);

		double sq3=sqrt(3.0);
		expol(0,0)=1.25+0.75*sq3;
		expol(0,1)=-0.25-0.25*sq3;
		expol(0,2)=-0.25+0.25*sq3;
		expol(0,3)=-0.25-0.25*sq3;
		expol(0,4)=-0.25-0.25*sq3;
		expol(0,5)=-0.25+0.25*sq3;
		expol(0,6)=1.25-0.75*sq3;
		expol(0,7)=-0.25+0.25*sq3;
		expol(1,1)=1.25+0.75*sq3;
		expol(1,2)=-0.25-0.25*sq3;
		expol(1,3)=-0.25+0.25*sq3;
		expol(1,4)=-0.25+0.25*sq3;
		expol(1,5)=-0.25-0.25*sq3;
		expol(1,6)=-0.25+0.25*sq3;
		expol(1,7)=1.25-0.75*sq3;
		expol(2,2)=1.25+0.75*sq3;
		expol(2,3)=-0.25-0.25*sq3;
		expol(2,4)=1.25-0.75*sq3;
		expol(2,5)=-0.25+0.25*sq3;
		expol(2,6)=-0.25-0.25*sq3;
		expol(2,7)=-0.25+0.25*sq3;
		expol(3,3)=1.25+0.75*sq3;
		expol(3,4)=-0.25+0.25*sq3;
		expol(3,5)=1.25-0.75*sq3;
		expol(3,6)=-0.25+0.25*sq3;
		expol(3,7)=-0.25-0.25*sq3;
		expol(4,4)=1.25+0.75*sq3;
		expol(4,5)=-0.25-0.25*sq3;
		expol(4,6)=-0.25+0.25*sq3;
		expol(4,7)=-0.25-0.25*sq3;
		expol(5,5)=1.25+0.75*sq3;
		expol(5,6)=-0.25-0.25*sq3;
		expol(5,7)=-0.25+0.25*sq3;
		expol(6,6)=1.25+0.75*sq3;
		expol(6,7)=-0.25-0.25*sq3;
		expol(7,7)=1.25+0.75*sq3;

		for (int k=0;k<8;++k)
		{
			for (int l=0;l<k;++l)
				expol(k,l)=expol(l,k);
		}
		//end calculation of expol

		//loop over all nodes
		for (int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
		{
			//get current node
			DRT::Node* actnode = discret_->lRowNode(lnodeid);
			//get all adjacent elements to this node
			DRT::Element** elements=actnode->Elements();

			//variable, that will become the new signed distance value of this node
			double phi_new = 0;
			int numadnodes = actnode->NumElement();

			//gradient of phi at the one relevant node
			Epetra_SerialDenseMatrix gradphi(3,1);
			for (int j=0; j<3; j++)
				gradphi(j,0)=2.0;

			//loop over all adjacent elements
			for (int j=0; j<numadnodes; ++j)
			{
				//get current element
				DRT::Element* ele = elements[j];

				//number of nodes in current element
				int iel = ele->NumNode();

				// create vector "ephinp" holding scalar phi values for this element
				Epetra_SerialDenseMatrix ephinp(iel,1);

				//temporal vector necessary just for function ExtractMyValues...
				//that is requiring vectors and not matrices
				vector<double> etemp(iel);

				// remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
				// arguments, but ele->NodeIds delivers an "int*" argument
				vector<int> lm(iel);

				// get vector of node GIDs for this element
				const int* nodeids = ele->NodeIds();
				for (int inode=0; inode < iel; inode++)
					lm[inode] = nodeids[inode];

				// get entries in "gfuncvalues" corresponding to node GIDs "lm" and store them in "ephinp"
				DRT::UTILS::ExtractMyValues(*phinp_, etemp, lm);

				for (int k=0; k<iel; k++)
					ephinp(k,0) = etemp[k];

				const DRT::Element::DiscretizationType distype = ele->Shape();

				//gradphi is calculated by extrapolating values of gausspoints
				//LINALG::Matrix<3,1> xsi_;
				//current node is number iquad of current element
				int iquad = 0;
				for (int i=0; i<iel; i++)
				{
					if (actnode->Id()==nodeids[i])
						iquad = i;
				}

				//note: the calculation of deriv_ is only implemented for hex8-elements
				//calculation of deriv_ as weighted value over the gausspoints with extrapolation

				//matrix with the derivatives of the form functions (in paramspace)
				Epetra_SerialDenseMatrix gradphi_gps(8,3);

				// the following is necessary for getting the gaussian points
				const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_8point;
				const DRT::UTILS::IntegrationPoints3D  intpoints(gaussrule);

				//calculation of gradphi_gps(8,3)
				//loop over all gaussian points
				for (int i=0; i<8; i++)
				{
					//coordinates of current gaussian point
					LINALG::Matrix<3,1> xsi_;
					for (int idim = 0; idim<3; idim++)
						xsi_(idim) = intpoints.qxg[i][idim];

					Epetra_SerialDenseMatrix deriv_(3,8);
					DRT::UTILS::shape_function_3D_deriv1(deriv_,xsi_(0),xsi_(1),xsi_(2),distype);

					//calculation of xji_ (inverse of transposed jacobian)

					//xyze are the positions of the nodes in the global coordinate system
				Epetra_SerialDenseMatrix xyze(3,iel);
				GEO::fillInitialPositionArray(ele, xyze);

				// get transposed of the jacobian matrix d x / d \xi
				Epetra_SerialDenseMatrix xjm(3,3);		//should be the "this"-object
				//computing: this = 0*this+1*deriv_*(xyze)T
				xjm.Multiply('N','T',1.0,deriv_,xyze,0.0);

				// inverse of jacobian
				//xji = xjm^-1
				Epetra_SerialDenseMatrix xji(xjm);
					//for inverting: works with LINALG-matrices, therefore using a temporal matrix
				//just for inverting
				LINALG::Matrix<3,3> xjm_temp;
				for (int i1=0; i1<3; i1++)
					for (int i2=0; i2<3; i2++)
						xjm_temp(i1,i2) = xjm(i1,i2);
				LINALG::Matrix<3,3> xji_temp;

				const double det = xji_temp.Invert(xjm_temp);
				if (det < 1e-16)
					dserror("zero or negative jacobian determinant");

				for (int i1=0; i1<3; i1++)
					for (int i2=0; i2<3; i2++)
						xji(i1,i2) = xji_temp(i1,i2);

				Epetra_SerialDenseMatrix derxy_(3,iel);
					//compute global derivatives
				derxy_.Multiply('N','N',1.0,xji,deriv_,0.0);

					//gradient of phi at gaussian point shall be calculated in the end
				Epetra_SerialDenseMatrix gradphi_gp(3,1);
				gradphi_gp.Multiply('N','N',1.0,derxy_,ephinp,0.0);

				//writing this result in gradphi_gps(NUMGPT_SOH8,3)
					for (int i1=0; i1<3; i1++)
						gradphi_gps(i,i1)=gradphi_gp(i1,0);

				}	//end loop over all gaussian points
				//end calculation of gradphi_gps

				//gradphi at all nodes
				Epetra_SerialDenseMatrix gradphi_nodes(8,3);

				//now, gradphi is computed at all nodes
				//room for improvement, as only the values at one node (iquad) are required...
				gradphi_nodes.Multiply('N','N',1.0,expol,gradphi_gps,0.0);



				//here, the highly important selection of gradphi in an upwind sense is done
				//the sign of the single derivatives is here not important
				for (int i1=0; i1<3; i1++)
				{
					if (abs(gradphi_nodes(iquad,i1))<abs(gradphi(i1,0)))
						gradphi(i1,0)=gradphi_nodes(iquad,i1);
				}
			}//end loop over all adjacent elements

			//What has been achieved by now:
			//gradient in current element has been computed (actually not more...)

			//now compute d_h = phi(x)/|gradphi(x)| as new signed distance value

			double norm_grad = sqrt(gradphi(0,0)*gradphi(0,0)+gradphi(1,0)*gradphi(1,0)+gradphi(2,0)*gradphi(2,0));
			vector<double> etemp(1);
			vector<int> lm(1);
			lm[0]=lnodeid;
			DRT::UTILS::ExtractMyValues(*phinp_, etemp, lm);

			if (norm_grad == 0.0)
				//no idea, whether this is useful. Actually, this should never happen...
				phi_new = etemp[0];
			else
				phi_new = etemp[0]/norm_grad;

			//What has been achieved by now:
			//the new signed distance value has been computed for the current node

			//now phi_new has to be written in the discretization
				phinpnew_->ReplaceMyValues(1,&phi_new,&lnodeid);

				//the following is not done, as huge mistakes occur, when some values are already changed
				//and these are used for calculating the other values
				//therefore, calculation of all new values always on the basis of the old values
				//and the replacement of phin_ and phinp_
				phinp_->ReplaceMyValues(1,&phi_new,&lnodeid);

		}//end loop over all nodes

	}//end of Reinitialization according to Parolini/Burman

	return;
} //SCATRA::ScaTraTimIntImpl::Reinitialize()
//end REINHARD


//REINHARD
/*----------------------------------------------------------------------*
 | calculate mass loss of Zalesaks disk
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalculateMassLoss()
{
	//basic idea: loop over all elements and calculate area of inner region
	//important: assume discretization with nice rectangular hex8-elements, one element as height
	//therefore a pseudo-3d problem...

	int numele = discret_->NumGlobalElements();

	//length of one hex8-element (here directly designed for Zalesaks disk!!)
	double h = 1/sqrt(numele);
	double area = 0;

	if (masscalc_ == INPAR::SCATRA::masscalc_squares)
	{
		//first alternative for calculating the mass loss
		//here, a loop over all nodes is done. For each node with negative phi-value, the area of a square around this
		//node is added to the inner area
		cout << "In Massenberechnung Squares" << endl;
		double h_h = h*h;
		//loop all nodes
		for (int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
		{
	      		// get pointer to current node
	      		DRT::Node* node = discret_->lRowNode(lnodeid);
	      		// build vector of all node GIDs of this element
	      		//lm[inode] = nodeid;
	      		const double* pos = node->X();
	      		//we are only interested in half of the nodes, that is for example
	      		//all nodes with negative z-value...
	      		if (pos[2]<=0)
	      		{
	    	  		//control, whether phi of this element is negative
	    	  		vector<double> etemp(1);
	    	  		vector<int> lm(1);
	    	  		lm[0] = lnodeid;
	    	  		DRT::UTILS::ExtractMyValues(*phinp_, etemp, lm);
	    	  		if (etemp[0]<0)
	    	  		{
	    			 	// get number of adjacent elements
	    			 	// 3 possibilities: 1,2,4
	    		 	 	// get GID of this node
	    		 	 	int numadnodes = node->NumElement();
	    		 	 	switch (numadnodes)
	    		 	 	{
	    	  				case 1:
	    	  				{
	    	  					area = area + 0.25*h_h;
	    	  					break;
	    	  				}
	    	  				case 2:
	    	  				{
	    	  					area = area + 0.5*h_h;
	    	  					break;
	    	  				}
	    	  				case 4:
	    	  				{
	    	  					area = area + h_h;
	    	  					break;
	    	  					}
	    	  				default:
	    	  					dserror ("wrong number of adjacent elements!");
	    		  		}

	    	 		}
	      		}

		}
	}//end first alternative

	else if (masscalc_ == INPAR::SCATRA::masscalc_interpolated)
	{
		//here the mass is calculated by calculating interface points and interpolating between these points
		//the interpolation procedure is depending on the number of nodes with negative phi-values
		//the loop is now element based
		cout << "In Massenberechnung Interpolated" << endl;
		const Epetra_Map* elemrowmap = discret_->ElementRowMap();

		// loop all elements
		for (int i=0; i<numele; ++i)
		{
			//get current element
			DRT::Element* actele = discret_->gElement(elemrowmap->GID(i));
			// extract phi-values of current element
			// remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
			// arguments, but ele->NodeIds delivers an "int*" argument
			vector<int> lm(actele->NumNode());

	    		// get vector of node GIDs for this element
	    		const int* nodeids = actele->NodeIds();
	    		for (unsigned inode=0; inode < lm.size(); inode++)
	    			lm[inode] = nodeids[inode];

	    		// create vector "ephinp" holding scalar phi values for this element
	    		vector<double> ephinp(actele->NumNode());

	    		//do extraction at last
	    		DRT::UTILS::ExtractMyValues(*phinp_, ephinp, lm);

	    		//compare the values of those 4 nodes, that have got a negative z-value
	    		//first of all get vector with these four relevant phi-values
	    		vector<double> phi_values(4);			//directly constructed for Zalesaks disk
	    		int curnum = 0;
	    		for (int j=0; j<actele->NumNode(); j++)
	    		{
	    			DRT::Node* node = discret_->lRowNode(lm[j]);
	    			const double* pos = node->X();
	    			if (pos[2] < 0)
	    			{
	    				if (curnum == 4)
	    					dserror("Mistake when calculating the mass loss");
	    				phi_values[curnum]=ephinp[j];
	    				curnum = curnum+1;
	    			}

	    		}

			//number of nonpositive values of the first 4 values of phi_values (negative means inside)
			int numneg = 0;
			for (int j=0; j<4; ++j)
			{
				if (phi_values[j] <= 0)
					numneg = numneg + 1;
			}

			switch (numneg)
			{
				case 0: //element is completely in outer region; no contribution to area
				{
					break;
				}
				case 1:
				{
					//get number of the negative node (as k)
					int k = 0;
					for (int j=0; j<4; ++j)
					{
						if (phi_values[j] <= 0)
							k=j;
					}

					//compare value phi_values[k] with neighbor values and compute lengths of triangle
					double l1 = 0;
					double l2 = 0;
					if (k==3)
						l1 = -h*phi_values[3]/(phi_values[0]-phi_values[3]);
					else
						l1 = -h*phi_values[k]/(phi_values[k+1]-phi_values[k]);

					if (k==0)
						l2 = -h*phi_values[0]/(phi_values[3]-phi_values[0]);
					else
						l2 = -h*phi_values[k]/(phi_values[k-1]-phi_values[k]);

					area = area + 0.5*l1*l2;

					break;
				}
				case 2:
				{
					//get number of the negative nodes (as k1 and k2)
					int k1 = 0;
					int k2 = 0;
					for (int j=0; j<4; ++j)
					{
						if (phi_values[j] <= 0)
						{
							k1 = j;
							break;
						}
					}
					for (int j=3; j>=0; --j)
					{
						if (phi_values[j] <= 0)
						{
							k2 = j;
							break;
						}
					}

					//two cases: nonpositive nodes are neighbors or not
					//but it is always ensured, that k1<k2
					if ((k1+1==k2)||((k1==0)&&(k2==3)))		//they are neighbors
					{
						double l1 = 0;
						double l2 = 0;
						if (k1+1 == k2)
						{
							if (k1 == 0)
								l1 = -h*phi_values[0]/(-phi_values[0]+phi_values[3]);
							else
								l1 = -h*phi_values[k1]/(-phi_values[k1]+phi_values[k1-1]);
							if (k2 == 3)
								l2 = -h*phi_values[3]/(-phi_values[3]+phi_values[0]);
							else
								l2 = -h*phi_values[k2]/(-phi_values[k2]+phi_values[k2+1]);
						}
						else		//k1==0 and k2==3
						{
							l1 = -h*phi_values[0]/(-phi_values[0]+phi_values[1]);
							l2 = -h*phi_values[3]/(-phi_values[3]+phi_values[2]);
						}

						area = area + 0.5*h*(l1+l2);		//trapezoid
					}
					else 			//they are not neighbors
					{
						//assumption: inside area is not interrupted
						//would be possible in coarser grid
						double l1 = 0;
						double l2 = 0;
						double l3 = 0;
						double l4 = 0;
						//only two possibilities k1==0 and k2==2 or k1==1 and k2==3
						if (k1==0)		//and k2==2
						{
							l1 = -h*phi_values[0]/(-phi_values[0]+phi_values[1]);
							l2 = -h*phi_values[0]/(-phi_values[0]+phi_values[3]);
							l3 = -h*phi_values[2]/(-phi_values[2]+phi_values[1]);
							l4 = -h*phi_values[2]/(-phi_values[2]+phi_values[3]);
						}
						else			//k1==1 and k2==3
						{
							l1 = -h*phi_values[1]/(-phi_values[1]+phi_values[2]);
							l2 = -h*phi_values[1]/(-phi_values[1]+phi_values[0]);
							l3 = -h*phi_values[3]/(-phi_values[3]+phi_values[2]);
							l4 = -h*phi_values[3]/(-phi_values[3]+phi_values[0]);
						}
						area = area + h*h - (0.5*(h-l1)*(h-l3)+0.5*(h-l2)*(h-l4));
					}
					break;
				}
				case 3:			//actually quite similar to case 1
				{
					//get number of the positive node (as k)
					int k = 0;
					for (int j=0; j<4; ++j)
					{
						if (phi_values[j] > 0)
							k=j;
					}

					//compare value ephinp[k] with neighbor values and compute lengths of triangle
					double l1 = 0;
					double l2 = 0;
					if (k==3)
						l1 = h*phi_values[3]/(-phi_values[0]+phi_values[3]);
					else
						l1 = h*phi_values[k]/(-phi_values[k+1]+phi_values[k]);

					if (k==0)
						l2 = h*phi_values[0]/(-phi_values[3]+phi_values[0]);
					else
						l2 = h*phi_values[k]/(-phi_values[k-1]+phi_values[k]);

					area = area + h*h - 0.5*l1*l2;

					break;
				}
				case 4:
				{
					area = area + h*h;
					break;
				}
				default:
					dserror("problem with number of nonpositive values");
			}


		}//end loop over all elements
	}//end second alternative

	double area_old = 0.05822070305889007;
	double mass_loss = (area_old-area)/area_old;

	//Output to screen of mass loss (in percent of the original mass)
	std::cout << "area_old: " << area_old << std::endl;
	std::cout << "area: " << area << std::endl;
	std::cout << "mass loss: " << mass_loss <<std::endl;

	return;
} // ScaTraTimIntImpl::CalculateMassLoss
//end REINHARD


//REINHARD
/*--------------------------------------------------------------------------*
 | update the velocity field during reinitialization (according to Sussman) |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetReinitVelocityField()
{
    const int numdim = 3; // the velocity field is always 3D

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
		//get current node
		DRT::Node* actnode = discret_->lRowNode(lnodeid);
		//get all adjacent elements to this node
		DRT::Element** elements=actnode->Elements();

	    //phi-value at current node
	    double phiatnode = 0.0;

		double norm_gradphi=0.0;

	    //gradient of phi at the one relevant node
	    Epetra_SerialDenseMatrix gradphi(3,1);

	    for (int i=0; i<3; i++)
	    	gradphi(i,0) = 0.0;

		int numadnodes = actnode->NumElement();

		//loop over all adjacent elements
	    for (int elecur=0; elecur<numadnodes; ++elecur)
	    {
	    	Epetra_SerialDenseMatrix gradphi_node(3,1);

	    	//get current element
	    	const DRT::Element* ele = elements[elecur];

		    //number of nodes in current element
		    int iel = ele->NumNode();

		    // create vector "ephinp" holding scalar phi values for this element
		    Epetra_SerialDenseMatrix ephinp(iel,1);

		    //temporal vector necessary just for function ExtractMyValues...
		    //that is requiring vectors and not matrices
		    vector<double> etemp(iel);

		    // remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
		    // arguments, but ele->NodeIds delivers an "int*" argument
	    	vector<int> lm(iel);

		    // get vector of node GIDs for this element
		    const int* nodeids = ele->NodeIds();
		    for (int inode=0; inode < iel; inode++)
		    	lm[inode] = nodeids[inode];

		    // get entries in "gfuncvalues" corresponding to node GIDs "lm" and store them in "ephinp"
		    DRT::UTILS::ExtractMyValues(*phinp_, etemp, lm);

		    for (int k=0; k<iel; k++)
		    	ephinp(k,0) = etemp[k];

		    const DRT::Element::DiscretizationType distype = ele->Shape();

		    //current node is number iquad of current element
		    int iquad = 0;
		    for (int k=0; k<iel; k++)
		    {
		    	if (actnode->Id()==nodeids[k])
		    		iquad = k;
		    }
		    phiatnode = ephinp(iquad,0);

		    //local coordinates of current node
		    LINALG::Matrix<3,1> xsi;
		    LINALG::SerialDenseMatrix eleCoordMatrix=DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
		    for (int k=0; k<3; k++)
		    	xsi(k,0)=eleCoordMatrix(k,iquad);

		    Epetra_SerialDenseMatrix deriv(3,iel);
	    	DRT::UTILS::shape_function_3D_deriv1(deriv,xsi(0),xsi(1),xsi(2),distype);

		    //xyze are the positions of the nodes in the global coordinate system
            Epetra_SerialDenseMatrix xyze(3,iel);			//Epetra_SerialDenseMatrix
            const DRT::Node*const* nodes = ele->Nodes();
            for (int inode=0; inode<iel; inode++)
            {
                const double* x = nodes[inode]->X();
                xyze(0,inode) = x[0];
                xyze(1,inode) = x[1];
                xyze(2,inode) = x[2];
            }

            //get transposed of the jacobian matrix
            Epetra_SerialDenseMatrix xjm(3,3);
            //computing: this = 0*this+1*deriv_*(xyze)T
            xjm.Multiply('N','T',1.0,deriv,xyze,0.0);

            // inverse of jacobian
            //xji = xjm^-1
            Epetra_SerialDenseMatrix xji(xjm);
		    //for inverting: works with LINALG-matrices, therefore using a temporal matrix
            //just for inverting
            LINALG::Matrix<3,3> xjm_temp;
            for (int i1=0; i1<3; i1++)
            	for (int i2=0; i2<3; i2++)
            		xjm_temp(i1,i2) = xjm(i1,i2);
            LINALG::Matrix<3,3> xji_temp;

            const double det = xji_temp.Invert(xjm_temp);
            if (det < 1e-16)
            	dserror("zero or negative jacobian determinant");

            for (int i1=0; i1<3; i1++)
            	for (int i2=0; i2<3; i2++)
            		xji(i1,i2) = xji_temp(i1,i2);

            Epetra_SerialDenseMatrix derxy(3,iel);
		    //compute global derivatives
            derxy.Multiply('N','N',1.0,xji,deriv,0.0);

            gradphi_node.Multiply('N','N',1.0,derxy,ephinp,0.0);

		    double norm = sqrt(gradphi_node(0,0)*gradphi_node(0,0)+gradphi_node(1,0)*gradphi_node(1,0)+gradphi_node(2,0)*gradphi_node(2,0));

//		    //the element is chosen, where the gradient is closest to 1
//    		if (abs(norm-1.0)<abs(norm_gradphi-1.0))
//    		{
//    			norm_gradphi=norm;
//    			gradphi=gradphi_node;
//    		}

		    //a middle value is calculated of all adjacent elements to get second order accuracy
		    norm_gradphi = norm_gradphi+norm;
		    for (int k=0; k<3; k++)
		    	gradphi(k,0) = gradphi(k,0)+gradphi_node(k,0);
		    if (elecur==(numadnodes-1))
		    {
		    	norm_gradphi=norm_gradphi/numadnodes;
		    	for (int k=0; k<3; k++)
		    		gradphi(k,0) = gradphi(k,0)/numadnodes;
		    }
	    }//end loop over all adjacent elements

	    if (norm_gradphi == 0.0)
	    {
	    	norm_gradphi = 1e-9;
	    	cout << "Here, the computed gradient is 0.0" << endl;
	    }

	    double signum = 0.0;
	    double epsilon = 0.015;		//vgl. Sussman1994 (1.5*h)
    	if (phiatnode<-epsilon)
    		signum = -1.0;
    	else if (phiatnode>epsilon)
    		signum = 1.0;
    	else
    		signum = phiatnode/epsilon + sin(PI*phiatnode/epsilon)/PI;

	    for (int index=0; index<numdim; index++)
	    {
	    	double value = signum*gradphi(index,0)/norm_gradphi;
	    	convel_->ReplaceMyValue(lnodeid,index,value);
	    }


	}//end loop over all nodes

  return;

} // ScaTraImplicitTimeInt::SetReinitVelocityField
//end REINHARD


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

        // get rigid body modes
        const vector<double>* mode = KSPcond[imode]->Get<vector<double> >("mode");

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

        if(*definition == "pointvalues")
        {
          dserror("option pointvalues not implemented");
        }
        else if(*definition == "integration")
        {
          ParameterList mode_params;

          // set parameters for elements
          mode_params.set("action","integrate_shape_functions");
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
              (mode_params           ,
              Teuchos::null      ,
              Teuchos::null      ,
              wi                 ,
              Teuchos::null      ,
              Teuchos::null      ,
              "KrylovSpaceProjection");

        }
        else
        {
          dserror("unknown definition of weight vector w for restriction of Krylov space");
        }

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
