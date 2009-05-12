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
#include "../drt_fluid/fluid_utils.H" // for conpotsplitter

// for the condition writer output
/*
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
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
  stepmax_  (params_->get<int>("max number timesteps")),
  maxtime_  (params_->get<double>("total time")),
  timealgo_ (params_->get<INPAR::SCATRA::TimeIntegrationScheme>("time int algo")),
  upres_    (params_->get<int>("write solution every")),
  uprestart_(params_->get<int>("write restart every")),
  writeflux_(params_->get<string>("write flux")),
  dta_      (params_->get<double>("time step size")),
  dtp_      (params_->get<double>("time step size")),
  cdvel_    (params_->get<int>("velocity field")),
  convform_ (params_->get<string>("form of convective term")),
  neumannin_(params_->get<string>("Neumann inflow")),
  fssgd_    (params_->get<string>("fs subgrid diffusivity")),
  frt_      (96485.3399/(8.314472 * params_->get<double>("TEMPERATURE",298.15))),
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
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_,conpotsplitter_);
    if (myrank_==0)
    {
      cout<<"\nSetup of conpotsplitter: numscal = "<<numscal_<<endl;
      cout<<"Temperature value T (Kelvin)     = "<<params_->get<double>("TEMPERATURE",298.0)<<endl;
      cout<<"Constant F/RT                    = "<<frt_<<endl;
    }
  }

  if (params_->get<int>("BLOCKPRECOND",0) )
  {
    // we need a block sparse matrix here
    if (prbtype_ != "elch") 
      dserror("Block-Preconditioning is only for ELCH problems");
    // initial guess for non-zeros per row: 27 neighboring nodes for hex8 times (numscal_+1) dofs
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(conpotsplitter_,conpotsplitter_,27*(numscal_+1),false,true));
    blocksysmat->SetNumdim(numscal_);
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
  // (get noderowmap of discretization for creating this multivector9
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  convel_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // (negative) fluid momentum residual (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector9
  fluidres_ = rcp(new Epetra_MultiVector(*noderowmap,3,true));

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
  // set initial density to 1.0:
  // - used throughout simulation for non-temperature case
  // - used as good initial guess for stationary temperature case
  // -------------------------------------------------------------------
  densnp_->PutScalar(1.0);

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
  const double& time,
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
  const double& time,
  const Teuchos::RCP<Epetra_Vector>& phinp,
  Teuchos::RCP<Epetra_Vector>& neumann_loads
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
 | export velocity to column map and add it to parameter list  gjb 03/09|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddVelocityToParameterList(Teuchos::ParameterList& p)
{
  //provide velocity field for usage on element level
  // -> export to column map is necessary for parallel evaluation
  //SetState cannot be used since this Multivector is nodebased and not dofbased!
  const Epetra_Map* nodecolmap = discret_->NodeColMap();
  RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
  LINALG::Export(*convel_,*tmp);
  p.set("velocity field",tmp);

  // for low-Mach-number flow, also provide (negative) fluid momentum residual for
  // obtaining subgrid-scale velocity field for usage on element level
  if (prbtype_ == "loma")
  {
    RefCountPtr<Epetra_MultiVector> tmp2 = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*fluidres_,*tmp2);
    p.set("fluid momentum residual",tmp2);
  }

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

  while (stopnonliniter==false)
  {
    itnum++;

    double tcpu;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix
    // -------------------------------------------------------------------
    {
      // get cpu time
      tcpu=ds_cputime();

      // zero out matrix entries
      sysmat_->Zero();

      // reset the residual vector 
      residual_->PutScalar(0.0);

      {
      // time measurement: element calls
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

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
      eleparams.set("frt",frt_);// ELCH specific factor F/RT

      //provide velocity field (export to column map necessary for parallel evaluation)
      AddVelocityToParameterList(eleparams);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("hist" ,hist_);

      // add element parameters and density state according to time-int. scheme
      AddSpecificTimeIntegrationParameters(eleparams);

      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,residual_);
        discret_->ClearState();
      }

      // finalize the complete matrix
      sysmat_->Complete();

      // end time measurement for element
      dtele_=ds_cputime()-tcpu;

      } // time measurement for element
    }


    // scaling to get true residual vector for all time integration schemes
    // we store this vector BEFORE contributions due to further b.c. are added
    // to the residual. This way, at Dirichlet-, Neumann- and ElektrodeKinetics 
    // boundaries the boundary flux values can be computed form the trueresidual
    trueresidual_->Update(ResidualScaling(),*residual_,0.0);

    // add Neumann b.c. scaled with a factor due to time discretization
    AddNeumannToResidual();

    // add potential Neumann inflow
    if (neumanninflow_) ComputeNeumannInflow(sysmat_,residual_);

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
      tcpu=ds_cputime();

      // time measurement: call linear solver
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + call linear solver");

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        solver_->AdaptTolerance(ittol,actresidual,adaptolbetter);
      }

      // print (DEBUGGING!)
      //LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
      //(SystemMatrix()->EpetraMatrix())->Print(cout);

      solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,itnum==1);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpu;
    }

    //------------------------------------------------ update solution vector
/*    if (itnum == 1)
        phinp_->Update(0.25,*increment_,1.0);
    else  */
    phinp_->Update(1.0,*increment_,1.0);

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
    Teuchos::RCP<Epetra_Vector> onlycon = conpotsplitter_.ExtractOtherVector(residual_);
    onlycon->Norm2(&conresnorm);

    conpotsplitter_.ExtractOtherVector(increment_,onlycon);
    onlycon->Norm2(&incconnorm_L2);

    conpotsplitter_.ExtractOtherVector(phinp_,onlycon);
    onlycon->Norm2(&connorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypot = conpotsplitter_.ExtractCondVector(residual_);
    onlypot->Norm2(&potresnorm);

    conpotsplitter_.ExtractCondVector(increment_,onlypot);
    onlypot->Norm2(&incpotnorm_L2);

    conpotsplitter_.ExtractCondVector(phinp_,onlypot);
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
 | contains the solver                                          vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Solve()
{
  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  PrintTimeStepInfo();

  double tcpu;

  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------
  {
    // time measurement: element calls
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");
    // get cpu time
    tcpu=ds_cputime();

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

    //provide velocity field (export to column map necessary for parallel evaluation)
    AddVelocityToParameterList(eleparams);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("hist",hist_);
    if (turbmodel_) discret_->SetState("subgrid diffusivity",subgrdiff_);

    // AVM3 separation
    if (incremental_ and fssgd_ != "No") AVM3Separation();

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
    dtele_=ds_cputime()-tcpu;
  }

  // scaling to get true residual vector for all time integration schemes
  // in incremental case: boundary flux values can be computed from trueresidual
  if (incremental_) trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  AddNeumannToResidual();

  // add potential Neumann inflow
  if (neumanninflow_) ComputeNeumannInflow(sysmat_,residual_);

  // Apply Dirichlet boundary conditions to system matrix and solve system in
  // incremental or non-incremental case
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
      tcpu=ds_cputime();

      solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,true);

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpu;
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
      tcpu=ds_cputime();

      solver_->Solve(sysmat_->EpetraOperator(),phinp_,residual_,true,true);

      // end time measurement for solver
      dtsolve_=ds_cputime()-tcpu;
    }
  }

  return;
} // ScaTraTimIntImpl::Solve


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
 | set the actual velocity fields and subgrid diffusivities   gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(
RCP<const Epetra_Vector> extvel,
RCP<const Epetra_Vector> extsubgrvisc,
RCP<const Epetra_Vector> extresidual,
RCP<DRT::Discretization> fluiddis)
{
  if (cdvel_ != 2)
    dserror("Wrong SetVelocity() called for velocity field type %d!",cdvel_);

  if (prbtype_ == "loma")
  {
    // store temperature and velocity of previous iteration for convergence check
    tempincnp_->Update(1.0,*phinp_,0.0);
    //velincnp_->Update(1.0,*convel_,0.0);
  }

  // check vector compatibility and determine space dimension
  int numdim =-1;
  if (extvel->MyLength()<= (4* convel_->MyLength()) and
      extvel->MyLength() > (3* convel_->MyLength()))
    numdim = 3;
  else if (extvel->MyLength()<= (3* convel_->MyLength()))
    numdim = 2;
  else
    dserror("fluid velocity vector too large");

  // get noderowmap of scatra discretization
  const Epetra_Map* scatranoderowmap = discret_->NodeRowMap();

  // get dofrowmap of fluid discretization
  const Epetra_Map* fluiddofrowmap = fluiddis->DofRowMap();

  vector<int>    Indices(1);
  vector<double> Values(1);

  // loop over local nodes of scatra discretization
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // first of all, assume the present node is not a slavenode
    bool slavenode=false;

    // get the processor-local scatra node
    DRT::Node*  scatralnode = discret_->lRowNode(lnodeid);

    // get the processor-local fluid node
    DRT::Node*  fluidlnode = fluiddis->lRowNode(lnodeid);

    // the set of degrees of freedom associated with the fluid node
    vector<int> fluidnodedofset = fluiddis->Dof(fluidlnode);

    // check whether we have a pbc condition on this scatra node
    vector<DRT::Condition*> mypbc;
    scatralnode->GetCondition("SurfacePeriodic",mypbc);

    // yes, we have a periodic boundary condition on this scatra node
    if (mypbc.size()>0)
    {
      // get master and list of all his slavenodes
      map<int, vector<int> >::iterator master = pbcmapmastertoslave_->find(scatralnode->Id());

      // check whether this is a slavenode
      if (master == pbcmapmastertoslave_->end())
      {
        // indeed a slavenode
        slavenode = true;
      }
      else
      {
        // we have a masternode: set values for all slavenodes
        vector<int>::iterator i;
        for(i=(master->second).begin();i!=(master->second).end();++i)
        {
          // global and processor-local scatra node ID for slavenode
          int globalslaveid = *i;
          int localslaveid  = scatranoderowmap->LID(globalslaveid);

          // get the processor-local fluid slavenode
          DRT::Node*  fluidlslavenode = fluiddis->lRowNode(localslaveid);

          // the set of degrees of freedom associated with the node
          vector<int> slavenodedofset = fluiddis->Dof(fluidlslavenode);

          // global and processor-local fluid dof ID
          int fgid = slavenodedofset[0];
          int flid = fluiddofrowmap->LID(fgid);

          // get velocity for this processor-local fluid dof
          double velocity = (*extvel)[flid];
          // insert velocity value in vector
          convel_->ReplaceMyValue(lnodeid, 0, velocity);

          // get (negative) fluid residual value for this processor-local fluid dof
          double residual = (*extresidual)[flid];
          // insert fluid residual value in vector
          fluidres_->ReplaceMyValue(lnodeid, 0, residual);

          // get subgrid viscosity for this processor-local fluid dof
          // and divide by turbulent Prandtl number to get diffusivity
          Indices[0] = localslaveid;
          Values[0]  = (*extsubgrvisc)[flid]/tpn_;
          subgrdiff_->ReplaceMyValues(1,&Values[0],&Indices[0]);

          for(int index=1;index<numdim;++index)
          {
            // global and processor-local fluid dof ID
            fgid = slavenodedofset[index];
            flid = fluiddofrowmap->LID(fgid);

            // get velocity for this processor-local fluid dof
            double velocity =(*extvel)[flid];
            // insert velocity value in vector
            convel_->ReplaceMyValue(localslaveid, index, velocity);

            // get (negative) fluid residual value for this processor-local fluid dof
            double residual = (*extresidual)[flid];
            // insert fluid residual value in vector
            fluidres_->ReplaceMyValue(lnodeid, index, residual);
          }
        }
      }
    }

    // do this for all nodes other than slavenodes
    if (slavenode == false)
    {
      // global and processor-local fluid dof ID
      int fgid = fluidnodedofset[0];
      int flid = fluiddofrowmap->LID(fgid);

      // get velocity for this processor-local fluid dof
      double velocity = (*extvel)[flid];
      // insert velocity value in vector
      convel_->ReplaceMyValue(lnodeid, 0, velocity);

      // get (negative) fluid residual value for this processor-local fluid dof
      double residual = (*extresidual)[flid];
      // insert fluid residual value in vector
      fluidres_->ReplaceMyValue(lnodeid, 0, residual);

      // get subgrid viscosity for this processor-local fluid dof
      // and divide by turbulent Prandtl number to get diffusivity
      Indices[0] = lnodeid;
      Values[0]  = (*extsubgrvisc)[flid]/tpn_;
      subgrdiff_->ReplaceMyValues(1,&Values[0],&Indices[0]);

      for(int index=1;index<numdim;++index)
      {
        // global and processor-local fluid dof ID
        fgid = fluidnodedofset[index];
        flid = fluiddofrowmap->LID(fgid);

        // get velocity for this processor-local fluid dof
        double velocity = (*extvel)[flid];
        // insert velocity value in vector
        convel_->ReplaceMyValue(lnodeid, index, velocity);

        // get (negative) fluid residual value for this processor-local fluid dof
        double residual = (*extresidual)[flid];
        // insert fluid residual value in vector
        fluidres_->ReplaceMyValue(lnodeid, index, residual);
      }
    }
  }

  // initial velocity field has now been set
  if (step_ == 0) initialvelset_ = true;

  return;

} // ScaTraTimIntImpl::SetVelocityField


/*----------------------------------------------------------------------*
 |  set initial field for phi                                 gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetInitialField(int init, int startfuncno)
{
  if (init == 0) // zero_field
  {
    phin_-> PutScalar(0);
    phinp_-> PutScalar(0);
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
      int err =0;

      // random noise is relative to difference of max-min values of initial profile
      double perc = params_->sublist("TURBULENCE PARAMETERS").get<double>("CHAN_AMPL_INIT_DIST",0.1);

      // out to screen
      if (myrank_==0)
      {
        cout << "Disturbed initial scalar profile:   max. " << perc*100 << "% random perturbation\n";
        cout << "\n\n";
      }

      double thisphi=0;
      double mymaxphi=0;
      double myminphi=10000000.0;
      double maxphi=0;
      double minphi=10000000.0;
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
          const int dofgid = nodedofset[0];
          int doflid = dofrowmap->LID(dofgid);

          thisphi=(*phinp_)[doflid];
          if (mymaxphi*mymaxphi < thisphi*thisphi) mymaxphi=thisphi;
          if (myminphi*myminphi > thisphi*thisphi) myminphi=thisphi;
        }
      }

      // get overall max and min values and range between min and max
      discret_->Comm().MaxAll(&mymaxphi,&maxphi,1);
      discret_->Comm().MaxAll(&myminphi,&minphi,1);
      double range = abs(maxphi - minphi);

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = discret_->Dof(lnode);

        // check whether we have a pbc condition on this node
        vector<DRT::Condition*> mypbc;

        lnode->GetCondition("SurfacePeriodic",mypbc);

        // check whether a periodic boundary condition is active on this node
        if (mypbc.size()>0)
        {
          // yes, we have one

          // get the list of all his slavenodes
          map<int, vector<int> >::iterator master = pbcmapmastertoslave_->find(lnode->Id());

          // slavenodes are ignored
          if(master == pbcmapmastertoslave_->end()) continue;
        }

        int numdofs = nodedofset.size();
        for (int k=0;k< numdofs;++k)
        {
          int dofgid = nodedofset[k];

          double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);

          double noise = perc * range * randomnumber;

          err += phinp_->SumIntoGlobalValues(1,&noise,&dofgid);
          err += phin_ ->SumIntoGlobalValues(1,&noise,&dofgid);
        }

        if(err!=0) dserror("dof not on proc");
      }
    }
  }
  else if (init==2) // field_by_condition
  {
    dserror("Initialfield by condition not finished yet;");
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
            // get initial value from condition
            double phi0 = 2.0;
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
  else
    dserror("unknown option for initial field: %d", init);

  return;
} // ScaTraTimIntImpl::SetInitialField


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                   gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl()
{
  return;
}


#endif /* CCADISCRET       */
