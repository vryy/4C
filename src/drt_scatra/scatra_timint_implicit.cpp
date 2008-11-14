/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_implicit.cpp
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme

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
#include "../drt_fluid/vm3_solver.H"
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
  time_(0.0),
  step_(0),
  prbtype_  (params_->get<string>("problem type")),
  stepmax_  (params_->get<int>("max number timesteps")),
  maxtime_  (params_->get<double>("total time")),
  timealgo_ (params_->get<INPUTPARAMS::ScaTraTimeIntegrationScheme>("time int algo")),
  upres_    (params_->get<int>("write solution every")),
  uprestart_(params_->get<int>("write restart every")),
  writeflux_(params_->get<string>("write flux")),
  dta_      (params_->get<double>("time step size")),
  dtp_      (params_->get<double>("time step size")),
  theta_    (params_->get<double>("theta")),
  cdvel_    (params_->get<int>("velocity field")),
  fssgd_    (params_->get<string>("fs subgrid diffusivity")),
  frt_      (96485.3399/(8.314472 * params_->get<double>("TEMPERATURE",298.0)))
{
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
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d mesh we have 27 adjacent
  // nodes with numdof DOF each.
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate
  const int numdof = discret_->NumDof(discret_->lRowNode(0));
  const int numscal = numdof - 1;
  if (prbtype_ == "elch")
  {
    // set up the condif-el.potential splitter
    FLD::UTILS::SetupFluidSplit(*discret_,numscal,conpotsplitter_);
    if (myrank_==0)
    {
      cout<<"\nSetup of conpotsplitter: numscal = "<<numscal<<endl;
      cout<<"Temperature value T (Kelvin)     = "<<params_->get<double>("TEMPERATURE",298.0)<<endl;
      cout<<"Constant F/RT                    = "<<frt_<<endl;
    }
  }

  if (params_->get<int>("BLOCKPRECOND",0) )
  {
    // we need a block sparse matrix here
    if (prbtype_ != "elch") dserror("Block-Preconditioning is only for ELCH problems");
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(conpotsplitter_,conpotsplitter_,27,false,true));
    blocksysmat->SetNumdim(numscal);
    sysmat_ = blocksysmat;
  }
  else
  {
    // we need the 'standard' sparse matrix
    if (fssgd_ == "No")
    {
      // initialize standard (stabilized) system matrix (and save its graph!)
      sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
    }
    else
    {
      // do not save the graph for this application
      sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));
    }
  }

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // solutions at time n+1 and n
  phinp_        = LINALG::CreateVector(*dofrowmap,true);
  phin_         = LINALG::CreateVector(*dofrowmap,true);

  // density at time n+1
  densnp_        = LINALG::CreateVector(*dofrowmap,true);

  if (prbtype_ == "loma")
  {
    // density at times n and n-1
    densn_        = LINALG::CreateVector(*dofrowmap,true);
    densnm_       = LINALG::CreateVector(*dofrowmap,true);

    // density increment at time n+1
    densincnp_    = LINALG::CreateVector(*dofrowmap,true);
  }

  // histvector --- a linear combination of phinm, phin (BDF)
  //                or phin, phidtn (One-Step-Theta)
  hist_         = LINALG::CreateVector(*dofrowmap,true);

  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  /// convective velocity (always three velocity components per node)
  convel_         = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

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

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // the residual vector --- more or less the rhs
  residual_     = LINALG::CreateVector(*dofrowmap,true);

  // incremental solution vector
  increment_     = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // initialize subgrid-diffusivity matrix + respective ouptput
  // -------------------------------------------------------------------
  if (fssgd_ != "No")
  {
    sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

    if (myrank_ == 0)
    {
      // Output
      cout << "Fine-scale subgrid-diffusivity approach based on AVM3: ";
      cout << &endl << &endl;
      cout << params_->get<string>("fs subgrid diffusivity");
      cout << &endl << &endl;
    }
  }

  // set initial field
  SetInitialField(params_->get<int>("scalar initial field"), params_->get<int>("scalar initial field func number"));

  // set initial density to 1.0:
  // used throughout simulation for non-temperature case
  // used as good initial guess for stationary temperature case
  densnp_->PutScalar(1.0);

  return;

} // ScaTraTimIntImpl::ScaTraTimIntImpl


/*----------------------------------------------------------------------*
| returns matching string for each time integration scheme   gjb 08/08 |
*----------------------------------------------------------------------*/
std::string SCATRA::ScaTraTimIntImpl::MapTimIntEnumToString
(
   const enum INPUTPARAMS::ScaTraTimeIntegrationScheme term
)
{
  // length of return string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPUTPARAMS::timeint_one_step_theta :
    return "One-Step-Theta";
    break;
  case INPUTPARAMS::timeint_bdf2 :
    return "    BDF2      ";
    break;
  case INPUTPARAMS::timeint_stationary :
    return "  Stationary  ";
    break;
  case INPUTPARAMS::timeint_gen_alpha :
    return "  Gen. Alpha  ";
    break;
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }
}


/*----------------------------------------------------------------------*
 | Start the time integration. Allows                                   |
 |                                                                      |
 |  o starting steps with different algorithms                          |
 |  o the "standard" time integration                                   |
 |                                                              vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Integrate()
{
  // solve the problem
  TimeLoop();

  // just beauty
  cout<<endl<<endl;

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // ScaTraTimIntImpl::Integrate


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
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    Solve();

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
  if (step_==0) PrepareFirstTimeStep();

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

  // apply DBCs
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

  // set needed parameters in parameter list
  ParameterList p;
  p.set("total time", time); // actual time t_{n+1}

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp);
  discret_->SetState("densnp",densnp_);
  // evaluate Neumann conditions at actual time t_{n+1}
  discret_->EvaluateNeumann(p,*neumann_loads);
  discret_->ClearState();

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

      // reset the residual vector and add actual Neumann loads
      // scaled with a factor resulting from time discretization
      residual_->Update(theta_*dta_,*neumann_loads_,0.0);

      if (prbtype_=="elch")
      { // evaluate electrode kinetics conditions

        // time measurement: evaluate condition 'ElectrodeKinetics'
        TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ElectrodeKinetics'");

        // create an parameter list
        ParameterList condparams;

        // action for elements
        condparams.set("action","calc_elch_electrode_kinetics");
        condparams.set("frt",frt_); // factor F/RT
        condparams.set("total time",time_);
        condparams.set("thsl",theta_*dta_);
        if (MethodName()==INPUTPARAMS::timeint_stationary)
          condparams.set("using stationary formulation",true);
        else
          condparams.set("using stationary formulation",false);

        // set vector values needed by elements
        discret_->ClearState();
        discret_->SetState("phinp",phinp_);

        std::string condstring("ElectrodeKinetics");
        discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
        discret_->ClearState();
      }

      {
      // time measurement: element calls
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

      // create the parameters for the discretization
      ParameterList eleparams;

      // action for elements
      eleparams.set("action","calc_condif_systemmat_and_residual");

      // other parameters that might be needed by the elements
      eleparams.set("total time",time_);
      eleparams.set("time-step length",dta_);
      eleparams.set("thsl",theta_*dta_);
      eleparams.set("problem type",prbtype_);
      eleparams.set("fs subgrid diffusivity",fssgd_);
      eleparams.set("frt",frt_);// ELCH specific factor F/RT
      if (MethodName()==INPUTPARAMS::timeint_stationary)
        eleparams.set("using stationary formulation",true);
      else
        eleparams.set("using stationary formulation",false);

      //provide velocity field (export to column map necessary for parallel evaluation)
      //SetState cannot be used since this Multivector is nodebased and not dofbased
      const Epetra_Map* nodecolmap = discret_->NodeColMap();
      RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
      LINALG::Export(*convel_,*tmp);
      eleparams.set("velocity field",tmp);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp",phinp_);
      discret_->SetState("densnp",densnp_);
      discret_->SetState("hist" ,hist_);

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

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    stopnonliniter = AbortNonlinIter(itnum,itemax,ittol,actresidual);
    if (stopnonliniter == true) break;

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual values are supposed to be zero at Dirichlet boundaries
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
  bool stopnonliniter = false;

  //----------------------------------------------------- compute norms
  double incconnorm_L2(0.0);
  double incpotnorm_L2(0.0);

  double connorm_L2(0.0);
  double potnorm_L2(0.0);

  double conresnorm(0.0);
  double potresnorm(0.0);

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
      - convergence check is not required (we solve at least once!)    */
  if (itnum == 1)
  {
    if (myrank_ == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
          itnum,itemax,ittol,conresnorm,potresnorm);
      printf(" (      --     ,te=%10.3E",dtele_);
      printf(")\n");
    }
  }
  /* ordinary case later iteration steps:
      - solution increment can be printed
      - convergence check should be done*/
  else
  {
    // this is the convergence check
    // We always require at least one solve. Otherwise the
    // perturbation at the FSI interface might get by unnoticed.
    if (conresnorm <= ittol and potresnorm <= ittol and
        incconnorm_L2/connorm_L2 <= ittol and incpotnorm_L2/potnorm_L2 <= ittol)
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
            itnum,itemax,ittol,conresnorm,potresnorm,
            incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
        printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
        printf(")\n");
        printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

        FILE* errfile = params_->get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
              itnum,itemax,ittol,conresnorm,potresnorm,
              incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
        }
      }
      //break;
    }
    else // if not yet converged
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
            itnum,itemax,ittol,conresnorm,potresnorm,
            incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
        printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
        printf(")\n");
      }
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if ((itnum == itemax))
  {
    stopnonliniter=true;
    if (myrank_ == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");

      FILE* errfile = params_->get<FILE*>("err file",NULL);
      if (errfile!=NULL)
      {
        fprintf(errfile,"elch unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
            itnum,itemax,ittol,conresnorm,potresnorm,
            incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
      }
    }
    //break;
  }

  // return the maximum residual value -> used for adaptivity of linear solver tolarance
  actresidual = max(conresnorm,potresnorm);
  actresidual = max(actresidual,incconnorm_L2/connorm_L2);
  actresidual = max(actresidual,incpotnorm_L2/potnorm_L2);

  return stopnonliniter;
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

    // reset the residual vector and add actual Neumann loads
    // scaled with a factor resulting from time discretization
    residual_->Update(theta_*dta_,*neumann_loads_,0.0);

    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_condif_systemmat_and_residual");

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("time-step length",dta_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("problem type",prbtype_);
    eleparams.set("fs subgrid diffusivity",fssgd_);
    if (MethodName()==INPUTPARAMS::timeint_stationary)
      eleparams.set("using stationary formulation",true);
    else
      eleparams.set("using stationary formulation",false);

    //provide velocity field (export to column map necessary for parallel evaluation)
    //SetState cannot be used since this Multivector is nodebased and not dofbased
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*convel_,*tmp);
    eleparams.set("velocity field",tmp);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    discret_->SetState("densnp",densnp_);
    discret_->SetState("hist" ,hist_ );

    // decide whether AVM3-based solution approach or standard approach
    if (fssgd_ != "No") AVM3Scaling(eleparams);
    else
    {
      // call standard loop over elements
      discret_->Evaluate(eleparams,sysmat_,residual_);
      discret_->ClearState();
    }

    // finalize the complete matrix
    sysmat_->Complete();

    // end time measurement for element
    dtele_=ds_cputime()-tcpu;
  }

  // Apply dirichlet boundary conditions to system matrix
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

  return;
} // ScaTraTimIntImpl::Solve


/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                          gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Output()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  if (step_%upres_==0)  //yes, output is desired
  {
    // write state vectors
    OutputState();

    // write domain decomposition for visualization (only once!)
    if (step_==upres_)
      output_->WriteElementData();

    //add restart data
    if (step_%uprestart_==0)
      OutputRestart();

    // write flux vector field
    if (writeflux_!="No")
      OutputFlux();
  }

  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((step_%uprestart_== 0) && (step_%upres_!=0))
  {
    OutputState();   // write state vectors
    OutputRestart(); // add restart data
    if (writeflux_!="No") // write flux vector field
      OutputFlux();
  }

  return;
} // ScaTraTimIntImpl::Output


/*----------------------------------------------------------------------*
 |  write current state to BINIO                             gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputState()
{
  output_->NewStep    (step_,time_);
  output_->WriteVector("phinp", phinp_);
  output_->WriteVector("convec_velocity", convel_,IO::DiscretizationWriter::nodevector);

  return;
}


/*----------------------------------------------------------------------*
 |  prepare AVM3-based scale separation                        vg 10/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AVM3Preparation()
{
  {// time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // create normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Zero();

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements, time factor and stationary flag
  eleparams.set("action","calc_subgrid_diffusivity_matrix");
  eleparams.set("thsl",theta_*dta_);
  if (MethodName()==INPUTPARAMS::timeint_stationary)
    eleparams.set("using stationary formulation",true);
  else
    eleparams.set("using stationary formulation",false);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_sd_,residual_);
  discret_->ClearState();

  // finalize the normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Complete();

  // apply DBC to normalized all-scale subgrid-diffusivity matrix
  LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,*(dbcmaps_->CondMap()));

  // extract the ML parameters
  ParameterList&  mllist = solver_->Params().sublist("ML Parameters");

  // call the VM3 constructor
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    vm3_solver_ = rcp(new FLD::VM3_Solver(sysmat_sd_,dirichtoggle,mllist,true,false) );
  }
  }// time measurement: avm3

  return;
}


/*----------------------------------------------------------------------*
 |  scaling of AVM3-based subgrid-diffusivity matrix           vg 10/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AVM3Scaling(ParameterList& eleparams)
{

  {// time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // subgrid-diffusivity-scaling vector
  subgrdiff_ = LINALG::CreateVector(*dofrowmap,true);
  }// time measurement: avm3

  // call loop over elements (one matrix + subgr.-visc.-scal. vector)
  discret_->Evaluate(eleparams,sysmat_,null,residual_,subgrdiff_);
  discret_->ClearState();

  {// time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // check whether VM3 solver exists
  if (vm3_solver_ == null) dserror("vm3_solver not allocated");

  // call the VM3 scaling:
  // scale precomputed matrix product by subgrid-viscosity-scaling vector
  Teuchos::RCP<LINALG::SparseMatrix> sysmat = SystemMatrix();
  vm3_solver_->Scale(sysmat_sd_,sysmat,zeros_,zeros_,subgrdiff_,zeros_,false );
  }// time measurement: avm3

  return;
}


/*----------------------------------------------------------------------*
 | update the velocity field                                  gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(int veltype, int velfuncno)
{
    if (veltype != cdvel_)
        dserror("velocity field type does not match: got %d, but expected %d!",veltype,cdvel_);

    if (veltype == 0) // zero
        convel_->PutScalar(0); // just to be sure!
    else if (veltype == 1)  // function
    {
    int numdim =3;
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value=DRT::UTILS::FunctionManager::Instance().Funct(velfuncno-1).Evaluate(index,lnode->X());
        convel_->ReplaceMyValue (lnodeid, index, value);
      }
    }
    }
    else
        dserror("error in setVelocityField");

    return;

} // ScaTraImplicitTimeInt::SetVelocityField


/*----------------------------------------------------------------------*
 | update the velocity field                                  gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel)
{
  if (veltype != cdvel_)
    dserror("velocity field type does not match: got %d, but expected %d!",veltype,cdvel_);

  // check vector compatibility and determine space dimension
  int numdim =-1;
  if (extvel->MyLength()== (2* convel_->MyLength()))
    numdim = 2;
  else if (extvel->MyLength()== (3* convel_->MyLength()))
    numdim = 3;
  else
    dserror("velocity vectors do not match in size");

  if ((numdim == 3) or (numdim == 2))
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      for(int index=0;index<numdim;++index)
      {
        double value = (*extvel)[lnodeid*numdim + index];
        //printf("myvelocityvalue[%d][%d] = %3.16lf\n",lnodeid,index,value);
        convel_->ReplaceMyValue(lnodeid, index, value);
      }
    }
  }

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
  else if (init == 1)  // field_by_function
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
        double initialval=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(k,lnode->X());
        phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        phinp_->ReplaceMyValues(1,&initialval,&doflid);
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
 | set velocity field for low-Mach-number flow                 vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetLomaVelocity(RCP<const Epetra_Vector> extvel,
    RCP<DRT::Discretization> fluiddis)
{
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
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // get dofrowmap of fluid discretization
  const Epetra_Map* dofrowmap = fluiddis->DofRowMap();

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
    vector<int> nodedofset = fluiddis->Dof(fluidlnode);

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
          int localslaveid  = noderowmap->LID(globalslaveid);

          // get the processor-local fluid slavenode
          DRT::Node*  fluidlslavenode = fluiddis->lRowNode(localslaveid);

          // the set of degrees of freedom associated with the node
          vector<int> slavenodedofset = fluiddis->Dof(fluidlslavenode);

          for(int index=0;index<numdim;++index)
          {
            // global and processor-local fluid dof ID
            int gid = slavenodedofset[index];
            int lid = dofrowmap->LID(gid);

            // get density for this processor-local scatra node
            double dens  = (*densnp_)[localslaveid];
            // get velocity for this processor-local fluid dof
            double velocity =(*extvel)[lid];
            // insert velocity*density-value in vector
            convel_->ReplaceMyValue(localslaveid, index, velocity*dens);
          }
        }
      }
    }

    // do this for all nodes other than slavenodes
    if (slavenode == false)
    {
      for(int index=0;index<numdim;++index)
      {
        // global and processor-local fluid dof ID
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        // get density for this processor-local scatra node
        double dens  = (*densnp_)[lnodeid];
        // get velocity for this processor-local fluid dof
        double velocity = (*extvel)[lid];
        // insert velocity*density-value in vector
        convel_->ReplaceMyValue(lnodeid, index, velocity*dens);
      }
    }
  }

  return;
} // ScaTraTimIntImpl::SetLomaVelocity


/*----------------------------------------------------------------------*
 | predict density for next time step for low-Mach-number flow vg 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PredictDensity()
{
  // predictor: linear extrapolation in time
  //densnp_->Update(2.0,*densn_,-1.0,*densnm_,0.0);

  // predictor: constant extrapolation in time
  densnp_->Update(1.0,*densn_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | compute density for low-Mach-number flow                    vg 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeDensity()
{
  // store density of previous iteration for convergence check
  densincnp_->Update(1.0,*densnp_,0.0);

  // set thermodynamic pressure p_therm (in N/m^2 = kg/(m*s^2) = J/m^3): 98100.0
  // constantly set to atmospheric pressure, for the time being -> dp_therm/dt=0
  // set specific gas constant R (in J/(kg*K)): 287.05
  // compute density based on equation of state: rho = (p_therm/R)*(1/T)
  double fac = 98100.0/287.05;
  densnp_->Reciprocal(*phinp_);
  densnp_->Scale(fac);

  //densnp_->PutScalar(1.0);

  return;
}


/*----------------------------------------------------------------------*
 | compute density for low-Mach-number flow                    vg 08/08 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::DensityConvergenceCheck(int          itnum,
                                                       int          itmax,
                                                       const double ittol)
{
  bool stopnonliniter = false;

  double densincnorm_L2;
  double densnorm_L2;

  densincnp_->Update(1.0,*densnp_,-1.0);
  densincnp_->Norm2(&densincnorm_L2);
  densnp_->Norm2(&densnorm_L2);

  // care for the case that there is (almost) zero density
  if (densnorm_L2 < 1e-6) densnorm_L2 = 1.0;

  if (myrank_==0)
  {
    cout<<"\n**********************\n OUTER ITERATION STEP \n**********************\n";
    printf("+------------+-------------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- dens-norm -|-- dens-inc --|\n");
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
             itnum,itmax,ittol,densnorm_L2,densincnorm_L2/densnorm_L2);
    printf("\n");
    printf("+------------+-------------------+--------------+--------------+\n");
  }

  if (densincnorm_L2/densnorm_L2 <= ittol) stopnonliniter=true;

  // warn if itemax is reached without convergence, but proceed to next timestep
  if ((itnum == itmax) and (densincnorm_L2/densnorm_L2 > ittol))
  {
    stopnonliniter=true;
    if (myrank_==0)
    {
      printf("|    >>>>>> not converged in itemax steps!                     |\n");
      printf("+--------------------------------------------------------------+\n");
    }
  }

  return stopnonliniter;
}


/*----------------------------------------------------------------------*
 | update density at n-1 and n for low-Mach-number flow        vg 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::UpdateDensity()
{
  densnm_->Update(1.0,*densn_ ,0.0);
  densn_->Update(1.0,*densnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputMeanTempAndDens()
{
  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  discret_->SetState("densnp",densnp_);
  // set action for elements
  ParameterList eleparams;
  eleparams.set("action","calc_temp_and_dens");

  // variables for integrals of temperature, density and domain
  double tempint = 0.0;
  double densint = 0.0;
  double domint  = 0.0;
  eleparams.set("temperature integral",tempint);
  eleparams.set("density integral",    densint);
  eleparams.set("domain integral",     domint);

  // evaluate integrals of temperature, density and domain
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // get integral values on this proc
  tempint = eleparams.get<double>("temperature integral");
  densint = eleparams.get<double>("density integral");
  domint  = eleparams.get<double>("domain integral");

  // get integral values in parallel case
  double partempint = 0.0;
  double pardensint = 0.0;
  double pardomint  = 0.0;
  discret_->Comm().SumAll(&tempint,&partempint,1);
  discret_->Comm().SumAll(&densint,&pardensint,1);
  discret_->Comm().SumAll(&domint,&pardomint,1);

  // print out mean values of temperature and density
  if (myrank_ == 0)
  {
    cout << "Mean temperature: " << tempint/domint << endl;
    cout << "Mean density: " << densint/domint << endl;
  }

  // clean up
  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputFlux()
{
  RCP<Epetra_MultiVector> flux = CalcFlux();
  int numscal = flux->GlobalLength()/discret_->NumGlobalNodes();
  if (prbtype_=="elch") numscal -= 1; // ELCH case

  // post_drt_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filter

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk = rcp(new Epetra_MultiVector(*noderowmap,3,true));
  for(int k=1;k<=numscal;++k)
  {
    ostringstream temp;
    temp << k;
    string name = "flux_phi_"+temp.str();
    for (int i = 0;i<fluxk->MyLength();++i)
    {
      DRT::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(actnode,k-1);
      fluxk->ReplaceMyValue(i,0,((*flux)[0])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,1,((*flux)[1])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,2,((*flux)[2])[(flux->Map()).LID(dofgid)]);
    }
    if (numscal==1)
      output_->WriteVector("flux", fluxk, IO::DiscretizationWriter::nodevector);
    else
      output_->WriteVector(name, fluxk, IO::DiscretizationWriter::nodevector);
  }
  // that's it
  return;
}


/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector                        gjb   04/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFlux()
{
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // empty vector for (normal) mass or heat flux vectors (3D)
  Teuchos::RCP<Epetra_MultiVector> flux = rcp(new Epetra_MultiVector(*dofrowmap,3,true));

  // we have only 1 dof per node, so we have to treat each spatial direction separately
  Teuchos::RCP<Epetra_Vector> fluxx = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxy = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxz = LINALG::CreateVector(*dofrowmap,true);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  discret_->SetState("densnp",densnp_);
  // set action for elements
  ParameterList eleparams;
  eleparams.set("action","calc_condif_flux");
  eleparams.set("problem type",prbtype_);
  eleparams.set("frt",frt_);

  //provide velocity field (export to column map necessary for parallel evaluation)
  const Epetra_Map* nodecolmap = discret_->NodeColMap();
  RefCountPtr<Epetra_MultiVector> vel = rcp(new Epetra_MultiVector(*nodecolmap,3));
  LINALG::Export(*convel_,*vel);
  eleparams.set("velocity field",vel);

  // set control parameters
  string fluxcomputation("nowhere"); // domain / boundary
  string fluxtype("noflux"); // noflux / totalflux / diffusiveflux
  if (writeflux_!="No")
  {
    size_t pos = writeflux_.find("_");    // find position of "_" in str
    fluxcomputation = writeflux_.substr(pos+1);   // get from "_" to the end
    fluxtype = writeflux_.substr(0,pos); // get from beginning to "_"
  }
  eleparams.set("fluxtype",fluxtype);

  // now compute the fluxes
  if (fluxcomputation=="domain")
  {
    // evaluate fluxes in the whole computational domain (e.g., for visualization of particle path-lines)
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz);
  }
  if (fluxcomputation=="boundary")
  {
    // evaluate fluxes on surface condition: only normal fluxes
    // normal flux value is stored at fluxx, fluxy and fluxz are set to 0.0
    string condstring("ElectrodeKinetics");

    // calculate integral of normal fluxes over indicated boundary
    // may be sum over several boundary parts given in input file
    double normfluxintegral = 0.0;
    eleparams.set("normfluxintegral",normfluxintegral);

    // do calculation on boundaries
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz,condstring);

    // get integral of normal flux on this proc
    normfluxintegral = eleparams.get<double>("normfluxintegral");

    // get integral of normal flux in parallel case
    double parnormfluxintegral = 0.0;
    discret_->Comm().SumAll(&normfluxintegral,&parnormfluxintegral,1);

    // print out integral of normal fluxes over indicated boundary
    if (myrank_ == 0)
      cout << "Integral of normal flux at indicated boundary: " << parnormfluxintegral << endl;
  }

  // clean up
  discret_->ClearState();

  // insert values into final flux vector for visualization
  for (int i = 0;i<flux->MyLength();++i)
  {
    flux->ReplaceMyValue(i,0,(*fluxx)[i]);
    flux->ReplaceMyValue(i,1,(*fluxy)[i]);
    flux->ReplaceMyValue(i,2,(*fluxz)[i]);
  }

  return flux;
} // ScaTraImplicitTimeInt::CalcNormalFlux



/*----------------------------------------------------------------------*
 |  calculate error compared to analytical solution            gjb 10/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol()
{
  int calcerr = params_->get<int>("CALCERROR");

  switch (calcerr)
  {
  case 0: // do nothing (the usual case)
    break;
  case 1:
  {
    //------------------------------------------------- Kwok et Wu,1995
    //   Reference:
    //   Kwok, Yue-Kuen and Wu, Charles C. K.
    //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
    //   Numerical Methods for Partial Differential Equations
    //   1995, Vol 11, 389-397

    // create the parameters for the discretization
    ParameterList p;

    // parameters for the elements
    p.set("action","calc_elch_kwok_error");
    p.set("total time",time_);
    p.set("frt",frt_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(p, errors);
    discret_->ClearState();

    // for the L2 norm, we need the square root
    double conerr1 = sqrt((*errors)[0]);
    double conerr2 = sqrt((*errors)[1]);
    double poterr  = sqrt((*errors)[2]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for Kwok et Wu:\n concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
             conerr1,conerr2,poterr);
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  }
  return;
} // end EvaluateErrorComparedToAnalyticalSol


/*----------------------------------------------------------------------*
 | construct toggle vector for Dirichlet dofs                  gjb 11/08|
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> SCATRA::ScaTraTimIntImpl::DirichletToggle()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
  return dirichtoggle;
}


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                   gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl()
{
  return;
}


#endif /* CCADISCRET       */
