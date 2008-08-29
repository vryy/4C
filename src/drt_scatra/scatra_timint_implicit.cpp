/*!----------------------------------------------------------------------
\file scatra_timint_implicit.cpp
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_timint_implicit.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_fluid/vm3_solver.H"
#include "../drt_lib/drt_function.H"

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
  fssgd_    (params_->get<string>("fs subgrid diffusivity","No"))
{
  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(discret_);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

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
  // nodes with 1 dof each.
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  // initialize standard (stabilized) system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // solutions at time n+1 and n
  phinp_        = LINALG::CreateVector(*dofrowmap,true);
  phin_         = LINALG::CreateVector(*dofrowmap,true);
  
  // state vector for solution at time n-1
  if (timealgo_==INPUTPARAMS::timeint_bdf2)
  {phinm_      = LINALG::CreateVector(*dofrowmap,true);}


  // histvector --- a linear combination of phinm, phin (BDF)
  //                or phin, phidtn (One-Step-Theta)
  hist_         = LINALG::CreateVector(*dofrowmap,true);

  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  /// convective velocity (always three velocity components per node)
  convel_         = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // toggle vector indicating which dofs have Dirichlet BCs
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  invtoggle_    = LINALG::CreateVector(*dofrowmap,false);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // the residual vector --- more or less the rhs
  residual_     = LINALG::CreateVector(*dofrowmap,true);

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
  // bound for the number of startsteps
  int    numstasteps         =params_->get<int>   ("number of start steps");

  if (timealgo_==INPUTPARAMS::timeint_stationary) // stationary case
    SolveStationaryProblem();
  else  // instationary case
  {
    // start procedure
    if (step_<numstasteps)
    {
      if (numstasteps>stepmax_)
      {
        dserror("more startsteps than steps");
      }

      dserror("no starting steps supported");
    }

    // continue with the final time integration
    TimeLoop();
  }

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // ScaTraTimIntImpl::Integrate


/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::TimeLoop()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  while (step_<stepmax_ and time_<maxtime_)
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
  //              initialization
  // -------------------------------------------------------------------
  if (step_==0) 
  {
    //ApplyDirichletBC(time_, phin_);
    //PrepareFirstTimeStep();
  }

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  IncrementTimeAndStep();

  // for bdf2 theta is set by the timestepsizes, 2/3 for const. dt
  if (timealgo_==INPUTPARAMS::timeint_bdf2)
  {
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  }

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  PrintTimeStepInfo();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  ApplyDirichletBC(time_,phinp_);
  ApplyNeumannBC(time_,phinp_,neumann_loads_);

  return;

} // ScaTraTimIntImpl::PrepareTimeStep


/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}           gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletBC
(
  const double& time,
  Teuchos::RCP<Epetra_Vector>& phinp
)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:      + apply dirich cond.");

  // apply DBCs
  // needed parameters
  ParameterList p;
  p.set("total time",time);  // actual time t_{n+1}

  // predicted Dirichlet values
  // \c  phinp then also holds prescribed new dirichlet values
  //     dirichtoggle_ is 1 for dirichlet dofs, 0 elsewhere
  discret_->ClearState();
  discret_->EvaluateDirichlet(p,phinp,Teuchos::null,Teuchos::null,dirichtoggle_);
  discret_->ClearState();

  // compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

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
  // evaluate Neumann conditions at actual time t_{n+1}
  discret_->EvaluateNeumann(p,*neumann_loads);
  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 | contains the solver                                          vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Solve(
  bool is_stat //if true, stationary formulations are used in the element
  )
{
  double tcpu;
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

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
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("fs subgrid diffusivity",fssgd_);
    eleparams.set("using stationary formulation",is_stat);

    //provide velocity field (export to column map necessary for parallel evaluation)
    //SetState cannot be used since this Multivector is nodebased and not dofbased
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*convel_,*tmp);
    eleparams.set("velocity field",tmp);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    discret_->SetState("hist" ,hist_ );

#if 0
    // AVMS solver: VM3 solution approach with extended matrix system
    // (may replace direct algebraic VM3 solution approach)
    // begin first encapsulation of AVMS solution approach
    if (fssgd_ != "No")
    {
      // create all-scale subgrid-diffusivity matrix
      sysmat_sd_->Zero();

      // call loop over elements (two matrices)
      discret_->Evaluate(eleparams,sysmat_,sysmat_sd_,residual_);
      discret_->ClearState();
    }
    // end first encapsulation of AVMS solution approach
#endif
    // direct algebraic VM3 solution approach
    // begin encapsulation of direct algebraic VM3 solution approach
    if (fssgd_ != "No")
    {
      {// time measurement: avm3
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

      // subgrid-viscosity-scaling vector
      sugrvisc_ = LINALG::CreateVector(*dofrowmap,true);

      }// time measurement: avm3

      if (step_ == 1)
      {
        {// time measurement: avm3
        TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

        // create normalized all-scale subgrid-diffusivity matrix
        sysmat_sd_->Zero();

        }// time measurement: avm3

        // call loop over elements (two matrices + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,sysmat_sd_,residual_,sugrvisc_);
        discret_->ClearState();

        {// time measurement: avm3
        TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

        // finalize the normalized all-scale subgrid-diffusivity matrix
        sysmat_sd_->Complete();

        // apply DBC to normalized all-scale subgrid-diffusivity matrix
        LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,dirichtoggle_);

        // extract the ML parameters
        ParameterList&  mllist = solver_->Params().sublist("ML Parameters");

        // call the VM3 constructor
        vm3_solver_ = rcp(new FLD::VM3_Solver(sysmat_sd_,dirichtoggle_,mllist,true,false) );

        }// time measurement: avm3
      }
      else
      {
        // call loop over elements (one matrix + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,null,residual_,sugrvisc_);
        discret_->ClearState();
      }

      // time measurement: avm3
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

      // check whether VM3 solver exists
      if (vm3_solver_ == null) dserror("vm3_solver not allocated");

      // call the VM3 scaling:
      // scale precomputed matrix product by subgrid-viscosity-scaling vector
      vm3_solver_->Scale(sysmat_sd_,sysmat_,zeros_,zeros_,sugrvisc_,zeros_,false );

    }
    // end encapsulation of direct algebraic VM3 solution approach
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

    LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,phinp_,dirichtoggle_);
  }

  //-------solve
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");
    // get cpu time
    tcpu=ds_cputime();

#if 0
    // AVMS solver: VM3 solution approach with extended matrix system
    // (may replace direct algebraic VM3 solution approach)
    // begin second encapsulation of AVMS solution approach
    if (fssgd_ != "No")
    {
      // add standard matrix to subgrid-diffusivity matrix: fine-scale matrix
      sysmat_sd_->Add(*sysmat_,false,1.0,1.0);

      // finalize fine-scale matrix
      sysmat_sd_->Complete();

      // apply DBC to fine-scale matrix
      LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,dirichtoggle_);

      // extract the ML parameters
      ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

      // call the AVMS constructor
      avms_solver_ = rcp(new AVMS_Solver(sysmat_sd_,sysmat_,dirichtoggle_,mllist) );

      // Apply the AVMS solver
      avms_solver_->Solve(*residual_,*phinp_,solver_.Params());
    }
    else
    // end second encapsulation of AVMS solution approach
#endif
      solver_->Solve(sysmat_->EpetraMatrix(),phinp_,residual_,true,true);

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
  //output_->WriteVector("residual", residual_);
  return;
}


/*----------------------------------------------------------------------*
 | solve stationary convection-diffusion problem               vg 05/07 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SolveStationaryProblem()
{
  // time measurement: time loop (stationary) 
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  if (myrank_==0) printf("Stationary Solver\n");
  step_=1;

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  ApplyDirichletBC(time_,phinp_);
  ApplyNeumannBC(time_,phinp_,neumann_loads_);

  // -------------------------------------------------------------------
  //                     solve nonlinear equation
  // -------------------------------------------------------------------
  Solve(true);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  return;
} // ScaTraTimIntImpl::SolveStationaryProblem


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
  { // just to be sure!
    phinp_->PutScalar(0);
    phin_-> PutScalar(0);
    if (timealgo_==INPUTPARAMS::timeint_bdf2)
    {
    phinm_->PutScalar(0);
    }
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

        phinp_->ReplaceMyValues(1,&initialval,&doflid);
        phin_->ReplaceMyValues(1,&initialval,&doflid);
        if (timealgo_==INPUTPARAMS::timeint_bdf2)
        {
        phinm_->ReplaceMyValues(1,&initialval,&doflid);
        }
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
            phinp_->ReplaceMyValues(1,&phi0,&doflid);
            phin_->ReplaceMyValues(1,&phi0,&doflid);
            if (timealgo_==INPUTPARAMS::timeint_bdf2)
            {
            phinm_->ReplaceMyValues(1,&phi0,&doflid);
            }
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
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputFlux()
{
  RCP<Epetra_MultiVector> flux = CalcFlux();
  int numdofpernode = flux->GlobalLength()/discret_->NumGlobalNodes();
  // post_drt_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filter

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk = rcp(new Epetra_MultiVector(*noderowmap,3,true));
  for(int k=0;k<numdofpernode;++k)
  {
    ostringstream temp;
    temp << k;
    string name = "flux_phi_"+temp.str();
    for (int i = 0;i<fluxk->MyLength();++i)
    {
      DRT::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(actnode,k);
      fluxk->ReplaceMyValue(i,0,((*flux)[0])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,1,((*flux)[1])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,2,((*flux)[2])[(flux->Map()).LID(dofgid)]);
    }
    if (numdofpernode==1)
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
  // set action for elements
  ParameterList eleparams;
  eleparams.set("action","calc_condif_flux");

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
  if (fluxcomputation=="condition")
  {
    // evaluate fluxes on surface condition only
    // if restriction to normal(!) fluxes is needed put it here
    string condstring("FluxCalculation");
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz,condstring);
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
 | Destructor dtor (public)                                   gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl()
{
  return;
}


#endif /* CCADISCRET       */
