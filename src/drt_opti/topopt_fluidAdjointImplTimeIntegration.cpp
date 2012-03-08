/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointImplTimeIntegration.cpp

\brief Control routine for fluid adjoint (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifdef CCADISCRET

// include Gmsh output
#define GMSHOUTPUT


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/time_integration_scheme.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_opti/topopt_optimizer.H"

#include "topopt_fluidAdjointImplTimeIntegration.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
TOPOPT::ADJOINT::ImplicitTimeInt::ImplicitTimeInt(RefCountPtr<DRT::Discretization> actdis,
                                                LINALG::Solver&       solver,
                                                ParameterList&        params,
                                                IO::DiscretizationWriter& output,
                                                bool alefluid)
  :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  step_(0),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1))
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint initialization");

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------

  // type of time-integration
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params_, "time int algo");

  // time-step size
  dt_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = params_.get<double>("total time");

  // set initial time = endtime so that it fits to the fluid parameter setting
  // potentially we do one step less here than in fluid
  // fluid criteria are:
  // o endtime <= numstep * dt
  // o endtime < maxtime + dt
  //
  // additionally evaluate the correct number of time steps since it is
  // required for evaluation of the fluid velocity at the correct step
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    time_ = dt_;
    stepmax_ = 1;
  }
  else
  {
    if (fabs(maxtime_-dt_*stepmax_)>1.0e-14)
      dserror("Fix total simulation time, time step size and number of time steps\n"
          "so that: sim_time = dt * num_timesteps");

    time_ = dt_*stepmax_; // stepmax_ remains
  }

  // parameter theta for time-integration schemes
  theta_      = params_.get<double>("theta");
  theta_pre_  = params_.get<double>("theta_pre");
  theta_div_  = params_.get<double>("theta_div");
  // compute or set 1.0 - theta for time-integration schemes
  omtheta_ = 1.0 - theta_;
  omtheta_pre_ = 1.0 - theta_pre_;
  omtheta_div_ = 1.0 - theta_div_;


  // -------------------------------------------------------------------
  // account for potential Neuman inflow terms if required
  // -------------------------------------------------------------------
  neumanninflow_ = false;
  if (params_.get<string>("Neumann inflow","no") == "yes") neumanninflow_ = true;

  // -------------------------------------------------------------------
  // care for periodic boundary conditions
  // -------------------------------------------------------------------
  pbcmapmastertoslave_ = params_.get<RCP<map<int,vector<int> > > >("periodic bc");
  discret_->ComputeNullSpaceIfNecessary(solver_.Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  numdim_ = params_.get<int>("number of velocity degrees of freedom");

  velpressplitter_ = rcp(new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*discret_,numdim_,*velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  if (not params_.get<int>("Simple Preconditioner",0) &&
      not params_.get<int>("AMG BS Preconditioner",0))
  {
    // initialize standard (stabilized) system matrix
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(*velpressplitter_,*velpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim_);
    sysmat_ = blocksysmat;
  }

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*dofrowmap,true);
  veln_  = LINALG::CreateVector(*dofrowmap,true);

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
  }

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------
  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_      = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize pseudo-porosity vector for topology optimization as null
  topopt_porosity_ = Teuchos::null;

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralAdjointParameter();
  SetElementTimeParameter();
} // ImplicitTimeInt::ImplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration. Allows                                   |
 |                                                                      |
 |  o starting steps with different algorithms                          |
 |  o the "standard" time integration                                   |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::Integrate()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

    cout << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                             " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";
    cout << "                             " << "Tau Type        = " << stabparams->get<string>("DEFINITION_TAU") <<"\n";
    cout << "                             " << "Evaluation Tau  = " << stabparams->get<string>("EVALUATION_TAU") <<"\n";
    cout << "\n";

    if(stabparams->get<string>("TDS") == "quasistatic")
    {
      if(stabparams->get<string>("TRANSIENT")=="yes_transient")
      {
        dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
      }
    }
    cout <<  "                             " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
    cout <<  "                             " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
    cout <<  "                             " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
    cout <<  "                             " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
    cout <<  "                             " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
    cout <<  "                             " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
    cout <<  "                             " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
    cout << endl;
    cout << "                             " << "Evaluation Mat  = " << stabparams->get<string>("EVALUATION_MAT") <<"\n";
    cout << "\n";
  }

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary) SolveStationaryProblem();
  else                                             TimeLoop();

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // ImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::TimeLoop()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint time loop");

  while (step_<stepmax_)
  {
    PrepareTimeStep();
    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch (timealgo_)
      {
      case INPAR::FLUID::timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta (%0.2f)   STEP = %4d/%4d \n",
              time_,maxtime_,dt_,theta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_afgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Af-Generalized-Alpha  STEP = %4d/%4d \n",
               time_,maxtime_,dt_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_npgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Np-Generalized-Alpha  STEP = %4d/%4d \n",
               time_,maxtime_,dt_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
               time_,maxtime_,dt_,step_,stepmax_);
        break;
      default:
        dserror("parameter out of range: IOP\n");
      } /* end of switch(timealgo) */
    }

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    NonLinearSolve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();


    // -------------------------------------------------------------------
    //  output
    // -------------------------------------------------------------------
//    Output(); TODO activate

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();
  } // end time loop
} // ImplicitTimeInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::PrepareTimeStep()
{

  // set time-dependent parameters
  IncrementTimeAndStep();

  // Set time parameter for element call
  SetElementTimeParameter();

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);

    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);

    discret_->ClearState();

// TODO make this work
//    // evaluate Neumann conditions
//    neumann_loads_->PutScalar(0.0);
//    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
//    discret_->ClearState();
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::NonLinearSolve()
{
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. adjoint solve");

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  int  itemax = 0;
  bool stopnonliniter = false;

  // REMARK:
  // commented reduced number of iterations out as it seems that more iterations
  // are necessary before sampling to obtain a converged result
//  // currently default for turbulent channel flow: only one iteration before sampling
//  if (special_flow_ == "channel_flow_of_height_2" && step_<samstart_ )
//       itemax  = 2;
//  else
  itemax  = params_.get<int>   ("max nonlin iter steps");

  dtsolve_  = 0.0;
  dtele_    = 0.0;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // time measurement: element
      TEUCHOS_FUNC_TIME_MONITOR("      + adjoint element calls");

      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      sysmat_->Zero();

      // create the parameters for the discretization
      ParameterList eleparams;

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      discret_->ClearState();

      // set action type
      eleparams.set("action","calc_adjoint_systemmat_and_residual");

      //set additional pseudo-porosity field for topology optimization
      eleparams.set("topopt_porosity",topopt_porosity_);

      // set fluid velocities of current and last time step
      RCP<Epetra_Vector> fluidveln = Teuchos::null;
      RCP<Epetra_Vector> fluidvelnp = Teuchos::null;
      if (timealgo_ == INPAR::FLUID::timeint_stationary)
      {
        fluidveln = fluid_vels_->find(1)->second;
        fluidvelnp = fluidveln;
      }
      else
      {
        // fluid starts at step 0 (=initial values), stops at stepmax
        // adjoint starts at step stepmax (=initial values), stops at step 0
        // evaluate fluid fields with according appropriate, reverse counted step
        fluidveln = fluid_vels_->find(stepmax_-step_)->second;
        fluidvelnp = fluid_vels_->find(stepmax_-step_+1)->second;
      }
      discret_->SetState("fluidveln",fluidveln);
      discret_->SetState("fluidvelnp",fluidvelnp);

      discret_->SetState("veln",veln_);
      discret_->SetState("velnp",velnp_);


      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          ||
          (params_.get<string>("CONVCHECK","L_2_norm")
              !=
                  "L_2_norm_without_residual_at_itemax"))
      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
        discret_->ClearState();

        //----------------------------------------------------------------------
        // account for potential Neumann inflow terms
        //----------------------------------------------------------------------
        if (neumanninflow_)
        {
          // create parameter list
          ParameterList condparams;

          // action for elements
          condparams.set("action","calc_adjoint_Neumann_inflow");

          // set vector values needed by elements
          discret_->ClearState();

          // set scheme-specific element parameters and vector values
          discret_->SetState("velnp",velnp_);

          std::string condstring("FluidNeumannInflow");
          discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
          discret_->ClearState();
        }

        // finalize the complete matrix
        sysmat_->Complete();
      }

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_->ExtractOtherVector(residual_);
    onlyvel->Norm2(&vresnorm_);

    velpressplitter_->ExtractOtherVector(incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2_);

    velpressplitter_->ExtractOtherVector(velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2_);

    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_->ExtractCondVector(residual_);
    onlypre->Norm2(&presnorm_);

    velpressplitter_->ExtractCondVector(incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2_);

    velpressplitter_->ExtractCondVector(velnp_,onlypre);
    onlypre->Norm2(&prenorm_L2_);

    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2_ < 1e-5) velnorm_L2_ = 1.0;
    if (prenorm_L2_ < 1e-5) prenorm_L2_ = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
            itnum,itemax,ittol,vresnorm_,presnorm_);
        printf(" (      --     ,te=%10.3E)\n",dtele_);
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
      if (vresnorm_ <= ittol and presnorm_ <= ittol and
          incvelnorm_L2_/velnorm_L2_ <= ittol and incprenorm_L2_/prenorm_L2_ <= ittol)
      {
        stopnonliniter=true;
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
              itnum,itemax,ittol,vresnorm_,presnorm_,
              incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

          FILE* errfile = params_.get<FILE*>("err file",NULL);
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                itnum,itemax,ittol,vresnorm_,presnorm_,
                incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
          }
        }
        break;
      }
      else // if not yet converged
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
              itnum,itemax,ittol,vresnorm_,presnorm_,
              incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm_ > ittol or presnorm_ > ittol or
        incvelnorm_L2_/velnorm_L2_ > ittol or
        incprenorm_L2_/prenorm_L2_ > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_.get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
              itnum,itemax,ittol,vresnorm_,presnorm_,
              incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
        }
      }
      break;
    }

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    {
      // time measurement: application of dbc
      TEUCHOS_FUNC_TIME_MONITOR("      + apply adjoint DBC");
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // time measurement: solver
      TEUCHOS_FUNC_TIME_MONITOR("      + adjoint solver calls");

      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = max(vresnorm_,presnorm_);
        currresidual = max(currresidual,incvelnorm_L2_/velnorm_L2_);
        currresidual = max(currresidual,incprenorm_L2_/prenorm_L2_);
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

      solver_.Solve(
          sysmat_->EpetraOperator(),
          incvel_,
          residual_,
          true,
          itnum==1,
          Teuchos::null,
          Teuchos::null,
          false);

      solver_.ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    velnp_->Update(1.0,*incvel_,1.0);
  } // end iteration loop
} // ImplicitTimeInt::LinearSolve



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                                      |
 | One-step-Theta: (step>1)                                             |
 |  velnm_ =veln_                                                       |
 |  veln_  =velnp_                                                      |
 |                                                      winklmaier 03/12|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::TimeUpdate()
{
  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  veln_ ->Update(1.0,*velnp_,0.0);

  return;
}// ImplicitTimeInt::TimeUpdate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::Output()
{

  //  ART_exp_timeInt_->Output();
  // output of solution
  if (step_%upres_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);
    output_.WriteVector("tract_resid",residual_);
    output_.WriteVector("neumann_loads",neumann_loads_);
    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_->ExtractCondVector(velnp_);
    output_.WriteVector("pressure", pressure);

#ifdef GMSHOUTPUT
    OutputToGmsh(step_, time_,false);
#endif

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      output_.WriteVector("veln", veln_);
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);

    // velocity/pressure vector at time n and n-1
    output_.WriteVector("veln", veln_);
    output_.WriteVector("neumann_loads",neumann_loads_);
  }

  if (topopt_porosity_!=Teuchos::null)
    optimizer_->ImportAdjointFluidData(velnp_,step_);

  return;
} // ImplicitTimeInt::Output


void TOPOPT::ADJOINT::ImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time,
    const bool inflow
    ) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  // 20 steps are kept
  std::string filename = "dummy";
  if (inflow)
  {
    filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_velpres_inflow", step, 20, screen_out, discret_->Comm().MyPID());
    //std::ofstream gmshfilecontent(filename.c_str());
  }
  else
  {
    filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_velpres", step, 20, screen_out, discret_->Comm().MyPID());
    //std::ofstream gmshfilecontent(filename.c_str());
  }
  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "velocity solution \" {" << endl;
    IO::GMSH::VelocityPressureFieldDofBasedToGmsh(discret_, velnp_ , "velocity", gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "pressure solution\" {" << endl;
    IO::GMSH::VelocityPressureFieldDofBasedToGmsh(discret_, velnp_, "pressure",gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << endl;

 return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::ReadRestart(int step)
{
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |set restart values (turbulent inflow only)             rasthofer 06/11|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::SetRestart(
  const int step,
  const double time,
  Teuchos::RCP<const Epetra_Vector> readvelnp,
  Teuchos::RCP<const Epetra_Vector> readveln,
  Teuchos::RCP<const Epetra_Vector> readvelnm,
  Teuchos::RCP<const Epetra_Vector> readaccnp,
  Teuchos::RCP<const Epetra_Vector> readaccn)
{
  time_ = time;
  step_ = step;

  velnp_->Update(1.0,*readvelnp,0.0);
  veln_->Update(1.0,*readveln,0.0);

  SetElementTimeParameter();
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set initial flow field for test cases                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = discret_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // initialize veln_ as well. That's what we actually want to do here!
    veln_->Update(1.0,*velnp_ ,0.0);

    // add random perturbation of certain percentage to function
    if (initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      int err =0;

      // random noise is perc percent of the initial profile
      double perc = params_.sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST",0.1);

      // out to screen
      if (myrank_==0)
      {
        cout << "Disturbed initial profile:   max. " << perc*100 << "% random perturbation\n";
        cout << "\n\n";
      }

      double bmvel=0;
      double mybmvel=0;
      double thisvel=0;
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = discret_->Dof(lnode);

        for(int index=0;index<numdim_;++index)
        {
          int gid = nodedofset[index];
          int lid = dofrowmap->LID(gid);

          thisvel=(*velnp_)[lid];
          if (mybmvel*mybmvel < thisvel*thisvel) mybmvel=thisvel;
        }
      }

      // the noise is proportional to the bulk mean velocity of the
      // undisturbed initial field (=2/3*maximum velocity)
      mybmvel=2*mybmvel/3;
      discret_->Comm().MaxAll(&mybmvel,&bmvel,1);

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

        // add random noise on initial function field
        for(int index=0;index<numdim_;++index)
        {
          int gid = nodedofset[index];

          double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);

          double noise = perc * bmvel * randomnumber;

          err += velnp_->SumIntoGlobalValues(1,&noise,&gid);
          err += veln_ ->SumIntoGlobalValues(1,&noise,&gid);
        }

        if(err!=0)
        {
          dserror("dof not on proc");
        }
      }
    }
  }
  else
  {
    dserror("Only initial fields by (un-)disturbed functions are available up to now!");
  }

  return;
} // end SetInitialFlowField



/*----------------------------------------------------------------------*
 | sent density field for topology optimization         winklmaier 12/11|
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::SetTopOptData(
    Teuchos::RCP<std::map<int,RCP<Epetra_Vector> > > fluidvelocities,
    RCP<Epetra_Vector> porosity,
    RCP<TOPOPT::Optimizer> optimizer
)
{
  fluid_vels_= fluidvelocities;
  topopt_porosity_ = porosity;
  optimizer_=optimizer;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
{

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params_,"calculate error");

  switch(calcerr)
  {
  case INPAR::FLUID::no_error_calculation:
    // do nothing --- no analytical solution available
    break;
  case INPAR::FLUID::beltrami_flow:
  case INPAR::FLUID::channel2D:
  case INPAR::FLUID::gravitation:
  case INPAR::FLUID::shear_flow:
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_adjoint_error");
    eleparams.set<int>("calculate error",calcerr);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",velnp_);

    // get (squared) error values
    // 0: vel_mag
    // 1: p
    // 2: u_mag,analytical
    // 3: p_analytic
    // (4: vel_x)
    // (5: vel_y)
    // (6: vel_z)
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(2+2));
    //  = Teuchos::rcp(new Epetra_SerialDenseVector(numdim_+2+2))

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    double velerr = 0.0;
    double preerr = 0.0;

    // integrated analytic solution in order to compute relative error
    double velint = 0.0;
    double pint = 0.0;

    // error in the single velocity components
    //double velerrx = 0.0;
    //double velerry = 0.0;
    //double velerrz = 0.0;

    // for the L2 norm, we need the square root
    velerr = sqrt((*errors)[0]);
    preerr = sqrt((*errors)[1]);

    // analytical vel_mag and p_mag
    velint= sqrt((*errors)[2]);
    pint = sqrt((*errors)[3]);

    if (myrank_ == 0)
    {
      {
        cout.precision(8);
        cout << endl << "----relative L_2 error norm for analytical solution Nr. " <<
          DRT::INPUT::get<INPAR::FLUID::CalcError>(params_,"calculate error") <<
          " ----------" << endl;
        cout << "| velocity:  " << velerr/velint << endl;
        cout << "| pressure:  " << preerr/pint << endl;
        cout << "--------------------------------------------------------------------" << endl << endl;
      }

      //velerrx = sqrt((*errors)[4]);
      //velerry = sqrt((*errors)[5]);
      //if (numdim_==3)
      //  velerrz = sqrt((*errors)[6]);

      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation+".relerror";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }

      ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.relerror";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());
        f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  }
  return;
} // end EvaluateErrorComparedToAnalyticalSol

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | solve stationary fluid problem                              gjb 10/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::SolveStationaryProblem()
{
  // time measurement: time loop (stationary) --- start TimeMonitor tm2
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint time loop");

  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  while (step_< stepmax_)
  {
    // -------------------------------------------------------------------
    //              set (pseudo-)time-dependent parameters
    // -------------------------------------------------------------------
    IncrementTimeAndStep();

   // -------------------------------------------------------------------
   //                         out to screen
   // -------------------------------------------------------------------
   if (myrank_==0)
   {
    printf("Stationary Fluid Adjoint Solver - STEP = %4d/%4d \n",step_,stepmax_);
   }

    SetElementTimeParameter();

    // -------------------------------------------------------------------
    //         evaluate Dirichlet and Neumann boundary conditions
    // -------------------------------------------------------------------
    {
      ParameterList eleparams;

      // other parameters needed by the elements
      eleparams.set("total time",time_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velaf",velnp_);
      // predicted dirichlet values
      // velnp then also holds prescribed new dirichlet values
      discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);

      discret_->ClearState();

      // TODO make this work
//      neumann_loads_->PutScalar(0.0);
//      discret_->EvaluateNeumann(eleparams,*neumann_loads_);
//      discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonLinearSolve();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
//    Output(); TODO activate
  } // end of time loop
} // ImplicitTimeInt::SolveStationaryProblem



/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
TOPOPT::ADJOINT::ImplicitTimeInt::~ImplicitTimeInt()
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::VelocityRowMap()
{
  return velpressplitter_->OtherMap();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::PressureRowMap()
{
  return velpressplitter_->CondMap();
}



double TOPOPT::ADJOINT::ImplicitTimeInt::ResidualScaling() const
{
  if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
    return 1.0;
  else if (TimIntScheme()==INPAR::FLUID::timeint_one_step_theta)
    return 1.0/(theta_*dt_);
  else
    dserror("time integration scheme not implemented");

  return 0.0;
}



void TOPOPT::ADJOINT::ImplicitTimeInt::SetVelocityField(
    Teuchos::RCP<const Epetra_Vector> setvelnp
)
{
  velnp_->Update(1.0,*setvelnp,0.0); return;
}



Teuchos::RCP<LINALG::SparseMatrix> TOPOPT::ADJOINT::ImplicitTimeInt::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}



Teuchos::RCP<LINALG::BlockSparseMatrixBase> TOPOPT::ADJOINT::ImplicitTimeInt::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}



// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementGeneralAdjointParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_general_adjoint_parameter");

  // set if objective contains dissipation
  eleparams.set<bool>("dissipation",params_.get<bool>("OBJECTIVE_DISSIPATION"));
  // set if objective contains inlet pressure
  eleparams.set<bool>("inletPres",params_.get<bool>("OBJECTIVE_INLET_PRESSURE"));
  // set if objective contains pressure drop
  eleparams.set<bool>("presDrop",params_.get<bool>("OBJECTIVE_PRESSURE_DROP"));
  // set objective's dissipation factor
  eleparams.set<double>("dissipationFac" ,params_.get<double>("DISSIPATION_FAC"));
  // set objective's inlet pressure factor
  eleparams.set<double>("inletPresFac" ,params_.get<double>("PRESSURE_INLET_FAC"));
  // set objective's pressure drop factor
  eleparams.set<double>("presDropFac" ,params_.get<double>("PRESSURE_DROP_FAC"));

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementTimeParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_adjoint_time_parameter");

  // set general element parameters
  eleparams.set("dt",dt_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",omtheta_);
  eleparams.set("theta_pre",theta_pre_);
  eleparams.set("omtheta_pre",omtheta_pre_);
  eleparams.set("theta_div",theta_div_);
  eleparams.set("omtheta_div",omtheta_div_);

  // set scheme-specific element parameters and vector values
  if ((timealgo_==INPAR::FLUID::timeint_stationary) ||
      (timealgo_==INPAR::FLUID::timeint_one_step_theta))
  {
    eleparams.set("total time",time_);
  }
  else
    dserror("time integration scheme not implemented");

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}
#endif  // #ifdef CCADISCRET
