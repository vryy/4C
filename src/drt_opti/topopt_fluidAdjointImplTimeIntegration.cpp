/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointImplTimeIntegration.cpp

\brief fluid adjoint implicit time integration for topology optimization


<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_fluidAdjointImplTimeIntegration.H"
#include "topopt_fluidAdjointResulttest.H"
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/time_integration_scheme.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_opti/topopt_optimizer.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                               winklmaier 03/12 |
 *----------------------------------------------------------------------*/
TOPOPT::ADJOINT::ImplicitTimeInt::ImplicitTimeInt(
    RCP<DRT::Discretization>      actdis,
    RCP<LINALG::Solver>           solver,
    RCP<ParameterList>            params,
    RCP<IO::DiscretizationWriter> output)
  : FluidAdjointTimeInt(actdis,solver,params,output)
{
  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint initialization");

  // parameter theta for time-integration schemes
  theta_      = params_->get<double>("theta");
  theta_pre_  = params_->get<double>("theta_pre");
  theta_div_  = params_->get<double>("theta_div");
  // compute or set 1.0 - theta for time-integration schemes
  omtheta_ = 1.0 - theta_;
  omtheta_pre_ = 1.0 - theta_pre_;
  omtheta_div_ = 1.0 - theta_div_;


  // -------------------------------------------------------------------
  // account for potential Neuman inflow terms if required
  // -------------------------------------------------------------------
  neumanninflow_ = false;
  if (params_->get<string>("Neumann inflow","no") == "yes") neumanninflow_ = true;

  // -------------------------------------------------------------------
  // care for periodic boundary conditions
  // -------------------------------------------------------------------
  pbcmapmastertoslave_ = params_->get<RCP<map<int,vector<int> > > >("periodic bc");
  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

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

  if (not params_->get<int>("Simple Preconditioner",0) &&
      not params_->get<int>("AMG BS Preconditioner",0))
  {
    // initialize standard (stabilized) system matrix
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  }
  else
  {
    dserror("understand this before using it :)");
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
  // adjoint velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*dofrowmap,true);
  veln_  = LINALG::CreateVector(*dofrowmap,true);
  velnm_ = LINALG::CreateVector(*dofrowmap,true);

  // adjoint velocity/pressure at time n+1, n and n-1
  fluidvelnp_ = LINALG::CreateVector(*dofrowmap,true);
  fluidveln_  = LINALG::CreateVector(*dofrowmap,true);

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

    // the dirichlet values of the standard fluid are contained in zeros_
    // -> clear them since they are unused
    zeros_->PutScalar(0.0);
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


/*----------------------------------------------------------------------*
 | Start the time integration                          winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::Integrate()
{
  // output of stabilization details
  if (myrank_==0)
  {
    cout << "\n\n\n";
    cout << "           ------------------------------------------------------------------------------           \n";
    cout << "          |                          solving adjoint equations                           |          \n";
    cout << "           ------------------------------------------------------------------------------           \n";
    cout << "\n" << endl;
  }

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary) SolveStationaryProblem();
  else                                             TimeLoop();

  return;
} // ImplicitTimeInt::Integrate



/*----------------------------------------------------------------------*
 | contains the time loop                              winklmaier 03/12 |
 *----------------------------------------------------------------------*/
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
    Output();

  } // end time loop
} // ImplicitTimeInt::TimeLoop


/*----------------------------------------------------------------------*
 | solve the stationary problem                        winklmaier 03/12 |
 *----------------------------------------------------------------------*/
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
    PrepareTimeStep();

    // output to screen
    if (myrank_==0)
    {
      printf("Stationary Fluid Adjoint Solver - STEP = %4d/%4d \n",step_,stepmax_);
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonLinearSolve();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();
  } // end of time loop
}



/*----------------------------------------------------------------------*
 | prepare one time step of instationary simulation    winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::PrepareTimeStep()
{
  // set (pseudo-)time-dependent parameters
  IncrementTimeAndStep();

  // Set time parameter for element call
  SetElementTimeParameter();

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    TEUCHOS_FUNC_TIME_MONITOR("      + evaluate adjoint boundary conditions");

    // we evaluate the Dirichlet boundary condition(s) without the
    // discretization since the dirichlet values are NOT based on
    // input functions (timecurve and function) instead the objective's
    // data is required
    EvaluateDirichlet();

    // we evaluate the Neumann boundary condition(s) with the evaluate-condition
    // function since the Neumann terms are different for the adjoints and
    // usually the evaluateNeumann enters the fluid Neumann conditions where
    // we want another condition
    ParameterList nbcparams;

    // set action for elements
    nbcparams.set<int>("action",FLD::ba_calc_adjoint_neumann);

    // set flag for test case
    nbcparams.set("special test case",params_->get<INPAR::TOPOPT::AdjointTestCases>("special test case"));

    // in instationary case: fluid starts at step 0 (=initial values), stops at
    // stepmax adjoint starts at step stepmax (=initial values), stops at step 0
    // evaluate fluid fields with according appropriate, reverse counted step
    if (timealgo_ == INPAR::FLUID::timeint_stationary)
    {
      fluidveln_ = fluid_vels_->find(1)->second;
      fluidvelnp_ = fluidveln_;
    }
    else
    {
      fluidveln_ = fluid_vels_->find(stepmax_-step_+1)->second;
      fluidvelnp_ = fluid_vels_->find(stepmax_-step_)->second;
    }

    discret_->SetState("fluidveln",fluidveln_);
    discret_->SetState("fluidvelnp",fluidvelnp_);

    discret_->SetState("veln",veln_);
    discret_->SetState("velnp",velnp_);

    neumann_loads_->PutScalar(0.0);

    if (DRT::Problem::Instance()->NDim()==2) // 2D -> 1D (line) neumann surface
    {
      discret_->EvaluateConditionUsingParentData(
          nbcparams            ,
          Teuchos::null        ,
          Teuchos::null        ,
          neumann_loads_       ,
          Teuchos::null        ,
          Teuchos::null        ,
          "LineNeumann");
    }
    else if (DRT::Problem::Instance()->NDim()==3) // 3D -> 2D (surface) neumann surface
    {
      discret_->EvaluateConditionUsingParentData(
          nbcparams            ,
          Teuchos::null        ,
          Teuchos::null        ,
          neumann_loads_       ,
          Teuchos::null        ,
          Teuchos::null        ,
          "SurfaceNeumann");
    }
    else
      dserror("Dimension not implemented");

    // clear state
    discret_->ClearState();
  }

  return;
}


/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop               winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::NonLinearSolve()
{
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. adjoint solve");

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_->get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_->get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

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
  itemax  = params_->get<int>   ("max nonlin iter steps");

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

      // reset sysmat and residual and add neumann terms
      sysmat_->Zero();

      residual_->Update(1.0,*neumann_loads_,0.0);

      // create the parameters for the discretization
      ParameterList eleparams;

      // set action type
      eleparams.set<int>("action",FLD::calc_adjoint_systemmat_and_residual);

      //set additional pseudo-porosity field for topology optimization
      eleparams.set("topopt_porosity",topopt_porosity_);

      discret_->ClearState();

      discret_->SetState("fluidveln",fluidveln_);
      discret_->SetState("fluidvelnp",fluidvelnp_);

      discret_->SetState("veln",veln_);
      discret_->SetState("velnp",velnp_);


      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          ||
          (params_->get<string>("CONVCHECK","L_2_norm")
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

          FILE* errfile = params_->get<FILE*>("err file",NULL);
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
      {
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
              itnum,itemax,ittol,vresnorm_,presnorm_,
              incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
        }
      }
//      dserror("stop here");
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

        FILE* errfile = params_->get<FILE*>("err file",NULL);
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
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

      solver_->Solve(
          sysmat_->EpetraOperator(),
          incvel_,
          residual_,
          true,
          itnum==1,
          Teuchos::null,
          Teuchos::null,
          false);

      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    velnp_->Update(1.0,*incvel_,1.0);
  } // end iteration loop
} // ImplicitTimeInt::NonLinearSolve



/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                                      |
 | One-step-Theta: (step>1)                                             |
 |  veln_  =velnp_                                                      |
 |                                                     winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::TimeUpdate()
{
  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

  return;
}// ImplicitTimeInt::TimeUpdate



/*----------------------------------------------------------------------*
 | output of solution vector to binio                   winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::Output() const
{
  // output of solution
  if (step_%upres_ == 0)
  {
    // step number and time
    output_->NewStep(step_,time_);
    // velocity/pressure vector
    output_->WriteVector("adjoint_velnp",velnp_);
    output_->WriteVector("adjoint_neumann_loads",neumann_loads_);
    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_->ExtractCondVector(velnp_);
    output_->WriteVector("adjoint_pressure", pressure);

    if (params_->get<bool>("GMSH_OUTPUT"))
      OutputToGmsh(step_, time_,false);

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_->WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      output_->WriteVector("adjoint_veln", veln_);
      output_->WriteVector("adjoint_velnm",velnm_);

    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_->NewStep(step_,time_);

    // velocity/pressure vector
    output_->WriteVector("adjoint_velnp",velnp_);

    // velocity/pressure vector at time n and n-1
    output_->WriteVector("adjoint_veln", veln_);
    output_->WriteVector("adjoint_velnm",velnm_);
    output_->WriteVector("adjoint_neumann_loads",neumann_loads_);
  }

  if (topopt_porosity_!=Teuchos::null)
  {
    if (timealgo_==INPAR::FLUID::timeint_stationary)
      optimizer_->ImportAdjointFluidData(velnp_,1);
    else
    {
      // we have to switch the step here since we start with step 0 at end-time,
      // but the step of adjoint fluid and fluid equations shall fit
      optimizer_->ImportAdjointFluidData(velnp_,stepmax_-step_);

      // initial solution (=v_t_max) is old solution at step 1
      if (step_==1)
        optimizer_->ImportAdjointFluidData(velnm_,stepmax_); // currently velnm contains veln because timeupdate was called before
    }
  }

  return;
} // ImplicitTimeInt::Output



/*----------------------------------------------------------------------*
 | output of solution vector to gmsh                   winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time,
    const bool inflow
) const
{
  // get Gmsh postprocessing filename
  std::ostringstream filename;
  const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());
  std::ostringstream pid_stream;
  pid_stream << ".p" << std::setw(2) << setfill('0') << myrank_;
  filename    << filebase << "_adjoint.solution_velpres_" << std::setw(5) << setfill('0') << step           << pid_stream.str() << ".pos";

  // delete old Gmsh postprocessing file
  int step_diff = 20; // stepdiff files are kept
  std::ostringstream filenamedel;
  filenamedel << filebase << "_adjoint_" << std::setw(5) << setfill('0') << step-step_diff << pid_stream.str() << ".pos";
  std::remove(filenamedel.str().c_str());

  std::ofstream gmshfilecontent(filename.str().c_str());
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

 return;
}


/*----------------------------------------------------------------------*
 | Read restart of instationary simulation             winklmaier 03/12 |
 -----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::ReadRestart(int step)
{
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();
}


/*----------------------------------------------------------------------*
 |  set initial flow field for test cases              winklmaier 03/12 |
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::SetInitialAdjointField(
  const INPAR::TOPOPT::InitialAdjointField initfield,
  const int startfuncno
  )
{
  // initial zero field (the default)
  switch (initfield)
  {
  case INPAR::TOPOPT::initadjointfield_zero_field:
  {
    veln_->PutScalar(0.0);
    velnp_->PutScalar(0.0);
    break;
  }
  case INPAR::TOPOPT::initadjointfield_field_by_function:// initial field by undisturbed
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
        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);
        veln_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // if we have no test case we set the new solution equal to the old although we have a linear system
    // -> smaller numerical errors around 10^-12
    // for test cases we start with zero field so that more terms are active and tested
//    if (params_->get<INPAR::TOPOPT::AdjointTestCases>("special test case") == INPAR::TOPOPT::adjointtest_no)
      *velnp_ = *veln_;

    break;
  }
  default:
    dserror("Type of initial field not available up to now!");
  }

  return;
} // end SetInitialAdjointField



/*----------------------------------------------------------------------*
 | sent density field for topology optimization         winklmaier 12/11|
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::SetTopOptData(
    Teuchos::RCP<const std::map<int,RCP<Epetra_Vector> > > fluidvelocities,
    RCP<const Epetra_Vector> porosity,
    RCP<TOPOPT::Optimizer>& optimizer
)
{
  fluid_vels_= fluidvelocities;
  topopt_porosity_ = porosity;
  optimizer_=optimizer;
}



/*----------------------------------------------------------------------*
 | evaluate dirichlet condition                         winklmaier 03/12|
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::EvaluateDirichlet()
{
  DRT::Node* node = NULL;


  if (params_->get<INPAR::TOPOPT::AdjointTestCases>("special test case") == INPAR::TOPOPT::adjointtest_no)
  {
    // TODO dbc due to objective function
  }
  else // special cases
  {
    INPAR::TOPOPT::AdjointTestCases testcase = params_->get<INPAR::TOPOPT::AdjointTestCases>("special test case");

    for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
    {
      node = discret_->lRowNode(inode);
      DRT::Condition* cond = node->GetCondition("Dirichlet");
      if (cond==NULL) continue;

      const int nsd = DRT::Problem::Instance()->NDim();
      double values[3] = {0.0}; // dimension is <= 3 so this is enough

      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      // get global coordinates of gauss point
      if (nsd == 2)
      {
        LINALG::Matrix<2,1> coords(node->X());
        x = coords(0);
        y = coords(1);
      }
      else if (nsd == 3)
      {
        LINALG::Matrix<3,1> coords(node->X());
        x = coords(0);
        y = coords(1);
        z = coords(2);
      }

      switch (testcase)
      {
      case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
      {
        values[0] = 1.0;
        values[1] = 0.0;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
      {
        values[0] = 5.0*x + 2.0*y;
        values[1] = 3.0*x + 7.0*y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
      {
        values[0] = x*x;
        values[1] = 2*y*y - 3*x*x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
      {
        values[0] = x*x + x*y;
        values[1] = 2*y*y - 3*x*x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_varying_theta:
      {
        double t = time_;
        values[0] = x-t;
        values[1] = 2*t-y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
      {
        double t = time_;
        values[0] = x*x + x*y*t;
        values[1] = 2*y*y*t*t - 3*x*x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
      {
        double t = time_;
        values[0] = x*y*t + x;
        values[1] = y*t - x*t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_primal:
        break;
      default:
        dserror("no dirichlet condition implemented for special test case");
      }

      vector<int> gdofs = discret_->Dof(node);
      velnp_->ReplaceGlobalValues(nsd,(double*)values,&gdofs[0]); // &dofs[0] gives pointer to dofs
    }
  }
}

Teuchos::RCP<DRT::ResultTest> TOPOPT::ADJOINT::ImplicitTimeInt::CreateFieldTest() const
{
  return Teuchos::rcp(new TOPOPT::ADJOINT::FluidAdjointResultTest(*this));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::VelocityRowMap() const
{
  return velpressplitter_->OtherMap();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::PressureRowMap() const
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



Teuchos::RCP<const LINALG::SparseMatrix> TOPOPT::ADJOINT::ImplicitTimeInt::SystemMatrix() const
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}



Teuchos::RCP<const LINALG::BlockSparseMatrixBase> TOPOPT::ADJOINT::ImplicitTimeInt::BlockSystemMatrix() const
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}



// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementGeneralAdjointParameter() const
{
  // get the optimization material
  const MAT::PAR::NewtonianFluid* mat = NULL;
  const int nummat = DRT::Problem::Instance()->Materials()->Num();
  for (int id = 1; id-1 < nummat; ++id)
  {
    Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

    if (imat == Teuchos::null)
      dserror("Could not find material Id %d", id);
    else
    {
      if (imat->Type() == INPAR::MAT::m_fluid)
      {
        const MAT::PAR::Parameter* matparam = imat->Parameter();
        mat = static_cast<const MAT::PAR::NewtonianFluid* >(matparam);
        break;
      }
    }
  }
  if (mat==NULL)
    dserror("optimization material not found");

  // and then the material parameters
  const double density = mat->density_;
  const double viscosity = mat->viscosity_;


  ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_general_adjoint_parameter);

  // set material parameters
  eleparams.set<double>("density" ,density);
  eleparams.set<double>("viscosity" ,viscosity);

  // set if objective contains dissipation
  eleparams.set<bool>("dissipation",params_->get<bool>("OBJECTIVE_DISSIPATION"));
  // set if objective contains pressure drop
  eleparams.set<bool>("presDrop",params_->get<bool>("OBJECTIVE_PRESSURE_DROP"));
  // set objective's dissipation factor
  eleparams.set<double>("dissipationFac" ,params_->get<double>("DISSIPATION_FAC"));
  // set objective's pressure drop factor
  eleparams.set<double>("presDropFac" ,params_->get<double>("PRESSURE_DROP_FAC"));

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set flag for test cases
  eleparams.set<INPAR::TOPOPT::AdjointTestCases>("special test case",params_->get<INPAR::TOPOPT::AdjointTestCases>("special test case"));

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementTimeParameter() const
{
  ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_adjoint_time_parameter);

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
