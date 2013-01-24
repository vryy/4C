
/*!----------------------------------------------------------------------
\file airwayimplicitintegration.cpp
\brief Control routine for reduced airway solvers,

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*----------------------------------------------------------------------*/
#include <stdio.h>

#include "airwayimplicitintegration.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

AIRWAY::RedAirwayImplicitTimeInt::RedAirwayImplicitTimeInt(RCP<DRT::Discretization>  actdis,
                                                           LINALG::Solver  &         solver,
                                                           Teuchos::ParameterList&   params,
                                                           IO::DiscretizationWriter& output)
  :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  time_(0.0),
  step_(0),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  coupledTo3D_(false)
{

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // time measurement: initialization
  if(!coupledTo3D_)
  {
    // time measurement: initialization
    TEUCHOS_FUNC_TIME_MONITOR(" + initialization");
  }

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = dtp_*double(stepmax_);
  // maximum iteration steps
  maxiter_     = params_.get<int> ("maximum iteration steps");
  // tolerance of nonlinear solution
  non_lin_tol_ = params_.get<double> ("tolerance");


  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap      = discret_->DofRowMap();
  const Epetra_Map* dofcolmap      = discret_->DofColMap();
  const Epetra_Map* elementcolmap  = discret_->ElementColMap();
  const Epetra_Map* elementrowmap  = discret_->ElementRowMap();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the volumetric flow rate dofs and for one vector which only
  // contains cross-sectional area degrees of freedom.
  // -------------------------------------------------------------------


  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Each node has 3 adjacent nodes (including itself), each
  // with 1 dofs. (3*1=3)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  // initialize standard (stabilized) system matrix
  sysmat_  = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,3,false,true));

  // Vectors passed to the element
  // -----------------------------

  // Pressures at time n+1, n and n-1
  pnp_          = LINALG::CreateVector(*dofrowmap,true);
  pn_           = LINALG::CreateVector(*dofrowmap,true);
  pnm_          = LINALG::CreateVector(*dofrowmap,true);

  p_nonlin_     = LINALG::CreateVector(*dofrowmap,true);
  sysmat_iad_   = LINALG::CreateVector(*dofrowmap,true);

  // Inlet volumetric flow rates at time n+1, n and n-1
  qin_np_       = LINALG::CreateVector(*elementcolmap,true);
  qin_n_        = LINALG::CreateVector(*elementcolmap,true);
  qin_nm_       = LINALG::CreateVector(*elementcolmap,true);

  // outlet volumetric flow rates at time n+1, n and n-1
  qout_np_      = LINALG::CreateVector(*elementcolmap,true);
  qout_n_       = LINALG::CreateVector(*elementcolmap,true);
  qout_nm_      = LINALG::CreateVector(*elementcolmap,true);

  // This vector will be used for exportation and restart reasons
  qexp_         = LINALG::CreateVector(*elementrowmap,true);
  pexp_         = LINALG::CreateVector(*dofrowmap,true);

  // This vector will be used to test convergence
  residual_     = LINALG::CreateVector(*dofrowmap,true);
  bc_residual_  = LINALG::CreateVector(*dofcolmap,true);

  // Volumetric flow rates at time n+1, n and n-1
  //  qcnp_          = LINALG::CreateVector(*elementrowmap,true);
  //  qcn_           = LINALG::CreateVector(*elementrowmap,true);
  //  qcnm_          = LINALG::CreateVector(*elementrowmap,true);

  // vectors for postprocessing, Element Node Ids, radii, generations, etc ...
  nodeIds_      = LINALG::CreateVector(*dofrowmap,true);
  radii_        = LINALG::CreateVector(*dofrowmap,true);
  generations_  = LINALG::CreateVector(*elementcolmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  // This part might be optimized later
  bcval_   = LINALG::CreateVector(*dofrowmap,true);
  dbctog_  = LINALG::CreateVector(*dofrowmap,true);

  acini_bc_             = LINALG::CreateVector(*elementcolmap,true);
  acini_e_volume0_      = LINALG::CreateVector(*elementcolmap,true);
  acini_e_volumenm_     = LINALG::CreateVector(*elementcolmap,true);
  acini_e_volumen_      = LINALG::CreateVector(*elementcolmap,true);
  acini_e_volumenp_     = LINALG::CreateVector(*elementcolmap,true);
  acini_e_volume_strain_= LINALG::CreateVector(*elementcolmap,true);

  num_of_inter_acinar_linkers_ = LINALG::CreateVector(*dofrowmap,true);
  // Vectors used for solution process
  // ---------------------------------


  // right hand side vector and right hand side corrector
  rhs_     = LINALG::CreateVector(*dofrowmap,true);

  // ---------------------------------------------------------------------------------------
  // Initialize all the arteries' cross-sectional areas to the initial crossectional area Ao
  // and the volumetric flow rate to 0
  // ---------------------------------------------------------------------------------------
  ParameterList eleparams;

  // loop all elements and initialize all of the values

  eleparams.set("p0np",pnp_);
  eleparams.set("p0n",pn_);
  eleparams.set("p0nm",pnm_);

  eleparams.set("radii",radii_);
  eleparams.set("generations",generations_);
  eleparams.set("acini_bc",acini_bc_);
  eleparams.set("acini_e_volume",acini_e_volumenp_);

  eleparams.set("action","get_initial_state");
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  acini_e_volumen_->Update(1.0,*acini_e_volumenp_,0.0);
  acini_e_volumenm_->Update(1.0,*acini_e_volumenp_,0.0);
  acini_e_volume0_->Update(1.0,*acini_e_volumenp_,0.0);

  // Fill the NodeId vector
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lColElement(nele);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmstride;
    //vector<int> lmowner;
    RCP<std::vector<int> > lmowner = Teuchos::rcp(new std::vector<int>);
    ele->LocationVector(*discret_,lm,*lmowner,lmstride);

    // loop all nodes of this element, add values to the global vectors

    if(myrank_ == (*lmowner)[0])
    {
      int    gid = lm[0];
      double val = gid;
      nodeIds_->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank_ == (*lmowner)[1])
    {
      int    gid = lm[1];
      double val = gid;
      nodeIds_->ReplaceGlobalValues(1,&val,&gid);
    }

#if 0
  // get number of inter-acinar linkers effecting an acinus which should be:
    // ((number_of_elements) - 1)/2
//    if (ele->Name()=="RedAcinusType")
    {
      for(int j=1;j<2;j++)
      {
        if(myrank_ == (*lmowner)[j])
        {
          DRT::Node ** nodes = ele->Nodes();
          DRT::Element ** elems = nodes[j]->Elements();
          for (int k=0;k<nodes[j]->NumElement();k++)
          {
            if (elems[k]->Owner()==myrank_)
            {
              int gid = lm[j];
              int lid = num_of_inter_acinar_linkers_->Map().LID(gid);
              double val = (*num_of_inter_acinar_linkers_)[lid]+1;
              num_of_inter_acinar_linkers_->ReplaceGlobalValues(1,&val,&gid);
            }
          }
        }
      }
    }
#endif
  }

#if 1
  const Epetra_Map* nodeColwmap      = discret_->NodeColMap();
  for (int nnode=0;nnode<discret_->NumMyColNodes();++nnode)
  {
    DRT::Node * node = discret_->lColNode(nnode);
    int dofgid = discret_->Dof(node,0) ;
    int doflid = num_of_inter_acinar_linkers_->Map().LID(dofgid);

    for (int k=0;k<node->NumElement();k++)
    {
      {
        double val = (*num_of_inter_acinar_linkers_)[doflid]+1;
        num_of_inter_acinar_linkers_->ReplaceGlobalValues(1,&val,&dofgid);
      }
    }
  }
#endif

} // RedAirwayImplicitTimeInt::RedAirwayImplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration.                                          |
 |                                                          ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Integrate()
{
  RCP<Teuchos::ParameterList> param;
  Integrate(false, param);
} //RedAirwayImplicitTimeInt::Integrate()


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration.                                          |
 |                                                          ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Integrate(bool CoupledTo3D,
                                                 RCP<Teuchos::ParameterList> CouplingParams)
{
  coupledTo3D_ = CoupledTo3D;
  if (CoupledTo3D && CouplingParams.get() == NULL)
  {
    dserror("Coupling parameter list is not allowed to be empty, If a 3-D/reduced-D coupling is defined\n");
  }

  TimeLoop(CoupledTo3D,CouplingParams);

  // print the results of time measurements
  if (!coupledTo3D_)
  {
    TimeMonitor::summarize();
  }

  return;
} // RedAirwayImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                   ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::TimeLoop(bool CoupledTo3D,
                                                RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  coupledTo3D_ = CoupledTo3D;
  // time measurement: time loop
  if(!coupledTo3D_)
  {
    TEUCHOS_FUNC_TIME_MONITOR(" + time loop");
  }

  while (step_<stepmax_ and time_<maxtime_)
  {
    this->OneStepTimeLoop(CoupledTo3D, CouplingTo3DParams);

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
    if (CoupledTo3D)
    {
      break;
    }
  }

} // RedAirwayImplicitTimeInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the one step time loop                          ismail 09/12|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::OneStepTimeLoop(bool CoupledTo3D,
                                                RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  coupledTo3D_ = CoupledTo3D;

  double time3D = time_;
  if(coupledTo3D_)
  {
    time3D  = CouplingTo3DParams->get<double>("time");
  }
  if(time3D!=time_ || !coupledTo3D_)
  {
    PrepareTimeStep();
  }
  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    if(!coupledTo3D_)
    {
      printf("TIME: %11.4E/%11.4E  DT = %11.4E   Solving Reduced Dimensional Airways    STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      }
      else
      {
        printf("SUBSCALE_TIME: %11.4E/%11.4E  SUBSCALE_DT = %11.4E   Solving Reduced Dimensional Airways    SUBSCALE_STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
      }
    }

  // -------------------------------------------------------------------
  //                            SolverParams
  // -------------------------------------------------------------------
  if(params_.get<string> ("solver type") == "Nonlinear")
  {
    NonLin_Solve(CouplingTo3DParams);
    if (!myrank_)
      cout<<endl;
  }
  else if (params_.get<string> ("solver type") == "Linear")
  {
    Solve(CouplingTo3DParams);
    if (!myrank_)
      cout<<endl<<endl;
  }
  else
  {
    dserror("[%s] is not a defined solver",(params_.get<string> ("solver type")).c_str() );
  }

  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  if (!CoupledTo3D)
  {
    TimeUpdate();
  }

  // -------------------------------------------------------------------
  //  lift'n'drag forces, statistics time sample and output of solution
  //  and statistics
  // -------------------------------------------------------------------
  if (!CoupledTo3D)
  {
    Output(CoupledTo3D,CouplingTo3DParams);
  }

  // -------------------------------------------------------------------
  //                       update time step sizes
  // -------------------------------------------------------------------
  dtp_ = dta_;

} // RedAirwayImplicitTimeInt::OneStepTimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the one step time loop                          ismail 09/12|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::IntegrateStep(RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E   Solving Reduced Dimensional Airways    STEP = %4d/%4d \n", time_,maxtime_,dta_,step_,stepmax_);

  }

  // -------------------------------------------------------------------
  //                            SolverParams
  // -------------------------------------------------------------------
  if(params_.get<string> ("solver type") == "Nonlinear")
  {
    NonLin_Solve(CouplingTo3DParams);
    if (!myrank_)
      cout<<endl;
  }
  else if (params_.get<string> ("solver type") == "Linear")
  {
    Solve(CouplingTo3DParams);
    if (!myrank_)
      cout<<endl<<endl;
  }
  else
  {
    dserror("[%s] is not a defined solver",(params_.get<string> ("solver type")).c_str() );
  }

} // RedAirwayImplicitTimeInt::IntegrateStep

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::PrepareTimeStep()
{
  rhs_->PutScalar(0.0);
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

} //RedAirwayImplicitTimeInt::PrepareTimeStep


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | the nonlinear iterative solver for reduced  ismail             01/11 |
 | dimensional airway                                                   |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*
Some detials!!
*/
void AIRWAY::RedAirwayImplicitTimeInt::NonLin_Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  //--------------------------------------------------------------------
  //
  //--------------------------------------------------------------------
  double error_norm1 = 1.e7;
  double error_norm2 = 1.e7;

  //--------------------------------------------------------------------
  // Evaluate Lung volume
  //--------------------------------------------------------------------
  double lung_volume_np = 0.0;
  bool err = this->SumAllColElemVal(acini_e_volumenp_,lung_volume_np);
  if(err)
  {
    dserror("Error by summing all acinar volumes");
  }


  // print out the total lung volume
  if(!myrank_)
  {
    cout<<"time: "<<time_-dta_<<" LungVolume: "<<lung_volume_np<<endl;
  }

  for (int i =1; i<=maxiter_; i++)
  {
    //------------------------------------------------------------------
    // update the pressures of the previous time step
    //------------------------------------------------------------------
    p_nonlin_->Update(1.0,*pnp_     ,0.0);

    //------------------------------------------------------------------
    // Solve the reduced dimensional model
    //------------------------------------------------------------------
    this->Solve(CouplingTo3DParams);

    //------------------------------------------------------------------
    // Find the change of pressure between the last two iteration steps
    //------------------------------------------------------------------
    p_nonlin_->Update(1.0,*pnp_     ,-1.0);

    //------------------------------------------------------------------
    // Evaluate the L2 norm of the difference
    //------------------------------------------------------------------
    p_nonlin_->Norm2 (&error_norm1);
    //    error_norm1 /= sqrt(double(p_nonlin_->GlobalLength()));

    this->EvalResidual(CouplingTo3DParams);
    residual_->Norm2 (&error_norm2);
    //    cout<<"Proc: "<<myrank_<<" Norm: "<<error_norm<<" length: "<<p_nonlin_->GlobalLength()<<endl;

    //------------------------------------------------------------------
    // if L2 norm is smaller then tolerance then proceed
    //------------------------------------------------------------------
    if (!myrank_)
    {
      printf("iteration step %4d/%4d ",i,maxiter_);
      printf(" | ||P{%d}-P{%d}||_L2 = %10.3E\t\t|Qresidual|_2 = %10.3E\n",i-1,i,error_norm1,error_norm2);
    }
    if(error_norm1 <= non_lin_tol_)
      break;

  }
  if (!myrank_)
    printf("\n");
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | the solver for reduced dimensional airway               ismail 01/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*
Some detials!!
*/
void AIRWAY::RedAirwayImplicitTimeInt::Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  // time measurement: Airways
  if(!coupledTo3D_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("   + solving reduced dimensional airways");
  }

  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------

  // get cpu time
  //  const double tcpuele = Teuchos::Time::wallTime();
  {
    // time measurement: element
    if(!coupledTo3D_)
    {
      TEUCHOS_FUNC_TIME_MONITOR("      + element calls");
    }

    // set both system matrix and rhs vector to zero
    sysmat_->Zero();
    rhs_->PutScalar(0.0);


    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_sys_matrix_rhs");
    eleparams.set("time step size",dta_);

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("pn" ,pn_ );
    discret_->SetState("pnm",pnm_);

    eleparams.set("qin_np",qin_np_);
    eleparams.set("qin_n" ,qin_n_);
    eleparams.set("qin_nm",qin_nm_);

    eleparams.set("acinar_vn" ,acini_e_volumen_);
    eleparams.set("acinar_vnp",acini_e_volumenp_);
    //    discret_->SetState("acinar_vn" ,acini_e_volumen_);
    //    discret_->SetState("acinar_vnp",acini_e_volumenp_);

    eleparams.set("qout_np",qout_np_);
    eleparams.set("qout_n" ,qout_n_ );
    eleparams.set("qout_nm",qout_nm_ );

    //------------------------------------------------------------------
    // Evaluate Lung volumes
    //------------------------------------------------------------------
    double lung_volume_np = 0.0;
    bool err = this->SumAllColElemVal(acini_e_volumenp_,lung_volume_np);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }

    double lung_volume_n  = 0.0;
    err = this->SumAllColElemVal(acini_e_volumen_,lung_volume_n);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }

    double lung_volume_nm = 0.0;
    err = this->SumAllColElemVal(acini_e_volumenm_,lung_volume_nm);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }

    eleparams.set("lungVolume_np",lung_volume_np);
    eleparams.set("lungVolume_n" ,lung_volume_n);
    eleparams.set("lungVolume_nm",lung_volume_nm);

    sysmat_iad_->PutScalar(0.0);
    eleparams.set("sysmat_iad",sysmat_iad_);
//    eleparams.set("num_of_inter_acinar_linkers",num_of_inter_acinar_linkers_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    //    sysmat_->Complete();
    discret_->ClearState();

#if 0  // Exporting some values for debugging purposes

    {
      cout<<"----------------------- My SYSMAT IS ("<<myrank_<<"-----------------------"<<endl;
      RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
      if (A_debug != Teuchos::null)
      {
        // print to screen
        (A_debug->EpetraMatrix())->Print(cout);
      }
      cout<<"Map is: ("<<myrank_<<")"<<endl<<*(discret_->DofRowMap())<<endl;
      cout<<"---------------------------------------("<<myrank_<<"------------------------"<<endl;
    }
#endif

  }
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_sys_matrix_rhs_iad");
    discret_->SetState("sysmat_iad",sysmat_iad_);
    discret_->SetState("num_of_inter_acinar_linkers",num_of_inter_acinar_linkers_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
    discret_->ClearState();

#if 0 // Exporting some values for debugging purposes

    {
      cout<<"----------------------- My SYSMAT IS ("<<myrank_<<"-----------------------"<<endl;
      RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
      if (A_debug != Teuchos::null)
      {
        // print to screen
        (A_debug->EpetraMatrix())->Print(cout);
      }
      cout<<"Map is: ("<<myrank_<<")"<<endl<<*(discret_->DofRowMap())<<endl;
      cout<<"---------------------------------------("<<myrank_<<"------------------------"<<endl;
    }
#endif

  }
  // end time measurement for element


  // -------------------------------------------------------------------
  // Solve the boundary conditions
  // -------------------------------------------------------------------
  bcval_->PutScalar(0.0);
  dbctog_->PutScalar(0.0);
  // Solve terminal BCs
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","set_bc");

    // set vecotr values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("pn" ,pn_ );
    discret_->SetState("pnm",pnm_);

    eleparams.set("acinar_vn" ,acini_e_volumen_);
    eleparams.set("acinar_vnp",acini_e_volumenp_);

    //    discret_->SetState("qcnp",qcnp_);
    //    discret_->SetState("qcn" ,qcn_ );
    //    discret_->SetState("qcnm",qcnm_);

    eleparams.set("qin_np",qin_np_);
    eleparams.set("qin_n",qin_n_);
    //    discret_->SetState("qin_n" ,qin_n_ );
    //    discret_->SetState("qin_nm",qin_nm_);

    eleparams.set("qout_np",qout_np_);
    //    discret_->SetState("qout_n" ,qout_n_ );
    //    discret_->SetState("qout_nm",qout_nm_);

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);
    eleparams.set("bcval",bcval_);
    eleparams.set("dbctog",dbctog_);

    // get lung volume
    double lung_volume_np = 0.0;
    bool err = this->SumAllColElemVal(acini_e_volumenp_,lung_volume_np);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }
    eleparams.set("lungVolume_np",lung_volume_np);

    //    eleparams.set("abc",abc_);
    //    eleparams.set("rhs",rhs_);

    // Add the parameters to solve terminal BCs coupled to 3D fluid boundary
    eleparams.set("coupling with 3D fluid params",CouplingTo3DParams);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();
  }


  double norm_bc_tog = 0.0;
  rhs_->Norm1(&norm_bc_tog);

  // -------------------------------------------------------------------
  // Apply the BCs to the system matrix and rhs
  // -------------------------------------------------------------------
  {
    // time measurement: application of dbc
    if(!coupledTo3D_)
    {
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
    }

    LINALG::ApplyDirichlettoSystem(sysmat_,pnp_,rhs_,bcval_,dbctog_);
  }

#if 0 // Exporting some values for debugging purposes

    {
      cout<<"----------------------- My SYSMAT IS ("<<myrank_<<"-----------------------"<<endl;
      RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
      if (A_debug != Teuchos::null)
      {
        // print to screen
        (A_debug->EpetraMatrix())->Print(cout);
      }
      cout<<"Map is: ("<<myrank_<<")"<<endl<<*(discret_->DofRowMap())<<endl;
      cout<<"---------------------------------------("<<myrank_<<"------------------------"<<endl;
    }
    cout<<"rhs: "<<*rhs_<<endl;

#endif

  //-------solve for total new velocities and pressures
  // get cpu time
  const double tcpusolve = Teuchos::Time::wallTime();
  {
    // time measurement: solver
    if(!coupledTo3D_)
    {
      TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");
    }

#if 0 // Exporting some values for debugging purposes

    RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
    if (A_debug != Teuchos::null)
    {
      // print to screen
      (A_debug->EpetraMatrix())->Print(cout);
    }
    //#else
    cout<<"DOF row map"<<*(discret_->DofRowMap())<<endl;
    cout<<"bcval: "<<*bcval_<<endl;
    cout<<"bctog: "<<*dbctog_<<endl;
    cout<<"pnp: "<<*pnp_<<endl;
    cout<<"rhs: "<<*rhs_<<endl;

#endif
    // call solver
    solver_.Solve(sysmat_->EpetraOperator(),pnp_,rhs_,true,true);
  }

  // end time measurement for solver
  dtsolve_ = Teuchos::Time::wallTime() - tcpusolve;

  if (myrank_ == 0)
  //  printf("te=%f, ts=%4.7f |",dtele_, dtsolve_);
  printf("ts=%4.7f |", dtsolve_);
  //  cout << "te=" << dtele_ << ", ts=" << dtsolve_ <<" | ";

  // find the flow rates
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_flow_rates");

    // set solution type
    eleparams.set("solver type", params_.get<string> ("solver type"));

    // set vecotr values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("pn" ,pn_ );
    discret_->SetState("pnm",pnm_);
    //    discret_->SetState("acinar_vn" ,acini_e_volumen_);
    eleparams.set("acinar_vn" ,acini_e_volumen_);

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);
    eleparams.set("qin_nm",qin_nm_);
    eleparams.set("qout_nm",qout_nm_);
    eleparams.set("qin_n",qin_n_);
    eleparams.set("qout_n",qout_n_);
    eleparams.set("qin_np",qin_np_);
    eleparams.set("qout_np",qout_np_);
    eleparams.set("sysmat_iad",sysmat_iad_);


    //    acini_e_volumenp_->PutScalar(0.0);
    //    acini_e_volume_strain_->PutScalar(0.0);

    eleparams.set("acinar_vnp_strain",acini_e_volume_strain_);
    eleparams.set("acinar_vnp",acini_e_volumenp_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,Teuchos::null ,rhs_,Teuchos::null,Teuchos::null);
    discret_->ClearState();
  }

  if(coupledTo3D_)
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","get_coupled_values");

    // set vecotr values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);

    eleparams.set("qin_np",qin_np_);
    eleparams.set("qin_n",qin_n_);

    eleparams.set("qout_np",qout_np_);
    eleparams.set("qout_n" ,qout_n_ );

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);

    // Add the parameters to solve terminal BCs coupled to 3D fluid boundary
    eleparams.set("coupling with 3D fluid params",CouplingTo3DParams);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();
  }

} // RedAirwayImplicitTimeInt::Solve




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble            |
 | this function will be kept empty untill further use     ismail 01/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::AssembleMatAndRHS()
{

  dtele_    = 0.0;
  dtfilter_ = 0.0;
  // time measurement: element
  if(!coupledTo3D_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("      + element calls");
  }

  // get cpu time
  //  const double tcpu=Teuchos::Time::wallTime();

} // RedAirwayImplicitTimeInt::AssembleMatAndRHS




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build system matrix and rhs                              ismail 01/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> qael)
{

}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                                      |
 |  pnm_  =  pn_                                                        |
 |  pn_   =  pnp_                                                       |
 |                                                                      |
 |                                                          ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::TimeUpdate()
{


  // Volumetric Flow rate/Cross-sectional area of this step become most recent
  pnm_->Update(1.0,*pn_ ,0.0);
  pn_ ->Update(1.0,*pnp_,0.0);

  qin_nm_->Update(1.0,*qin_n_ ,0.0);
  qin_n_ ->Update(1.0,*qin_np_,0.0);

  qout_nm_->Update(1.0,*qout_n_ ,0.0);
  qout_n_ ->Update(1.0,*qout_np_,0.0);

  acini_e_volumenm_->Update(1.0,*acini_e_volumen_,0.0);
  acini_e_volumen_->Update(1.0,*acini_e_volumenp_,0.0);

  return;
}// RedAirwayImplicitTimeInt::TimeUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                       ismail 07/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Output(bool               CoupledTo3D,
                                              RCP<Teuchos::ParameterList> CouplingParams)
{
  int step      = 0;
  int upres     = 0;
  int uprestart = 0;
  double time_backup = 0.0;
  // -------------------------------------------------------------------
  // if coupled to 3D problem, then get the export information from
  // the 3D problem
  // -------------------------------------------------------------------
  if (CoupledTo3D)
  {
    step        = step_;
    upres       = upres_;
    uprestart   = uprestart_;
    time_backup = time_;
    step_      = CouplingParams->get<int>("step");
    upres_     = CouplingParams->get<int>("upres");
    uprestart_ = CouplingParams->get<int>("uprestart");
    time_      = CouplingParams->get<double>("time");
  }

  if (step_%upres_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // "pressure" vectors
    output_.WriteVector("pnm",pnm_);
    output_.WriteVector("pn",pn_);
    output_.WriteVector("pnp",pnp_);

    // "flow" vectors for capacitances
    //    output_.WriteVector("qcnm",qcnm_);
    //    output_.WriteVector("qcn",qcn_);
    //    output_.WriteVector("qcnp",qcnp_);

    // write domain decomposition for visualization
    //    output_.WriteElementData();

    // write the flow values
    LINALG::Export(*qin_nm_,*qexp_);
    output_.WriteVector("qin_nm",qexp_);
    LINALG::Export(*qin_n_ ,*qexp_);
    output_.WriteVector("qin_n" ,qexp_ );
    LINALG::Export(*qin_np_,*qexp_);
    output_.WriteVector("qin_np",qexp_);

    LINALG::Export(*qout_nm_,*qexp_);
    output_.WriteVector("qout_nm",qexp_);
    LINALG::Export(*qout_n_ ,*qexp_);
    output_.WriteVector("qout_n" ,qexp_);
    LINALG::Export(*qout_np_,*qexp_);
    output_.WriteVector("qout_np",qexp_);

#if 0
    // write the acinar values
    LINALG::Export(*acini_e_volumenm_,*qexp_);
    output_.WriteVector("acini_vnm",qexp_);
    LINALG::Export(*acini_e_volumen_,*qexp_);
    output_.WriteVector("acini_vn",qexp_);
    LINALG::Export(*acini_e_volumenp_,*qexp_);
    output_.WriteVector("acini_vnp",qexp_);
    LINALG::Export(*acini_e_volume_strain_,*qexp_);
    output_.WriteVector("acini_volumetric_strain",qexp_);
    LINALG::Export(*acini_e_volume0_,*qexp_);
    output_.WriteVector("acini_v0",qexp_);
#else
    {
      Epetra_Export exporter(acini_e_volumenm_->Map(),qexp_->Map());
      int err = qexp_->Export(*acini_e_volumenm_,exporter,Zero);
      if (err) dserror("Export using exporter returned err=%d",err);
      output_.WriteVector("acini_vnm",qexp_);
    }
    {
      Epetra_Export exporter(acini_e_volumen_->Map(),qexp_->Map());
      int err = qexp_->Export(*acini_e_volumen_,exporter,Zero);
      if (err) dserror("Export using exporter returned err=%d",err);
      output_.WriteVector("acini_vn",qexp_);
    }
    {
      Epetra_Export exporter(acini_e_volumenp_->Map(),qexp_->Map());
      int err = qexp_->Export(*acini_e_volumenp_,exporter,Zero);
      if (err) dserror("Export using exporter returned err=%d",err);
      output_.WriteVector("acini_vnp",qexp_);
    }
    {
      Epetra_Export exporter(acini_e_volume_strain_->Map(),qexp_->Map());
      int err = qexp_->Export(*acini_e_volume_strain_,exporter,Zero);
      if (err) dserror("Export using exporter returned err=%d",err);
      output_.WriteVector("acini_volumetric_strain",qexp_);
    }
    {
      Epetra_Export exporter(acini_e_volume0_->Map(),qexp_->Map());
      int err = qexp_->Export(*acini_e_volume0_,exporter,Zero);
      if (err) dserror("Export using exporter returned err=%d",err);
      output_.WriteVector("acini_v0",qexp_);
    }
#endif

    if (step_==upres_)
    {
      output_.WriteVector("NodeIDs",nodeIds_);
      output_.WriteVector("radii",radii_);
      LINALG::Export(*generations_,*qexp_);
      output_.WriteVector("generations",qexp_);
      LINALG::Export(*acini_bc_,*qexp_);
      output_.WriteVector("acin_bc",qexp_);
      output_.WriteElementData();
    }

    // write mesh in each restart step --- the elements are required since
    // they contain history variables (the time dependent subscales)
    //    output_.WriteMesh(step_,time_);

    if (CoupledTo3D)
    {
      output_.WriteDouble("Actual_RedD_step", step);
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // "pressure" vectors
    output_.WriteVector("pnm",pnm_);
    output_.WriteVector("pn",pn_);
    output_.WriteVector("pnp",pnp_);


    // write the flow values
    LINALG::Export(*qin_nm_,*qexp_);
    output_.WriteVector("qin_nm",qexp_);
    LINALG::Export(*qin_n_ ,*qexp_);
    output_.WriteVector("qin_n" ,qexp_ );
    LINALG::Export(*qin_np_,*qexp_);
    output_.WriteVector("qin_np",qexp_);
    //
    LINALG::Export(*qout_nm_,*qexp_);
    output_.WriteVector("qout_nm",qexp_);
    LINALG::Export(*qout_n_ ,*qexp_);
    output_.WriteVector("qout_n" ,qexp_);
    LINALG::Export(*qout_np_,*qexp_);
    output_.WriteVector("qout_np",qexp_);
    //
    LINALG::Export(*acini_e_volumenm_,*qexp_);
    output_.WriteVector("acini_vnm",qexp_);
    LINALG::Export(*acini_e_volumen_,*qexp_);
    output_.WriteVector("acini_vn",qexp_);
    LINALG::Export(*acini_e_volumenp_,*qexp_);
    output_.WriteVector("acini_vnp",qexp_);
    LINALG::Export(*acini_e_volume_strain_,*qexp_);
    output_.WriteVector("acini_volumetric_strain",qexp_);
    LINALG::Export(*acini_e_volume0_,*qexp_);
    output_.WriteVector("acini_v0",qexp_);



    // write mesh in each restart step --- the elements are required since
    // they contain history variables (the time dependent subscales)
    output_.WriteMesh(step_,time_);

    if (CoupledTo3D)
    {
      output_.WriteDouble("Actual_RedD_step", step);
    }
  }
  // -------------------------------------------------------------------
  // if coupled to 3D problem, then retrieve the old information of the
  // the reduced model problem
  // -------------------------------------------------------------------
  if (CoupledTo3D)
  {
    step_     = step;
    upres_    = upres;
    uprestart_= uprestart;
    time_     = time_backup;
  }
  return;
} // RedAirwayImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | ReadRestart (public)                                     ismail 01/10|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::ReadRestart(int step)
{

  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");

  if (coupledTo3D_)
  {
    step_ = reader.ReadInt("Actual_RedD_step");
  }
  else
  {
    step_ = reader.ReadInt("step");
  }

  reader.ReadVector(pnp_,"pnp");
  reader.ReadVector(pn_,"pn");
  reader.ReadVector(pnm_,"pnm");

  reader.ReadVector(qexp_,"acini_vnm");
  LINALG::Export(*qexp_,*acini_e_volumenm_);
  reader.ReadVector(qexp_,"acini_vn");
  LINALG::Export(*qexp_,*acini_e_volumen_);
  reader.ReadVector(qexp_,"acini_vnp");
  LINALG::Export(*qexp_,*acini_e_volumenp_);
  reader.ReadVector(qexp_,"acini_volumetric_strain");
  LINALG::Export(*qexp_,*acini_e_volume_strain_);
  reader.ReadVector(qexp_,"acini_v0");
  LINALG::Export(*qexp_,*acini_e_volume0_);

  reader.ReadVector(qexp_, "qin_nm");
  LINALG::Export(*qexp_,*qin_nm_);
  reader.ReadVector(qexp_ , "qin_n" );
  LINALG::Export(*qexp_,*qin_n_);
  reader.ReadVector(qexp_, "qin_np");
  LINALG::Export(*qexp_,*qin_np_);

  reader.ReadVector(qexp_, "qout_nm");
  LINALG::Export(*qexp_,*qout_nm_);
  reader.ReadVector(qexp_ , "qout_n" );
  LINALG::Export(*qexp_,*qout_n_);
  reader.ReadVector(qexp_, "qout_np");
  LINALG::Export(*qexp_,*qout_np_);

  // read the previously written elements including the history data
  //reader.ReadMesh(step_);

}//RedAirwayImplicitTimeInt::ReadRestart


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#if 0
Teuchos::RCP<DRT::ResultTest>
AIRWAY::RedAirwayImplicitTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new AIRWAY::RedAirwayResultTest(*this));
}
#endif

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayImplicitTimeInt::EvalResidual( Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  residual_->PutScalar(0.0);
  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------
  {
    // set both system matrix and rhs vector to zero
    sysmat_->Zero();
    rhs_->PutScalar(0.0);

    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_sys_matrix_rhs");
    eleparams.set("time step size",dta_);

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("pn" ,pn_ );
    discret_->SetState("pnm",pnm_);


    eleparams.set("acinar_vn" ,acini_e_volumen_);
    eleparams.set("acinar_vnp",acini_e_volumenp_);
    //    discret_->SetState("acinar_vn" ,acini_e_volumen_);
    //    discret_->SetState("acinar_vnp",acini_e_volumenp_);

    eleparams.set("qin_np",qin_np_);
    eleparams.set("qin_n" ,qin_n_);
    eleparams.set("qin_nm",qin_nm_);


    eleparams.set("qout_np",qout_np_);
    eleparams.set("qout_n" ,qout_n_ );
    eleparams.set("qout_nm",qout_nm_ );

    // get lung volume
    double lung_volume_np = 0.0;
    bool err = this->SumAllColElemVal(acini_e_volumenp_,lung_volume_np);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }

    double lung_volume_n  = 0.0;
    err = this->SumAllColElemVal(acini_e_volumen_,lung_volume_n);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }

    double lung_volume_nm = 0.0;
    err = this->SumAllColElemVal(acini_e_volumenm_,lung_volume_nm);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }



    eleparams.set("lungVolume_np",lung_volume_np);
    eleparams.set("lungVolume_n" ,lung_volume_n);
    eleparams.set("lungVolume_nm",lung_volume_nm);

    eleparams.set("sysmat_iad",sysmat_iad_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    discret_->ClearState();
  }
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_sys_matrix_rhs_iad");
    discret_->SetState("sysmat_iad",sysmat_iad_);
    discret_->SetState("num_of_inter_acinar_linkers",num_of_inter_acinar_linkers_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
    discret_->ClearState();

#if 0  // Exporting some values for debugging purposes

    {
      cout<<"----------------------- My SYSMAT IS ("<<myrank_<<"-----------------------"<<endl;
      RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
      if (A_debug != Teuchos::null)
      {
        // print to screen
        (A_debug->EpetraMatrix())->Print(cout);
      }
      cout<<"Map is: ("<<myrank_<<")"<<endl<<*(discret_->DofRowMap())<<endl;
      cout<<"---------------------------------------("<<myrank_<<"------------------------"<<endl;
    }
#endif

  }
  // -------------------------------------------------------------------
  // Solve the boundary conditions
  // -------------------------------------------------------------------
  bcval_->PutScalar(0.0);
  dbctog_->PutScalar(0.0);
  // Solve terminal BCs
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","set_bc");

    // set vecotr values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("pn" ,pn_ );
    discret_->SetState("pnm",pnm_);
    //    discret_->SetState("qcnp",qcnp_);
    //    discret_->SetState("qcn" ,qcn_ );
    //    discret_->SetState("qcnm",qcnm_);

    eleparams.set("acinar_vn" ,acini_e_volumen_);
    eleparams.set("acinar_vnp",acini_e_volumenp_);

    eleparams.set("qin_np",qin_np_);
    eleparams.set("qin_n",qin_n_);

    eleparams.set("qout_np",qout_np_);

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);
    eleparams.set("bcval",bcval_);
    eleparams.set("dbctog",dbctog_);

    // Add the parameters to solve terminal BCs coupled to 3D fluid boundary
    eleparams.set("coupling with 3D fluid params",CouplingTo3DParams);

    // get lung volume
    double lung_volume_np = 0.0;
    bool err = this->SumAllColElemVal(acini_e_volumenp_,lung_volume_np);
    if(err)
    {
      dserror("Error by summing all acinar volumes");
    }
    eleparams.set("lungVolume_np",lung_volume_np);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();
  }

  // -------------------------------------------------------------------
  // Apply the BCs to the system matrix and rhs
  // -------------------------------------------------------------------
  {
    LINALG::ApplyDirichlettoSystem(sysmat_,pnp_,rhs_,bcval_,dbctog_);
  }

  // -------------------------------------------------------------------
  // Evaluate Residual
  // -------------------------------------------------------------------
  sysmat_->Multiply(false, *pnp_, *residual_);
  residual_-> Update(-1.0,*rhs_ ,1.0);
}//EvalResidual

/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                 ismail 01/10|
 *----------------------------------------------------------------------*/
AIRWAY::RedAirwayImplicitTimeInt::~RedAirwayImplicitTimeInt()
{
  return;
}


void AIRWAY::RedAirwayImplicitTimeInt::SetAirwayFluxFromTissue(Teuchos::RCP<Epetra_Vector> coupflux)
{
  const Epetra_BlockMap& condmap = coupflux->Map();

  for (int i=0; i<condmap.NumMyElements(); ++i)
  {
    int condID = condmap.GID(i);
    DRT::Condition* cond = coupcond_[condID];
    std::vector<double> newval(1,0.0);
    newval[0] = (*coupflux)[i];
    cond->Add("val",newval);
  }
}


void AIRWAY::RedAirwayImplicitTimeInt::SetupForCoupling()
{
  std::vector<DRT::Condition*> nodecond;
  discret_->GetCondition("RedAirwayPrescribedCond",nodecond);
  unsigned int numnodecond = nodecond.size();
  if (numnodecond == 0) dserror("no redairway prescribed conditions");

  std::vector<int> tmp;
  for (unsigned int i = 0; i < numnodecond; ++i)
  {
    DRT::Condition* actcond = nodecond[i];
    if (actcond->Type() == DRT::Condition::RedAirwayNodeTissue)
    {
      int condID = actcond->GetInt("coupling id");
      coupcond_[condID] = actcond;
      tmp.push_back(condID);
      pres_[condID] = 0.0;
    }
  }
  unsigned int numcond = tmp.size();
  if (numcond == 0) dserror("no coupling conditions found");
  coupmap_ = Teuchos::rcp(new Epetra_Map(tmp.size(),tmp.size(),&tmp[0],0,discret_->Comm()));
}


void AIRWAY::RedAirwayImplicitTimeInt::ExtractPressure(Teuchos::RCP<Epetra_Vector> couppres)
{

  for (int i=0; i<coupmap_->NumMyElements(); i++)
  {
    int condgid = coupmap_->GID(i);
    DRT::Condition * cond = coupcond_[condgid];
    const std::vector<int>* nodes = cond->Nodes();
    if (nodes->size()!=1)
      dserror("Too many nodes on coupling with tissue condition ID=[%d]\n",condgid);

    int gid = (*nodes)[0];
    double pressure = 0.0;
    if (discret_->HaveGlobalNode(gid))
    {
      DRT::Node* node = discret_->gNode(gid);
      if (myrank_==node->Owner())
      {
        int giddof = discret_->Dof(node, 0);
        int liddof = pnp_->Map().LID(giddof);
        pressure = (*pnp_)[liddof];
      }
    }
    double parpres = 0.;
    discret_->Comm().SumAll(&pressure,&parpres,1);
    (*couppres)[i] = parpres;
  }
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Sum all ColElement values                                            |
 |                                                          ismail 11/12|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
bool AIRWAY::RedAirwayImplicitTimeInt::SumAllColElemVal(Teuchos::RCP<Epetra_Vector> vec, double & sum)
{

  // Check if the vector is a ColElement vector
  const Epetra_Map* elementcolmap  = discret_->ElementColMap();
  if (!vec->Map().SameAs(*elementcolmap))
  {
    return true;
  }

  // define the volume of acini on the current processor
  double local_sum = 0.0;

  // Since the acinar_volume vector is a ColMap, we first need to export
  // it to a RowMap and eliminate the ghosted values
  {
    // define epetra exporter
    Epetra_Export exporter(vec->Map(),qexp_->Map());
    // export from ColMap to RowMap
    int err = qexp_->Export(*vec,exporter,Zero);
    if (err) dserror("Export using exporter returned err=%d",err);
  }
  // Get the mean acinar volume on the current processor
  qexp_->MeanValue(&local_sum);

  // Multiply the mean by the size of the vector to get the total
  // acinar volume on the current processor
  local_sum *= double(qexp_->MyLength());

  // Get the total volume of Acini on all processors
  sum = 0.0;
  discret_->Comm().SumAll(&local_sum,&sum,1);

  // return all is fine
  return false;
}

