
/*!----------------------------------------------------------------------
\file airwayimplicitintegration.cpp
\brief Control routine for artery solvers,

     including solvers based on

     o two-step Taylor-Galerking scheme


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "airwayimplicitintegration.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 01/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

AIRWAY::RedAirwayImplicitTimeInt::RedAirwayImplicitTimeInt(RCP<DRT::Discretization>  actdis,
                                                           LINALG::Solver  &         solver,
                                                           ParameterList&            params,
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
  upres_(params.get("write solution every", -1))
{

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = dtp_*double(stepmax_);


  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap  = discret_->DofRowMap();

  //  const Epetra_Map* dofcolmap  = discret_->DofColMap();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global node numbering
  // -------------------------------------------------------------------
  // const Epetra_Map* noderowmap = discret_->NodeRowMap();


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
  // Volumetric flow rate at time n+1, n and n-1
  qnp_          = LINALG::CreateVector(*dofrowmap,true);
  qn_           = LINALG::CreateVector(*dofrowmap,true);
  qnm_          = LINALG::CreateVector(*dofrowmap,true);

  // Pressures at time n+1, n and n-1
  pnp_          = LINALG::CreateVector(*dofrowmap,true);
  pn_           = LINALG::CreateVector(*dofrowmap,true);
  pnm_          = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  // This part might be optimized later
  bcval_   = LINALG::CreateVector(*dofrowmap,true);
  dbctog_  = LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------


  // right hand side vector and right hand side corrector
  rhs_     = LINALG::CreateVector(*dofrowmap,true);

  // ---------------------------------------------------------------------------------------
  // Initialize all the arteries' cross-sectional areas to the initial crossectional area Ao
  // and the volumetric flow rate to 0
  // ---------------------------------------------------------------------------------------
  ParameterList eleparams;
  discret_->ClearState();
  discret_->SetState("qnp",qnp_);
  discret_->SetState("pnp",pnp_);

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lColElement(nele);

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    //vector<int> lmowner;
    RCP<vector<int> > lmowner = rcp(new vector<int>);
    ele->LocationVector(*discret_,lm,*lmowner);

    // loop all nodes of this element, add values to the global vectors
    eleparams.set("p0",pnp_);
    eleparams.set("q0",qnp_);
    eleparams.set("lmowner",lmowner);
    eleparams.set("action","get_initial_state");
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  }


#if 0
  cout<<"|**************************************************************************|"<<endl;
  cout<<"|******************** The Initialize Vector qanp is ***********************|"<<endl;
  cout<<"|**************************************************************************|"<<endl;
  cout<<*pnp_<<endl;
  cout<<"|**************************************************************************|"<<endl;
#endif

} // ArtNetExplicitTimeInt::ArtNetExplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration.                                          |
 |                                                          ismail 12/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Integrate()
{
  RCP<ParameterList> param;
  Integrate(false, param);
} // ArtNetExplicitTimeInt::Integrate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration.                                          |
 |                                                          ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::Integrate(bool CoupledTo3D,
                                                 RCP<ParameterList> CouplingParams)
{
  if (CoupledTo3D && CouplingParams.get() == NULL)
  {
    dserror("Coupling parameter list is not allowed to be empty, If a 3-D/reduced-D coupling is defined\n");
  }

  TimeLoop(CoupledTo3D,CouplingParams);

  // print the results of time measurements
    TimeMonitor::summarize();

  return;
} // RedAirwayImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                   ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::TimeLoop(bool CoupledTo3D,
                                                RCP<ParameterList> CouplingTo3DParams)
{

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");
  
  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();
    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      printf("TIME: %11.4E/%11.4E  DT = %11.4E   Solving Artery    STEP = %4d/%4d \n",
              time_,maxtime_,dta_,step_,stepmax_);
    }

    Solve(CouplingTo3DParams);

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();

    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    Output();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

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
 | setup the variables to do a new time step                ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | the solver for artery                                   ismail 06/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*
Some detials!!
*/
void AIRWAY::RedAirwayImplicitTimeInt::Solve(Teuchos::RCP<ParameterList> CouplingTo3DParams)
{
  // time measurement: Artery
  TEUCHOS_FUNC_TIME_MONITOR("   + solving artery");


  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------

  // get cpu time
  //  const double tcpuele = ds_cputime();

  {
    // time measurement: element
    TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

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
    discret_->SetState("qnp",qnp_);
    discret_->SetState("pnp",pnp_);


    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
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
    discret_->SetState("qnp",qnp_);

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);
    eleparams.set("bcval",bcval_);
    eleparams.set("dbctog",dbctog_);
    eleparams.set("rhs",rhs_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
  }


  // -------------------------------------------------------------------
  // Apply the BCs to the system matrix and rhs
  // -------------------------------------------------------------------
  {
    // time measurement: application of dbc
    TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");

#if 0
    cout<<"Boundary values are: "<<endl<<*bcval_<<endl;
    cout<<"Boundary toggels are: "<<endl<<*dbctog_<<endl;
#endif

    LINALG::ApplyDirichlettoSystem(sysmat_,pnp_,rhs_,bcval_,dbctog_);
  }

  //-------solve for total new velocities and pressures
  // get cpu time
  const double tcpusolve = ds_cputime();
  {
    // time measurement: solver
    TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");

#if 0  // Exporting some values for debugging purposes

    RCP<LINALG::SparseMatrix> A_debug = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
    if (A_debug != Teuchos::null)
    {
      // print to screen
      (A_debug->EpetraMatrix())->Print(cout);
    }

    cout<<"pnp: "<<*pnp_<<endl;
    cout<<"rhs: "<<*rhs_<<endl;

#endif 

    // call solver
    solver_.Solve(sysmat_->EpetraOperator(),pnp_,rhs_,true,true);
  }


  // end time measurement for solver
  dtsolve_ = ds_cputime() - tcpusolve;

  if (myrank_ == 0)
    cout << "te=" << dtele_ << ", ts=" << dtsolve_ << "\n\n" ;

  // find the flow rate qnp
  #if 0
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calculate_flow_rate");

    // set vecotr values needed by elements
    discret_->ClearState();
    discret_->SetState("pnp",pnp_);
    discret_->SetState("qnp",qnp_);

    eleparams.set("time step size",dta_);
    eleparams.set("total time",time_);

    // call standard loop over all elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
  }
  #endif

} // ArtNetExplicitTimeInt:Solve




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble            |
 | this function will be kept empty untill further use     ismail 07/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::AssembleMatAndRHS()
{

  dtele_    = 0.0;
  dtfilter_ = 0.0;
  // time measurement: element
  TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

  // get cpu time
  //  const double tcpu=ds_cputime();

} // RedAirwayImplicitTimeInt::AssembleMatAndRHS




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build system matrix and rhs                              ismail 06/09|
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
 |  qnm_   =  qn_                                                       |
 |  arean_ = areap_                                                     |
 |                                                                      |
 |                                                          ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void AIRWAY::RedAirwayImplicitTimeInt::TimeUpdate()
{


  // Volumetric Flow rate/Cross-sectional area of this step become most recent
  qnm_->Update(1.0,*qn_ ,0.0);
  qn_ ->Update(1.0,*qnp_,0.0);

  pnm_->Update(1.0,*pn_ ,0.0);
  pn_ ->Update(1.0,*pnp_,0.0);

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
void AIRWAY::RedAirwayImplicitTimeInt::Output()
{

  if (step_%upres_ == 0)
  {
    // step number and time
    //    cout<<"My output is:"<<output_<<endl;
    output_.NewStep(step_,time_);

    // "volumetric flow rate/cross-sectional area" vector
    output_.WriteVector("qnp",qnp_);
    output_.WriteVector("pnp",pnp_);

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0)
    {
      // also write impedance bc information if required
      // Note: this method acts only if there is an impedance BC
      // impedancebc_->WriteRestart(output_);
    }
    //#endif

  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // "volumetric flow rate/cross-sectional area" vector
    output_.WriteVector("qnp",qnp_);
    output_.WriteVector("pnp",pnp_);
        
    // also write impedance bc information if required
    // Note: this method acts only if there is an impedance BC
    // impedancebc_->WriteRestart(output_);

  }

  return;
} // RedAirwayImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | ReadRestart (public)                                     ismail 07/09|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>// 
void AIRWAY::RedAirwayImplicitTimeInt::ReadRestart(int step)
{

  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(qnp_,"qnp");
  reader.ReadVector(pnp_,"pnp");

  // also read impedance bc information if required
  // Note: this method acts only if there is an impedance BC
  // impedancebc_->ReadRestart(reader);

}


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                 ismail 01/09|
 *----------------------------------------------------------------------*/
AIRWAY::RedAirwayImplicitTimeInt::~RedAirwayImplicitTimeInt()
{
  return;
}

#endif

