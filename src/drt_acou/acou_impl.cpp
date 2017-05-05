/*!----------------------------------------------------------------------
\file acou_impl.cpp
\brief Main control routine for acoustic simulations

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "acou_impl.H"
#include "acou_ele_action.H"
#include "acou_ele.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_io/io.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouImplicitTimeInt::AcouImplicitTimeInt(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  AcouTimeInt(actdis,solver,params,output),
  sourcefuncno_   ((params_->get<int>("SOURCETERMFUNCNO"))-1),
  dtele_          (0.0),
  dtsolve_        (0.0),
  acouopt_        (params_->get<bool>("acouopt")),
  outputcount_    (0),
  errormaps_      (DRT::INPUT::IntegralValue<bool>(*params_,"ERRORMAPS")),
  padapttol_      (params_->get<double>("P_ADAPT_TOL")),
  calcerr_        (false),
  allelesequal_   (DRT::INPUT::IntegralValue<bool>(*params_,"ALLELESEQUAL"))
{
  // some security checks
  if(dtp_ == 0.0)
    dserror("Can't work with time step size == 0.0");
  if(padaptivity_==true && errormaps_==false)
    dserror("If you want to do p-adaptivity, you also have to set the flag ERRORMAPS to Yes");
  if(padaptivity_==true && discret_->Comm().NumProc()>1)
    dserror("p-adaptivity does not yet work in parallel!"); // TODO
  if(acouopt_ == true && uprestart_<=0)
    dserror("for optimization with acoustical parameters, RESTARTEVRY must be >=1");

  // service
  if(params_->get<int>("CALCERRORFUNCNO")>0)
    calcerr_ = true;
  if(acouopt_&&adjoint_)
    outputcount_ = params_->get<int>("outputcount");

  // get the dof map
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector containing element based error values
  if(errormaps_)
    error_ = LINALG::CreateVector(*(discret_->ElementRowMap()),true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*dofrowmap,true);
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise

  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);
  }

  // print user information which might not be known by everyone
  if (errormaps_ && !myrank_ )
    std::cout<<"Local postprocessing is only effective when temporal accuracy is of order k+2. Did you choose your time integrator accordingly?"<<std::endl;

  // create system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  sysmat_->Zero();

  // Vector used for solution process
  residual_      = LINALG::CreateVector(*dofrowmap,true);

  output_->WriteMesh(0,0.0);
} // AcouImplicitTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouImplicitTimeInt::~AcouImplicitTimeInt()
{}

/*----------------------------------------------------------------------*
 |  Initialization of algorithm to zero (public)         schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::SetInitialZeroField()
{
  velnp_->PutScalar(0.0);

  ACOU::AcouTimeInt::SetInitialZeroField();

  return;
} // SetInitialZeroField()

/*----------------------------------------------------------------------*
 |  Initialization of algorithm by given function (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::SetInitialField(int startfuncno)
{
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;

  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::project_field);
  initParams.set<int>("funct",startfuncno);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  initParams.set<bool>("padaptivity",padaptivity_);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);

  DRT::Element::LocationArray la(2);
  int err = 0;
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    elevec1.Scale(0.0);elevec2.Scale(0.0);
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);

    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1,ele))
      elevec2.Shape(discret_->NumDof(1,ele), 1);

    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);
    // now fill the ele vector into the discretization
    for (unsigned int i=0; i<la[0].lm_.size(); ++i)
      la[0].lm_[i] = discret_->DofRowMap(0)->LID(la[0].lm_[i]);

    err += velnp_->ReplaceMyValues(la[0].lm_.size(),elevec1.A(),&la[0].lm_[0]);
  }
  //if(err!=0) dserror("Could not replace all values");

  return;
} // SetInitialField

/*----------------------------------------------------------------------*
 |  Time loop (public)                                   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::Integrate(Teuchos::RCP<Epetra_MultiVector> history)
{
  // time measurement: integration
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::Integrate");

  // if necessary, write a monitor file
  InitMonitorFile();

  // output of initial field (given by function for purely acoustic simulation or given by optics for PAT simulation)
  Output(history);

  // evaluate error
  EvaluateErrorComparedToAnalyticalSol();

  // intermediate integrate
  IntermediateIntegrate();

  // call elements to calculate system matrix/rhs and assemble
  AssembleMatAndRHS();

  // apply Dirichlet boundary conditions to system of equations
  ApplyDirichletToSystem();

  // time loop
  while (step_<stepmax_ and time_<maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    // output to screen
    OutputToScreen();

    // solve
    Solve();

    // p-adaptivity
    UpdatePolyAndState();

    // output of solution
    Output(history);

    // evaluate error
    EvaluateErrorComparedToAnalyticalSol();

    // intermediate integrate
    IntermediateIntegrate();
  } // while (step_<stepmax_ and time_<maxtime_)

  if (!myrank_) printf("\n");

  return;
} // Integrate

/*----------------------------------------------------------------------*
 |  Solve the system for trace and then interior field   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::Solve()
{
  // solve linear equation and timing
  const double tcpusolve=Teuchos::Time::wallTime();
  solver_->Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,false,Teuchos::null);
  dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;

  // update interior variables
  UpdateInteriorVariablesAndAssemebleRHS();
  ApplyDirichletToSystem();

  return;
} // Solve

/*----------------------------------------------------------------------*
 |  Dirichlet function (public)                    schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::ApplyDirichletToSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
  Teuchos::ParameterList params;
  params.set<double>("total time",time_);
  discret_->EvaluateDirichlet(params,zeros_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));
  return;
} // ApplyDirichletToSystem

/*----------------------------------------------------------------------*
 |  Calculate system matrix (public)                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::AssembleMatAndRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::AssembleMatAndRHS");

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual and sysmat
  residual_->Scale(0.0);
  sysmat_->Zero();

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret_->ClearState(true);
  if(!padaptivity_)
    discret_->SetState("trace",velnp_);

  // set time step size
  eleparams.set<double>("dt",dtp_);

  // call standard loop over elements
  bool resonly = false;// !(!bool(step_-1) || !bool(step_-restart_-1));

  // set information needed by the elements
  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<bool>("resonly",resonly);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<int>("useacouoptvecs",-1);
  eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<double>("time",time_);
  eleparams.set<double>("timep",time_+dtp_);
  eleparams.set<int>("step",step_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
  discret_->ClearState(true);

  if(!resonly)
  {
    // absorbing boundary conditions
    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_abc);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }
  if(adjoint_ && phys_==INPAR::ACOU::acou_lossless) // only needed for fluid, since the source term for the solid is calculated directly in the update routine
  {
    std::string condname = "PressureMonitor";
    std::vector<DRT::Condition*> pressuremon;
    discret_->GetCondition(condname,pressuremon);
    if(pressuremon.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_pressuremon);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }
  sysmat_->Complete();

  // residual_->Print(std::cout);
  // Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_)->EpetraMatrix()->Print(std::cout);

  return;
} // AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | Update interior field and calculate residual (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::UpdateInteriorVariablesAndAssemebleRHS()
{
  dtele_ = 0.0;

  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::UpdateInteriorVariablesAndAssemebleRHS");

  // get cpu time
  const double tcpu=Teuchos::Time::wallTime();

  Teuchos::RCP<std::vector<double> > elevals;
   if(errormaps_)
     elevals = Teuchos::rcp(new std::vector<double>(discret_->NumGlobalElements(),0.0));

  // create parameterlist
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<double>("dt",dtp_);
  eleparams.set<double>("time",time_);
  eleparams.set<double>("timep",time_+dtp_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<double>("padaptivitytol",padapttol_);
  eleparams.set<int>("useacouoptvecs",-1);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  eleparams.set<bool>("allelesequal",allelesequal_);
  eleparams.set<Teuchos::RCP<std::vector<double> > >("elevals",elevals);
  eleparams.set<int>("action",ACOU::update_secondary_solution_and_calc_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  bool calcgrad = acouopt_&&adjoint_;
  if(step_==stepmax_+1||step_==int(maxtime_/dtp_)+1) calcgrad=false;
  eleparams.set<bool>("calculategradient",calcgrad);
  int offset = uprestart_-step_%uprestart_;
  if(offset==uprestart_) offset = 0;
  eleparams.set<int>("calcgradoffset",offset);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);
  eleparams.set<bool>("resonly",true);

  residual_->Scale(0.0);
  discret_->SetState("trace",velnp_);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // update the error vector
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceMyValue(el,0,localvals[error_->Map().GID(el)]);
  }

  discret_->ClearState(true);

  // calculate source term for adjoint simulation
  if(adjoint_ && phys_==INPAR::ACOU::acou_lossless) // only needed for fluid, since the source term for the solid is calculated directly in the update routine
  {
    std::string condname = "PressureMonitor";
    std::vector<DRT::Condition*> pressuremon;
    discret_->GetCondition(condname,pressuremon);
    if(pressuremon.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_pressuremon);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpu;

  return;
} // UpdateInteriorVariablesAndAssemebleRHS


/*----------------------------------------------------------------------*
 | P-Adaptivity                                          schoeder 07/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::UpdatePolyAndState()
{
  /* This function serves to supply all required steps for p-adaptivity. What we need
   * is the following:
   * 1.) Do the local postprocessing, calculate delta_k (this is the amount the polynomial shape function needs to change)
   * 2.) Several things
   *     - Update the degree
   *     - Map / project the values
   *     - Rebuild the distributed vectors
   *     - Fill the distributed vectors
   * 3.) Do the next time step
   * UpdateInteriorVariables and ComputeResidual have to be separated (or not, if the element are samrt)
   */

  // 1.) can be fully done by the elements: if p-adaptivity is desired, then the elements
  //     do not store the error in the error vector, but the delta_k!
  if(!padaptivity_) return;

  for(int i=0; i<discret_->NumMyColElements(); ++i)
    dynamic_cast<DRT::ELEMENTS::Acou*>(discret_->lColElement(i))->SetDegree(int(error_->operator [](i)));

  // actually we don't want an entire FillComplete call, since nodes and elements remain the same
  // we only want the face and internal dofs and the faces are rebuild
  discret_->BuildFaces();
  discret_->BuildFaceRowMap();
  discret_->BuildFaceColMap();
  discret_->AssignDegreesOfFreedom(0);

  // update maps for global vectors
  velnp_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  residual_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  sysmat_ = Teuchos::null;
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),108,false,true));

  // now we have to call the calculation of the residual, because we skipped it in
  // UpdateInteriorVariablesAndComputeResidual
  AssembleMatAndRHS();

  return;
}
/*----------------------------------------------------------------------*
 |  Intermediate forward integration in adjoint run      schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::IntermediateIntegrate()
{
  // if this is no adjoint run or no acouopt_, return
  if(acouopt_==false || adjoint_ == false)
   return;

  // if this is no "starting step", return
  if(step_%uprestart_!=0)
    return;

  // save the trace field of the adjoint integration
  Teuchos::RCP<Epetra_Vector> adjointvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(0))));
  adjointvelnp->Update(1.0,*velnp_,0.0);

  // read the restart
  double time = 0.0;
  int restartstep = -1;
  {
    // determine "start" for the intermediate integration
    restartstep = (stepmax_<(maxtime_/dtp_)) ? stepmax_ : (maxtime_/dtp_);
    restartstep -= step_ + uprestart_;
    if(restartstep < 0) return; // happens in the last adjoint step

    // prepare and do the restart read
    std::stringstream name;
    name<< params_->get<std::string>("name");
    name<<"_invforward_acou_run_"<<outputcount_;
    Teuchos::RCP<IO::InputControl> inputfile = Teuchos::rcp(new IO::InputControl(name.str(), discret_->Comm()));
    IO::DiscretizationReader reader(discret_,inputfile,restartstep);
    time = reader.ReadDouble("time");

    // read the trace field
    reader.ReadVector(velnp_,"velnps");

    // read the internal field
    Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
    reader.ReadVector(intvelnp,"intvelnp");

    // bring the internal field to the elements BUT DO NOT OVERWRITE THEIR FIELD VECTORS
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",ACOU::ele_init_from_acouoptrestart);
    eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
    eleparams.set<bool>("padaptivity",padaptivity_);
    eleparams.set<int>("length",uprestart_+1); // +1 is for the last step
    discret_->SetState(1,"intvelnp",intvelnp);
    discret_->SetState("trace",velnp_);
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
    discret_->ClearState(true);
  }

  // perform the time steps!
  {
    // first step with setup of system matrix
    // reset residual
    residual_->Scale(0.0);

    // delete all entries of the system matrix
    sysmat_->Zero();

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;
    eleparams.set<double>("dt",dtp_);
    eleparams.set<int>("sourcefuncno",sourcefuncno_);
    eleparams.set<bool>("resonly",false);
    eleparams.set<bool>("padaptivity",padaptivity_);
    eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
    eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
    eleparams.set<bool>("adjoint",false); // this is forward integration!
    eleparams.set<int>("useacouoptvecs",0); // elements have to use the different vectors!!!
    eleparams.set<double>("time",time);
    eleparams.set<double>("timep",time+dtp_);
    eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

    // evaluate system matrix and residual
    discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
    discret_->ClearState(true);

    // add contribution to sysmat from absorbing boundary conditions
    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_abc);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }

    // sysmat is ready to go
    sysmat_->Complete();

    // apply the Dirichlet conditions (no function call due to different time)
    {
      Teuchos::ParameterList params;
      params.set<double>("total time",time);
      discret_->EvaluateDirichlet(params,zeros_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
      LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));
      zeros_->PutScalar(0.0);
    }

    // do the "time loop"
    int maxstep = restartstep + uprestart_;
    int minstep = restartstep;
    while (restartstep<maxstep)
    {
      // increment time and step
      restartstep++;
      time += dtp_;

      // output to screen
      if(!myrank_) printf(">");

      // solve
      solver_->Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,false,Teuchos::null);

      // update interior variables and compute residual
      {
        // parameters
        Teuchos::ParameterList eleparams;
        eleparams.set<int>("sourcefuncno",sourcefuncno_);
        eleparams.set<double>("dt",dtp_);
        eleparams.set<double>("time",time);
        eleparams.set<double>("timep",time+dtp_);
        eleparams.set<bool>("adjoint",false);
        eleparams.set<int>("useacouoptvecs",restartstep-minstep-1);
        eleparams.set<bool>("errormaps",errormaps_);
        eleparams.set<bool>("padaptivity",padaptivity_);
        eleparams.set<double>("padaptivitytol",padapttol_);
        eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
        eleparams.set<bool>("allelesequal",allelesequal_);
        eleparams.set<int>("action",ACOU::update_secondary_solution_and_calc_residual);
        eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
        //eleparams.set<int>("step",restartstep);
        eleparams.set<bool>("resonly",true);
        eleparams.set<bool>("calculategradient",false);

        // reset residual
        residual_->Scale(0.0);

        // evaluate
        discret_->SetState("trace",velnp_);
        discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
      }

      // apply dirichlet
      {
        Teuchos::ParameterList params;
        params.set<double>("total time",time);
        discret_->EvaluateDirichlet(params,zeros_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
        LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));
        zeros_->PutScalar(0.0);
      }
    } // while (restartstep<restartstep+uprestart_)
  }

  // now the elements hold the correct state vectors
  // act if nothing ever happened!
  velnp_->Update(1.0,*adjointvelnp,0.0);
  AssembleMatAndRHS();
  ApplyDirichletToSystem();

 return;
}

/*----------------------------------------------------------------------*
 |  Output time step information (public)                schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::OutputToScreen()
{
  // output to screen
  if (!myrank_)
  {
    if(invana_)
    {
      if(adjoint_)
        printf("<");
      else
        printf(">");
    }
    else
      printf("TIME: %11.4E/%11.4E  DT = %11.4E %s STEP = %4d/%4d, ts=%10.3E, te=%10.3E \n",time_,maxtime_,dtp_,Name().c_str(),step_,stepmax_,dtsolve_,dtele_);
  }

  return;
} // OutputToScreen

/*----------------------------------------------------------------------*
 |  Calculate node based values (public)                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec)
{
  // output of solution
  Teuchos::RCP<Epetra_Vector> pressure;

  // get the node vector
  {
    const Epetra_Map* nodemap = discret_->NodeRowMap();
    pressure.reset(new Epetra_Vector(*nodemap));
  }

  // call element routine for interpolate HDG to elements
  Teuchos::ParameterList params;
  params.set<int>("action",ACOU::interpolate_hdg_to_node);
  params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  params.set<bool>("padaptivity",false);
  params.set<int>("useacouoptvecs",-1);
  params.set<bool>("writestress",writestress_);

  DRT::Element::LocationArray la(2);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  Epetra_SerialDenseVector interpolVec;
  std::vector<unsigned char> touchCount(pressure->MyLength());
  pressure->PutScalar(0.);

  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);
    if (interpolVec.M() == 0)
    {
      if(phys_==INPAR::ACOU::acou_solid)
        interpolVec.Resize(ele->NumNode()*(2*numdim_+2+6)+2);
      else
        interpolVec.Resize(ele->NumNode()*(numdim_+2)+1);
    }

    ele->Evaluate(params,*discret_,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

    // sum values on nodes into vectors and record the touch count (build average of values)
    for (int i=0; i<ele->NumNode(); ++i)
    {
      DRT::Node* node = ele->Nodes()[i];
      const int localIndex = pressure->Map().LID(node->Id());

      if (localIndex < 0)
        continue;

      touchCount[localIndex]++;
      (*pressure)[localIndex] += interpolVec(i+numdim_*ele->NumNode());
    }
  }

  for (int i=0; i<pressure->MyLength(); ++i)
    (*pressure)[i] /= touchCount[i];

  for(int i=0; i<pressure->MyLength(); ++i)
    outvec->ReplaceMyValue(i,0,pressure->operator [](i));

  return;
} // NodalPressurField

/*----------------------------------------------------------------------*
 |  Return discretization (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
{
  if(calcerr_)
  {
    // call element routine
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::calc_acou_error);
    params.set<double>("time",time_);
    params.set<bool>("padaptivity",padaptivity_);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
    params.set<int>("funct",params_->get<int>("CALCERRORFUNCNO"));
    params.set<int>("useacouoptvecs",-1);
    params.set<double>("dt",dtp_);

    discret_->SetState(0,"trace",velnp_);

    Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(6));

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(params, errors);
    discret_->ClearState(true);

    // std::vector containing
    // [0]: L2 pressure error
    // [1]: L2 pressure norm
    // [2]: L2 velocity error
    // [3]: L2 velocity norm
    // [4]: L2 velocity gradient error
    // [5]: L2 velocity gradient norm

    Teuchos::RCP<std::vector<double> > relerror = Teuchos::rcp(new std::vector<double>(3));

    if ( (*errors)[1] != 0.0 )
      (*relerror)[0] = sqrt((*errors)[0])/sqrt((*errors)[1]);
    else if ((*errors)[0] != 0.0)
      (*relerror)[0] = 1.0;
    else
      (*relerror)[0] = 0.0;

    if ( (*errors)[3] != 0.0 )
      (*relerror)[1] = sqrt((*errors)[2])/sqrt((*errors)[3]);
    else if ((*errors)[2] != 0.0)
      (*relerror)[1] = 1.0;
    else
      (*relerror)[1] = 0.0;

    if ( (*errors)[5] != 0.0 )
      (*relerror)[2] = sqrt((*errors)[4])/sqrt((*errors)[5]);
    else if ((*errors)[4] != 0.0)
      (*relerror)[2] = 1.0;
    else
      (*relerror)[2] = 0.0;

    if(!myrank_)
    {
      //std::cout<<"time "<<time_<<" relative L2 pressure error "<<(*relerror)[0]<<" absolute L2 pressure error "<<sqrt((*errors)[0])<< " L2 pressure norm "<<sqrt((*errors)[1])<<std::endl;
      std::cout<<"time "<<time_<<" absolute L2 pressure error "<<sqrt((*errors)[0])<<" absolute L2 postprocessed pressure error "<<sqrt((*errors)[2])<<" L2 pressure norm "<<sqrt((*errors)[1])<<std::endl;
      if(phys_==INPAR::ACOU::acou_solid)
      {
        std::cout<<"time "<<time_<<" relative L2 velocity error "<<(*relerror)[1]<<" absolute L2 velocity error "<<sqrt((*errors)[2])<< " L2 velocity norm "<<sqrt((*errors)[3])<<std::endl;
        std::cout<<"time "<<time_<<" relative L2 velgradi error "<<(*relerror)[2]<<" absolute L2 velgradi error "<<sqrt((*errors)[4])<< " L2 velgradi norm "<<sqrt((*errors)[5])<<std::endl;
      }
    }
  }
  return;
}

/****************************************************************************************/
// ADER TRI TET MIT BACI
/****************************************************************************************/
ACOU::AcouTimeIntAderTriTet::AcouTimeIntAderTriTet(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  AcouImplicitTimeInt(actdis,solver,params,output)
{
  tempsrc_ = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
}

void ACOU::AcouTimeIntAderTriTet::Integrate(Teuchos::RCP<Epetra_MultiVector> history)
{

  // output of initial field (given by function for purely acoustic simulation or given by optics for PAT simulation)
  Output(history);

  // evaluate error
  //EvaluateErrorComparedToAnalyticalSol();

  AssembleMatAndRHS();

  // time loop
  while (step_<stepmax_ and time_<maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    // output to screen
    OutputToScreen();

    // solve
    Solve();

    // output of solution
    Output(history);

    // evaluate error
    EvaluateErrorComparedToAnalyticalSol();
  }
}

void ACOU::AcouTimeIntAderTriTet::Solve()
{
  // for postprocessing-correction:
  {
    // first step: evaluate right hand side for lambda system
    Teuchos::ParameterList eleparams;
    eleparams.set<double>("dt",dtp_);
    eleparams.set<bool>("padaptivity",padaptivity_);
    eleparams.set<int>("useacouoptvecs",-1);
    eleparams.set<int>("action",ACOU::prepare_ader_postprocessing);
    eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
    eleparams.set<bool>("adjoint",adjoint_);
    eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
    eleparams.set<bool>("allelesequal",false);

    // evaluate
    residual_->Scale(0.0);
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
    ApplyDirichletToSystem();

    // solve for lambda
    Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    solver_->Solve(sysmat_->EpetraOperator(),lambda,residual_,true,false,Teuchos::null);

    // second step: go back to the elements and store improved gradient and improved divergence
    eleparams.set<int>("action",ACOU::ader_postpro_gradanddiv);
    discret_->SetState(0,"trace",lambda);
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
    discret_->ClearState(true);
  }


  AssembleMatAndRHS();
  ApplyDirichletToSystem();

  solver_->Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,false,Teuchos::null);
  UpdateInteriorVariablesAndAssemebleRHS();

  return;
}

void ACOU::AcouTimeIntAderTriTet::AssembleMatAndRHS()
{
  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual and sysmat
  residual_->Scale(0.0);
  sysmat_->Zero();

  // set general vector values needed by elements
  discret_->ClearState(true);

  // set element parameters
  eleparams.set<double>("dt",dtp_);
  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<bool>("resonly",false);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<int>("useacouoptvecs",-1);
  eleparams.set<int>("action",ACOU::calc_ader_sysmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<double>("time",time_);
  eleparams.set<double>("timep",time_+dtp_);
  eleparams.set<int>("step",step_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  // set storage for intermediate values
  tempsrc_->PutScalar(0.0);
  discret_->SetState(1,"tempsrc",tempsrc_);

  // evaluate
  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  tempsrc_->Update(1.0,*(discret_->GetState(1,"tempsrc")),0.0);
  discret_->ClearState(true);

  // complete
  sysmat_->Complete();

  return;
}

void ACOU::AcouTimeIntAderTriTet::UpdateInteriorVariablesAndAssemebleRHS()
{
  // create parameterlist
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<double>("dt",dtp_);
  eleparams.set<double>("time",time_);
  eleparams.set<double>("timep",time_+dtp_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<double>("padaptivitytol",padapttol_);
  eleparams.set<int>("useacouoptvecs",-1);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  eleparams.set<bool>("allelesequal",allelesequal_);
  eleparams.set<int>("action",ACOU::update_ader_solution);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("calculategradient",false);
  eleparams.set<int>("step",step_);

  // evaluate
  residual_->Scale(0.0);
  discret_->SetState("trace",velnp_);
  discret_->SetState(1,"tempsrc",tempsrc_);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // clear
  discret_->ClearState(true);

  return;
}

void ACOU::AcouTimeIntAderTriTet::EvaluateErrorComparedToAnalyticalSol()
{
  if(calcerr_)
  {
    // call element routine
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::calc_acou_error);
    params.set<double>("time",time_);
    params.set<bool>("padaptivity",padaptivity_);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
    params.set<int>("funct",params_->get<int>("CALCERRORFUNCNO"));
    params.set<int>("useacouoptvecs",-1);
    params.set<double>("dt",dtp_);
    params.set<bool>("adjoint",false);
    params.set<bool>("allelesequal",allelesequal_);
    params.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
    params.set<bool>("calculategradient",false);
    params.set<int>("step",step_);

    discret_->SetState(0,"trace",velnp_);

    Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(6));

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(params, errors);
    discret_->ClearState(true);

    if(!myrank_)
      std::cout<<"time "<<time_<<" absolute L2 pressure error "<<sqrt((*errors)[0])<<" L2 pressure norm "<<sqrt((*errors)[1])<<std::endl;

    // postprocess solution!
    if(1)
    {
      params.set<int>("action",ACOU::prepare_ader_postprocessing);

      // evaluate
      residual_->Scale(0.0);
      discret_->Evaluate(params,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
      ApplyDirichletToSystem();

      // solve for lambda
      Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
      solver_->Solve(sysmat_->EpetraOperator(),lambda,residual_,true,false,Teuchos::null);

      // calculate improved gradient and divergence
      params.set<int>("action",ACOU::ader_postprocessing);
      double error_post_val;
      double* error_post = &error_post_val;
      *error_post = 0.0;
      params.set<double*>("error_post",error_post);

      // evaluate
      discret_->SetState(0,"trace",lambda);
      discret_->Evaluate(params,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
      discret_->ClearState(true);

      std::cout<<"time "<<time_<<" absolute L2 postprocessed pressure error "<<sqrt(*error_post)<<std::endl;

    }

  }
  return;
}
