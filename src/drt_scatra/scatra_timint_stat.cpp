/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_stat.cpp
\brief solution algorithm for stationary problems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_stat.H"
#include "scatra_utils.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntStationary::TimIntStationary(
  Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::Init()
{
  // initialize base class
  ScaTraTimIntImpl::Init();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // fine-scale vector
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    fsphinp_ = LINALG::CreateVector(*dofrowmap,true);
  if (turbmodel_ != INPAR::FLUID::no_model)
    dserror("Turbulence is not stationary problem!");

  // -------------------------------------------------------------------
  // set element parameters
  // -------------------------------------------------------------------
  // note: - this has to be done before element routines are called
  //       - order is important here: for savety checks in SetElementGeneralScaTraParameter(),
  //         we have to konw the time-integration parameters
  SetElementTimeParameter();
  SetElementGeneralScaTraParameter();
  SetElementTurbulenceParameter();

  // initialize time-dependent electrode kinetics variables (galvanostatic mode)
  if (IsElch(scatratype_))
    ComputeTimeDerivPot0(true);

  // Important: this adds the required ConditionID's to the single conditions.
  // It is necessary to do this BEFORE ReadRestart() is called!
  // Output to screen and file is suppressed
  OutputElectrodeInfo(false,false);

  // setup krylov
  PrepareKrylovProjection();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntStationary::~TimIntStationary()
{
  return;
}


/*----------------------------------------------------------------------*
 | set time parameter for element evaluation (usual call)    ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",true);
  eleparams.set<bool>("incremental solver",incremental_);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",time_);
  eleparams.set<double>("time factor",1.0);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation                ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetElementTimeParameterForForcedIncrementalSolve()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",true);
  // this is important to have here and the only difference compared to SetElementTimeParameter()
  eleparams.set<bool>("incremental solver",true);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",time_);
  eleparams.set<double>("time factor",1.0);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | set time parameter for element evaluation                 ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetElementTimeParameterInitial()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",true);
  eleparams.set<bool>("incremental solver",incremental_);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",0.0); // only different
  eleparams.set<double>("time factor",1.0);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetOldPartOfRighthandside()
{
  hist_->PutScalar(0.0);

  return;
}


/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntStationary::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElectrodeKinetics",cond);
  int numcond = cond.size();

  for(int icond = 0; icond < numcond; icond++)
  {
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      if(dlcap!=0.0)
        dserror("Double layer charging and galvanostatic mode are not implemented for BDF2! You have to use one-step-theta time integration scheme");
        //dlcapexists_=true;

      if(DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC")==true)
          dserror("Double layer charging and galvanostatic mode are not implemented for BDF2! You have to use one-step-theta time integration scheme");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetTimeForNeumannEvaluation(
  Teuchos::ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                    vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AddNeumannToResidual()
{
  residual_->Update(1.0,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false,*phinp_,*fsphinp_);

  // set fine-scale vector
  discret_->SetState("fsphinp",fsphinp_);

  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AddTimeIntegrationSpecificVectors()
{
  //SetElementTimeParameter();

  discret_->SetState("hist",hist_);
  discret_->SetState("phinp",phinp_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 09/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  if (myrank_==0)
    std::cout<<"Reading ScaTra restart data (time="<<time_<<" ; step="<<step_<<")"<<std::endl;

  // read state vectors that are needed for restart
  reader.ReadVector(phinp_, "phinp");

  // for elch problems with moving boundary
  // if(isale_)
  //  reader.ReadVector(trueresidual_, "trueresidual");

  // Initialize Nernst-BC
  InitNernstBC();

  // restart for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      std::vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");
      std::vector<DRT::Condition*>::iterator fool;
      bool read_pot=false;

      // read desired values from the .control file and add/set the value to
      // the electrode kinetics boundary condition representing the cathode
      for (fool=cond.begin(); fool!=cond.end(); ++fool)
      {
        DRT::Condition* mycond = (*(fool));
        const int condid = mycond->GetInt("ConditionID");
        if (condid_cathode==condid)
        {
          double pot = reader.ReadDouble("pot");
          mycond->Add("pot",pot);
          read_pot=true;
          if (myrank_==0)
            std::cout<<"Successfully read restart data for galvanostatic mode (condid "<<condid<<")"<<std::endl;
        }
      }
      if (!read_pot)
        dserror("Reading of electrode potential for restart not successful.");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | update of solution at end of time step                     gjb 12/10 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::Update(const int num)
{
  // for the stationary scheme there is nothing to do except this:

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (writeflux_!=INPAR::SCATRA::flux_no)
  {
    if (DoOutput() or DoBoundaryFluxStatistics())
      flux_ = CalcFlux(true, num);
  }
  return;
};

/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 10/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::OutputRestart()
{
  // This feature enables starting a time-dependent simulation from
  // a non-trivial steady-state solution that was calculated before.
  output_->WriteVector("phin", phinp_);  // for OST and BDF2
  output_->WriteVector("phinm", phinp_); // for BDF2
  output_->WriteVector("phidtn", zeros_); // for OST

  // for elch problems with moving boundary
  //if (isale_)
  //  output_->WriteVector("trueresidual", trueresidual_);

  // write additional restart data for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      std::vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");

      std::vector<DRT::Condition*>::iterator fool;
      // loop through conditions and find the cathode
      for (fool=cond.begin(); fool!=cond.end(); ++fool)
      {
        DRT::Condition* mycond = (*(fool));
        const int condid = mycond->GetInt("ConditionID");
        if (condid_cathode==condid)
        {
          // electrode potential of the adjusted electrode kinetics BC at time n+1
          double pot = mycond->GetDouble("pot");
          output_->WriteDouble("pot",pot);
        }
      }
    }
  }

  return;
}

