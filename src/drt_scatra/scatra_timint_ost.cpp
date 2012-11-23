/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time-integration scheme

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_ost.H"
#include "scatra_ele_action.H"
#include "scatra_utils.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_inpar/inpar_elch.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fluid/dyn_smag.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
  Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  theta_(params_->get<double>("THETA"))
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // temporal solution derivative at time n
  phidtn_ = LINALG::CreateVector(*dofrowmap,true);

  //solution at time n-1, for level set problems
  if (scatratype_ == INPAR::SCATRA::scatratype_levelset)
     phinm_  = LINALG::CreateVector(*dofrowmap,true);

  // ELCH with natural convection
  if (extraparams_->isSublist("ELCH CONTROL"))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"NATURAL_CONVECTION") == true)
    {
      // density at time n
      elchdensn_ = LINALG::CreateVector(*dofrowmap,true);
      elchdensn_->PutScalar(1.0);
    }
  }

  // fine-scale vector at time n+1
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    fsphinp_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize time-dependent electrode kinetics variables (galvanostatic mode)
  ElectrodeKineticsTimeUpdate(true);

  // Important: this adds the required ConditionID's to the single conditions.
  // It is necessary to do this BEFORE ReadRestart() is called!
  // Output to screen and file is suppressed
  OutputElectrodeInfo(false,false);

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    gjb 08/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta()
{
  return;
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen                   |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    IO::cout << "TIME: "
             << std::setw(11) << std::setprecision(4) << std::scientific << time_ << "/"
             << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_ << "  DT = "
             << std::setw(11) << std::setprecision(4) << std::scientific << dta_ << "  "
             << MethodTitle() << " (theta = "
             << std::setw(3)  << std::setprecision(2) << theta_ << ") STEP = "
             << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_ << IO::endl;
  }
  return;
}


/*----------------------------------------------------------------------*
| verify if given coefficients are in admissable range        gjb 04/10 |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::VerifyCoeff()
{
  if ( (theta_ <= 0.0) or (theta_ > 1.0) )
    dserror("theta out of range (0.0,1.0]");
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  // hist_ = phin_ + dt*(1-Theta)*phidtn_
  hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);

  // for electrochemical applications
  ElectrodeKineticsSetOldPartOfRHS();

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dta_, *phidtn_,1.0);

  // for the electric potential we just use the 'old' value of last time step
  Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(phin_);
  splitter_->InsertCondVector(onlypot, phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-derivative predictor (currently not used)
  //thermpressnp_ += dta_*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetTimeForNeumannEvaluation(
  Teuchos::ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddNeumannToResidual()
{
  residual_->Update(theta_*dta_,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AVM3Separation()
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
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCs()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(convel_,phinp_,thermpressnp_,dirichtoggle,*extraparams_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddSpecificTimeIntegrationParameters(
  Teuchos::ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);
  params.set("time factor",theta_*dta_);
  params.set("alpha_F",1.0);

  if (scatratype_==INPAR::SCATRA::scatratype_loma)
  {
    params.set("thermodynamic pressure",thermpressnp_);
    params.set("time derivative of thermodynamic pressure",thermpressdtnp_);
  }

  discret_->SetState("hist",hist_);
  discret_->SetState("phinp",phinp_);
  if (scatratype_==INPAR::SCATRA::scatratype_levelset)
  {
//    discret_->SetState("phin",phin_);
  }

  if(reinitswitch_) params.set("theta_reinit", theta_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeThermPressure()
{
  // compute "history" part
  const double hist = thermpressn_ + (1.0-theta_)*dta_*thermpressdtn_;

  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams,INPAR::SCATRA::flux_diffusive_domain);

  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // provide velocity field (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"velocity field",vel_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set action for elements
  eleparams.set<int>("action",SCATRA::calc_domain_and_bodyforce);
  eleparams.set("total time",time_);
  eleparams.set<int>("scatratype",scatratype_);

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate domain and bodyforce integral
  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint  = (*scalars)[0];
  double parbofint  = (*scalars)[1];

  // set action for elements
  eleparams.set<int>("action",SCATRA::bd_calc_loma_therm_press);

  // variables for integrals of normal velocity and diffusive flux
  double normvelint      = 0.0;
  double normdifffluxint = 0.0;
  eleparams.set("normal velocity integral",normvelint);
  eleparams.set("normal diffusive flux integral",normdifffluxint);

  // evaluate velocity-divergence and diffusive (minus sign!) flux on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  vector<std::string> condnames;
  condnames.push_back("ScaTraFluxCalc");
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condnames[i]);
  }

  // get integral values on this proc
  normvelint      = eleparams.get<double>("normal velocity integral");
  normdifffluxint = eleparams.get<double>("normal diffusive flux integral");

  // get integral values in parallel case
  double parnormvelint      = 0.0;
  double parnormdifffluxint = 0.0;
  discret_->Comm().SumAll(&normvelint,&parnormvelint,1);
  discret_->Comm().SumAll(&normdifffluxint,&parnormdifffluxint,1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  const double lhs = theta_*dta_*shr*parnormvelint/pardomint;
  const double rhs = theta_*dta_*(shr-1.0)*(-parnormdifffluxint+parbofint)/pardomint;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    cout << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
    cout << "Data output for instationary thermodynamic pressure:" << endl;
    cout << "Velocity in-/outflow at indicated boundary: " << parnormvelint << endl;
    cout << "Diffusive flux at indicated boundary: "       << parnormdifffluxint << endl;
    cout << "Thermodynamic pressure: "                     << thermpressnp_ << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative                                     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeTimeDerivative()
{
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0/(theta_*dta_);
  const double fact2 = 1.0 - (1.0/theta_);
  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

  // We know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  // However, we do not want to break the linear relationship
  // as stated above. We do not want to set Dirichlet values for
  // dependent values like phidtnp_. This turned out to be inconsistent.
  // ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative of thermodynamic pressure           vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeThermPressureTimeDerivative()
{
  // time derivative of thermodynamic pressure:
  // tpdt(n+1) = (tp(n+1)-tp(n))/(theta*dt)+((theta-1)/theta)*tpdt(n)
  double fact1 = 1.0/(theta_*dta_);
  double fact2 = (theta_-1.0)/theta_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update()
{
  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (writeflux_!=INPAR::SCATRA::flux_no)
  {
    if (DoOutput() or DoBoundaryFluxStatistics())
      flux_ = CalcFlux(true);
  }

  // after the next command (time shift of solutions) do NOT call
  // ComputeTimeDerivative() anymore within the current time step!!!

  // phinm is needed for restart of level set problems
  if (scatratype_ == INPAR::SCATRA::scatratype_levelset)
  {
     phinm_ ->Update(1.0,*phin_,0.0);
  }

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  // perform update of time-dependent electrode variables
  ElectrodeKineticsTimeUpdate();

  // potential time update of time-dependent materials
  ElementMaterialTimeUpdate();

  return;
}

/*----------------------------------------------------------------------*
 | update level set after reinitialization              rasthofer 02/10 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateReinit()
{
  //phinm is needed for restart of level set problems
  phinm_ ->Update(1.0,*phin_,0.0);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // compute time derivative at time n (and n+1)
  CalcInitialPhidt();
}

/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n for ELCH natural convection            gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateDensityElch()
{
  elchdensn_->Update(1.0,*elchdensnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  // phinm is needed to reconstruct the interface
  if (scatratype_ == INPAR::SCATRA::scatratype_levelset)
  {
     output_->WriteVector("phinm", phinm_);
  }

  // for elch problems with moving boundary
  if (isale_)
    output_->WriteVector("trueresidual", trueresidual_);

  // write additional restart data for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");

       vector<DRT::Condition*>::iterator fool;
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

           // electrode potential of the adjusted electrode kinetics BC at time n
           double potn = mycond->GetDouble("potn");
           output_->WriteDouble("potn",potn);

           // electrode potential time derivative of the adjusted electrode kinetics BC at time n
           double potdtn = mycond->GetDouble("potdtn");
           output_->WriteDouble("potdtn",potdtn);

           // history of electrode potential of the adjusted electrode kinetics BC
           double pothist = mycond->GetDouble("pothist");
           output_->WriteDouble("pothist",pothist);
         }
       }
    }
  }

  // write additional restart data for loma
  // required for restart of closed systems
  if (scatratype_ == INPAR::SCATRA::scatratype_loma)
  {
    // thermodynamic pressure at time n+1
    output_->WriteDouble("thermpressnp",thermpressnp_);
    // thermodynamic pressure at time n
    output_->WriteDouble("thermpressn",thermpressn_);
    // time derivative of thermodynamic pressure at time n+1
    output_->WriteDouble("thermpressdtnp",thermpressdtnp_);
    // time derivative of thermodynamic pressure at time n
    output_->WriteDouble("thermpressdtn",thermpressdtn_);
    // as well as initial mass
    output_->WriteDouble("initialmass",initialmass_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  if (myrank_==0)
    cout<<"Reading ScaTra restart data (time="<<time_<<" ; step="<<step_<<")"<<endl;

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_, "phinp");
  reader.ReadVector(phin_,  "phin");
  reader.ReadVector(phidtn_,"phidtn");

  // for elch problems with moving boundary
  if(isale_)
    reader.ReadVector(trueresidual_, "trueresidual");

  // phinm is needed for restart of level set problems
  if (scatratype_ == INPAR::SCATRA::scatratype_levelset)
  {
     reader.ReadVector(phinm_,  "phinm");
  }

  // restart for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");
      vector<DRT::Condition*>::iterator fool;
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
          double potn = reader.ReadDouble("potn");
          mycond->Add("potn",potn);
          double pothist = reader.ReadDouble("pothist");
          mycond->Add("pothist",pothist);
          double potdtn = reader.ReadDouble("potdtn");
          mycond->Add("potdtn",potdtn);
          read_pot=true;
          if (myrank_==0)
            cout<<"Successfully read restart data for galvanostatic mode (condid "<<condid<<")"<<endl;
        }
      }
      if (!read_pot)
        dserror("Reading of electrode potential for restart not successful.");
    }
  }

  // restart data of loma problems
  // required for restart of closed systems
  if (scatratype_ == INPAR::SCATRA::scatratype_loma)
  {
    // thermodynamic pressure at time n+1
    thermpressnp_ = reader.ReadDouble("thermpressnp");
    // thermodynamic pressure at time n
    thermpressn_ = reader.ReadDouble("thermpressn");
    // time derivative of thermodynamic pressure at time n+1
    thermpressdtnp_ = reader.ReadDouble("thermpressdtnp");
    // time derivative of thermodynamic pressure at time n
    thermpressdtn_ = reader.ReadDouble("thermpressdtn");
    // as well as initial mass
    initialmass_ = reader.ReadDouble("initialmass");
  }

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    AVM3Preparation();

  return;
}


/*--------------------------------------------------------------------------------------------*
 | Redistribute the scatra discretization and vectors according to nodegraph  rasthofer 07/11 |
 |                                                                            DA wichmann     |
 *--------------------------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  // let the base class do the basic redistribution and transfer of the base class members
  ScaTraTimIntImpl::Redistribute(nodegraph);

  // now do all the ost specfic steps
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> old;

  if (fsphinp_ != Teuchos::null)
  {
    old = fsphinp_;
    fsphinp_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *fsphinp_);
  }

  if (phinm_ != Teuchos::null)
  {
    old = phinm_;
    phinm_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *phinm_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step        gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrepareFirstTimeStep()
{
  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  //ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  // compute initial field for electric potential (ELCH)
  CalcInitialPotentialField();

  // compute time derivative of phi at time t=0
  CalcInitialPhidt();

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ElectrodeKineticsTimeUpdate(const bool init)
{
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);
      for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
      {
        // compute time derivative and perform time-shift of variables
        {
          double potnp = cond[i]->GetDouble("pot");
          if (init) // create and initialize additional b.c. entries if desired
          {
            cond[i]->Add("potn",potnp);
            cond[i]->Add("potdtn",0.0);
            cond[i]->Add("pothist",0.0);
          }
          double potn = cond[i]->GetDouble("potn");
          double potdtn = cond[i]->GetDouble("potdtn");
          // update the time derivative
          double potdtn_new = (1/(theta_*dta_))*(potnp- potn) + (1 - (1/theta_))*potdtn;
          cond[i]->Add("potdtn",potdtn_new);
          // time-shift of status variables
          cond[i]->Add("potn",potnp);
        }

        // prepare "old part of rhs" for galvanostatic equation (to be used at next time step)
#if 0
        {
          // re-read values (just to be really sure no mix-up occurs)
          double potn = cond[i]->GetDouble("potn");
          double potdtn = cond[i]->GetDouble("potdtn");
          // prepare old part of rhs for galvanostatic mode
          double pothist = potn + (1.0-theta_)*dta_*potdtn;
          cond[i]->Add("pothist",pothist);
        }
#endif
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | set old part of RHS for galvanostatic equation             gjb 04/10 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ElectrodeKineticsSetOldPartOfRHS()
{
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);
      for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
      {
        // prepare "old part of rhs" for galvanostatic equation (to be used at this time step)
        {
          // re-read values (just to be really sure no mix-up occurs)
          double potn = cond[i]->GetDouble("potn");
          double potdtn = cond[i]->GetDouble("potdtn");
          // prepare old part of rhs for galvanostatic mode
          double pothist = potn + (1.0-theta_)*dta_*potdtn;
          cond[i]->Add("pothist",pothist);
        }
      }
    }
  }
  return;
}



