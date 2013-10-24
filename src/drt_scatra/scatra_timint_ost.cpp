/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time-integration scheme

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_ost.H"
#include "scatra_ele_action.H"
#include "scatra_utils.H"
#include "turbulence_hit_scalar_forcing.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_inpar/inpar_elch.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"


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

  // initialize time-dependent electrode kinetics variables (galvanostatic mode or double layer contribution)
  if (IsElch(scatratype_))
    ComputeTimeDerivPot0(true);

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
 | dynamic Smagorinsky model                           krank  09/13     |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCv()
{
  if (turbmodel_==INPAR::FLUID::dynamic_vreman)
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    Vrem_->ApplyFilterForDynamicComputationOfDt(convel_,phinp_,thermpressnp_,dirichtoggle,*extraparams_);
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
  std::vector<std::string> condnames;
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
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "Data output for instationary thermodynamic pressure:" << std::endl;
    std::cout << "Velocity in-/outflow at indicated boundary: " << parnormvelint << std::endl;
    std::cout << "Diffusive flux at indicated boundary: "       << parnormdifffluxint << std::endl;
    std::cout << "Thermodynamic pressure: "                     << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
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

/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElectrodeKinetics",cond);
  int numcond = cond.size();

  for(int icond = 0; icond < numcond; icond++)
  {
    double pot0np =  cond[icond]->GetDouble("pot");
    const int curvenum =   cond[icond]->GetInt("curve");
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      // create and initialize additional b.c. entries for galvanostatic simulations or
      // simulations including a double layer
      cond[icond]->Add("pot0n",0.0);
      cond[icond]->Add("pot0dtnp",0.0);
      cond[icond]->Add("pot0dtn",0.0);
      cond[icond]->Add("pot0hist",0.0);

      if(dlcap!=0.0)
        dlcapexists_=true;
    }
    else
    {
      // compute time derivative of applied potential
      if (curvenum>=0)
      {
        const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time_);
        // adjust potential at metal side accordingly
        pot0np *= curvefac;
      }
      // compute time derivative of applied potential pot0
      // pot0dt(n+1) = (pot0(n+1)-pot0(n)) / (theta*dt) + (1-(1/theta))*pot0dt(n)
      double pot0n = cond[icond]->GetDouble("pot0n");
      double pot0dtn = cond[icond]->GetDouble("pot0dtn");
      double pot0dtnp =(pot0np-pot0n)/(dta_*theta_) + (1-(1/theta_))*pot0dtn;
      // add time derivative of applied potential pot0dtnp to BC
      cond[icond]->Add("pot0dtnp", pot0dtnp);
    }
  }

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
void SCATRA::TimIntOneStepTheta::Update(const int num)
{
  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (writeflux_!=INPAR::SCATRA::flux_no)
  {
    if (DoOutput() or DoBoundaryFluxStatistics())
      flux_ = CalcFlux(true, num);
  }

  // after the next command (time shift of solutions) do NOT call
  // ComputeTimeDerivative() anymore within the current time step!!!

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  // call time update of forcing routine
  if (homisoturb_forcing_ != Teuchos::null)
    homisoturb_forcing_->TimeUpdateForcing();

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

  // for elch problems with moving boundary
  if (isale_)
    output_->WriteVector("trueresidual", trueresidual_);

  // write additional restart data for galvanostatic applications or simulations including a double layer formulation
  if (IsElch(scatratype_))
  {
    if((DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC")) or
        dlcapexists_==true)
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
         // galvanostatic mode: only applied potential of cathode is adapted
         if (condid_cathode==condid or dlcapexists_==true)
         {
           std::stringstream temp;
           temp << condid;

           // electrode potential of the adjusted electrode kinetics BC at time n+1
           double pot = mycond->GetDouble("pot");
           output_->WriteDouble("pot_"+temp.str(),pot);

           // electrode potential of the adjusted electrode kinetics BC at time n
           double pot0n = mycond->GetDouble("pot0n");
           output_->WriteDouble("pot0n_"+temp.str(),pot0n);

           // electrode potential time derivative of the adjusted electrode kinetics BC at time n
           double pot0dtn = mycond->GetDouble("pot0dtn");
           output_->WriteDouble("pot0dtn_"+temp.str(),pot0dtn);

           // history of electrode potential of the adjusted electrode kinetics BC
           double pothist = mycond->GetDouble("pot0hist");
           output_->WriteDouble("pot0hist_"+temp.str(),pothist);
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

  if (scatratype_ == INPAR::SCATRA::scatratype_cardio_monodomain)
  {
   output_->WriteMesh(step_,time_); // add info to control file for reading all variables in restart
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
    std::cout<<"Reading ScaTra restart data (time="<<time_<<" ; step="<<step_<<")"<<std::endl;

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_, "phinp");
  reader.ReadVector(phin_,  "phin");
  reader.ReadVector(phidtn_,"phidtn");

  // for elch problems with moving boundary
  if(isale_)
    reader.ReadVector(trueresidual_, "trueresidual");

  // restart for galvanostatic applications
  if (IsElch(scatratype_))
  {
    // Initialize Nernst-BC
    InitNernstBC();

    if((DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC")) or
        dlcapexists_==true)
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
        // galvanostatic mode: only applied potential of cathode is adapted
        if (condid_cathode==condid or dlcapexists_==true)
        {
          std::stringstream temp;
          temp << condid;

          double pot = reader.ReadDouble("pot_"+temp.str());
          mycond->Add("pot",pot);
          double pot0n = reader.ReadDouble("pot0n_"+temp.str());
          mycond->Add("pot0n",pot0n);
          double pot0hist = reader.ReadDouble("pot0hist_"+temp.str());
          mycond->Add("pot0hist",pot0hist);
          double pot0dtn = reader.ReadDouble("pot0dtn_"+temp.str());
          mycond->Add("pot0dtn",pot0dtn);
          read_pot=true;
          if (myrank_==0)
            std::cout<<"Successfully read restart data for galvanostatic mode (condid "<<condid<<")"<<std::endl;
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

  if (scatratype_ == INPAR::SCATRA::scatratype_cardio_monodomain)
  {
    reader.ReadVector(activation_time_np_, "activation_time_np");
    reader.ReadMesh(step); // Read all saved data in nodes and elements und call nodal and element Unpacking each global variable has to be read

/*    Teuchos::ParameterList params; // THIS HAS TO BE EREASED LATER
    params.set<int>("scatratype", scatratype_);
    params.set<int>("action", SCATRA::set_material_internal_state);
    for(int k = 0; k < material_internal_state_np_->NumVectors(); ++k)
      {
        std::ostringstream temp;
        temp << k+1;
        reader.ReadVector(material_internal_state_np_component_, "mat_int_state_"+temp.str());
        params.set< int >("k", k);
        params.set< Teuchos::RCP<Epetra_Vector> >("material_internal_state_component", material_internal_state_np_component_);     // Probably do it once at the beginning
        discret_->Evaluate(params);
        //(*material_internal_state_np_)(k) = *material_internal_state_np_component_;
      }
      */
  }

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
void SCATRA::TimIntOneStepTheta::ElectrodeKineticsTimeUpdate()
{
  if (IsElch(scatratype_))
  {
    if((DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC")) or
        dlcapexists_==true)
    {
      ComputeTimeDerivPot0(false);

      std::vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);
      for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
      {
        {
          double pot0np = cond[i]->GetDouble("pot");
          cond[i]->Add("pot0n",pot0np);

          double pot0dtnp = cond[i]->GetDouble("pot0dtnp");
          cond[i]->Add("pot0dtn",pot0dtnp);
        }
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
    if((DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC")) or
        dlcapexists_==true)
    {
      std::vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);
      for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
      {
        // prepare "old part of rhs" for galvanostatic equation (to be used at this time step)
        {
          // re-read values (just to be really sure no mix-up occurs)
          double pot0n = cond[i]->GetDouble("pot0n");
          double pot0dtn = cond[i]->GetDouble("pot0dtn");
          // prepare old part of rhs for galvanostatic mode
          double pothist = pot0n + (1.0-theta_)*dta_*pot0dtn;
          cond[i]->Add("pot0hist",pothist);
        }
      }
    }
  }
  return;
}



