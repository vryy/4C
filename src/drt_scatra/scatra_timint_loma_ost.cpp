/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_loma_ost.cpp
\brief One-step-theta time-integration scheme with extensions for
       loma problems

\level 2

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_loma_ost.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "turbulence_hit_scalar_forcing.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       vg 11/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntLomaOST::TimIntLomaOST(
  Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   sctratimintparams,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
  : ScaTraTimIntImpl(actdis,solver,sctratimintparams,extraparams,output),
    ScaTraTimIntLoma(actdis,solver,params,sctratimintparams,extraparams,output),
    TimIntOneStepTheta(actdis,solver,sctratimintparams,extraparams,output)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  ScaTraTimIntLoma::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     vg 11/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntLomaOST::~TimIntLomaOST()
{
  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-derivative predictor (currently not used)
  //thermpressnp_ += dta_*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::DynamicComputationOfCs()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(phinp_,thermpressnp_,dirichtoggle,*extraparams_,nds_vel_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           krank  09/13     |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::DynamicComputationOfCv()
{
  if (turbmodel_==INPAR::FLUID::dynamic_vreman)
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    Vrem_->ApplyFilterForDynamicComputationOfDt(phinp_,thermpressnp_,dirichtoggle,*extraparams_,nds_vel_);
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | add thermodynamic pressure to parameter list for element evaluation rasthofer 12/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::AddThermPressToParameterList(
  Teuchos::ParameterList& params //!< parameter list
)
{
  params.set("thermodynamic pressure",thermpressnp_);
  params.set("time derivative of thermodynamic pressure",thermpressdtnp_);
  return;
}

/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::ComputeThermPressure()
{
  // compute "history" part
  const double hist = thermpressn_ + (1.0-theta_)*dta_*thermpressdtn_;

  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams);

  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // provide numbers of dofsets associated with velocity and displacement dofs
  eleparams.set<int>("ndsvel",nds_vel_);
  if (isale_)
    eleparams.set<int>("ndsdisp",nds_disp_);

  // set action for elements
  eleparams.set<int>("action",SCATRA::calc_domain_and_bodyforce);
  SetElementTimeParameter();

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
 | compute time derivative of thermodynamic pressure           vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::ComputeThermPressureTimeDerivative()
{
  // time derivative of thermodynamic pressure:
  // tpdt(n+1) = (tp(n+1)-tp(n))/(theta*dt)+((theta-1)/theta)*tpdt(n)
  double fact1 = 1.0/(theta_*dta_);
  double fact2 = (theta_-1.0)/theta_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::OutputRestart()
{
  // do standard output
  TimIntOneStepTheta::OutputRestart();

  // write additional restart data for loma
  // required for restart of closed systems

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

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntLomaOST::ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input)
{
  // do standard output
  TimIntOneStepTheta::ReadRestart(step,input);

  // restart data of loma problems
  // required for restart of closed systems

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if(input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_,step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_,input,step));

  // thermodynamic pressure at time n+1
  thermpressnp_ = reader->ReadDouble("thermpressnp");
  // thermodynamic pressure at time n
  thermpressn_ = reader->ReadDouble("thermpressn");
  // time derivative of thermodynamic pressure at time n+1
  thermpressdtnp_ = reader->ReadDouble("thermpressdtnp");
  // time derivative of thermodynamic pressure at time n
  thermpressdtn_ = reader->ReadDouble("thermpressdtn");
  // as well as initial mass
  initialmass_ = reader->ReadDouble("initialmass");

  return;
}
