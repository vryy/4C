/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_loma.cpp

  \brief scatra time integration for loma

<pre>
Maintainer: Ursula Rasthofer / Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/45
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "scatra_timint_loma.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/sutherland.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_fluid/fluid_utils.H" // for splitter


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 12/13 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntLoma::ScaTraTimIntLoma(
        Teuchos::RCP<DRT::Discretization>        dis,
        Teuchos::RCP<LINALG::Solver>             solver,
        Teuchos::RCP<Teuchos::ParameterList>     params,
        Teuchos::RCP<Teuchos::ParameterList>     sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList>     extraparams,
        Teuchos::RCP<IO::DiscretizationWriter>   output)
  : ScaTraTimIntImpl(dis,solver,sctratimintparams,extraparams,output),
    lomaparams_(params),
    initialmass_(0.0),
    thermpressn_(0.0),
    thermpressnp_(0.0),
    thermpressdtn_(0.0),
    thermpressdtnp_(0.0)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                 rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::Init()
{
  // set up a species-temperature splitter (if more than one scalar)
  if (numscal_ > 1)
  {
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_-1,*splitter_);
  }

  // safety check
  if (DRT::INPUT::IntegralValue<int>(*lomaparams_,"SGS_MATERIAL_UPDATE"))
    dserror("Material update using subgrid-scale temperature currently not supported for loMa problems. Read remark in file 'scatra_ele_calc_loma.H'!");

  return;
}


/*----------------------------------------------------------------------*
 | set initial thermodynamic pressure                          vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::SetInitialThermPressure()
{
  // get thermodynamic pressure from material parameters
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_sutherland);
  if (id !=-1) // i.e., Sutherland material found
  {
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::Sutherland* actmat = static_cast<const MAT::PAR::Sutherland*>(mat);

    thermpressn_ = actmat->thermpress_;
  }
  else thermpressn_ = 0.0;


  // initialize also value at n+1
  // (computed if not constant, otherwise prescribed value remaining)
  thermpressnp_ = thermpressn_;

  // initialize time derivative of thermodynamic pressure at n+1 and n
  // (computed if not constant, otherwise remaining zero)
  thermpressdtnp_ = 0.0;
  thermpressdtn_  = 0.0;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  // -> For constant thermodynamic pressure, this is done here once and
  // for all simulation time.
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntLoma::SetInitialThermPressure


/*----------------------------------------------------------------------*
 | compute initial time derivative of thermodynamic pressure   vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeInitialThermPressureDeriv()
{
  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams,INPAR::SCATRA::flux_diffusive_domain);

  // set scalar vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phin_);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  discret_->AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
  discret_->AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
  discret_->AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);

  // provide displacement field in case of ALE
  if (isale_) discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set parameters for element evaluation
  eleparams.set<int>("action",SCATRA::calc_domain_and_bodyforce);
  eleparams.set<int>("scatratype",scatratype_);

  // the time = 0.0, since this function is called BEFORE the first IncrementTimeAndStep() in InitialCalculations()
  // therefore, the standard SetElementTimeParameter() can be used for this method

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

  // compute initial time derivative of thermodynamic pressure
  // (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  thermpressdtn_ = (-shr*thermpressn_*parnormvelint
                    + (shr-1.0)*(-parnormdifffluxint+parbofint))/pardomint;

  // set time derivative of thermodynamic pressure at n+1 equal to the one at n
  // for following evaluation of intermediate values
  thermpressdtnp_ = thermpressdtn_;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntLoma::ComputeInitialThermPressureDeriv


/*----------------------------------------------------------------------*
 | compute initial total mass in domain                        vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeInitialMass()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phin_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mean_scalars);
  eleparams.set<int>("scatratype",scatratype_);
  // inverted scalar values are required here
  eleparams.set("inverting",true);

  //provide displacement field in case of ALE
  if (isale_) discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();   // clean up

  // compute initial mass times gas constant: R*M_0 = int(1/T_0)*tp
  initialmass_ = (*scalars)[0]*thermpressn_;

  // print out initial total mass
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "Initial total mass in domain (times gas constant): " << initialmass_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
  }

  return;
} // SCATRA::ScaTraTimIntLoma::ComputeInitialMass


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure from mass conservation       vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeThermPressureFromMassCons()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mean_scalars);
  eleparams.set<int>("scatratype",scatratype_);
  // inverted scalar values are required here
  eleparams.set("inverting",true);

  //provide displacement field in case of ALE
  if (isale_) discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();   // clean up

  // compute thermodynamic pressure: tp = R*M_0/int(1/T)
  thermpressnp_ = initialmass_/(*scalars)[0];

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "Thermodynamic pressure from mass conservation: " << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  ComputeThermPressureTimeDerivative();

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntLoma::ComputeThermPressureFromMassCons


/*----------------------------------------------------------------------*
 | convergence check (only for low-Mach-number flow)           vg 09/11 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntLoma::ConvergenceCheck(int          itnum,
                                                int          itmax,
                                                const double ittol)
{
  bool stopnonliniter = false;

  // define L2-norm of residual, incremental scalar and scalar
  double resnorm_L2(0.0);
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);

  // for the time being, only one scalar considered for low-Mach-number flow
  /*if (numscal_>1)
  {
    Teuchos::RCP<Epetra_Vector> onlyphi = splitter_->ExtractCondVector(increment_);
    onlyphi->Norm2(&phiincnorm_L2);

    splitter_->ExtractCondVector(phinp_,onlyphi);
    onlyphi->Norm2(&phinorm_L2);
  }
  else*/
  residual_ ->Norm2(&resnorm_L2);
  increment_->Norm2(&phiincnorm_L2);
  phinp_    ->Norm2(&phinorm_L2);

  // check for any INF's and NaN's
  if (std::isnan(resnorm_L2) or
      std::isnan(phiincnorm_L2) or
      std::isnan(phinorm_L2))
    dserror("At least one of the calculated vector norms is NaN.");

  if (abs(std::isinf(resnorm_L2)) or
      abs(std::isinf(phiincnorm_L2)) or
      abs(std::isinf(phinorm_L2)))
    dserror("At least one of the calculated vector norms is INF.");

  // for scalar norm being (close to) zero, set to one
  if (phinorm_L2 < 1e-5) phinorm_L2 = 1.0;

  if (myrank_==0)
  {
    printf("+------------+-------------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|- residual   -|- scalar-inc -|\n");
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
         itnum,itmax,ittol,resnorm_L2,phiincnorm_L2/phinorm_L2);
    printf("\n");
    printf("+------------+-------------------+--------------+--------------+\n");
  }

  if ((resnorm_L2 <= ittol) and
      (phiincnorm_L2/phinorm_L2 <= ittol)) stopnonliniter=true;

  // warn if itemax is reached without convergence, but proceed to next timestep
  if ((itnum == itmax) and
      ((resnorm_L2 > ittol) or (phiincnorm_L2/phinorm_L2 > ittol)))
  {
    stopnonliniter=true;
    if (myrank_==0)
    {
      printf("|            >>>>>> not converged in itemax steps!             |\n");
      printf("+--------------------------------------------------------------+\n");
    }
  }

  return stopnonliniter;
} // SCATRA::ScaTraTimIntLoma::ConvergenceCheck


/*----------------------------------------------------------------------*
 | add parameters depending on the problem              rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::AddProblemSpecificParametersAndVectors(
  Teuchos::ParameterList& params //!< parameter list
)
{
  AddThermPressToParameterList(params);
  return;
}


/*--------------------------------------------------------------------------*
 | add parameters depending on the problem for inital phidt rasthofer 12/13 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::AddProblemSpecificParametersAndVectorsForCalcInitialPhiDt(
  Teuchos::ParameterList& params //!< parameter list
)
{
  params.set("thermodynamic pressure",thermpressn_);
  params.set("time derivative of thermodynamic pressure",thermpressdtn_);
  return;
}
