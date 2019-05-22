/*!-----------------------------------------------------------------------------------------------*
\brief scatra time integration for loma
\level 2
\maintainer Anh-Tu Vuong
 *------------------------------------------------------------------------------------------------*/
#include "../drt_fluid/fluid_utils.H"  // for splitter

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/tempdepwater.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"

#include "scatra_timint_loma.H"


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 12/13 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntLoma::ScaTraTimIntLoma(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
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
 | initialize algorithm                                     rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::Init()
{
  // safety check
  if (DRT::INPUT::IntegralValue<int>(*lomaparams_, "SGS_MATERIAL_UPDATE"))
    dserror(
        "Material update using subgrid-scale temperature currently not supported for loMa "
        "problems. Read remark in file 'scatra_ele_calc_loma.H'!");

  return;
}


/*----------------------------------------------------------------------*
 | setup algorithm                                          rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::Setup()
{
  SetupSplitter();
  return;
}

/*----------------------------------------------------------------------*
 | setup splitter                                          deanda 11/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::SetupSplitter()
{
  // set up a species-temperature splitter (if more than one scalar)
  if (NumScal() > 1)
  {
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_, NumScal() - 1, *splitter_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | set initial thermodynamic pressure                          vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::SetInitialThermPressure()
{
  // get thermodynamic pressure from material parameters
  int id = problem_->Materials()->FirstIdByType(INPAR::MAT::m_sutherland);
  if (id != -1)  // i.e., Sutherland material found
  {
    const MAT::PAR::Parameter* mat = problem_->Materials()->ParameterById(id);
    const MAT::PAR::Sutherland* actmat = static_cast<const MAT::PAR::Sutherland*>(mat);

    thermpressn_ = actmat->thermpress_;
  }
  else
  {
    // No Sutherland material found -> now check for temperature-dependent water,
    // which is allowed to be used in TFSI
    int id = problem_->Materials()->FirstIdByType(INPAR::MAT::m_tempdepwater);
    if (id != -1)  // i.e., temperature-dependent water found
    {
      // set thermodynamic pressure to zero once and for all
      thermpressn_ = 0.0;
    }
    else
      dserror(
          "Neiter Sutherland material nor temperature-dependent water found for initial setting of "
          "thermodynamic pressure!");
  }

  // initialize also value at n+1
  // (computed if not constant, otherwise prescribed value remaining)
  thermpressnp_ = thermpressn_;

  // initialize time derivative of thermodynamic pressure at n+1 and n
  // (computed if not constant, otherwise remaining zero)
  thermpressdtnp_ = 0.0;
  thermpressdtn_ = 0.0;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  // -> For constant thermodynamic pressure, this is done here once and
  // for all simulation time.
  ComputeThermPressureIntermediateValues();

  return;
}  // SCATRA::ScaTraTimIntLoma::SetInitialThermPressure


/*----------------------------------------------------------------------*
 | compute initial time derivative of thermodynamic pressure   vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeInitialThermPressureDeriv()
{
  // check for temperature-dependent water, which is allowed to be used in TFSI
  int id = problem_->Materials()->FirstIdByType(INPAR::MAT::m_tempdepwater);
  if (id != -1)
    dserror(
        "Temperature-dependent water found for initial computation of derivative of thermodynamic "
        "pressure -> Set 'CONSTHERMPRES' to 'YES' in FS3I input section!");

  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams);

  // set scalar vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp", phin_);

  // provide numbers of dofsets associated with velocity and displacement dofs
  eleparams.set<int>("ndsvel", nds_vel_);
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // set parameters for element evaluation
  eleparams.set<int>("action", SCATRA::calc_domain_and_bodyforce);

  // the time = 0.0, since this function is called BEFORE the first IncrementTimeAndStep() in
  // InitialCalculations() therefore, the standard SetElementTimeParameter() can be used for this
  // method

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate domain and bodyforce integral
  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint = (*scalars)[0];
  double parbofint = (*scalars)[1];

  // set action for elements
  eleparams.set<int>("action", SCATRA::bd_calc_loma_therm_press);

  // variables for integrals of normal velocity and diffusive flux
  double normvelint = 0.0;
  double normdifffluxint = 0.0;
  eleparams.set("normal velocity integral", normvelint);
  eleparams.set("normal diffusive flux integral", normdifffluxint);

  // evaluate velocity-divergence and diffusive (minus sign!) flux on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  std::vector<std::string> condnames;
  condnames.push_back("ScaTraFluxCalc");
  for (unsigned int i = 0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams, Teuchos::null, Teuchos::null, Teuchos::null,
        Teuchos::null, Teuchos::null, condnames[i]);
  }

  // get integral values on this proc
  normvelint = eleparams.get<double>("normal velocity integral");
  normdifffluxint = eleparams.get<double>("normal diffusive flux integral");

  // get integral values in parallel case
  double parnormvelint = 0.0;
  double parnormdifffluxint = 0.0;
  discret_->Comm().SumAll(&normvelint, &parnormvelint, 1);
  discret_->Comm().SumAll(&normdifffluxint, &parnormdifffluxint, 1);

  // clean up
  discret_->ClearState();

  // compute initial time derivative of thermodynamic pressure
  // (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  thermpressdtn_ =
      (-shr * thermpressn_ * parnormvelint + (shr - 1.0) * (-parnormdifffluxint + parbofint)) /
      pardomint;

  // set time derivative of thermodynamic pressure at n+1 equal to the one at n
  // for following evaluation of intermediate values
  thermpressdtnp_ = thermpressdtn_;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
}  // SCATRA::ScaTraTimIntLoma::ComputeInitialThermPressureDeriv


/*----------------------------------------------------------------------*
 | compute initial total mass in domain                        vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeInitialMass()
{
  // check for temperature-dependent water, which is allowed to be used in TFSI
  int id = problem_->Materials()->FirstIdByType(INPAR::MAT::m_tempdepwater);
  if (id != -1)
    dserror(
        "Temperature-dependent water found for initial computation of mass -> Set 'CONSTHERMPRES' "
        "to 'YES' in FS3I input section!");

  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp", phin_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  // inverted scalar values are required here
  eleparams.set("inverting", true);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // compute initial mass times gas constant: R*M_0 = int(1/T_0)*tp
  initialmass_ = (*scalars)[0] * thermpressn_;

  // print out initial total mass
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
    std::cout << "Initial total mass in domain (times gas constant): " << initialmass_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
  }

  return;
}  // SCATRA::ScaTraTimIntLoma::ComputeInitialMass


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure from mass conservation       vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeThermPressureFromMassCons()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp", phinp_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  // inverted scalar values are required here
  eleparams.set("inverting", true);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // compute thermodynamic pressure: tp = R*M_0/int(1/T)
  thermpressnp_ = initialmass_ / (*scalars)[0];

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
    std::cout << "Thermodynamic pressure from mass conservation: " << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  ComputeThermPressureTimeDerivative();

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
}  // SCATRA::ScaTraTimIntLoma::ComputeThermPressureFromMassCons


/*----------------------------------------------------------------------*
 | add parameters depending on the problem              rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  AddThermPressToParameterList(params);
  return;
}
