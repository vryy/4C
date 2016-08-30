/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_std.cpp

\brief singleton class holding all static parameters required for the evaluation of a standard scalar transport element

This singleton class holds all static parameters required for the evaluation of a standard scalar transport element,
e.g., stabilization parameters and finite difference check parameters. All parameters are usually set only once at
the beginning of a simulation, namely during initialization of the global time integrator, and then never touched again
throughout the simulation. Enhanced scalar transport problems, such as electrochemistry and levelset problems, instantiate
additional, problem specific singleton classes holding additional static parameters required for element evaluation. These
additional singleton classes are not meant to be derived from, but rather to coexist with this general class.

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                             thon/vuong 07/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd* DRT::ELEMENTS::ScaTraEleParameterStd::Instance(
    const std::string&           disname,   //!< name of discretization
    const ScaTraEleParameterStd* delete_me  //!< creation/destruction indication
    )
{
  // each discretization is associated with exactly one instance of this class according to a static map
  static std::map<std::string,ScaTraEleParameterStd*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if not
  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterStd(disname);
  }

  // destruct instance given to the destructor
  else
  {
    for(std::map<std::string,ScaTraEleParameterStd*>::iterator i=instances.begin(); i!=instances.end(); ++i)
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterStd::Done()
{
  // delete singleton
  Instance("",this);

  return;
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd::ScaTraEleParameterStd(
    const std::string& disname   //!< name of discretization
    ) :
    is_ale_(false),
    is_conservative_(false),
    sphericalcoords_(false),
    calcflux_domain_(INPAR::SCATRA::flux_none),
    writefluxids_(Teuchos::null),
    fdcheck_(INPAR::SCATRA::fdcheck_none),
    fdcheckeps_(0.),
    fdchecktol_(0.),
    stabtype_(INPAR::SCATRA::stabtype_no_stabilization),
    whichtau_(INPAR::SCATRA::tau_zero),
    charelelength_(INPAR::SCATRA::streamlength),
    diffreastafac_(0.0),
    sgvel_(false),
    assgd_(false),
    whichassgd_(INPAR::SCATRA::assgd_artificial),
    tau_gp_(false),
    mat_gp_(false),
    tau_value_(0.),
    probnum_(0),
    semiimplicit_(false),
    // we have to know the time parameters here to check for illegal combinations
    scatraparatimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname))
{
  return;
}


/*----------------------------------------------------------------------*
 | set parameters                                            ehrl 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterStd::SetParameters(
    Teuchos::ParameterList& parameters   //!< parameter list
    )
{
  // set ale case
  is_ale_ = parameters.get<bool>("isale",false);

  // set flag for conservative form
  const INPAR::SCATRA::ConvForm convform =
    DRT::INPUT::get<INPAR::SCATRA::ConvForm>(parameters, "convform");

  is_conservative_ = false;
  if (convform ==INPAR::SCATRA::convform_conservative) is_conservative_ = true;

  // flag for writing the flux vector fields
  calcflux_domain_ =  DRT::INPUT::get<INPAR::SCATRA::FluxType>(parameters, "calcflux_domain");

  //! vector containing ids of scalars for which flux vectors are calculated
  if (calcflux_domain_ != INPAR::SCATRA::flux_none)
    writefluxids_ =  parameters.get<Teuchos::RCP<std::vector<int> > >("writeflux_ids");

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = parameters.sublist("stabilization");

  // get definition for stabilization parameter tau
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(stablist,"DEFINITION_TAU");

  // set correct stationary definition for stabilization parameter automatically
  // and ensure that exact stabilization parameter is only used in stationary case
  if (scatraparatimint_->IsStationary())
  {
    if (whichtau_ == INPAR::SCATRA::tau_taylor_hughes_zarins)
      whichtau_ = INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_franca_valentin)
      whichtau_ = INPAR::SCATRA::tau_franca_valentin_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_shakib_hughes_codina)
      whichtau_ = INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_codina)
      whichtau_ = INPAR::SCATRA::tau_codina_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_franca_madureira_valentin)
      whichtau_ = INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt;
  }
  else
  {
    if (whichtau_ == INPAR::SCATRA::tau_exact_1d)
      dserror("exact stabilization parameter only available for stationary case");
  }

  if (whichtau_ == INPAR::SCATRA::tau_numerical_value)
    tau_value_ = parameters.sublist("stabilization").get<double>("TAU_VALUE");

  // get characteristic element length for stabilization parameter definition
  charelelength_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLength>(stablist,"CHARELELENGTH");

  // set (sign) factor for diffusive and reactive stabilization terms
  // (factor is zero for SUPG) and overwrite tau definition when there
  // is no stabilization
  stabtype_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(stablist,"STABTYPE");
  switch(stabtype_)
  {
  case INPAR::SCATRA::stabtype_no_stabilization:
    whichtau_ = INPAR::SCATRA::tau_zero;
    break;
  case INPAR::SCATRA::stabtype_SUPG:
    diffreastafac_ = 0.0;
    break;
  case INPAR::SCATRA::stabtype_GLS:
    diffreastafac_ = 1.0;
    break;
  case INPAR::SCATRA::stabtype_USFEM:
    diffreastafac_ = -1.0;
    break;
  case INPAR::SCATRA::stabtype_hdg_centered:
  case INPAR::SCATRA::stabtype_hdg_upwind:
    if(whichtau_ != INPAR::SCATRA::tau_numerical_value or tau_value_ <= 0.0)
      dserror("Wrong definition for tau for hdg stabilization, only tau_numerical_value is allowed with tau>0");
    break;
  default:
    dserror("unknown definition for stabilization parameter");
    break;
  }

  // set flags for subgrid-scale velocity and all-scale subgrid-diffusivity term
  // (default: "false" for both flags)
  sgvel_ = DRT::INPUT::IntegralValue<int>(stablist,"SUGRVEL");
  assgd_ = DRT::INPUT::IntegralValue<int>(stablist,"ASSUGRDIFF");

  // select type of all-scale subgrid diffusivity if included
  whichassgd_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::AssgdType>(stablist,"DEFINITION_ASSGD");

  // set flags for potential evaluation of tau and material law at int. point
  const INPAR::SCATRA::EvalTau tauloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(stablist,"EVALUATION_TAU");
  tau_gp_ = (tauloc == INPAR::SCATRA::evaltau_integration_point); // set true/false

  const INPAR::SCATRA::EvalMat matloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist,"EVALUATION_MAT");
  mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point); // set true/false

  // check for illegal combinations
  if (sgvel_ or assgd_)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
     dserror("Evaluation of material and stabilization parameters need to be done at the integration points if subgrid-scale velocity is included!");
  }

  // get quantities for finite difference check
  fdcheck_ = DRT::INPUT::get<INPAR::SCATRA::FDCheck>(parameters,"fdcheck");
  fdcheckeps_ = parameters.get<double>("fdcheckeps");
  fdchecktol_ = parameters.get<double>("fdchecktol");

  // set flag for use of spherical coordinates
  sphericalcoords_ = parameters.get<bool>("sphericalcoords",false);

  // set problem number
  probnum_ = parameters.get<int>("probnum",0);

  // set evaluation type
  semiimplicit_ = parameters.get<bool>("semiimplicit",false);

  return;
}
