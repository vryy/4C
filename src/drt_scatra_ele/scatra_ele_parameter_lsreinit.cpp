/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_lsreinit.cpp

\brief Setting of scatra parameter for element evaluation of reinitialization equation

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_parameter_lsreinit.H"
#include "../drt_lib/drt_dserror.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterLsReinit* DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance( bool create )
{
  static ScaTraEleParameterLsReinit* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleParameterLsReinit();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterLsReinit::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterLsReinit::ScaTraEleParameterLsReinit()
  : DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter(),
    signtype_(INPAR::SCATRA::signtype_nonsmoothed),
    charelelengthreinit_(INPAR::SCATRA::root_of_volume_reinit),
    interfacethicknessfac_(1.0),
    useprojectedreinitvel_(false),
    linform_(INPAR::SCATRA::fixed_point),
    artdiff_(INPAR::SCATRA::artdiff_none)
{
}


//----------------------------------------------------------------------*
//  set reinitialization specific parameters            rasthofer 12/13 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterLsReinit::SetElementLsReinitScaTraParameter(
  Teuchos::ParameterList& params)
{
  // get reinitialization parametre list
  Teuchos::ParameterList& reinitlist = params.sublist("REINITIALIZATION");
  // get signum function
  signtype_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::SmoothedSignType>(reinitlist, "SMOOTHED_SIGN_TYPE");
  // characteristic element length for signum function
  charelelengthreinit_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLengthReinit>(reinitlist, "CHARELELENGTHREINIT");
  // interface thickness for signum function
  interfacethicknessfac_ = reinitlist.get<double>("INTERFACE_THICKNESS");
  // form of linearization for nonlinear terms
  linform_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::LinReinit>(reinitlist, "LINEARIZATIONREINIT");

  // set form of velocity evaluation
  INPAR::SCATRA::VelReinit velreinit = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelReinit>(reinitlist, "VELREINIT");
  if (velreinit == INPAR::SCATRA::vel_reinit_node_based) useprojectedreinitvel_ = true;

  // get definition for stabilization parameter tau (get definition of tau before stabtype is set,
  // otherwise tau_zero in case of stabtype_no_stabilization (calc initial phidt) will be overwritten)
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(reinitlist,"DEFINITION_TAU_REINIT");

  // overwrite stabilization type
  stabtype_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(reinitlist,"STABTYPEREINIT");
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
  default:
    dserror("unknown definition for stabilization parameter");
    break;
  }

  // set flags for subgrid-scale velocity and artificial diffusion term
  // (default: "false" for both flags)
  sgvel_ = false;
  artdiff_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::ArtDiff>(reinitlist,"ARTDIFFREINIT");

  // select type of artificial diffusion if included
  whichassgd_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::AssgdType>(reinitlist,"DEFINITION_ARTDIFFREINIT");

  // check for illegal combinations
  if (artdiff_ != INPAR::SCATRA::artdiff_none)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
     dserror("Evaluation of material and stabilization parameters need to be done at the integration points for reinitialization");
    // due to artificial diff
  }

  return;
}

