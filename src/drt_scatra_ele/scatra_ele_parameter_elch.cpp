/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.cpp

\brief Setting of elch scatra parameter for element evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_parameter_elch.H"

#include "../drt_lib/drt_dserror.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch* DRT::ELEMENTS::ScaTraEleParameterElch::Instance( bool create )
{
  static ScaTraEleParameterElch* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleParameterElch();
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
void DRT::ELEMENTS::ScaTraEleParameterElch::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch()
  : DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter(),
  elchtype_(INPAR::ELCH::elchtype_undefined),
  equpot_(INPAR::ELCH::equpot_undefined),
  frt_(0.0),
  cursolvar_(false),
  diffusioncoefbased_(true),
  newmanconsta_(2.0),
  newmanconstb_(-2.0),
  newmanconstc_(-1.0)
{
  return;
}

//----------------------------------------------------------------------*
//  set elch-specific parameters                             ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::SetElementElchScaTraParameter(
  Teuchos::ParameterList& params,
  int myrank )
{
  // type of elch problem
  elchtype_ = DRT::INPUT::get<INPAR::ELCH::ElchType>(params, "elchtype");

  // type of closing equation for electric potential
  equpot_ = DRT::INPUT::get<INPAR::ELCH::EquPot>(params, "equpot");

  // get parameter F/RT
  frt_ = params.get<double>("frt");

  // safety check - only stabilization of SUPG-type available
  if ((stabtype_ !=INPAR::SCATRA::stabtype_no_stabilization) and (stabtype_ !=INPAR::SCATRA::stabtype_SUPG))
    dserror("Only SUPG-type stabilization available for ELCH.");

  // safety check - only stabilization of SUPG-type available
  if ((mat_gp_ == false) and (elchtype_ == INPAR::ELCH::elchtype_diffcond))
    dserror("Since most of the materials of the Diffusion-conduction formulation depend on the concentration,\n"
            "an evaluation of the material at the element center is disabled.");

  return;
}


//-----------------------------------------------------------------------------*
//  set elch diffusion-conduction formulation specific parameters   ehrl 04/10 |
//-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::SetElementElchDiffCondScaTraParameter(
  Teuchos::ParameterList& params,
  int myrank )
{
  // set diffusion-conduction sepecific parameter only when diffusion-conduction formulation is used
  if(elchtype_ == INPAR::ELCH::elchtype_diffcond)
  {
    Teuchos::ParameterList& diffcondparams = params.sublist("DIFFCOND");

    // flag if current is used as a solution variable
    cursolvar_ = DRT::INPUT::IntegralValue<int>(diffcondparams,"CURRENT_SOLUTION_VAR");

    // mat_diffcond: flag if diffusion potential is based on diffusion coefficients or transference number
    diffusioncoefbased_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams,"MAT_DIFFCOND_DIFFBASED");

    // switch for dilute and concentrated solution theory (diffusion potential in current equation):
    //    A          B
    //   |--|  |----------|
    //   z_1 + (z_2 - z_1) t_1
    // ------------------------ (RT/F kappa 1/c_k grad c_k)
    //      z_1 z_2
    //     |________|
    //         C
    newmanconsta_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_A");
    newmanconstb_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_B");
    newmanconstc_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_C");

    // safety checks - no stabilization for diffusion-conduction formulation
    if(stabtype_ !=INPAR::SCATRA::stabtype_no_stabilization)
      dserror("No stabilization available for the diffusion-conduction formulation \n"
              "since we had no problems so far.");

    if(whichtau_ != INPAR::SCATRA::tau_zero)
      dserror("No stabilization available for the diffusion-conduction formulation \n"
              "since we had no problems so far.");

    if (mat_gp_ == false and tau_gp_==false)
      dserror("Since most of the materials of the Diffusion-conduction formulation depend on the concentration,\n"
              "an evaluation of the material (and the potential stabilization parameter) at the element center is disabled.");
  } // if(elchtype_ == INPAR::ELCH::elchtype_diffcond)

  return;
}


