/*----------------------------------------------------------------------------*/
/*! \file

\brief data container holding all beam to beam contact input parameters

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_beaminteraction_beam_to_beam_contact_params.hpp"

#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamContactParams::BeamToBeamContactParams()
    : isinit_(false),
      issetup_(false),
      strategy_(Inpar::BEAMCONTACT::bstr_none),
      penalty_law_(Inpar::BEAMCONTACT::pl_lp),
      btb_penalty_law_regularization_g0_(-1.0),
      btb_penalty_law_regularization_f0_(-1.0),
      btb_penalty_law_regularization_c0_(-1.0),
      gap_shift_(0.0),
      btb_point_penalty_param_(-1.0),
      btb_line_penalty_param_(-1.0),
      btb_perp_shifting_angle1_(-1.0),
      btb_perp_shifting_angle2_(-1.0),
      btb_parallel_shifting_angle1_(-1.0),
      btb_parallel_shifting_angle2_(-1.0),
      segangle_(-1.0),
      num_integration_intervals_(0),
      btb_basicstiff_gap_(-1.0),
      btb_endpoint_penalty_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactParams::Init()
{
  issetup_ = false;

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_contact_params_list =
      Global::Problem::Instance()->beam_contact_params();

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  strategy_ = Core::UTILS::IntegralValue<Inpar::BEAMCONTACT::Strategy>(
      beam_contact_params_list, "BEAMS_STRATEGY");

  if (strategy_ != Inpar::BEAMCONTACT::bstr_penalty)
    FOUR_C_THROW(
        "currently only a penalty strategy is supported for beam contact"
        " if not using the 'old' beam contact manager!");

  /****************************************************************************/
  penalty_law_ = Core::UTILS::IntegralValue<Inpar::BEAMCONTACT::PenaltyLaw>(
      beam_contact_params_list, "BEAMS_PENALTYLAW");

  /****************************************************************************/
  btb_penalty_law_regularization_g0_ = beam_contact_params_list.get<double>("BEAMS_PENREGPARAM_G0");
  btb_penalty_law_regularization_f0_ = beam_contact_params_list.get<double>("BEAMS_PENREGPARAM_F0");
  btb_penalty_law_regularization_c0_ = beam_contact_params_list.get<double>("BEAMS_PENREGPARAM_C0");

  // Todo check and refine these safety checks
  if (penalty_law_ != Inpar::BEAMCONTACT::pl_lp and penalty_law_ != Inpar::BEAMCONTACT::pl_qp)
  {
    if (btb_penalty_law_regularization_g0_ == -1.0 or btb_penalty_law_regularization_f0_ == -1.0 or
        btb_penalty_law_regularization_c0_ == -1.0)
      FOUR_C_THROW(
          "Regularized penalty law chosen, but not all regularization parameters are set!");
  }

  /****************************************************************************/
  // Todo check this parameter
  gap_shift_ = beam_contact_params_list.get<double>("BEAMS_GAPSHIFTPARAM");

  if (gap_shift_ != 0.0 and penalty_law_ != Inpar::BEAMCONTACT::pl_lpqp)
    FOUR_C_THROW("BEAMS_GAPSHIFTPARAM only possible for penalty law LinPosQuadPen!");

  /****************************************************************************/
  btb_point_penalty_param_ = beam_contact_params_list.get<double>("BEAMS_BTBPENALTYPARAM");

  if (btb_point_penalty_param_ < 0.0)
    FOUR_C_THROW("beam-to-beam point penalty parameter must not be negative!");


  // input parameters required for all-angle-beam contact formulation ...
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_SEGCON"))
  {
    /****************************************************************************/
    btb_line_penalty_param_ = beam_contact_params_list.get<double>("BEAMS_BTBLINEPENALTYPARAM");

    if (btb_line_penalty_param_ < 0.0)
      FOUR_C_THROW(
          "You chose all-angle-beam contact algorithm: thus, beam-to-beam line"
          " penalty parameter must not be negative!");

    /****************************************************************************/
    // Todo find more verbose and expressive naming
    // note: conversion from degrees (input parameter) to radians (class variable) done here!
    btb_perp_shifting_angle1_ =
        beam_contact_params_list.get<double>("BEAMS_PERPSHIFTANGLE1") / 180.0 * M_PI;
    btb_perp_shifting_angle2_ =
        beam_contact_params_list.get<double>("BEAMS_PERPSHIFTANGLE2") / 180.0 * M_PI;

    btb_parallel_shifting_angle1_ =
        beam_contact_params_list.get<double>("BEAMS_PARSHIFTANGLE1") / 180.0 * M_PI;
    btb_parallel_shifting_angle2_ =
        beam_contact_params_list.get<double>("BEAMS_PARSHIFTANGLE2") / 180.0 * M_PI;

    if (btb_perp_shifting_angle1_ < 0.0 or btb_perp_shifting_angle2_ < 0.0 or
        btb_parallel_shifting_angle1_ < 0.0 or btb_parallel_shifting_angle2_ < 0.0)
      FOUR_C_THROW(
          "You chose all-angle-beam contact algorithm: thus, shifting angles for"
          " beam-to-beam contact fade must be >= 0 degrees");

    if (btb_perp_shifting_angle1_ > 0.5 * M_PI or btb_perp_shifting_angle2_ > 0.5 * M_PI or
        btb_parallel_shifting_angle1_ > 0.5 * M_PI or btb_parallel_shifting_angle2_ > 0.5 * M_PI)
      FOUR_C_THROW(
          "You chose all-angle-beam contact algorithm: thus, Shifting angles for"
          " beam-to-beam contact fade must be <= 90 degrees");

    if (btb_parallel_shifting_angle2_ <= btb_perp_shifting_angle1_)
      FOUR_C_THROW("No angle overlap between large-angle and small-angle contact!");

    /****************************************************************************/
    // note: conversion from degrees (input parameter) to radians (class variable) done here!
    segangle_ = beam_contact_params_list.get<double>("BEAMS_SEGANGLE") / 180.0 * M_PI;

    if (segangle_ <= 0.0) FOUR_C_THROW("Segmentation angle must be greater than zero!");

    /****************************************************************************/
    num_integration_intervals_ = beam_contact_params_list.get<int>("BEAMS_NUMINTEGRATIONINTERVAL");

    if (num_integration_intervals_ <= 0)
      FOUR_C_THROW("Number of integration intervals must be greater than zero!");
  }

  /****************************************************************************/
  // Todo check need and usage of this parameter
  btb_basicstiff_gap_ = beam_contact_params_list.get<double>("BEAMS_BASICSTIFFGAP");

  /****************************************************************************/
  btb_endpoint_penalty_ =
      Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_ENDPOINTPENALTY");

  /****************************************************************************/
  // safety checks for currently unsupported parameter settings
  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_NEWGAP"))
    FOUR_C_THROW("BEAMS_NEWGAP currently not supported!");

  /****************************************************************************/
  // for the time being only allow all-angle-beam contact formulation ...
  if (not Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_SEGCON"))
    FOUR_C_THROW(
        "only all-angle-beam contact (BEAMS_SEGCON) formulation tested yet"
        " in new beam interaction framework!");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_DEBUG"))
    FOUR_C_THROW("get rid of this nasty BEAMS_DEBUG flag");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_INACTIVESTIFF"))
    FOUR_C_THROW("get rid of BEAMS_INACTIVESTIFF flag; no longer supported!");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_BTSOL") or
      beam_contact_params_list.get<double>("BEAMS_BTSPENALTYPARAM") != 0.0)
    FOUR_C_THROW("currently only beam-to-(BEAM/SPHERE) contact supported!");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_SMOOTHING") !=
      Inpar::BEAMCONTACT::bsm_none)
    FOUR_C_THROW("BEAMS_SMOOTHING currently not supported!");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_DAMPING") !=
          Inpar::BEAMCONTACT::bd_no or
      beam_contact_params_list.get<double>("BEAMS_DAMPINGPARAM") != -1000.0 or
      beam_contact_params_list.get<double>("BEAMS_DAMPREGPARAM1") != -1000.0 or
      beam_contact_params_list.get<double>("BEAMS_DAMPREGPARAM2") != -1000.0)
    FOUR_C_THROW("BEAMS_DAMPING currently not supported!");

  /****************************************************************************/
  if (beam_contact_params_list.get<double>("BEAMS_MAXDISISCALEFAC") != -1.0 or
      beam_contact_params_list.get<double>("BEAMS_MAXDELTADISSCALEFAC") != -1.0)
    FOUR_C_THROW("BEAMS_MAXDISISCALEFAC and BEAMS_MAXDELTADISSCALEFAC currently not supported!");

  /****************************************************************************/
  if (btb_basicstiff_gap_ != -1.0) FOUR_C_THROW("BEAMS_BASICSTIFFGAP currently not supported!");

  /****************************************************************************/
  if (Core::UTILS::IntegralValue<Inpar::BEAMCONTACT::OctreeType>(
          beam_contact_params_list, "BEAMS_OCTREE") != Inpar::BEAMCONTACT::boct_none or
      Core::UTILS::IntegralValue<int>(beam_contact_params_list, "BEAMS_ADDITEXT") != true or
      beam_contact_params_list.get<int>("BEAMS_TREEDEPTH") != 6 or
      beam_contact_params_list.get<int>("BEAMS_BOXESINOCT") != 8)
    FOUR_C_THROW(
        "you seem to have set a search-related parameter in the beam contact section! "
        "this is not applicable in case of binning!");

  // Todo BEAMS_EXTVAL is missing here

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactParams::Setup()
{
  check_init();

  // empty for now

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
