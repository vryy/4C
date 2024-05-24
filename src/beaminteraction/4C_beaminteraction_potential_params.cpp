/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_beaminteraction_potential_params.hpp"

#include "4C_beaminteraction_potential_runtime_visualization_output_params.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamPotentialParams::BeamPotentialParams()
    : isinit_(false),
      issetup_(false),
      pot_law_exponents_(Teuchos::null),
      pot_law_prefactors_(Teuchos::null),
      potential_type_(INPAR::BEAMPOTENTIAL::beampot_vague),
      strategy_(INPAR::BEAMPOTENTIAL::strategy_vague),
      cutoff_radius_(0.0),
      regularization_type_(INPAR::BEAMPOTENTIAL::regularization_none),
      regularization_separation_(0.0),
      num_integration_segments_(-1),
      num_gp_s_(-1),
      use_fad_(false),
      choice_master_slave_(INPAR::BEAMPOTENTIAL::MasterSlaveChoice::choice_master_slave_vague),
      visualization_output_(false),
      params_runtime_visualization_output_btb_potential_(Teuchos::null)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::Init(const double restart_time)
{
  issetup_ = false;

  // Teuchos parameter list for beam potential-based interactions
  const Teuchos::ParameterList& beam_potential_params_list =
      GLOBAL::Problem::Instance()->beam_potential_params();

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/

  pot_law_prefactors_ = Teuchos::rcp(new std::vector<double>);
  pot_law_exponents_ = Teuchos::rcp(new std::vector<double>);
  pot_law_prefactors_->clear();
  pot_law_exponents_->clear();
  // read potential law parameters from input and check
  {
    std::istringstream PL(
        Teuchos::getNumericStringParameter(beam_potential_params_list, "POT_LAW_EXPONENT"));
    std::string word;
    char* input;
    while (PL >> word) pot_law_exponents_->push_back(std::strtod(word.c_str(), &input));
  }
  {
    std::istringstream PL(
        Teuchos::getNumericStringParameter(beam_potential_params_list, "POT_LAW_PREFACTOR"));
    std::string word;
    char* input;
    while (PL >> word) pot_law_prefactors_->push_back(std::strtod(word.c_str(), &input));
  }
  if (!pot_law_prefactors_->empty())
  {
    if (pot_law_prefactors_->size() != pot_law_exponents_->size())
      FOUR_C_THROW(
          "number of potential law prefactors does not match number of potential law exponents."
          " Check your input file!");

    for (unsigned int i = 0; i < pot_law_exponents_->size(); ++i)
      if (pot_law_exponents_->at(i) <= 0)
        FOUR_C_THROW(
            "only positive values are allowed for potential law exponent."
            " Check your input file");
  }

  /****************************************************************************/
  strategy_ = CORE::UTILS::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialStrategy>(
      beam_potential_params_list, "STRATEGY");

  if (strategy_ == INPAR::BEAMPOTENTIAL::strategy_vague)
    FOUR_C_THROW("You must specify a strategy to be used to evaluate beam interaction potential!");

  /****************************************************************************/
  potential_type_ = CORE::UTILS::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialType>(
      beam_potential_params_list, "BEAMPOTENTIAL_TYPE");

  if (potential_type_ == INPAR::BEAMPOTENTIAL::beampot_vague)
    FOUR_C_THROW("You must specify the type of the specified beam interaction potential!");

  if (potential_type_ == INPAR::BEAMPOTENTIAL::beampot_surf and
      strategy_ != INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox)
  {
    FOUR_C_THROW("Surface interaction is not implemented for this strategy yet!");
  }

  /****************************************************************************/
  cutoff_radius_ = beam_potential_params_list.get<double>("CUTOFF_RADIUS");

  if (cutoff_radius_ != -1.0 and cutoff_radius_ <= 0.0)
    FOUR_C_THROW("Invalid cutoff radius! Must be positive value or -1 to deactivate.");

  /****************************************************************************/
  regularization_type_ =
      CORE::UTILS::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialRegularizationType>(
          beam_potential_params_list, "REGULARIZATION_TYPE");

  if ((regularization_type_ != INPAR::BEAMPOTENTIAL::regularization_none and
          strategy_ == INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox) or
      (regularization_type_ == INPAR::BEAMPOTENTIAL::regularization_constant and
          strategy_ == INPAR::BEAMPOTENTIAL::strategy_singlelengthspec_smallsepapprox))
  {
    FOUR_C_THROW(
        "This kind of regularization of the force law is not implemented for this strategy yet!");
  }

  /****************************************************************************/
  regularization_separation_ = beam_potential_params_list.get<double>("REGULARIZATION_SEPARATION");

  if (regularization_type_ != INPAR::BEAMPOTENTIAL::regularization_none and
      regularization_separation_ <= 0.0)
  {
    FOUR_C_THROW(
        "Invalid regularization separation! Must be a positive value since force law "
        "is not defined for separations <= 0!");
  }

  /****************************************************************************/
  num_integration_segments_ = beam_potential_params_list.get<int>("NUM_INTEGRATION_SEGMENTS");

  if (num_integration_segments_ <= 0)
    FOUR_C_THROW("Invalid number of integration segments per element!");

  /****************************************************************************/
  num_gp_s_ = beam_potential_params_list.get<int>("NUM_GAUSSPOINTS");

  if (num_gp_s_ <= 0) FOUR_C_THROW("Invalid number of Gauss points per integration segment!");

  /****************************************************************************/
  use_fad_ =
      CORE::UTILS::IntegralValue<int>(beam_potential_params_list, "AUTOMATIC_DIFFERENTIATION");

  /****************************************************************************/
  choice_master_slave_ = Teuchos::getIntegralValue<INPAR::BEAMPOTENTIAL::MasterSlaveChoice>(
      beam_potential_params_list, "CHOICE_MASTER_SLAVE");

  if (choice_master_slave_ == INPAR::BEAMPOTENTIAL::MasterSlaveChoice::choice_master_slave_vague)
  {
    FOUR_C_THROW("Invalid choice of master and slave!");
  }

  /****************************************************************************/
  // check for vtk output which is to be handled by an own writer object
  visualization_output_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_potential_params_list.sublist("RUNTIME VTK OUTPUT"), "VTK_OUTPUT_BEAM_POTENTIAL");

  // create and initialize parameter container object for runtime output
  if (visualization_output_)
  {
    params_runtime_visualization_output_btb_potential_ =
        Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams(restart_time));

    params_runtime_visualization_output_btb_potential_->Init(
        beam_potential_params_list.sublist("RUNTIME VTK OUTPUT"));
    params_runtime_visualization_output_btb_potential_->Setup();
  }


  /****************************************************************************/
  // safety checks for currently unsupported parameter settings
  /****************************************************************************/

  // outdated: octtree for search of potential-based interaction pairs
  if (CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(
          beam_potential_params_list, "BEAMPOT_OCTREE") != INPAR::BEAMCONTACT::boct_none)
  {
    FOUR_C_THROW("Octree-based search for potential-based beam interactions is deprecated!");
  }

  // outdated: flags to indicate, if beam-to-solid or beam-to-sphere potential-based interaction is
  // applied
  if (CORE::UTILS::IntegralValue<int>(beam_potential_params_list, "BEAMPOT_BTSOL") != 0)
  {
    FOUR_C_THROW(
        "The flag BEAMPOT_BTSOL is outdated! remove them as soon"
        "as old beamcontact_manager is gone!");
  }

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::Setup()
{
  throw_error_if_not_init();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::throw_error_if_not_init_and_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::throw_error_if_not_init() const
{
  if (!is_init()) FOUR_C_THROW("Init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
