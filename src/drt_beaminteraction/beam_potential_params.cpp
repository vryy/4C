/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_potential_params.H"
#include "beam_potential_runtime_vtk_output_params.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

// Todo get rid of this as soon as historic dependency on INPAR::BEAMCONTACT::OctreeType is fully
// gone
#include "../drt_inpar/inpar_beamcontact.H"

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
      num_GPs_(-1),
      useFAD_(false),
      choice_master_slave_(INPAR::BEAMPOTENTIAL::MasterSlaveChoice::choice_master_slave_vague),
      vtk_output_(false),
      params_runtime_vtk_BTB_potential_(Teuchos::null)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::Init()
{
  issetup_ = false;

  // Teuchos parameter list for beam potential-based interactions
  const Teuchos::ParameterList& beam_potential_params_list =
      DRT::Problem::Instance()->BeamPotentialParams();

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
      dserror(
          "number of potential law prefactors does not match number of potential law exponents."
          " Check your input file!");

    for (unsigned int i = 0; i < pot_law_exponents_->size(); ++i)
      if (pot_law_exponents_->at(i) <= 0)
        dserror(
            "only positive values are allowed for potential law exponent."
            " Check your input file");
  }

  /****************************************************************************/
  strategy_ = DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialStrategy>(
      beam_potential_params_list, "STRATEGY");

  if (strategy_ == INPAR::BEAMPOTENTIAL::strategy_vague)
    dserror("You must specify a strategy to be used to evaluate beam interaction potential!");

  /****************************************************************************/
  potential_type_ = DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialType>(
      beam_potential_params_list, "BEAMPOTENTIAL_TYPE");

  if (potential_type_ == INPAR::BEAMPOTENTIAL::beampot_vague)
    dserror("You must specify the type of the specified beam interaction potential!");

  if (potential_type_ == INPAR::BEAMPOTENTIAL::beampot_surf and
      strategy_ != INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox)
  {
    dserror("Surface interaction is not implemented for this strategy yet!");
  }

  /****************************************************************************/
  cutoff_radius_ = beam_potential_params_list.get<double>("CUTOFF_RADIUS");

  if (cutoff_radius_ != -1.0 and cutoff_radius_ <= 0.0)
    dserror("Invalid cutoff radius! Must be positive value or -1 to deactivate.");

  /****************************************************************************/
  regularization_type_ =
      DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialRegularizationType>(
          beam_potential_params_list, "REGULARIZATION_TYPE");

  if ((regularization_type_ != INPAR::BEAMPOTENTIAL::regularization_none and
          strategy_ == INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox) or
      (regularization_type_ == INPAR::BEAMPOTENTIAL::regularization_constant and
          strategy_ == INPAR::BEAMPOTENTIAL::strategy_singlelengthspec_smallsepapprox))
  {
    dserror(
        "This kind of regularization of the force law is not implemented for this strategy yet!");
  }

  /****************************************************************************/
  regularization_separation_ = beam_potential_params_list.get<double>("REGULARIZATION_SEPARATION");

  if (regularization_type_ != INPAR::BEAMPOTENTIAL::regularization_none and
      regularization_separation_ <= 0.0)
  {
    dserror(
        "Invalid regularization separation! Must be a positive value since force law "
        "is not defined for separations <= 0!");
  }

  /****************************************************************************/
  num_integration_segments_ = beam_potential_params_list.get<int>("NUM_INTEGRATION_SEGMENTS");

  if (num_integration_segments_ <= 0)
    dserror("Invalid number of integration segments per element!");

  /****************************************************************************/
  num_GPs_ = beam_potential_params_list.get<int>("NUM_GAUSSPOINTS");

  if (num_GPs_ <= 0) dserror("Invalid number of Gauss points per integration segment!");

  /****************************************************************************/
  useFAD_ = DRT::INPUT::IntegralValue<int>(beam_potential_params_list, "AUTOMATIC_DIFFERENTIATION");

  /****************************************************************************/
  choice_master_slave_ = Teuchos::getIntegralValue<INPAR::BEAMPOTENTIAL::MasterSlaveChoice>(
      beam_potential_params_list, "CHOICE_MASTER_SLAVE");

  if (choice_master_slave_ == INPAR::BEAMPOTENTIAL::MasterSlaveChoice::choice_master_slave_vague)
  {
    dserror("Invalid choice of master and slave!");
  }

  /****************************************************************************/
  // check for vtk output which is to be handled by an own writer object
  vtk_output_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_potential_params_list.sublist("RUNTIME VTK OUTPUT"), "VTK_OUTPUT_BEAM_POTENTIAL");

  // create and initialize parameter container object for runtime vtk output
  if (vtk_output_)
  {
    params_runtime_vtk_BTB_potential_ =
        Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams);

    params_runtime_vtk_BTB_potential_->Init(
        beam_potential_params_list.sublist("RUNTIME VTK OUTPUT"));
    params_runtime_vtk_BTB_potential_->Setup();
  }


  /****************************************************************************/
  // safety checks for currently unsupported parameter settings
  /****************************************************************************/

  // outdated: octtree for search of potential-based interaction pairs
  if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(
          beam_potential_params_list, "BEAMPOT_OCTREE") != INPAR::BEAMCONTACT::boct_none)
  {
    dserror("Octree-based search for potential-based beam interactions is deprecated!");
  }

  // outdated: flags to indicate, if beam-to-solid or beam-to-sphere potential-based interaction is
  // applied
  if (DRT::INPUT::IntegralValue<int>(beam_potential_params_list, "BEAMPOT_BTSOL") != 0)
  {
    dserror(
        "The flag BEAMPOT_BTSOL is outdated! remove them as soon"
        "as old beamcontact_manager is gone!");
  }

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::Setup()
{
  ThrowErrorIfNotInit();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::ThrowErrorIfNotInitAndSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::ThrowErrorIfNotInit() const
{
  if (!IsInit()) dserror("Init() has not been called, yet!");
}
