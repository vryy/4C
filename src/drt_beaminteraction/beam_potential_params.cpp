/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_potential_params.cpp

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_potential_params.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

// Todo get rid of this as soon as historic dependency on INPAR::BEAMCONTACT::OctreeType is fully gone
#include "../drt_inpar/inpar_beamcontact.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamPotentialParams::BeamPotentialParams()
: isinit_(false),
  issetup_(false),
  pot_law_exponents_(Teuchos::null),
  pot_law_prefactors_(Teuchos::null),
  potential_type_(INPAR::BEAMPOTENTIAL::beampot_vague)
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
        Teuchos::getNumericStringParameter(beam_potential_params_list,"POT_LAW_EXPONENT") );
    std::string word;
    char* input;
    while (PL >> word)
      pot_law_exponents_->push_back(std::strtod(word.c_str(), &input));
  }
  {
    std::istringstream PL(
        Teuchos::getNumericStringParameter(beam_potential_params_list,"POT_LAW_PREFACTOR"));
    std::string word;
    char* input;
    while (PL >> word)
      pot_law_prefactors_->push_back(std::strtod(word.c_str(), &input));
  }
  if (!pot_law_prefactors_->empty())
  {
    if (pot_law_prefactors_->size() != pot_law_exponents_->size())
      dserror("number of potential law prefactors does not match number of potential law exponents."
          " Check your input file!");

    for (unsigned int i=0; i<pot_law_exponents_->size(); ++i)
      if (pot_law_exponents_->at(i) <= 0)
        dserror("only positive values are allowed for potential law exponent."
            " Check your input file");
  }

  /****************************************************************************/
  potential_type_ = DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialType>(
      beam_potential_params_list,"BEAMPOTENTIAL_TYPE");

  if (potential_type_ == INPAR::BEAMPOTENTIAL::beampot_vague)
    dserror("You must specify the type of the specified beam interaction potential!");



  /****************************************************************************/
  // safety checks for currently unsupported parameter settings
  /****************************************************************************/
  // read cutoff radius for search of potential-based interaction pairs
  if ( beam_potential_params_list.get<double>("CUTOFFRADIUS") != -1.0 )
    dserror("The parameter CUTOFFRADIUS in beam potential input section is deprecated!"
        " Choose your cutoff globally in MESHFREE section and remove this parameter as"
        " soon as old code is completely gone");

  // outdated: octtree for search of potential-based interaction pairs
  if ( DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(
      beam_potential_params_list,"BEAMPOT_OCTREE") != INPAR::BEAMCONTACT::boct_none )
  {
    dserror("Octree-based search for potential-based beam interactions is deprecated!");
  }

  // outdated: flags to indicate, if beam-to-solid or beam-to-sphere potential-based interaction is applied
  if ( DRT::INPUT::IntegralValue<int>(beam_potential_params_list,"BEAMPOT_BTSOL") != 0)
  {
    dserror("The flag BEAMPOT_BTSOL is outdated! remove them as soon"
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
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialParams::ThrowErrorIfNotInit() const
{
  if (!IsInit())
    dserror("Init() has not been called, yet!");
}
