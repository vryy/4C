/*----------------------------------------------------------------------------*/
/*! \file
\brief data container holding all crosslinking input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------------*/

#include "../drt_beaminteraction/crosslinking_params.H"

#include "../drt_mat/crosslinkermat.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::CrosslinkingParams::CrosslinkingParams()
    : isinit_(false), issetup_(false), viscosity_(0.0), kt_(0.0), deltatime_(0.0), init_box_(true)
{
  maxnum_init_crosslinker_pertype_.clear();
  numcrosslinkerpertype_.clear();
  matcrosslinkerpertype_.clear();
  linkertypes_.clear();
  max_num_bonds_per_filament_bspot_.clear();
  filamentbspotintervalglobal_.clear();
  filamentbspotintervallocal_.clear();
  filamentbspotrangeglobal_.clear();
  filamentbspotrangelocal_.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::CrosslinkingParams::Init(STR::TIMINT::BaseDataGlobalState const& gstate)
{
  issetup_ = false;

  const Teuchos::ParameterList& crosslinking_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("CROSSLINKING");

  // viscosity
  viscosity_ = crosslinking_params_list.get<double>("VISCOSITY");

  // thermal energy
  kt_ = crosslinking_params_list.get<double>("KT");

  // time step for stochastic events concering crosslinking
  deltatime_ = crosslinking_params_list.get<double>("TIMESTEP");

  init_box_.PutScalar(1.0e12);
  std::istringstream init_box_stream(
      Teuchos::getNumericStringParameter(crosslinking_params_list, "INIT_LINKER_BOUNDINGBOX"));
  for (int col = 0; col < 2; ++col)
  {
    for (int row = 0; row < 3; ++row)
    {
      double value = 1.0e12;
      if (init_box_stream >> value)
        init_box_(row, col) = value;
      else
        dserror(
            " Specify six values for bounding box in three dimensional problem."
            " Fix your input file.");
    }
  }

  bool feasibleboxinput = true;
  for (int col = 0; col < 2; ++col)
    for (int row = 0; row < 3; ++row)
      if (init_box_(row, col) > 1.0e11) feasibleboxinput = false;

  if (not feasibleboxinput)
  {
    std::istringstream pbb_stream(Teuchos::getNumericStringParameter(
        DRT::Problem::Instance()->BinningStrategyParams(), "BOUNDINGBOX"));
    for (int col = 0; col < 2; ++col)
    {
      for (int row = 0; row < 3; ++row)
      {
        double value = 1.0e12;
        if (pbb_stream >> value)
          init_box_(row, col) = value;
        else
          dserror(
              " Specify six values for bounding box in three dimensional problem."
              " Fix your input file.");
      }
    }
  }

  // safety check
  // todo: maybe make input of time step obligatory
  if (deltatime_ < 0.0)
  {
    deltatime_ = (*gstate.GetDeltaTime())[0];
    if (gstate.GetMyRank() == 0)
      std::cout << " Time step " << (*gstate.GetDeltaTime())[0]
                << " form Structural Dynamic section "
                   "used for crosslinking.\n"
                << std::endl;
  }

  // number of linker in simulation volume
  {
    numcrosslinkerpertype_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list, "NUMCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word) numcrosslinkerpertype_.push_back(std::strtod(word.c_str(), &input));

    // safety check
    for (int i = 0; i < static_cast<int>(numcrosslinkerpertype_.size()); ++i)
      if (numcrosslinkerpertype_[i] < 0)
        dserror(" negative number of crosslinker does not make sense.");
  }

  // material numbers for crosslinker types
  {
    matcrosslinkerpertype_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list, "MATCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word) matcrosslinkerpertype_.push_back(std::strtod(word.c_str(), &input));

    // safety check
    for (int i = 0; i < static_cast<int>(matcrosslinkerpertype_.size()); ++i)
      if (matcrosslinkerpertype_[i] < 0) dserror(" negative material number does not make sense.");
  }

  // safety check
  if (numcrosslinkerpertype_.size() != matcrosslinkerpertype_.size())
    dserror("number of crosslinker types does not fit number of assigned materials");

  // compute number of linker types
  linkertypes_.clear();
  for (unsigned int type_i = 0; type_i < matcrosslinkerpertype_.size(); ++type_i)
  {
    if (not(std::find(linkertypes_.begin(), linkertypes_.end(), matcrosslinkerpertype_[type_i]) !=
            linkertypes_.end()))
      linkertypes_.push_back(Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(
          MAT::Material::Factory(matcrosslinkerpertype_[type_i]))
                                 ->LinkerType());
  }

  // number of initially set linker
  {
    maxnum_init_crosslinker_pertype_.clear();
    std::vector<int> maxnuminitcrosslinkerpertype;
    std::istringstream PL(Teuchos::getNumericStringParameter(
        crosslinking_params_list, "MAXNUMINITCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word) maxnuminitcrosslinkerpertype.push_back(std::strtod(word.c_str(), &input));

    if (maxnuminitcrosslinkerpertype.size() > 1 or maxnuminitcrosslinkerpertype[0] != 0)
    {
      // safety checks
      for (int i = 0; i < static_cast<int>(maxnuminitcrosslinkerpertype.size()); ++i)
        if (maxnuminitcrosslinkerpertype[i] < 0)
          dserror(" negative number of crosslinker does not make sense.");
      if (maxnuminitcrosslinkerpertype.size() != numcrosslinkerpertype_.size())
        dserror("number of initial set crosslinker types does not fit number of crosslinker types");

      for (int i = 0; i < static_cast<int>(maxnuminitcrosslinkerpertype.size()); ++i)
        maxnum_init_crosslinker_pertype_[matcrosslinkerpertype_[i]] =
            maxnuminitcrosslinkerpertype[i];
    }
  }

  // maximal number of bonds per filament binding spot
  {
    max_num_bonds_per_filament_bspot_.clear();
    std::istringstream PL(Teuchos::getNumericStringParameter(
        crosslinking_params_list, "MAXNUMBONDSPERFILAMENTBSPOT"));
    std::string word;
    char* input;
    int count = 0;
    while (PL >> word)
    {
      max_num_bonds_per_filament_bspot_[linkertypes_[count]] = std::strtod(word.c_str(), &input);
      if (max_num_bonds_per_filament_bspot_.at(linkertypes_[count]) < 0)
        dserror(" Choose a number of bonds per filament binding spot type >= 0. ");
      ++count;
    }

    if (max_num_bonds_per_filament_bspot_.size() != linkertypes_.size())
      dserror(" Num linker types %i does not match num input for MAXNUMBONDSPERFILAMENTBSPOT %i. ",
          linkertypes_.size(), max_num_bonds_per_filament_bspot_.size());
  }

  // distance between the two binding spots on each filament the same
  {
    filamentbspotintervalglobal_.clear();
    std::istringstream PL(Teuchos::getNumericStringParameter(
        crosslinking_params_list, "FILAMENTBSPOTINTERVALGLOBAL"));
    std::string word;
    char* input;
    int count = 0;
    while (PL >> word)
    {
      filamentbspotintervalglobal_[linkertypes_[count]] = std::strtod(word.c_str(), &input);
      ++count;
    }
  }

  // distance between the two binding spots on a filament as percentage of current filament
  // reference length
  {
    filamentbspotintervallocal_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list, "FILAMENTBSPOTINTERVALLOCAL"));
    std::string word;
    char* input;
    int count = 0;
    while (PL >> word)
    {
      filamentbspotintervallocal_[linkertypes_[count]] = std::strtod(word.c_str(), &input);
      ++count;
    }
  }

  if (linkertypes_.size() != filamentbspotintervalglobal_.size() and
      linkertypes_.size() != filamentbspotintervallocal_.size())
    dserror("You need to specify filament binding spots for all your linker types");

  // safety checks for feasibility of input
  if (filamentbspotintervalglobal_.size() == filamentbspotintervallocal_.size())
  {
    for (auto const& iter : filamentbspotintervalglobal_)
    {
      // safety feasibility checks
      if (iter.second <= 0.0 and not(filamentbspotintervallocal_.at(iter.first) > 0.0 and
                                     filamentbspotintervallocal_.at(iter.first) <= 1.0))
        dserror(
            " Choose realistic value for FILAMENTBSPOTINTERVAL (i.e. distance between "
            "two binding spots on a filament) in input file. ");
      if (iter.second > 0.0 and filamentbspotintervallocal_.at(iter.first) > 0.0)
        dserror(" You can only set either a global or a local filament binding spot interval");
    }
  }

  // start and end arc parameter for binding spots on a filament
  {
    filamentbspotrangeglobal_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list, "FILAMENTBSPOTRANGEGLOBAL"));
    std::string word;
    char* input;
    int count = 0;
    while (PL >> word)
    {
      std::pair<double, double> pair;
      pair.first = std::strtod(word.c_str(), &input);
      if (PL >> word)
        pair.second = std::strtod(word.c_str(), &input);
      else
        dserror("Filament binding spot range needs to be specified via two values");
      filamentbspotrangeglobal_[linkertypes_[count]] = pair;

      if (pair.first > 0.0 and pair.second > 0.0 and (pair.first > pair.second))
        dserror(" lower bound > upper bound, fix FILAMENTBSPOTRANGEGLOBAL in input file ");

      ++count;
    }
  }

  // start and end arc parameter for binding spots on a filament
  {
    filamentbspotrangelocal_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list, "FILAMENTBSPOTRANGELOCAL"));
    std::string word;
    char* input;
    int count = 0;
    while (PL >> word)
    {
      std::pair<double, double> pair;
      pair.first = std::strtod(word.c_str(), &input);
      if (PL >> word)
        pair.second = std::strtod(word.c_str(), &input);
      else
        dserror("Filament binding spot range needs to be specified via two values");
      filamentbspotrangelocal_[linkertypes_[count]] = pair;

      if (pair.first > 0.0 and pair.second > 0.0 and (pair.first > pair.second))
        dserror(" lower bound > upper bound, fix FILAMENTBSPOTRANGEGLOCAL in input file ");
      if (pair.first > 1.0 or pair.second > 1.0)
        dserror("values > 1.0 do not make sense for local filament binding spot range");

      ++count;
    }
  }

  if (linkertypes_.size() != filamentbspotrangeglobal_.size() and
      linkertypes_.size() != filamentbspotrangelocal_.size())
    dserror("You need to specify filament binding spots for all your linker types");

  // safety checks for feasibility of input
  if (filamentbspotrangeglobal_.size() == filamentbspotrangelocal_.size())
  {
    for (auto const& iter : filamentbspotrangeglobal_)
    {
      if (filamentbspotrangelocal_.at(iter.first).first > 0.0 and
          filamentbspotrangelocal_.at(iter.first).second > 0.0 and
          (iter.second.first > 0.0 or iter.second.second > 0.0))
        dserror("either local or global binding spot range can be specified");
    }
  }


  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::CrosslinkingParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}
