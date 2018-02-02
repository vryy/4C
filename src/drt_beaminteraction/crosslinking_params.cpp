/*----------------------------------------------------------------------------*/
/*!
\file crosslinking_params.cpp

\brief data container holding all crosslinking input parameters

\level 3

\maintainer Jonas Eichinger, Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "crosslinking_params.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../drt_lib/drt_globalproblem.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::CrosslinkingParams::CrosslinkingParams()
  : isinit_(false),
    issetup_(false),
    viscosity_(0.0),
    kt_(0.0),
    deltatime_(0.0),
    filamentbspotinterval_(-1.0),
    max_num_bonds_per_filament_bspot_(-1.0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::CrosslinkingParams::Init( STR::TIMINT::BaseDataGlobalState const& gstate )
{
  issetup_ = false;

  const Teuchos::ParameterList& crosslinking_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("CROSSLINKING");


  {
    // number of crosslinkers in the simulated volume
    numcrosslinkerpertype_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list,"NUMCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word)
      numcrosslinkerpertype_.push_back( std::strtod( word.c_str(), &input ) );

    // safety check
    for( int i = 0; i < static_cast<int>( numcrosslinkerpertype_.size() ); ++i )
      if ( numcrosslinkerpertype_[i] < 0 )
        dserror(" negative number of crosslinker does not make sense.");
  }

  {
    // material numbers for crosslinker types
    matcrosslinkerpertype_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list,"MATCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word)
      matcrosslinkerpertype_.push_back( std::strtod( word.c_str(), &input ) );

    // safety check
    for( int i = 0; i < static_cast<int>( matcrosslinkerpertype_.size() ); ++i )
      if ( matcrosslinkerpertype_[i] < 0 )
        dserror(" negative material number does not make sense.");
   }

  // safety check
  if ( numcrosslinkerpertype_.size() != matcrosslinkerpertype_.size() )
    dserror("number of crosslinker types does not fit number of assigned materials");

  {
    // number of crosslinkers in the simulated volume
    maxnum_init_crosslinker_pertype_.clear();
    std::vector<int> maxnuminitcrosslinkerpertype;
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list,"MAXNUMINITCROSSLINKERPERTYPE"));
    std::string word;
    char* input;
    while (PL >> word)
      maxnuminitcrosslinkerpertype.push_back( std::strtod( word.c_str(), &input ) );

    if ( maxnuminitcrosslinkerpertype.size() > 1 or maxnuminitcrosslinkerpertype[0] != 0 )
    {
      // safety checks
      for( int i = 0; i < static_cast<int>( maxnuminitcrosslinkerpertype.size() ); ++i )
        if ( maxnuminitcrosslinkerpertype[i] < 0 )
          dserror(" negative number of crosslinker does not make sense.");
      if ( maxnuminitcrosslinkerpertype.size() != numcrosslinkerpertype_.size() )
        dserror("number of initial set crosslinker types does not fit number of crosslinker types");

      for( int i = 0; i < static_cast<int>( maxnuminitcrosslinkerpertype.size() ); ++i )
        maxnum_init_crosslinker_pertype_[matcrosslinkerpertype_[i]] = maxnuminitcrosslinkerpertype[i];
    }
  }



  // viscosity
  viscosity_ =  crosslinking_params_list.get<double> ("VISCOSITY");

  // thermal energy
  kt_ = crosslinking_params_list.get<double> ("KT");

  // time step for stochastic events concering crosslinking
  deltatime_ = crosslinking_params_list.get<double> ("TIMESTEP");

  // safety check
  // todo: maybe make input of time step obligatory
  if ( deltatime_ < 0.0 )
  {
    deltatime_ = (*gstate.GetDeltaTime())[0];
    if ( gstate.GetMyRank() == 0 )
      std::cout << " Time step " << (*gstate.GetDeltaTime())[0] << " form Structural Dynamic section "
          "used for crosslinking.\n" << std::endl;
  }

  // distance between the two binding spots on a filament
  filamentbspotinterval_ = crosslinking_params_list.get<double>("FILAMENTBSPOTINTERVAL");
  if ( not (filamentbspotinterval_ > 0.0) )
    dserror(" Choose realistic value for FILAMENTBSPOTINTERVAL (i.e. distance between "
        "two binding spots on a filament) in input file. ");

  max_num_bonds_per_filament_bspot_ = crosslinking_params_list.get<int>("MAXNUMBONDSPERFILAMENTBSPOT");
  if ( max_num_bonds_per_filament_bspot_ < 0 )
    dserror( " Choose a number of bonds per filament binding spot >= 0. " );


  {
    // start and end arc parameter for binding spots on a filament
    filamentbspotrange_.clear();
    std::istringstream PL(
        Teuchos::getNumericStringParameter(crosslinking_params_list,"FILAMENTBSPOTRANGE"));
    std::string word;
    char* input;
    while (PL >> word)
      filamentbspotrange_.push_back( std::strtod( word.c_str(), &input ) );

    // sanity check
    if ( static_cast<int>(filamentbspotrange_.size()) != 2 )
      dserror(" Enter two values for FILAMENTBSPOTRANGE in input file. ");
    if ( filamentbspotrange_[0] > 0 and filamentbspotrange_[1] > 0 and (filamentbspotrange_[0] > filamentbspotrange_[1]) )
      dserror(" lower bound > upper bound, fix FILAMENTBSPOTRANGE in input file ");
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
