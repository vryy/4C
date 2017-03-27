/*----------------------------------------------------------------------------*/
/*!
\file crosslinking_params.cpp

\brief data container holding all crosslinking input parameters

\level 3

\maintainer Jonas Eichinger, Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "crosslinking_params.H"

#include "../drt_lib/drt_globalproblem.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::CrosslinkingParams::CrosslinkingParams()
  : isinit_(false),
    issetup_(false),
    viscosity_(0.0),
    kt_(0.0),
    filamentbspotinterval_(0.0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::CrosslinkingParams::Init()
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
      if( numcrosslinkerpertype_[i] < 0 )
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
      if( matcrosslinkerpertype_[i] < 0 )
        dserror(" negative material number does not make sense.");
   }

  // safety check
  if ( numcrosslinkerpertype_.size() != matcrosslinkerpertype_.size() )
    dserror("number of crosslinker types does not fit number of assigned materials");


  // viscosity
  viscosity_ =  crosslinking_params_list.get<double> ("VISCOSITY");

  // thermal energy
  kt_ = crosslinking_params_list.get<double> ("KT");

  // distance between the two binding spots on a filament
  filamentbspotinterval_ = crosslinking_params_list.get<double>("FILAMENTBSPOTINTERVAL");
  if( not (filamentbspotinterval_ > 0.0) )
    dserror(" Choose realistic value for FILAMENTBSPOTINTERVAL (i.e. distance between "
        "two binding spots on a filament) in input file. ");

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
