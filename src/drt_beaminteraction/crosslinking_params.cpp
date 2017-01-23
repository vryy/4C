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
    numcrosslink_(0),
    viscosity_(0.0),
    kt_(0.0),
    linkinglength_(0.0),
    linkinglengthtol_(0.0),
    alink_(0.0),
    ilink_(0.0),
    iplink_(0.0),
    linkingangle_(0.0),
    linkingangletol_(0.0),
    k_on_(0.0),
    k_off_(0.0),
    deltabelleq_(0.0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::CrosslinkingParams::Init()
{
  issetup_ = false;

  const Teuchos::ParameterList& crosslinking_params_list =
      DRT::Problem::Instance()->CrosslinkingParams();

  // number of crosslinker
  numcrosslink_ = crosslinking_params_list.get<int> ("NUMCROSSLINK");
  // viscosity
  viscosity_ =  crosslinking_params_list.get<double> ("VISCOSITY");
  // thermal energy
  kt_ = crosslinking_params_list.get<double> ("KT");
  // distance between the two binding domains of a linkers
  linkinglength_ = crosslinking_params_list.get<double>("LINKINGLENGTH");
  // tolerance for linking length
  linkinglengthtol_ = crosslinking_params_list.get<double>("LINKINGLENGTHTOL");
  // crossection of linker
  alink_ = crosslinking_params_list.get<double>("ALINK");
  // moment of inertia of linker
  ilink_ = crosslinking_params_list.get<double>("ILINK");
  // 2* ilink_
  iplink_ = crosslinking_params_list.get<double>("IPLINK");
  // linking angle
  linkingangle_ = crosslinking_params_list.get<double> ("LINKINGANGLE");
  // tolerance for linking angle
  linkingangletol_ = crosslinking_params_list.get<double>("LINKINGANGLETOL");
  // on rate
  k_on_ = crosslinking_params_list.get<double>("K_ON");
  // off rate
  k_off_ = crosslinking_params_list.get<double> ("K_OFF");
    // off rate
  deltabelleq_ = crosslinking_params_list.get<double> ("DELTABELLEQ");


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
