/*!----------------------------------------------------------------------
\file red_airway_input.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_RED_AIRWAYS
#ifdef CCADISCRET

#include "red_airway.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

extern struct _GENPROB     genprob;

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  if (genprob.ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional AIRWAY element.",genprob.ndim);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractString("TYPE",elemType_);
  if (elemType_ == "PoiseuilleResistive" || elemType_ == "InductoResistive" || elemType_ == "ComplientResistive" || elemType_ == "RLC" || elemType_ == "SUKI")
  {
    linedef->ExtractDouble("WallCompliance",Ew_);
    linedef->ExtractDouble("AirCompliance",Ea_);
    linedef->ExtractDouble("WallThickness",tw_);
    linedef->ExtractDouble("Area",A_);
  }
  else
  {
    dserror("Reading type of RED_AIRWAY element failed: ComplientResistive/PoiseuilleResistive/InductoResistive/RLC/SUKI");
    exit(1);
  }
  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_RED_AIRWAY
