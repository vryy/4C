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
    double Ew, Ea, tw, A;
    linedef->ExtractDouble("WallCompliance",Ew);
    linedef->ExtractDouble("AirCompliance",Ea); 
    linedef->ExtractDouble("WallThickness",tw);
    linedef->ExtractDouble("Area",A);
    
    elemParams_["WallCompliance"] = Ew;
    elemParams_["AirCompliance"]  = Ea; 
    elemParams_["WallThickness"]  = tw;
    elemParams_["Area"]           = A;

  }
  else
  {
    dserror("Reading type of RED_AIRWAY element failed: ComplientResistive/PoiseuilleResistive/InductoResistive/RLC/SUKI");
    exit(1);
  }

  // build the elemParams Vector (This will be ommited later)
  elemVars_["flow_in"]  = 0.0;
  elemVars_["flow_out"] = 0.0;
  if (elemType_ ==  "PoiseuilleResistive")
  {
    
  }
  else if (elemType_ == "InductoResistive")
  {
    elemVars_["inductor_pressure"]  = 0.0;
    elemVars_["inductor_flow"]      = 0.0;
  }
  else if (elemType_ == "ComplientResistive")
  {
    elemVars_["capacitor_pressure"] = 0.0;
    elemVars_["capacitor_flow"]     = 0.0;
  }
  else if (elemType_ == "RLC")
  {
    elemVars_["capacitor_pressure"] = 0.0;
    elemVars_["capacitor_flow"]     = 0.0;
    elemVars_["inductor_pressure"]  = 0.0;
    elemVars_["inductor_flow"]      = 0.0;
  }
  else if (elemType_ == "SUKI")
  {
  }
  else
  {
  }
  

  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_RED_AIRWAY
