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

#include "red_airway.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional AIRWAY element.",ndim);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractString("TYPE",elemType_);
  if (elemType_ == "Resistive" || elemType_ == "InductoResistive" || elemType_ == "ComplientResistive" || elemType_ == "RLC" || elemType_ == "SUKI")
  {
    linedef->ExtractString("Resistance",resistance_);
    double Ew, Ea, tw, A;
    int generation;
    linedef->ExtractDouble("WallCompliance",Ew);
    linedef->ExtractDouble("AirCompliance",Ea); 
    linedef->ExtractDouble("WallThickness",tw);
    linedef->ExtractDouble("Area",A);
    linedef->ExtractInt("Generation",generation);
    
    elemParams_["WallCompliance"] = Ew;
    elemParams_["AirCompliance"]  = Ea; 
    elemParams_["WallThickness"]  = tw;
    elemParams_["Area"]           = A;
    generation_                   = generation;

  }
  else
  {
    dserror("Reading type of RED_AIRWAY element failed: ComplientResistive/PoiseuilleResistive/TurbulentPoiseuilleResistive/InductoResistive/RLC/SUKI");
    exit(1);
  }

  return true;
}

#endif  // #ifdef D_RED_AIRWAY
