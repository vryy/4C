/*-----------------------------------------------------------*/
/*!
\file inpar_beaminteraction.cpp

\brief input parameter for beaminteraction

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/


#include "inpar_beaminteraction.H"

#include "drt_validparameters.H"


void INPAR::BEAMINTERACTION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAM INTERACTION",false,"");


  /*----------------------------------------------------------------------*/
  /* parameters for crosslinking submodel */

  Teuchos::ParameterList& crosslinking = beaminteraction.sublist("CROSSLINKING",false,"");

  // remove this some day
  setStringToIntegralParameter<int>("CROSSLINKER","No", "Crosslinker in problem", yesnotuple,yesnovalue,&crosslinking);

  // number of crosslinker of certain type
  setNumericStringParameter("NUMCROSSLINKERPERTYPE","0","number of crosslinker of certain type ",&crosslinking);
  // material number characterizing crosslinker type
  setNumericStringParameter("MATCROSSLINKERPERTYPE","-1","material number characterizing crosslinker type ",&crosslinking);
  //Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY",0.0,"viscosity",&crosslinking);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&crosslinking);
  // distance between two binding spots on a filament
  DoubleParameter("FILAMENTBSPOTINTERVAL",-1.0 ,"distance between two binding spots on a filament",&crosslinking);
  // start and end for bspots on a filament
  setNumericStringParameter("FILAMENTBSPOTRANGE","-1.0 -1.0", "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);


  /*----------------------------------------------------------------------*/
  /* parameters for contractile cell submodel */

  Teuchos::ParameterList& contraccells = beaminteraction.sublist("CONTRACTILE CELLS",false,"");


  setStringToIntegralParameter<int>("CONTRACTILECELLS","No",
                                 "Contractile cells in problem",
                                 yesnotuple,yesnovalue,&contraccells);

  //number of overall crosslink molecules in the boundary volume
  IntParameter("NUMCELLS",0,"number of contracting cells in simulation volume",&contraccells);



}
