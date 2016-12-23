/*-----------------------------------------------------------*/
/*!
\file inpar_crosslinking.cpp

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_crosslinking.H"


void INPAR::CROSSLINKING::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& crosslinking = list->sublist("CROSSLINKING",false,"");


  setStringToIntegralParameter<int>("CROSSLINKER","No",
                                 "Crosslinker in problem",
                                 yesnotuple,yesnovalue,&crosslinking);



  //Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY",0.0,"viscosity",&crosslinking);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&crosslinking);

  //number of overall crosslink molecules in the boundary volume
  IntParameter("NUMCROSSLINK",0,"number of crosslinkers for switching on- and off-rates; if molecule diffusion model is used: number of crosslink molecules",&crosslinking);
  //Reading double parameter for crosslinker protein mean length
  DoubleParameter("LINKINGLENGTH",0.0,"Mean distance between two nodes connected by a crosslinker",&crosslinking);
  //Absolute value of difference between maximal/minimal and mean cross linker length
  DoubleParameter("LINKINGLENGTHTOL",0.0,"Absolute value of difference between maximal/minimal and mean cross linker length",&crosslinking);
  //angle between filament axes at crosslinked points with zero potential energy
  DoubleParameter("LINKINGANGLE",0.0,"equilibrium angle between crosslinker axis and filament at each binding site",&crosslinking);
  //only angles in the range PHIZERO +/- PHIODEV are admitted at all for the angle PHI between filament axes at crosslinked points;
  //the default value for this parameter is 2*pi so that by default any value is admitted
  DoubleParameter("LINKINGANGLETOL",6.28,"only angles in the range PHIZERO +/- PHIODEV",&crosslinking);
  //Moment of inertia of area of crosslinkers
  DoubleParameter("ILINK",0.0,"Moment of inertia of area of crosslinkers",&crosslinking);
  //Polar moment of inertia of area of crosslinkers
  DoubleParameter("IPLINK",0.0,"Polar moment of inertia of area of crosslinkers",&crosslinking);
  //Cross section of crosslinkers
  DoubleParameter("ALINK",0.0,"Cross section of crosslinkers",&crosslinking);
  //Reading double parameter for crosslinker on-rate
  DoubleParameter("K_ON",0.0,"crosslinker on-rate ",&crosslinking);
  //Reading double parameter for crosslinker off-rate
  DoubleParameter("K_OFF",0.0,"crosslinker off-rate ",&crosslinking);



}
