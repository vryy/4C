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

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAMINTERACTION",false,"");

  setStringToIntegralParameter<int>("BROWNDYNPROB","No",
                                 "Statistical mechanics problem",
                                 yesnotuple,yesnovalue,&beaminteraction);



}
