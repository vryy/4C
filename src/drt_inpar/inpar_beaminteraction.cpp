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
  /* parameters for one-step-theta structural integrator */

  Teuchos::ParameterList& contraccells = beaminteraction.sublist("CONTRACTILE CELLS",false,"");


  setStringToIntegralParameter<int>("CONTRACTILECELLS","No",
                                 "Contractile cells in problem",
                                 yesnotuple,yesnovalue,&contraccells);

  //number of overall crosslink molecules in the boundary volume
  IntParameter("NUMCELLS",0,"number of contracting cells in simulation volume",&contraccells);



}
