/*-----------------------------------------------------------*/
/*!
\file inpar_browniandyn.cpp

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/


#include "inpar_browniandyn.H"

#include "drt_validparameters.H"


void INPAR::BROWNIANDYN::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& browniandyn = list->sublist("BROWNIAN DYNAMICS",false,"");

  setStringToIntegralParameter<int>("BROWNDYNPROB","No",
                                 "Statistical mechanics problem",
                                 yesnotuple,yesnovalue,&browniandyn);



  //Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY",0.0,"viscosity",&browniandyn);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&browniandyn);
  //cutoff for random forces, which determines the maximal value
  DoubleParameter("MAXRANDFORCE",-1.0,"Any random force beyond MAXRANDFORCE*(standard dev.) will be omitted and redrawn. -1.0 means no bounds.'",&browniandyn);
  //time interval in which random numbers are constant
  DoubleParameter("TIMEINTCONSTRANDFORCES",-1.0,"Within this time interval the random numbers remain constant. -1.0 means no prescribed time interval.'",&browniandyn);


}
