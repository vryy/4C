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
#include "drt_validparameters.H"
#include "../drt_lib/drt_conditiondefinition.H"

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

  // number of initial (are set right in the beginning) crosslinker of certain type
  setNumericStringParameter("MAXNUMINITCROSSLINKERPERTYPE","0","number of initial crosslinker of certain type (additional to NUMCROSSLINKERPERTYPE) ",&crosslinking);
  // number of crosslinker of certain type
  setNumericStringParameter("NUMCROSSLINKERPERTYPE","0","number of crosslinker of certain type ",&crosslinking);
  // material number characterizing crosslinker type
  setNumericStringParameter("MATCROSSLINKERPERTYPE","-1","material number characterizing crosslinker type ",&crosslinking);
  // time step for stochastic events concerning crosslinking
  DoubleParameter("TIMESTEP",-1.0, "time step for stochastic events concerning crosslinking (e.g. diffusion, p_link, p_unlink) ", &crosslinking);
  //Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY",0.0,"viscosity",&crosslinking);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&crosslinking);
  // distance between two binding spots on a filament
  DoubleParameter("FILAMENTBSPOTINTERVAL",-1.0 ,"distance between two binding spots on a filament",&crosslinking);
  // maximal number of binding partner per filament binding spot
  IntParameter("MAXNUMBONDSPERFILAMENTBSPOT", 1 ,"maximal number of bonds per filament binding spot",&crosslinking);
  // start and end for bspots on a filament
  setNumericStringParameter("FILAMENTBSPOTRANGE","-1.0 -1.0", "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);


  /*----------------------------------------------------------------------*/
  /* parameters for contractile cell submodel */

  Teuchos::ParameterList& contraccells = beaminteraction.sublist("CONTRACTILE CELLS",false,"");

  setStringToIntegralParameter<int>( "CONTRACTILECELLS", "No",
                                 "Contractile cells in problem",
                                 yesnotuple, yesnovalue, &contraccells );

  //number of overall crosslink molecules in the boundary volume
  IntParameter( "NUMCELLS", 0, "number of contracting cells in simulation volume", &contraccells );


  /*----------------------------------------------------------------------*/
  /* parameters for beam to ? contact submodel*/
  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  /* parameters for beam to beam contact */

  Teuchos::ParameterList& beamtobeamcontact = beaminteraction.sublist("BEAM TO BEAM CONTACT",false,"");

  setStringToIntegralParameter<int>("STRATEGY","None","Type of employed solving strategy",
        tuple<std::string>("None","none",
                           "Penalty", "penalty"),
        tuple<int>(
                bstr_none, bstr_none,
                bstr_penalty, bstr_penalty),
        &beamtobeamcontact);

  // ...

  /*----------------------------------------------------------------------*/
  /* parameters for beam to sphere contact */

  Teuchos::ParameterList& beamtospherecontact = beaminteraction.sublist("BEAM TO SPHERE CONTACT",false,"");

  setStringToIntegralParameter<int>("STRATEGY","None","Type of employed solving strategy",
        tuple<std::string>("None","none",
                           "Penalty", "penalty"),
        tuple<int>(
                bstr_none, bstr_none,
                bstr_penalty, bstr_penalty),
        &beamtospherecontact);

  // ...

  /*----------------------------------------------------------------------*/
  /* parameters for beam to solid contact */

  Teuchos::ParameterList& beamtosolidcontact = beaminteraction.sublist("BEAM TO SOLID CONTACT",false,"");

  setStringToIntegralParameter<int>("STRATEGY","None","Type of employed solving strategy",
        tuple<std::string>("None","none",
                           "Penalty", "penalty"),
        tuple<int>(
                bstr_none, bstr_none,
                bstr_penalty, bstr_penalty),
        &beamtosolidcontact);

  // ...

}

void INPAR::BEAMINTERACTION::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*-------------------------------------------------------------------*/
  // beam potential interaction: atom/charge density per unit length on LINE
  Teuchos::RCP<ConditionDefinition> beam_filament_condition =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE BEAM FILAMENT CONDITIONS",
                                         "BeamLineFilamentCondition",
                                         "Beam_Line_Filament_Condition",
                                         DRT::Condition::FilamentBeamLineCondition,
                                         false,
                                         DRT::Condition::Line));

  beam_filament_condition->AddComponent( Teuchos::rcp( new SeparatorConditionComponent( "ID" ) ) );
  beam_filament_condition->AddComponent( Teuchos::rcp( new IntConditionComponent( "FilamentId", false, false ) ) );
  beam_filament_condition->AddComponent( Teuchos::rcp( new SeparatorConditionComponent( "TYPE", true ) ) );
  beam_filament_condition->AddComponent( Teuchos::rcp( new StringConditionComponent( "Type",
        "Arbitrary",
        Teuchos::tuple<std::string>( "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen" ),
        Teuchos::tuple<std::string>( "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen" ), true ) ) );

  condlist.push_back(beam_filament_condition);


}
