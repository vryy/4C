/*----------------------------------------------------------------------*/
/*!
\file inpar_windkessel.cpp

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_windkessel.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::WINDKESSEL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& windkstruct = list->sublist("WINDKESSEL-STRUCTURE COUPLING",false,"");

  DoubleParameter("TOLWINDKESSEL",1.0E-08,
                  "tolerance in the windkessel error norm for the Newton iteration",
                  &windkstruct);
  DoubleParameter("TOLWINDKESSELDOFINCR",1.0E-08,
                  "tolerance in the windkessel dof increment error norm for the Newton iteration",
                  &windkstruct);
  DoubleParameter("TIMINT_THETA",0.5,
                  "theta for one-step-theta time-integration scheme of Windkessel",
                  &windkstruct);
  setStringToIntegralParameter<int>("RESTART_WITH_WINDKESSEL","No","Must be chosen if a non-windkessel simulation is to be restarted as windkessel-structural coupled problem.",
                                 yesnotuple,yesnovalue,&windkstruct);
  setStringToIntegralParameter<int>("ENHANCED_OUTPUT","No","Set to yes for enhanced output (like e.g. derivative information)",
                                 yesnotuple,yesnovalue,&windkstruct);

  // linear solver id used for monolithic windkessel-structural problems
  IntParameter("LINEAR_WINDK_STRUCT_SOLVER",-1,"number of linear solver used for windkessel-structural problems",&windkstruct);

  setStringToIntegralParameter<int>("SOLALGORITHM","direct","",
                                 tuple<std::string>(
                                   "simple",
                                   "direct"),
                                 tuple<int>(
                                   INPAR::WINDKESSEL::windksolve_simple,
                                   INPAR::WINDKESSEL::windksolve_direct),
                                 &windkstruct);
}



void INPAR::WINDKESSEL::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;


  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a four-element Windkessel - mhv 11/13

  Teuchos::RCP<ConditionDefinition> windkesselcondition =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF WINDKESSEL CONDITIONS",
                                         "WindkesselStdStructureCond",
                                         "Surface Windkessel",
                                         DRT::Condition::WindkesselStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselcondition,"id");
  AddNamedReal(windkesselcondition,"C");
  AddNamedReal(windkesselcondition,"R_p");
  AddNamedReal(windkesselcondition,"Z_c");
  AddNamedReal(windkesselcondition,"L");
  AddNamedReal(windkesselcondition,"p_ref");
  AddNamedReal(windkesselcondition,"p_0");

  condlist.push_back(windkesselcondition);


  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a special heart valve arterial four-element Windkessel - mhv 02/14

  Teuchos::RCP<ConditionDefinition> windkesselheartvalvearterialcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF HEART VALVE ARTERIAL WINDKESSEL CONDITIONS",
                                         "WindkesselHeartValveArterialStructureCond",
                                         "Surface heart valve arterial Windkessel",
                                         DRT::Condition::WindkesselHeartValveArterialStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselheartvalvearterialcond,"id");
  AddNamedReal(windkesselheartvalvearterialcond,"R_arvalve_max");
  AddNamedReal(windkesselheartvalvearterialcond,"R_arvalve_min");
  AddNamedReal(windkesselheartvalvearterialcond,"R_atvalve_max");
  AddNamedReal(windkesselheartvalvearterialcond,"R_atvalve_min");
  AddNamedReal(windkesselheartvalvearterialcond,"k_p");
  AddNamedReal(windkesselheartvalvearterialcond,"C");
  AddNamedReal(windkesselheartvalvearterialcond,"R_p");
  AddNamedReal(windkesselheartvalvearterialcond,"Z_c");
  AddNamedReal(windkesselheartvalvearterialcond,"L");
  AddNamedReal(windkesselheartvalvearterialcond,"p_ref");
  AddNamedReal(windkesselheartvalvearterialcond,"p_ar_0");
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("P_AT")));
  AddNamedReal(windkesselheartvalvearterialcond,"fac");
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("crv")));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("VALVE")));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new StringConditionComponent("valvelaw", "smooth",
                                                                                       Teuchos::tuple<std::string>("smooth","pwlin"),
                                                                                       Teuchos::tuple<std::string>("smooth","pwlin"),
                                                                                       true)));

  condlist.push_back(windkesselheartvalvearterialcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a heart valve arterial Windkessel accounting for proximal and distal arterial pressure
  // formulation proposed by Cristobal Bertoglio - mhv 03/14

  Teuchos::RCP<ConditionDefinition> windkesselheartvalvearterialproxdistcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF HEART VALVE ARTERIAL PROX DIST WINDKESSEL CONDITIONS",
                                         "WindkesselHeartValveArterialProxDistStructureCond",
                                         "Surface heart valve arterial proximal and distal Windkessel",
                                         DRT::Condition::WindkesselHeartValveArterialProxDistStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselheartvalvearterialproxdistcond,"id");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_arvalve_max");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_arvalve_min");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_atvalve_max");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_atvalve_min");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"k_p");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"L_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"C_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"C_ard");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_ard");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_ref");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_arp_0");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"y_arp_0");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_ard_0");
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("P_AT")));
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"fac");
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("crv")));
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));

  condlist.push_back(windkesselheartvalvearterialproxdistcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a heart valve cardiovascular full Windkessel model (closed-loop circulatory system model)
  // mhv 02/15

  Teuchos::RCP<ConditionDefinition> windkesselheartvalvecardiovascularfullcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF HEART VALVE CARDIOVASCULAR FULL WINDKESSEL CONDITIONS",
                                         "WindkesselHeartValveCardiovascularFullStructureCond",
                                         "Surface heart valve cardiovascular full Windkessel",
                                         DRT::Condition::WindkesselHeartValveCardiovascularFullStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselheartvalvecardiovascularfullcond,"id");
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_arvalve_max"); // maximal arterial valve resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_arvalve_min"); // minimal arterial valve resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_atvalve_max"); // maximal atrial valve resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_atvalve_min"); // minimal atrial valve resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"E_at_max"); // maximum atrial elastance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"E_at_min"); // baseline atrial elastance
  windkesselheartvalvecardiovascularfullcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("at_act_curve"))); // atrial activation curve
  windkesselheartvalvecardiovascularfullcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true))); // atrial activation curve
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"C_ar"); // arterial compliance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_ar"); // arterial resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"L_ar"); // arterial inertance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"Z_ar"); // arterial characteristic impedance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"C_ven"); // venous compliance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"R_ven"); // venous resistance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"L_ven"); // venous inertance
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"p_at_0"); // initial atrial pressure
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"p_ar_0"); // initial arterial pressure
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"p_ven_0"); // initial venous pressure
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"q_ar_0"); // initial arterial flow
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"q_ven_0"); // initial venous flow
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"V_at_0"); // initial volume of atrium
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"V_ar_0"); // initial volume of arteries and capillaries
  AddNamedReal(windkesselheartvalvecardiovascularfullcond,"V_ven_0"); // initial volume of veins

  condlist.push_back(windkesselheartvalvecardiovascularfullcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and Windkessel: Neumann coupling surface - mhv 11/13

  Teuchos::RCP<ConditionDefinition> windkesselstructurecouplingcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF WINDKESSEL STRUCTURE COUPLING CONDITIONS",
                                         "SurfaceNeumann",
                                         "structure windkessel coupling surface condition",
                                         DRT::Condition::WindkesselStructureCoupling,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselstructurecouplingcond,"coupling_id");

  condlist.push_back(windkesselstructurecouplingcond);

}

