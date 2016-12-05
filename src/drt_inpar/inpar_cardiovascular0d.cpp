/*----------------------------------------------------------------------*/
/*!
\file inpar_cardiovascular0d.cpp

\brief Input parameters for 0d cardiovascular-structure coupling

\level 2

\maintainer Marc Hirschvogel
*/

/*----------------------------------------------------------------------*/



#include "inpar_cardiovascular0d.H"

#include "drt_validparameters.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::CARDIOVASCULAR0D::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& cardvasc0dstruct = list->sublist("CARDIOVASCULAR 0D-STRUCTURE COUPLING",false,"");

  DoubleParameter("TOL_CARDVASC0D_RES",1.0E-08,
                  "tolerance in the cardiovascular0d error norm for the Newton iteration",
                  &cardvasc0dstruct);
  DoubleParameter("TOL_CARDVASC0D_DOFINCR",1.0E-08,
                  "tolerance in the cardiovascular0d dof increment error norm for the Newton iteration",
                  &cardvasc0dstruct);
  DoubleParameter("TIMINT_THETA",0.5,
                  "theta for one-step-theta time-integration scheme of Cardiovascular0D",
                  &cardvasc0dstruct);
  setStringToIntegralParameter<int>("RESTART_WITH_CARDVASC0D","No","Must be chosen if a non-cardiovascular0d simulation is to be restarted as cardiovascular0d-structural coupled problem.",
                                 yesnotuple,yesnovalue,&cardvasc0dstruct);
  setStringToIntegralParameter<int>("ENHANCED_OUTPUT","No","Set to yes for enhanced output (like e.g. derivative information)",
                                 yesnotuple,yesnovalue,&cardvasc0dstruct);

  // linear solver id used for monolithic 0D cardiovascular-structural problems
  IntParameter("LINEAR_COUPLED_SOLVER",-1,"number of linear solver used for cardiovascular 0D-structural problems",&cardvasc0dstruct);

  setStringToIntegralParameter<int>("SOLALGORITHM","direct","",
                                 tuple<std::string>(
                                   "simple",
                                   "direct"),
                                 tuple<int>(
                                   INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_simple,
                                   INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_direct),
                                 &cardvasc0dstruct);



  Teuchos::ParameterList& cardvasc0dartvensyspulcoupled = cardvasc0dstruct.sublist("CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED PARAMETERS",false,"");

  DoubleParameter("R_arvalve_max_l",0.0,"maximal left arterial (semilunar) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_arvalve_min_l",0.0,"minimal left arterial (semilunar) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_atvalve_max_l",0.0,"maximal left atrial (atrioventricular) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_atvalve_min_l",0.0,"minimal left atrial (atrioventricular) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_arvalve_max_r",0.0,"maximal right arterial (semilunar) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_arvalve_min_r",0.0,"minimal right arterial (semilunar) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_atvalve_max_r",0.0,"maximal right atrial (atrioventricular) valve resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_atvalve_min_r",0.0,"minimal right atrial (atrioventricular) valve resistance",&cardvasc0dartvensyspulcoupled);

  setStringToIntegralParameter<int>("ATRIUM_MODEL","0D","",
                                 tuple<std::string>(
                                   "0D",
                                   "3D"),
                                 tuple<int>(
                                   INPAR::CARDIOVASCULAR0D::atr_elastance_0d,
                                   INPAR::CARDIOVASCULAR0D::atr_structure_3d),
                                 &cardvasc0dartvensyspulcoupled);
  IntParameter("Atrium_act_curve_l",-1,"0D left atrium activation curve",&cardvasc0dartvensyspulcoupled);
  IntParameter("Atrium_act_curve_r",-1,"0D right atrium activation curve",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_at_max_l",0.0,"0D maximum left atrial elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_at_min_l",0.0,"0D baseline left atrial elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_at_max_r",0.0,"0D maximum right atrial elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_at_min_r",0.0,"0D baseline right atrial elastance",&cardvasc0dartvensyspulcoupled);

  setStringToIntegralParameter<int>("VENTRICLE_MODEL","3D","",
                                 tuple<std::string>(
                                   "3D",
                                   "0D"),
                                 tuple<int>(
                                   INPAR::CARDIOVASCULAR0D::ventr_structure_3d,
                                   INPAR::CARDIOVASCULAR0D::ventr_elastance_0d),
                                 &cardvasc0dartvensyspulcoupled);
  IntParameter("Ventricle_act_curve_l",-1,"0D left ventricular activation curve",&cardvasc0dartvensyspulcoupled);
  IntParameter("Ventricle_act_curve_r",-1,"0D right ventricular activation curve",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_v_max_l",0.0,"0D maximum left ventricular elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_v_min_l",0.0,"0D baseline left ventricular elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_v_max_r",0.0,"0D maximum right ventricular elastance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("E_v_min_r",0.0,"0D baseline right ventricular elastance",&cardvasc0dartvensyspulcoupled);

  DoubleParameter("C_ar_sys",0.0,"systemic arterial compliance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_ar_sys",0.0,"systemic arterial resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("L_ar_sys",0.0,"systemic arterial inertance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("Z_ar_sys",0.0,"systemic arterial impedance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("C_ar_pul",0.0,"pulmonary arterial compliance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_ar_pul",0.0,"pulmonary arterial resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("L_ar_pul",0.0,"pulmonary arterial inertance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("Z_ar_pul",0.0,"pulmonary arterial impedance",&cardvasc0dartvensyspulcoupled);

  DoubleParameter("C_ven_sys",0.0,"systemic venous compliance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_ven_sys",0.0,"systemic venous resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("L_ven_sys",0.0,"systemic venous inertance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("C_ven_pul",0.0,"pulmonary venous compliance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("R_ven_pul",0.0,"pulmonary venous resistance",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("L_ven_pul",0.0,"pulmonary venous inertance",&cardvasc0dartvensyspulcoupled);

  // intital conditions
  DoubleParameter("q_v_in_l_0",0.0,"initial left ventricular in-flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_at_l_0",0.0,"initial left atrial pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_v_out_l_0",0.0,"initial left ventricular out-flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_v_l_0",0.0,"initial left ventricular pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_ar_sys_0",0.0,"initial systemic arterial pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_ar_sys_0",0.0,"initial systemic arterial flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_ven_sys_0",0.0,"initial systemic venous pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_ven_sys_0",0.0,"initial systemic venous flux",&cardvasc0dartvensyspulcoupled);

  DoubleParameter("q_v_in_r_0",0.0,"initial right ventricular in-flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_at_r_0",0.0,"initial right atrial pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_v_out_r_0",0.0,"initial right ventricular out-flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_v_r_0",0.0,"initial right ventricular pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_ar_pul_0",0.0,"initial pulmonary arterial pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_ar_pul_0",0.0,"initial pulmonary arterial flux",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("p_ven_pul_0",0.0,"initial pulmonary venous pressure",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("q_ven_pul_0",0.0,"initial pulmonary venous flux",&cardvasc0dartvensyspulcoupled);
  // volumes - only for postprocessing matters!
  DoubleParameter("V_at_l_0",0.0,"initial volume of left 0D atrium - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_v_l_0",0.0,"initial volume of left 0D ventricle - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_ar_sys_0",0.0,"initial volume of systemic arteries and capillaries - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_ven_sys_0",0.0,"initial volume of systemic veins - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_at_r_0",0.0,"initial volume of right 0D atrium - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_v_r_0",0.0,"initial volume of right 0D ventricle - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_ar_pul_0",0.0,"initial volume of pulmonary arteries and capillaries - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);
  DoubleParameter("V_ven_pul_0",0.0,"initial volume of pulmonary veins - only for postprocessing matters!",&cardvasc0dartvensyspulcoupled);

}



void INPAR::CARDIOVASCULAR0D::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;


  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a four-element windkessel - mhv 11/13

  Teuchos::RCP<ConditionDefinition> cardiovascular0dwindkesselonlycondition =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF CARDIOVASCULAR 0D WINDKESSEL ONLY CONDITIONS",
                                         "Cardiovascular0DWindkesselOnlyStructureCond",
                                         "Surface Cardiovascular0D",
                                         DRT::Condition::Cardiovascular0DWindkesselOnly_Structure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(cardiovascular0dwindkesselonlycondition,"id");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"C");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"R_p");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"Z_c");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"L");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"p_ref");
  AddNamedReal(cardiovascular0dwindkesselonlycondition,"p_0");

  condlist.push_back(cardiovascular0dwindkesselonlycondition);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and an arterial cardiovascular 0D flow model accounting for proximal and distal arterial pressure
  // formulation proposed by Cristobal Bertoglio - mhv 03/14

  Teuchos::RCP<ConditionDefinition> cardiovascular0darterialproxdistcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF CARDIOVASCULAR 0D ARTERIAL PROX DIST CONDITIONS",
                                         "Cardiovascular0DArterialProxDistStructureCond",
                                         "Surface 0D cardiovascular arterial proximal and distal",
                                         DRT::Condition::Cardiovascular0DArterialProxDist_Structure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(cardiovascular0darterialproxdistcond,"id");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_arvalve_max");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_arvalve_min");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_atvalve_max");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_atvalve_min");
  AddNamedReal(cardiovascular0darterialproxdistcond,"k_p");
  AddNamedReal(cardiovascular0darterialproxdistcond,"L_arp");
  AddNamedReal(cardiovascular0darterialproxdistcond,"C_arp");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_arp");
  AddNamedReal(cardiovascular0darterialproxdistcond,"C_ard");
  AddNamedReal(cardiovascular0darterialproxdistcond,"R_ard");
  AddNamedReal(cardiovascular0darterialproxdistcond,"p_ref");
  AddNamedReal(cardiovascular0darterialproxdistcond,"p_v_0");
  AddNamedReal(cardiovascular0darterialproxdistcond,"p_arp_0");
  AddNamedReal(cardiovascular0darterialproxdistcond,"y_arp_0");
  AddNamedReal(cardiovascular0darterialproxdistcond,"p_ard_0");
  cardiovascular0darterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("P_AT")));
  AddNamedReal(cardiovascular0darterialproxdistcond,"fac");
  cardiovascular0darterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("crv")));
  cardiovascular0darterialproxdistcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));

  condlist.push_back(cardiovascular0darterialproxdistcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a full closed-loop 0D cardiovascular flow model (closed-loop circulatory system model)
  // mhv 02/15

  Teuchos::RCP<ConditionDefinition> cardiovascular0darterialvenoussyspulcoupledcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED CONDITIONS",
                                         "Cardiovascular0DArterialVenousSysPulCoupledStructureCond",
                                         "Surface cardiovascular 0D arterial venous sys pul coupled",
                                         DRT::Condition::Cardiovascular0DArterialVenousSysPulCoupled_Structure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(cardiovascular0darterialvenoussyspulcoupledcond,"id");
  cardiovascular0darterialvenoussyspulcoupledcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("TYPE")));
  cardiovascular0darterialvenoussyspulcoupledcond->AddComponent(Teuchos::rcp(new StringConditionComponent("type", "ventricle_left",
                                                                                       Teuchos::tuple<std::string>("ventricle_left","ventricle_right","atrium_left","atrium_right","dummy"),
                                                                                       Teuchos::tuple<std::string>("ventricle_left","ventricle_right","atrium_left","atrium_right","dummy"),
                                                                                       false)));

  condlist.push_back(cardiovascular0darterialvenoussyspulcoupledcond);


  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and 0D cardiovascular flow models: Neumann coupling surface - mhv 11/13

  Teuchos::RCP<ConditionDefinition> cardiovascular0dstructurecouplingcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF CARDIOVASCULAR 0D-STRUCTURE COUPLING CONDITIONS",
                                         "SurfaceNeumannCardiovascular0D",
                                         "structure 0d cardiovascular coupling surface condition",
                                         DRT::Condition::Cardiovascular0DStructureCoupling,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(cardiovascular0dstructurecouplingcond,"coupling_id");

  condlist.push_back(cardiovascular0dstructurecouplingcond);



}

