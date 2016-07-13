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
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_arvalve_max"); // maximal arterial valve resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_arvalve_min"); // minimal arterial valve resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_atvalve_max"); // maximal atrial valve resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_atvalve_min"); // minimal atrial valve resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"E_at_max"); // maximum atrial elastance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"E_at_min"); // baseline atrial elastance
  cardiovascular0darterialvenoussyspulcoupledcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("at_act_curve"))); // atrial activation curve
  cardiovascular0darterialvenoussyspulcoupledcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true))); // atrial activation curve
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"C_ar"); // arterial compliance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_ar"); // arterial resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"L_ar"); // arterial inertance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"Z_ar"); // arterial characteristic impedance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"C_ven"); // venous compliance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"R_ven"); // venous resistance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"L_ven"); // venous inertance
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"p_at_0"); // initial atrial pressure
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"q_v_in_0"); // initial ventricular in-flux
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"q_v_out_0"); // initial ventricular out-flux
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"p_v_0"); // initial ventricular pressure
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"p_ar_0"); // initial arterial pressure
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"q_ar_0"); // initial arterial flux
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"p_ven_0"); // initial venous pressure
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"q_ven_0"); // initial venous flux
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"V_at_0"); // initial volume of atrium - only for postprocessing matters!
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"V_ar_0"); // initial volume of arteries and capillaries - only for postprocessing matters!
  AddNamedReal(cardiovascular0darterialvenoussyspulcoupledcond,"V_ven_0"); // initial volume of veins - only for postprocessing matters!

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

