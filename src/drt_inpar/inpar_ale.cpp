/*----------------------------------------------------------------------*/
/*!
\file inpar_ale.cpp

\maintainer Matthias Mayr

\brief Input parameters for ALE mesh motion

\level 2

</pre>
*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_ale.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::ALE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  DoubleParameter("TIMESTEP",0.1,"time step size",&adyn);
  IntParameter("NUMSTEP",41,"max number of time steps",&adyn);
  DoubleParameter("MAXTIME",4.0,"max simulation time",&adyn);

  setStringToIntegralParameter<int>("ALE_TYPE", "solid",
      "ale mesh movement algorithm",
      tuple<std::string>(
          "solid",
          "solid_linear",
          "laplace_material",
          "laplace_spatial",
          "springs_material",
          "springs_spatial"
          ),
      tuple<int>(
          solid,
          solid_linear,
          laplace_material,
          laplace_spatial,
          springs_material,
          springs_spatial
          ),
      &adyn);

  BoolParameter("ASSESSMESHQUALITY", "no",
      "Evaluate element quality measure according to [Oddy et al. 1988]",
      &adyn);

  BoolParameter("UPDATEMATRIX", "no",
      "Update stiffness matrix in every time step (only for linear/material strategies)",
      &adyn);

  IntParameter("MAXITER",1,"Maximum number of newton iterations.",&adyn);
  DoubleParameter("TOLRES",1.0e-06,"Absolute tolerance for length scaled L2 residual norm ",&adyn);
  DoubleParameter("TOLDISP",1.0e-06,"Absolute tolerance for length scaled L2 increment norm ",&adyn);

  IntParameter("NUM_INITSTEP",0,"",&adyn);
  IntParameter("RESTARTEVRY",1,"write restart data every RESTARTEVRY steps",&adyn);
  IntParameter("RESULTSEVRY",0,"write results every RESULTSTEVRY steps",&adyn);
  setStringToIntegralParameter<int>("DIVERCONT", "continue",
                                    "What to do if nonlinear solver does not converge?",
                                    tuple<std::string>(
                                        "stop",
                                        "continue"),
                                    tuple<int>(
                                        divcont_stop,
                                        divcont_continue),
                                    &adyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying and mesh sliding algorithm",
                                  tuple<std::string>(
                                    "no",
                                    "meshtying",
                                    "meshsliding"),
                                  tuple<int>(
                                      no_meshtying,
                                      meshtying,
                                      meshsliding),
                                  &adyn);

  // linear solver id used for scalar ale problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for ale problems...",&adyn);

}



void INPAR::ALE::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // Ale update boundary condition

  std::vector<Teuchos::RCP<ConditionComponent> > aleupdatecomponents;

  aleupdatecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPLING")));
  aleupdatecomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "coupling","lagrange",
        Teuchos::tuple<std::string>("lagrange","heightfunction","sphereHeightFunction","meantangentialvelocity","meantangentialvelocityscaled"),
        Teuchos::tuple<std::string>("lagrange","heightfunction","sphereHeightFunction","meantangentialvelocity","meantangentialvelocityscaled"),
        true)));

  aleupdatecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  aleupdatecomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 1)));

  aleupdatecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NODENORMALFUNCT")));
  aleupdatecomponents.push_back(Teuchos::rcp(new IntConditionComponent("nodenormalfunct")));

  Teuchos::RCP<ConditionDefinition> linealeupdate =
    Teuchos::rcp(new ConditionDefinition("DESIGN ALE UPDATE LINE CONDITIONS",
                                         "ALEUPDATECoupling",
                                         "ALEUPDATE Coupling",
                                         DRT::Condition::ALEUPDATECoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfaleupdate =
    Teuchos::rcp(new ConditionDefinition("DESIGN ALE UPDATE SURF CONDITIONS",
                                         "ALEUPDATECoupling",
                                         "ALEUPDATE Coupling",
                                         DRT::Condition::ALEUPDATECoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<aleupdatecomponents.size(); ++i)
  {
    linealeupdate->AddComponent(aleupdatecomponents[i]);
    surfaleupdate->AddComponent(aleupdatecomponents[i]);
  }

  condlist.push_back(linealeupdate);
  condlist.push_back(surfaleupdate);
}


