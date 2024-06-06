/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for immersed

\level 1


*/

/*----------------------------------------------------------------------*/


#include "4C_inpar_immersed.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::Immersed::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  Teuchos::ParameterList& immersedmethod =
      list->sublist("IMMERSED METHOD", false, "General parameters for any immersed problem");

  Teuchos::Tuple<std::string, 3> coupname;
  Teuchos::Tuple<int, 3> couplabel;

  coupname[0] = "basic_sequ_stagg";
  couplabel[0] = cell_basic_sequ_stagg;
  coupname[1] = "iter_stagg_fixed_rel_param";
  couplabel[1] = cell_iter_stagg_fixed_rel_param;
  coupname[2] = "iter_stagg_AITKEN_rel_param";
  couplabel[2] = cell_iter_stagg_AITKEN_rel_param;



  setStringToIntegralParameter<int>("COUPALGO", "partitioned",
      "Coupling strategies for immersed method.", tuple<std::string>("partitioned", "monolithic"),
      tuple<int>(partitioned, monolithic), &immersedmethod);

  setStringToIntegralParameter<int>("SCHEME", "dirichletneumann",
      "Coupling schemes for partitioned immersed method.",
      tuple<std::string>("neumannneumann", "dirichletneumann"),
      tuple<int>(neumannneumann, dirichletneumann), &immersedmethod);

  setStringToIntegralParameter<int>("DIVERCONT", "stop", "What to do after maxiter is reached.",
      tuple<std::string>("stop", "continue"), tuple<int>(nlnsolver_stop, nlnsolver_continue),
      &immersedmethod);

  setStringToIntegralParameter<int>("OUTPUT_EVRY_NLNITER", "no",
      "write output after every solution step of the nonlin. part. iter. scheme",
      tuple<std::string>("yes", "no"), tuple<int>(1, 0), &immersedmethod);

  setStringToIntegralParameter<int>("CORRECT_BOUNDARY_VELOCITIES", "no",
      "correct velocities in fluid elements cut by surface of immersed structure",
      tuple<std::string>("yes", "no"), tuple<int>(1, 0), &immersedmethod);

  setStringToIntegralParameter<int>("DEFORM_BACKGROUND_MESH", "no",
      "switch between immersed with fixed or deformable background mesh",
      tuple<std::string>("yes", "no"), tuple<int>(1, 0), &immersedmethod);

  setStringToIntegralParameter<int>("TIMESTATS", "everyiter",
      "summarize time monitor every nln iteration", tuple<std::string>("everyiter", "endofsim"),
      tuple<int>(1, 0), &immersedmethod);

  Core::UTILS::DoubleParameter(
      "FLD_SRCHRADIUS_FAC", 1.0, "fac times fluid ele. diag. length", &immersedmethod);
  Core::UTILS::DoubleParameter(
      "STRCT_SRCHRADIUS_FAC", 0.5, "fac times structure bounding box diagonal", &immersedmethod);
  Core::UTILS::IntParameter("NUM_GP_FLUID_BOUND", 8,
      "number of gp in fluid elements cut by surface of immersed structure (higher number yields "
      "better mass conservation)",
      &immersedmethod);

  /*----------------------------------------------------------------------*/
  /* parameters for paritioned immersed solvers */
  Teuchos::ParameterList& immersedpart = immersedmethod.sublist("PARTITIONED SOLVER", false, "");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &immersedpart);

  setStringToIntegralParameter<int>("PREDICTOR", "d(n)", "Predictor for interface displacements",
      tuple<std::string>(
          "d(n)", "d(n)+dt*(1.5*v(n)-0.5*v(n-1))", "d(n)+dt*v(n)", "d(n)+dt*v(n)+0.5*dt^2*a(n)"),
      tuple<int>(1, 2, 3, 4), &immersedpart);

  setStringToIntegralParameter<int>("COUPVARIABLE", "Displacement",
      "Coupling variable at the fsi interface", tuple<std::string>("Displacement", "Force"),
      tuple<int>(0, 1), &immersedpart);

  Core::UTILS::DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &immersedpart);
  Core::UTILS::DoubleParameter(
      "RELAX", 1.0, "fixed relaxation parameter for partitioned FSI solvers", &immersedpart);
  Core::UTILS::DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &immersedpart);
  Core::UTILS::IntParameter(
      "ITEMAX", 100, "Maximum number of iterations over fields", &immersedpart);
}



void Inpar::Immersed::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // IMMERSED FSI

  Teuchos::RCP<Core::Conditions::ConditionDefinition> immersedsearchbox =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOLUME IMMERSED SEARCHBOX",
          "ImmersedSearchbox", "Immersed Searchbox", Core::Conditions::ImmersedSearchbox, true,
          Core::Conditions::geometry_type_volume));

  condlist.push_back(immersedsearchbox);

  /*--------------------------------------------------------------------*/
  // IMMERSED COUPLING

  std::vector<Teuchos::RCP<Input::LineComponent>> immersedcomponents;

  immersedcomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineimmersed =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN IMMERSED COUPLING LINE CONDITIONS", "IMMERSEDCoupling", "IMMERSED Coupling",
          Core::Conditions::IMMERSEDCoupling, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfimmersed =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN IMMERSED COUPLING SURF CONDITIONS", "IMMERSEDCoupling", "IMMERSED Coupling",
          Core::Conditions::IMMERSEDCoupling, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < immersedcomponents.size(); ++i)
  {
    lineimmersed->AddComponent(immersedcomponents[i]);
    surfimmersed->AddComponent(immersedcomponents[i]);
  }

  condlist.push_back(lineimmersed);
  condlist.push_back(surfimmersed);
}

FOUR_C_NAMESPACE_CLOSE
