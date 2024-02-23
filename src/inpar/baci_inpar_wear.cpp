/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for wear

\level 1

*/

/*----------------------------------------------------------------------*/



#include "baci_inpar_wear.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN



void INPAR::WEAR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for wear */
  Teuchos::ParameterList& wear = list->sublist("WEAR", false, "");

  setStringToIntegralParameter<int>("WEARLAW", "None", "Type of wear law",
      tuple<std::string>("None", "none", "Archard", "archard"),
      tuple<int>(wear_none, wear_none, wear_archard, wear_archard), &wear);

  CORE::UTILS::BoolParameter("MATCHINGGRID", "Yes", "is matching grid", &wear);

  setStringToIntegralParameter<int>("WEARCOEFF_CONF", "material",
      "configuration in which wcoeff is defined",
      tuple<std::string>("material", "mat", "spatial", "sp"),
      tuple<int>(wear_coeff_mat, wear_coeff_mat, wear_coeff_sp, wear_coeff_sp), &wear);

  setStringToIntegralParameter<int>("WEAR_SHAPE_EVO", "material",
      "configuration for shape evolution step",
      tuple<std::string>("material", "mat", "spatial", "sp"),
      tuple<int>(wear_se_mat, wear_se_mat, wear_se_sp, wear_se_sp), &wear);

  setStringToIntegralParameter<int>("WEAR_SHAPEFCN", "std",
      "Type of employed set of shape functions for wear",
      tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
      tuple<int>(wear_shape_dual, wear_shape_dual, wear_shape_standard, wear_shape_standard,
          wear_shape_standard),
      &wear);

  CORE::UTILS::DoubleParameter("WEARCOEFF", 0.0, "Wear coefficient for slave surface", &wear);
  CORE::UTILS::DoubleParameter(
      "WEARCOEFF_MASTER", 0.0, "Wear coefficient for master surface", &wear);
  CORE::UTILS::DoubleParameter(
      "WEAR_TIMERATIO", 1.0, "Time step ratio between wear and spatial time scale", &wear);
  CORE::UTILS::DoubleParameter("SSSLIP", 1.0, "Fixed slip for steady state wear", &wear);

  CORE::UTILS::BoolParameter("SSWEAR", "No", "flag for steady state wear", &wear);

  CORE::UTILS::BoolParameter(
      "VOLMASS_OUTPUT", "No", "flag for output of mass/volume in ref,mat and cur. conf.", &wear);

  setStringToIntegralParameter<int>("WEAR_SIDE", "slave", "Definition of wear side",
      tuple<std::string>("s", "slave", "Slave", "both", "slave_master", "sm"),
      tuple<int>(wear_slave, wear_slave, wear_slave, wear_both, wear_both, wear_both), &wear);

  setStringToIntegralParameter<int>("WEARTYPE", "internal_state", "Definition of wear algorithm",
      tuple<std::string>("intstate", "is", "internal_state", "primvar", "pv", "primary_variable"),
      tuple<int>(
          wear_intstate, wear_intstate, wear_intstate, wear_primvar, wear_primvar, wear_primvar),
      &wear);

  setStringToIntegralParameter<int>("WEARTIMINT", "explicit", "Definition of wear time integration",
      tuple<std::string>("explicit", "e", "expl", "implicit", "i", "impl"),
      tuple<int>(wear_expl, wear_expl, wear_expl, wear_impl, wear_impl, wear_impl), &wear);

  setStringToIntegralParameter<int>("WEAR_COUPALGO", "stagg",
      "Definition of wear (ALE) coupling algorithm",
      tuple<std::string>("stagg", "s", "iterstagg", "is", "monolithic", "mono"),
      tuple<int>(
          wear_stagg, wear_stagg, wear_iterstagg, wear_iterstagg, wear_monolithic, wear_monolithic),
      &wear);

  setStringToIntegralParameter<int>("WEAR_TIMESCALE", "equal",
      "Definition wear time scale compares to std. time scale",
      tuple<std::string>("equal", "e", "different", "d"),
      tuple<int>(wear_time_equal, wear_time_equal, wear_time_different, wear_time_different),
      &wear);
}

BACI_NAMESPACE_CLOSE
