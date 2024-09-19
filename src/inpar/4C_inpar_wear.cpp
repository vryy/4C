/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for wear

\level 1

*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_wear.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Wear::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for wear */
  Teuchos::ParameterList& wear = list->sublist("WEAR", false, "");

  setStringToIntegralParameter<WearLaw>("WEARLAW", "None", "Type of wear law",
      tuple<std::string>("None", "none", "Archard", "archard"),
      tuple<WearLaw>(wear_none, wear_none, wear_archard, wear_archard), &wear);

  Core::UTILS::bool_parameter("MATCHINGGRID", "Yes", "is matching grid", &wear);

  setStringToIntegralParameter<WearShape>("WEAR_SHAPEFCN", "std",
      "Type of employed set of shape functions for wear",
      tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
      tuple<WearShape>(wear_shape_dual, wear_shape_dual, wear_shape_standard, wear_shape_standard,
          wear_shape_standard),
      &wear);

  Core::UTILS::double_parameter("WEARCOEFF", 0.0, "Wear coefficient for slave surface", &wear);
  Core::UTILS::double_parameter(
      "WEARCOEFF_MASTER", 0.0, "Wear coefficient for master surface", &wear);
  Core::UTILS::double_parameter(
      "WEAR_TIMERATIO", 1.0, "Time step ratio between wear and spatial time scale", &wear);
  Core::UTILS::double_parameter("SSSLIP", 1.0, "Fixed slip for steady state wear", &wear);

  Core::UTILS::bool_parameter("SSWEAR", "No", "flag for steady state wear", &wear);

  Core::UTILS::bool_parameter(
      "VOLMASS_OUTPUT", "No", "flag for output of mass/volume in ref,mat and cur. conf.", &wear);

  setStringToIntegralParameter<WearSide>("WEAR_SIDE", "slave", "Definition of wear side",
      tuple<std::string>("s", "slave", "Slave", "both", "slave_master", "sm"),
      tuple<WearSide>(wear_slave, wear_slave, wear_slave, wear_both, wear_both, wear_both), &wear);

  setStringToIntegralParameter<WearType>("WEARTYPE", "internal_state",
      "Definition of wear algorithm",
      tuple<std::string>("intstate", "is", "internal_state", "primvar", "pv", "primary_variable"),
      tuple<WearType>(
          wear_intstate, wear_intstate, wear_intstate, wear_primvar, wear_primvar, wear_primvar),
      &wear);

  setStringToIntegralParameter<WearTimInt>("WEARTIMINT", "explicit",
      "Definition of wear time integration",
      tuple<std::string>("explicit", "e", "expl", "implicit", "i", "impl"),
      tuple<WearTimInt>(wear_expl, wear_expl, wear_expl, wear_impl, wear_impl, wear_impl), &wear);

  setStringToIntegralParameter<WearTimeScale>("WEAR_TIMESCALE", "equal",
      "Definition wear time scale compares to std. time scale",
      tuple<std::string>("equal", "e", "different", "d"),
      tuple<WearTimeScale>(
          wear_time_equal, wear_time_equal, wear_time_different, wear_time_different),
      &wear);
}

FOUR_C_NAMESPACE_CLOSE
