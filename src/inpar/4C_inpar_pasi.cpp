/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_inpar_pasi.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | set valid parameters for pasi                                             |
 *---------------------------------------------------------------------------*/
void Inpar::PaSI::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& pasidyn = list->sublist("PASI DYNAMIC", false,
      "general control parameters for particle structure interaction problems");

  // time loop control
  Core::UTILS::int_parameter("RESULTSEVRY", 1, "Increment for writing solution", &pasidyn);
  Core::UTILS::int_parameter("RESTARTEVRY", 1, "Increment for writing restart", &pasidyn);
  Core::UTILS::double_parameter("TIMESTEP", 0.01, "Time increment dt", &pasidyn);
  Core::UTILS::int_parameter("NUMSTEP", 100, "Total number of Timesteps", &pasidyn);
  Core::UTILS::double_parameter("MAXTIME", 1.0, "Total simulation time", &pasidyn);

  // type of partitioned coupling
  setStringToIntegralParameter<int>("COUPLING", "partitioned_onewaycoup",
      "partitioned coupling strategies for particle structure interaction",
      tuple<std::string>("partitioned_onewaycoup", "partitioned_twowaycoup",
          "partitioned_twowaycoup_disprelax", "partitioned_twowaycoup_disprelaxaitken"),
      tuple<int>(partitioned_onewaycoup, partitioned_twowaycoup, partitioned_twowaycoup_disprelax,
          partitioned_twowaycoup_disprelaxaitken),
      &pasidyn);

  // partitioned iteration dependent parameters
  Core::UTILS::int_parameter(
      "ITEMAX", 10, "maximum number of partitioned iterations over fields", &pasidyn);

  Core::UTILS::double_parameter("CONVTOLSCALEDDISP", -1.0,
      "tolerance of dof and dt scaled interface displacement increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::double_parameter("CONVTOLRELATIVEDISP", -1.0,
      "tolerance of relative interface displacement increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::double_parameter("CONVTOLSCALEDFORCE", -1.0,
      "tolerance of dof and dt scaled interface force increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::double_parameter("CONVTOLRELATIVEFORCE", -1.0,
      "tolerance of relative interface force increments in partitioned iterations", &pasidyn);

  Core::UTILS::bool_parameter(
      "IGNORE_CONV_CHECK", "no", "ignore convergence check and proceed simulation", &pasidyn);

  // parameters for relaxation
  Core::UTILS::double_parameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &pasidyn);
  Core::UTILS::double_parameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &pasidyn);
  Core::UTILS::double_parameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &pasidyn);
}

FOUR_C_NAMESPACE_CLOSE
