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
void Inpar::PaSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& pasidyn = list->sublist("PASI DYNAMIC", false,
      "general control parameters for particle structure interaction problems");

  // time loop control
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &pasidyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &pasidyn);
  Core::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &pasidyn);
  Core::UTILS::IntParameter("NUMSTEP", 100, "Total number of Timesteps", &pasidyn);
  Core::UTILS::DoubleParameter("MAXTIME", 1.0, "Total simulation time", &pasidyn);

  // type of partitioned coupling
  setStringToIntegralParameter<int>("COUPLING", "partitioned_onewaycoup",
      "partitioned coupling strategies for particle structure interaction",
      tuple<std::string>("partitioned_onewaycoup", "partitioned_twowaycoup",
          "partitioned_twowaycoup_disprelax", "partitioned_twowaycoup_disprelaxaitken"),
      tuple<int>(partitioned_onewaycoup, partitioned_twowaycoup, partitioned_twowaycoup_disprelax,
          partitioned_twowaycoup_disprelaxaitken),
      &pasidyn);

  // partitioned iteration dependent parameters
  Core::UTILS::IntParameter(
      "ITEMAX", 10, "maximum number of partitioned iterations over fields", &pasidyn);

  Core::UTILS::DoubleParameter("CONVTOLSCALEDDISP", -1.0,
      "tolerance of dof and dt scaled interface displacement increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::DoubleParameter("CONVTOLRELATIVEDISP", -1.0,
      "tolerance of relative interface displacement increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::DoubleParameter("CONVTOLSCALEDFORCE", -1.0,
      "tolerance of dof and dt scaled interface force increments in partitioned iterations",
      &pasidyn);

  Core::UTILS::DoubleParameter("CONVTOLRELATIVEFORCE", -1.0,
      "tolerance of relative interface force increments in partitioned iterations", &pasidyn);

  Core::UTILS::BoolParameter(
      "IGNORE_CONV_CHECK", "no", "ignore convergence check and proceed simulation", &pasidyn);

  // parameters for relaxation
  Core::UTILS::DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &pasidyn);
  Core::UTILS::DoubleParameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &pasidyn);
  Core::UTILS::DoubleParameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &pasidyn);
}

FOUR_C_NAMESPACE_CLOSE
