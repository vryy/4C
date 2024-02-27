/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_pasi.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | set valid parameters for pasi                                             |
 *---------------------------------------------------------------------------*/
void INPAR::PASI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& pasidyn = list->sublist("PASI DYNAMIC", false,
      "general control parameters for particle structure interaction problems");

  // time loop control
  CORE::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &pasidyn);
  CORE::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &pasidyn);
  CORE::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &pasidyn);
  CORE::UTILS::IntParameter("NUMSTEP", 100, "Total number of Timesteps", &pasidyn);
  CORE::UTILS::DoubleParameter("MAXTIME", 1.0, "Total simulation time", &pasidyn);

  // type of partitioned coupling
  setStringToIntegralParameter<int>("COUPLING", "partitioned_onewaycoup",
      "partitioned coupling strategies for particle structure interaction",
      tuple<std::string>("partitioned_onewaycoup", "partitioned_twowaycoup",
          "partitioned_twowaycoup_disprelax", "partitioned_twowaycoup_disprelaxaitken"),
      tuple<int>(partitioned_onewaycoup, partitioned_twowaycoup, partitioned_twowaycoup_disprelax,
          partitioned_twowaycoup_disprelaxaitken),
      &pasidyn);

  // partitioned iteration dependent parameters
  CORE::UTILS::IntParameter(
      "ITEMAX", 10, "maximum number of partitioned iterations over fields", &pasidyn);

  CORE::UTILS::DoubleParameter("CONVTOLSCALEDDISP", -1.0,
      "tolerance of dof and dt scaled interface displacement increments in partitioned iterations",
      &pasidyn);

  CORE::UTILS::DoubleParameter("CONVTOLRELATIVEDISP", -1.0,
      "tolerance of relative interface displacement increments in partitioned iterations",
      &pasidyn);

  CORE::UTILS::DoubleParameter("CONVTOLSCALEDFORCE", -1.0,
      "tolerance of dof and dt scaled interface force increments in partitioned iterations",
      &pasidyn);

  CORE::UTILS::DoubleParameter("CONVTOLRELATIVEFORCE", -1.0,
      "tolerance of relative interface force increments in partitioned iterations", &pasidyn);

  CORE::UTILS::BoolParameter(
      "IGNORE_CONV_CHECK", "no", "ignore convergence check and proceed simulation", &pasidyn);

  // parameters for relaxation
  CORE::UTILS::DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &pasidyn);
  CORE::UTILS::DoubleParameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &pasidyn);
  CORE::UTILS::DoubleParameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &pasidyn);
}

BACI_NAMESPACE_CLOSE
