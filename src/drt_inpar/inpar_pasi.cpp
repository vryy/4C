/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "inpar_pasi.H"
#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | set valid parameters for pasi                                             |
 *---------------------------------------------------------------------------*/
void INPAR::PASI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& pasidyn = list->sublist("PASI DYNAMIC", false,
      "general control parameters for particle structure interaction problems");

  // time loop control
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &pasidyn);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &pasidyn);
  DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &pasidyn);
  IntParameter("NUMSTEP", 100, "Total number of Timesteps", &pasidyn);
  DoubleParameter("MAXTIME", 1.0, "Total simulation time", &pasidyn);

  // type of partitioned coupling
  setStringToIntegralParameter<int>("COUPLING", "partitioned_onewaycoup",
      "partitioned coupling strategies for particle structure interaction",
      tuple<std::string>("partitioned_onewaycoup", "partitioned_twowaycoup",
          "partitioned_twowaycoup_disprelax", "partitioned_twowaycoup_disprelaxaitken"),
      tuple<int>(partitioned_onewaycoup, partitioned_twowaycoup, partitioned_twowaycoup_disprelax,
          partitioned_twowaycoup_disprelaxaitken),
      &pasidyn);

  // partitioned iteration dependent parameters
  IntParameter("ITEMAX", 10, "maximum number of partitioned iterations over fields", &pasidyn);

  DoubleParameter("CONVTOLSCALEDDISP", -1.0,
      "tolerance of dof and dt scaled interface displacement increments in partitioned iterations",
      &pasidyn);

  DoubleParameter("CONVTOLRELATIVEDISP", -1.0,
      "tolerance of relative interface displacement increments in partitioned iterations",
      &pasidyn);

  DoubleParameter("CONVTOLSCALEDFORCE", -1.0,
      "tolerance of dof and dt scaled interface force increments in partitioned iterations",
      &pasidyn);

  DoubleParameter("CONVTOLRELATIVEFORCE", -1.0,
      "tolerance of relative interface force increments in partitioned iterations", &pasidyn);

  BoolParameter(
      "IGNORE_CONV_CHECK", "no", "ignore convergence check and proceed simulation", &pasidyn);

  // parameters for relaxation
  DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &pasidyn);
  DoubleParameter("MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &pasidyn);
  DoubleParameter("MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &pasidyn);
}
