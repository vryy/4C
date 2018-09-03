/*!----------------------------------------------------------------------
\file inpar_pasi.cpp

\brief Input parameters for particle structure interaction problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "inpar_pasi.H"
#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_conditiondefinition.H"

/*----------------------------------------------------------------------*
 | set valid parameters for pasi                         sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void INPAR::PASI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& pasidyn = list->sublist("PASI DYNAMIC", false,
      "general control parameters for particle structure interaction problems");

  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &pasidyn);

  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &pasidyn);

  DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &pasidyn);

  IntParameter("NUMSTEP", 100, "Total number of Timesteps", &pasidyn);

  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &pasidyn);

  setStringToIntegralParameter<int>("COUPALGO", "partitioned_onewaycoup",
      "Coupling strategies for particle structure interaction",
      tuple<std::string>("partitioned_onewaycoup", "partitioned_twowaycoup",
          "partitioned_twowaycoup_forcerelax", "partitioned_twowaycoup_forcerelaxaitken"),
      tuple<int>(partitioned_onewaycoup, partitioned_twowaycoup, partitioned_twowaycoup_forcerelax,
          partitioned_twowaycoup_forcerelaxaitken),
      &pasidyn);

  IntParameter("ITEMAX", 10, "Maximum number of iterations over fields", &pasidyn);

  // paramters for partitioned PASI
  Teuchos::ParameterList& pasidynpart = pasidyn.sublist("PARTITIONED", false,
      "Partitioned Particle Structure Interaction\n"
      "Control section for partitioned PASI");

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration within partitioned PASI", &pasidynpart);

  // ignore convergence check and proceed simulation
  BoolParameter(
      "IGNORE_CONV_CHECK", "no", "ignore convergence check and proceed simulation", &pasidynpart);

  // parameters for relaxation of partitioned PASI
  DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &pasidynpart);
  DoubleParameter("MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &pasidynpart);
  DoubleParameter("MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &pasidynpart);

}  // INPAR::PASI::SetValidParameters

/*----------------------------------------------------------------------*
 | set valid conditions for pasi                         sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void INPAR::PASI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  return;
}  // INPAR::PASI::SetValidConditions
