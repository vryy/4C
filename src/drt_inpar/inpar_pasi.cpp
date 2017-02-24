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
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& stidyn = list->sublist(
      "PASI DYNAMIC",
      false,
      "general control parameters for particle structure interaction problems"
      );

  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&stidyn);

  IntParameter("RESTARTEVRY",1, "Increment for writing restart", &stidyn);

  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&stidyn);

  IntParameter("NUMSTEP",100, "Total number of Timesteps", &stidyn);

  DoubleParameter("MAXTIME",1000.0, "Total simulation time", &stidyn);

  setStringToIntegralParameter<int>(
                               "COUPALGO","partitioned_onewaycoup",
                               "Coupling strategies for particle structure interaction",
                               tuple<std::string>(
                                 "partitioned_onewaycoup"),
                                 tuple<int>(
                                 partitioned_onewaycoup),
                                 &stidyn);

} // INPAR::PASI::SetValidParameters


/*----------------------------------------------------------------------*
 | set valid conditions for pasi                         sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
void INPAR::PASI::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  return;
} // INPAR::PASI::SetValidConditions
