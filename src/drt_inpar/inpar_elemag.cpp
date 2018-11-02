/*----------------------------------------------------------------------*/
/*!
\file inpar_elemag.cpp

\brief Input parameters for electromagnetic simulations

<pre>
\level 3

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_elemag.H"
#include "../drt_lib/drt_conditiondefinition.H"

void INPAR::ELEMAG::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& electromagneticdyn = list->sublist(
      "ELECTROMAGNETIC DYNAMIC", false, "control parameters for electromagnetic problems\n");

  // general settings for time-integration scheme
  DoubleParameter("TIMESTEP", 0.01, "Time-step length dt", &electromagneticdyn);
  IntParameter("NUMSTEP", 100, "Number of time steps", &electromagneticdyn);
  DoubleParameter("MAXTIME", 1.0, "Total simulation time", &electromagneticdyn);

  // additional parameters
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &electromagneticdyn);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &electromagneticdyn);
  IntParameter("LINEAR_SOLVER", -1, "Number of linear solver used for electromagnetic problem",
      &electromagneticdyn);
  IntParameter("STARTFUNCNO", -1, "Function for initial field", &electromagneticdyn);
  IntParameter("SOURCEFUNCNO", -1, "Function for source term in volume", &electromagneticdyn);
  // BoolParameter("DOUBLEORFLOAT","Yes","Yes, if evaluation with double, no if with
  // float",&electromagneticdyn); BoolParameter("ALLELESEQUAL","No","Yes, if all elements have same
  // shape and material",&electromagneticdyn);

  {
    // time integration

    Teuchos::Tuple<std::string, 5> name;
    Teuchos::Tuple<int, 5> label;
    name[0] = "One_Step_Theta";
    label[0] = elemag_ost;
    name[1] = "Implicit_Euler";
    label[1] = elemag_implicit_euler;
    name[2] = "Explicit_Euler";
    label[2] = elemag_explicit_euler;
    name[3] = "Runge_Kutta";
    label[3] = elemag_rk;
    name[4] = "Crank_Nicolson";
    label[4] = elemag_cn;

    setStringToIntegralParameter<int>("TIMEINT", "One_Step_Theta",
        "Type of time integration scheme", name, label, &electromagneticdyn);
  }

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 2> name;
    Teuchos::Tuple<int, 2> label;
    name[0] = "zero_field";
    label[0] = initfield_zero_field;
    name[1] = "field_by_function";
    label[1] = initfield_field_by_function;

    setStringToIntegralParameter<int>("INITIALFIELD", "zero_field", "Initial field for ele problem",
        name, label, &electromagneticdyn);

    // Error calculation
    BoolParameter("CALCERR", "No", "Calc the error wrt ERRORFUNCNO?", &electromagneticdyn);
  }

  IntParameter("ERRORFUNCNO", -1, "Function for error calculation", &electromagneticdyn);

  // PML
  // StringParameter("PML_DEFINITION_FILE","none.txt","Filename of file containing the pml
  // definition",&electromagneticdyn);
}

/// set specific electromagnetic conditions
void INPAR::ELEMAG::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;
}
