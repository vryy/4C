/*----------------------------------------------------------------------*/
/*!
\file inpar_elemag.cpp

\brief Input parameters for electromagnetic simulations

<pre>
\level 3

\maintainer Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_elemag.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::ELEMAG::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& electromagneticdyn = list->sublist("ELECTROMAGNETIC DYNAMIC",false,"control parameters for electromagnetic problems\n");

  // general settings for time-integration scheme
  DoubleParameter("TIMESTEP",0.01,"Time-step length dt",&electromagneticdyn);
  IntParameter("NUMSTEP",100,"Number of time steps",&electromagneticdyn);
  DoubleParameter("MAXTIME",1.0,"Total simulation time",&electromagneticdyn);

  // additional parameters
  IntParameter("CALCERRORFUNCNO",-1,"Function for error calculation",&electromagneticdyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&electromagneticdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&electromagneticdyn);
  IntParameter("LINEAR_SOLVER",-1,"Number of linear solver used for electromagnetic problem",&electromagneticdyn);
  IntParameter("STARTFUNCNO",-1,"Function for initial field",&electromagneticdyn);
  //IntParameter("SOURCETERMFUNCNO",-1,"Function for source term in volume",&electromagneticdyn);
  //BoolParameter("DOUBLEORFLOAT","Yes","Yes, if evaluation with double, no if with float",&electromagneticdyn);
  //BoolParameter("ALLELESEQUAL","No","Yes, if all elements have same shape and material",&electromagneticdyn);

  // time integration
  setStringToIntegralParameter<int>("TIMEINT","One_Step_Theta",
                    "Type of time integration scheme",
                    tuple<std::string>(
                    "One_Step_Theta"),
                    tuple<int>(
                    elemag_ost),
                    &electromagneticdyn);

  // PML
  //StringParameter("PML_DEFINITION_FILE","none.txt","Filename of file containing the pml definition",&electromagneticdyn);
}


/// set specific electromagnetic conditions
void INPAR::ELEMAG::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

}

