/*-----------------------------------------------------------*/
/*!
\file inpar_contact_xcontact.cpp

\brief Input parameters for the eXtended contact approach

\maintainer Michael Hiermeier

\date Jul 14, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "inpar_contact_xcontact.H"

#include "drt_validparameters.H"

void INPAR::XCONTACT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  // get the mutable xcontact sub-list
  Teuchos::ParameterList& p_xcontact = list->sublist("XCONTACT DYNAMIC", false,
      "Configuration and control of the eXtended contact dynamics simulation.");

  IntParameter("NUMSTEP",10,"Number of Time Steps",&p_xcontact);
  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&p_xcontact);
  DoubleParameter("MAXTIME",0.0,"Total simulation time",&p_xcontact);
  DoubleParameter("CONVTOL",1E-6,"Tolerance for convergence check",&p_xcontact);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&p_xcontact);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&p_xcontact);

  setStringToIntegralParameter<int>("COUPL_ALGO","None",
      "Define the algorithm to couple the level-set and the structural contact approaches.",
      tuple<std::string>("None","none",
                         "Monolithic","monolithic",
                         "Partitioned","partitioned"),
      tuple<int>(
                 field_coupl_algo_undefined, field_coupl_algo_undefined,
                 field_coupl_algo_monolithic,field_coupl_algo_monolithic,
                 field_coupl_algo_partitioned,field_coupl_algo_partitioned),
      &p_xcontact);

}
