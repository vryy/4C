/*----------------------------------------------------------------------*/
/*!
\file inpar_cut.cpp

\brief Input parameters for cut library

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_cut.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::CUT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  //  Teuchos::Array<std::string> yesnotuple =
  //      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  //  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& cut_general = list->sublist("CUT GENERAL", false, "");

  // Intersection precision (double or cln)
  setStringToIntegralParameter<int>("KERNEL_INTERSECTION_FLOATTYPE", "double",
      "The floattype of the cut surface-edge intersection", tuple<std::string>("cln", "double"),
      tuple<int>(INPAR::CUT::floattype_cln, INPAR::CUT::floattype_double), &cut_general);

  // Computing disctance surface to point precision (double or cln)
  setStringToIntegralParameter<int>("KERNEL_DISTANCE_FLOATTYPE", "double",
      "The floattype of the cut distance computation", tuple<std::string>("cln", "double"),
      tuple<int>(INPAR::CUT::floattype_cln, INPAR::CUT::floattype_double), &cut_general);

  // A general floattype for GEO::CUT::Position for Embedded Elements (ComputeDistance)
  // If specified this floattype is used for all computations of GEO::CUT::Position with embedded
  // elements
  setStringToIntegralParameter<int>("GENERAL_POSITON_DISTANCE_FLOATTYPE", "none",
      "A general floattype for GEO::CUT::Position for Embedded Elements (ComputeDistance)",
      tuple<std::string>("none", "cln", "double"),
      tuple<int>(
          INPAR::CUT::floattype_none, INPAR::CUT::floattype_cln, INPAR::CUT::floattype_double),
      &cut_general);

  // A general floattype for GEO::CUT::Position for Elements (ComputePosition)
  // If specified this floattype is used for all computations of GEO::CUT::Position
  setStringToIntegralParameter<int>("GENERAL_POSITON_POSITION_FLOATTYPE", "none",
      "A general floattype for GEO::CUT::Position Elements (ComputePosition)",
      tuple<std::string>("none", "cln", "double"),
      tuple<int>(
          INPAR::CUT::floattype_none, INPAR::CUT::floattype_cln, INPAR::CUT::floattype_double),
      &cut_general);
}
