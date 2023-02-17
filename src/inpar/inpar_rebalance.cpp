/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for rebalancing the discretization

\level 2

*/
/*-----------------------------------------------------------*/

#include "inpar_rebalance.H"
#include "inpar_parameterlist_utils.H"

void INPAR::REBALANCE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  using namespace DRT::INPUT;

  Teuchos::ParameterList& meshpartitioning = list->sublist("MESH PARTITIONING", false, "");

  setStringToIntegralParameter<RebalanceType>("METHOD", "hypergraph",
      "Type of rebalance/partition algorithm to be used for decomposing the entire mesh into "
      "subdomains for parallel computing.",
      tuple<std::string>("none", "hypergraph", "recursive_coordinate_bisection"),
      tuple<RebalanceType>(RebalanceType::none, RebalanceType::hypergraph,
          RebalanceType::recursive_coordinate_bisection),
      &meshpartitioning);

  DoubleParameter("IMBALANCE_TOL", 1.1,
      "Tolerance for relative imbalance of subdomain sizes for graph partitioning of unstructured "
      "meshes read from input files.",
      &meshpartitioning);
}
