/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of functions related to partitioning and parallel distribution

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_REBALANCE_HPP
#define FOUR_C_REBALANCE_HPP

#include "baci_config.hpp"

#include "baci_discretization_geometric_search_params.hpp"

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::REBALANCE
{
  /*!
  \brief Compute rebalanced node maps for a given discretization considering weights to model
  costs

  Starting from a given discretization, a new node row and column map is computed, such that
  they both are better balanced across a given number partitions.

  \note This just computes the new node row/column maps, but does not perform any redistribution
  of data among ranks.

  @param[in] initialGraph Initial graph used for rebalancing
  @param[in] rebalanceParams Parameter list with rebalancing options
  @param[in] initialNodeWeights Initial weights of the graph nodes
  @param[in] initialEdgeWeights Initial weights of the graph edges
  @param[in] initialNodeCoordinates Coordinates of the discretization

  @return Node row map and node column map after rebalancing with weights
  */
  std::pair<Teuchos::RCP<Epetra_Map>, Teuchos::RCP<Epetra_Map>> RebalanceNodeMaps(
      Teuchos::RCP<const Epetra_CrsGraph> initialGraph,
      const Teuchos::ParameterList& rebalanceParams,
      const Teuchos::RCP<Epetra_Vector>& initialNodeWeights = Teuchos::null,
      const Teuchos::RCP<Epetra_CrsMatrix>& initialEdgeWeights = Teuchos::null,
      const Teuchos::RCP<Epetra_MultiVector>& initialNodeCoordinates = Teuchos::null);

  /*!
  \brief Rebalance graph using node and edge weights based on the initial graph

  The partitioning will be done based on the method given in the parameter list. This method
  only makes sense with graph or hypergraph partitioning.

  @note Use Isorropia package to access Zoltan. By default, Isorropia will use Zoltan hypergraph
  partitioning, treating the graph columns as hyper-edges and the graph rows as vertices. The
   rebalanced graph will be FillComplete().

  @param[in] initialGraph Initial graph to be rebalanced
  @param[in] rebalanceParams Parameter list with rebalancing options
  @param[in] initialNodeWeights Initial weights of the graph nodes
  @param[in] initialEdgeWeights Initial weights of the graph edges
  @param[in] initialNodeCoordinates Coordinates of the discretization

  @return Rebalanced graph
  */
  Teuchos::RCP<Epetra_CrsGraph> RebalanceGraph(const Epetra_CrsGraph& initialGraph,
      const Teuchos::ParameterList& rebalanceParams,
      const Teuchos::RCP<Epetra_Vector>& initialNodeWeights = Teuchos::null,
      const Teuchos::RCP<Epetra_CrsMatrix>& initialEdgeWeights = Teuchos::null,
      const Teuchos::RCP<Epetra_MultiVector>& initialNodeCoordinates = Teuchos::null);

  /*!
  \brief Rebalance coordinates using weights based on the initial coordinates

  This method only makes sense with geometric partitioning methods like RCB. The partitioner
  will figure things out in the background.

  @param[in] initialCoordinates Initial coordinates to be rebalanced
  @param[in] initialWeights Initial weights of the coordinates
  @param[in] rebalanceParams Parameter list with rebalancing options

  @return Rebalanced coordinates
  */
  std::pair<Teuchos::RCP<Epetra_MultiVector>, Teuchos::RCP<Epetra_MultiVector>>
  RebalanceCoordinates(const Epetra_MultiVector& initialCoordinates,
      const Teuchos::ParameterList& rebalanceParams, const Epetra_MultiVector& initialWeights);

  /*!
  \brief Create node and edge weights based on element connectivity

  @param[in] dis Discretization used to build the weights

  @return Node and edge weights to be used for repartitioning
  */
  std::pair<Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_CrsMatrix>> BuildWeights(
      const DRT::Discretization& dis);

  /*!
  \brief Build node graph of a given  discretization

  \pre Discretization does NOT have to be FillComplete().

  @param[in] dis Discretization whose node graph will be build
  @param[in] roweles Element row map of this discretization

  @return Uncompleted node graph of input discretization
  */
  Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(
      Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<const Epetra_Map> roweles);

  /*!
  \brief Build monolithic node graph of a given discretization

  The monolithic graph is build by using a global collision search on the reference configuration,
  to obtain information about close elements. Based on this information, additional edges are
  build into the graph.

  \pre Discretization has to be FillComplete()!

  @param[in] dis Discretization whose monolithic node graph will be build

  @return Completed monolithic node graph of input discretization
  */
  Teuchos::RCP<const Epetra_CrsGraph> BuildMonolithicNodeGraph(
      const DRT::Discretization& dis, const CORE::GEOMETRICSEARCH::GeometricSearchParams& params);

}  // namespace CORE::REBALANCE

BACI_NAMESPACE_CLOSE

#endif
