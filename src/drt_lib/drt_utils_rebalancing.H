/*---------------------------------------------------------------------*/
/*!

\brief A collection of helper methods related to partitioning and parallel distribution

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_UTILS_REBALANCING_H
#define DRT_UTILS_REBALANCING_H

#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// forward declarations
namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace Isorropia
{
  namespace Epetra
  {
    class CostDescriber;
  }  // namespace Epetra
}  // namespace Isorropia

namespace DRT
{
  namespace UTILS
  {
    //! @name Rebalancing of nodal maps
    //! @{

    /*!
    \brief Compute rebalanced node maps for a given discretization

    Starting from a given discretization, a new node row and column map is computed, such that
    they both are better balanced across a given number partitions.

    \note This just computes the new node row/column maps, but does not perform any redistribution
    of data among ranks.

    @param[in] discretization Discretization whose graph shall be redistributed
    @param[in] elementRowMap Element row map of the discretization prior to the rebalancing
    @param[out] nodeRowMap Node row map after redistribution
    @param[out] nodeColumnMap Node column map after redistribution
    @param[in] comm Communicator (might be different from the one inside discretization)
    @param[in] outflag Flag to switch screen output on/off
    @param[in] numPartitions Desired number of partitions
    @param[in] imbalanceTol Max. relative imbalance of subdomain size
    */
    void ComputeRebalancedNodeMaps(Teuchos::RCP<DRT::Discretization> discretization,
        Teuchos::RCP<const Epetra_Map> elementRowMap, Teuchos::RCP<Epetra_Map>& nodeRowMap,
        Teuchos::RCP<Epetra_Map>& nodeColumnMap, Teuchos::RCP<const Epetra_Comm> comm,
        const bool outflag, const int numPartitions, const double imbalanceTol = 1.1);

    /*!
    \brief Compute rebalanced node maps for a given discretization

    Starting from a given discretization, a new node row and column map is computed, such that
    they both are better balanced across a given number partitions.

    \note This just computes the new node row/column maps, but does not perform any redistribution
    of data among ranks.

    \pre Input discretization has to be FillComplete().

    @param[in] dis Discretization whose node maps are to be rebalanced
    @param[out] rownodes Node row map after rebalancing
    @param[out] colnodes Node column map after rebalancing
    @param[in] outflag Flag to switch screen output on/off
    */
    void ComputeRebalancedNodeMaps(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes, const bool outflag);

    /*!
    \brief Compute rebalanced node maps for a given discretization considering weights to model
    costs

    Starting from a given discretization, a new node row and column map is computed, such that
    they both are better balanced across a given number partitions.

    \note This just computes the new node row/column maps, but does not perform any redistribution
    of data among ranks.

    \pre Input discretization has to be FillComplete().

    @param[in] dis Discretization whose node maps are to be rebalanced
    @param[out] rownodes Node row map after rebalancing
    @param[out] colnodes Node column map after rebalancing
    @param[in] outflag Flag to switch screen output on/off
    */
    void ComputeRebalancedNodeMapsUsingWeights(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes, const bool outflag);

    //! @}

    //! @name Redistribute discretizations
    //! @{

    /*!
    \brief Export and FillComplete a discretization based on given nodal row and column maps

    Given a nodal row and column map, this routines creates the corresponding element maps. We
    then export nodes and elements of the discretization to their new maps. Finally, the
    discretization has to perform a FillComplete().

    \post Discretization is FillComplete().

    @param[in/out] discretization Discretization to be redistributed
    @param[in] noderowmap Given target node row map
    @param[in] nodecolmap Given target node column map
    @param[in] assigndegreesoffreedom Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] initelements Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] doboundaryconditions Flag forwarded to DRT::Discretization::FillComplete()
    */
    void ExportAndFillCompleteDiscretization(DRT::Discretization& discretization,
        const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
        const bool assigndegreesoffreedom, const bool initelements,
        const bool doboundaryconditions);

    /*!
    \brief Redistribute and FillComplete() a discretization to a new partitioning using weights to
    model costs

    1. Compute node row and column maps with a better and more balanced partitioning
    2. Redistriibute and FillComplete() the discretization.

    @param[in/out] discretization Discretization to be repartitioned
    @param[in] assigndegreesoffreedom Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] initelements Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] doboundaryconditions Flag forwarded to DRT::Discretization::FillComplete()

    \sa DRT::UTILS::ComputeRebalancedNodeMapsUsingWeights()
    \sa DRT::Discretization::FillComplete()
    */
    void RedistributeAndFillCompleteDiscretizationUsingWeights(
        Teuchos::RCP<DRT::Discretization> discretization, const bool assigndegreesoffreedom,
        const bool initelements, const bool doboundaryconditions);

    //! @}

    //! @name Perform actual rebalancing via Isorropia/Zoltan
    //! @{

    /*!
    \brief Rebalance graph

    Use Isorropia package to access Zoltan. By default, Isorropia will use Zoltan hypergraph
    partitioning, treating the graph columns as hyper-edges and the graph rows as vertices.

    \note The rebalanced graph will be FillComplete().

    @param[in] initialGraph Initial graph to be rebalanced
    @param[in] rebalanceParams Parameter list with rebalancing options

    @return Rebalanced graph
    */
    Teuchos::RCP<Epetra_CrsGraph> RebalanceGraph(Teuchos::RCP<const Epetra_CrsGraph> initialGraph,
        const Teuchos::ParameterList& rebalanceParams);

    /*!
    \brief Rebalance graph using a cost model

    Use Isorropia package to access Zoltan. By default, Isorropia will use Zoltan hypergraph
    partitioning, treating the graph columns as hyper-edges and the graph rows as vertices.

    \note The rebalanced graph will be FillComplete().

    @param[in] initialGraph Initial graph to be rebalanced
    @param[in] costs Isorropia cost describer
    @param[in] rebalanceParams Parameter list with rebalancing options

    @return Rebalanced graph

    \sa DRT::UTILS::SetupCostDescriber()
    */
    Teuchos::RCP<Epetra_CrsGraph> RebalanceGraph(Teuchos::RCP<const Epetra_CrsGraph> initialGraph,
        Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs,
        const Teuchos::ParameterList& rebalanceParams);

    //! @}

    /*!
    \brief Create parameter list with rebalancing options for Isorropia/Zoltan

    @param[in] numPartitions Desired number of partitions
    @param[in] Max. relative imbalance of subdomain size
    */
    Teuchos::RCP<Teuchos::ParameterList> CreateRebalancingParameterList(
        const int numPartitions, const double imbalanceTol);

    /*!
    \brief Create Isorropia cost describer based on element connectivity

    @param[in] discretization Discretization used to build the cost model

    @return Cost describer ready to be used for repartitioning
    */
    Teuchos::RCP<Isorropia::Epetra::CostDescriber> SetupCostDescriber(
        const DRT::Discretization& discretization);

    /*!
    \brief Build graph of a given  discretization

    \pre Discretization does not have to be FillComplete().

    @param[in] dis Discretization whose graph will be build
    @param[in] roweles Element row map of this discretization
    @param[out] rownodes Node row map
    @param[in] outflag Flag to switch screen output on/off

    @return Uncompleted graph of input discretization
    */
    Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<const Epetra_Map> roweles, Teuchos::RCP<Epetra_Map>& rownodes,
        Teuchos::RCP<const Epetra_Comm> comm, const bool outflag);

  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifndef DRT_UTILS_REBALANCING_H
