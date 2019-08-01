/*---------------------------------------------------------------------*/
/*!

\brief A collection of helper methods related to partitioning and parallel distribution

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef UTILS_PARMETIS_H
#define UTILS_PARMETIS_H

#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace UTILS
  {
    /*!
    \brief Redistribute the underlying graph of a discretization to a new and more balanced
    partitioning

    @param[in] discretization Discretization whose graph shall be redistributed
    @param[out] nodeRowMap Node row map after redistribution
    @param[out] nodeColumnMap Node column map after redistribution
    @param[in] outflag Flag to switch screen output on/off
    @param[in] numPartitions Desired number of partitions
    @param[in] imbalanceTol Max. relative imbalance of subdomain size
    */
    void RedistributeGraphOfDiscretization(Teuchos::RCP<DRT::Discretization> discretization,
        Teuchos::RCP<const Epetra_Map> elementRowMap, Teuchos::RCP<Epetra_Map>& nodeRowMap,
        Teuchos::RCP<Epetra_Map>& nodeColumnMap, Teuchos::RCP<const Epetra_Comm> comm,
        const bool outflag, const int numPartitions, const double imbalanceTol = 1.1);

    /*!
    \brief Redistribute the underlying graph of a discretization to a new and more balanced
    partitioning while considering weight to model costs

    \pre Input discretization has to be FillComplete().

    @param[in/out] dis Discretization to be partitioned
    @param[out] rownodes Node row map after repartitioning
    @param[out] colnodes Node column map after repartitioning
    @param[in] outflag Flag to switch screen output on/off
    */
    void RedistributeGraphOfDiscretizationUsingWeights(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes, const bool outflag);

    /*!
    \brief Redistribute and FillComplete() a discretization to a new partitioning using weights to
    model costs

    First, call DRT::UTILS::RedistributeDiscretizationUsingWeights() to compute the new partitioning
    and to perform the redistribution. Then, call FillComplete() on the discretization.

    @param[in/out] dis Discretization to be repartitioned
    @param[in] assigndegreesoffreedom Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] initelements Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] doboundaryconditions Flag forwarded to DRT::Discretization::FillComplete()

    \sa DRT::UTILS::RedistributeGraphOfDiscretizationUsingWeights()
    \sa DRT::Discretization::FillComplete()
    */
    void RedistributeAndFillCompleteDiscretizationUsingWeights(
        Teuchos::RCP<DRT::Discretization> dis, const bool assigndegreesoffreedom,
        const bool initelements, const bool doboundaryconditions);

    /*!
    \brief Build graph of a given  discretization

    @param[in] dis Discretization whose graph will be build
    @param[in] roweles Element row map of this discretization
    @param[out] rownodes Node row map
    @param[in] outflag Flag to switch screen output on/off

    @return Uncompleted graph of input discretization
    */
    Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<const Epetra_Map> roweles, Teuchos::RCP<Epetra_Map>& rownodes,
        Teuchos::RCP<const Epetra_Comm> comm, const bool outflag);

    /*!
    \brief Rebalance graph

    Use Isorropia package to access Zoltan. By default, Isorropia will use Zoltan hypergraph
    partitioning, treating the graph columns as hyper-edges and the graph rows as vertices.

    \note The rebalanced graph will NOT be FillComplete().

    @param[in] initialGraph Initial graph to be rebalanced
    @param[in] rebalanceParams Parameter list with rebalancing options

    @return Rebalanced graph
    */
    Teuchos::RCP<Epetra_CrsGraph> RebalanceGraph(Teuchos::RCP<const Epetra_CrsGraph> initialGraph,
        const Teuchos::ParameterList& rebalanceParams);

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H
