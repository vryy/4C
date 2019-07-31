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

#include "Teuchos_RCP.hpp"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace UTILS
  {
    /*!
    \brief Redistribute a discretization using Zoltan

    Actually, only the \p rownodes and \p colnodes maps are created.

    \param[in] dis Discretization to be redistributed
    \param[in] roweles Map of row elements
    \param[in] rownodes Map of row nodes
    \param[in] colnodes Map of column nodes
    \param[in] comm Communicator object
    \param[in] outflag
    \param[in] parts Number of subdomains to be created
    \param[in] imbalance Max. relative imbalance of subdomain size

    \note For historic reasons, this routine is still named after the partitioning tool ParMetis.
    \todo ToDo Rename this routine to reflect that we don't use ParMetis anymore.
    */
    void PartUsingParMetis(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map> roweles,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes,
        Teuchos::RCP<Epetra_Comm> comm, const bool outflag, const int parts = -1,
        const double imbalance = -1.0);

    /*!
    \brief Redistribute a discretization to a new partitioning using weights to model costs

    Use Zoltan to compute a new partitioning that takes a cost model based on the graph's edge and
    vertex weights into account. Then, redistribute the discretization to the new partitioning.

    @param[in/out] dis Discretization to be partitioned
    @param[out] rownodes Node row map after repartitioning
    @param[out] colnodes Node column map after repartitioning
    @param[in] outflag Flag to switch screen output on/off

    \note We do not perform the partitioning ourselves, but call Zoltan through Trilinos's Isorropia
    package.
    */
    void RedistributeDiscretizationUsingWeights(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes, const bool outflag);

    Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map> roweles, Teuchos::RCP<Epetra_Map>& rownodes,
        Teuchos::RCP<Epetra_Comm> comm, bool outflag);

    /*!
    \brief Redistribute and FillComplete() a discretization to a new partitioning using weights to
    model costs

    First, call DRT::UTILS::RedistributeDiscretizationUsingWeights() to compute the new partitioning
    and to perform the redistribution. Then, call FillComplete() on the discretization.

    @param[in/out] dis Discretization to be repartitioned
    @param[in] assigndegreesoffreedom Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] initelements Flag forwarded to DRT::Discretization::FillComplete()
    @param[in] doboundaryconditions Flag forwarded to DRT::Discretization::FillComplete()

    \sa DRT::UTILS::RedistributeDiscretizationUsingWeights()
    \sa DRT::Discretization::FillComplete()
    */
    void RedistributeAndFillCompleteDiscretizationUsingWeights(
        Teuchos::RCP<DRT::Discretization> dis, const bool assigndegreesoffreedom,
        const bool initelements, const bool doboundaryconditions);

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H
