/*!----------------------------------------------------------------------
\file drt_utils_parmetis.H

\brief A collection of helper methods related to partitioning and parallel distribution

\level 0

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
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

    void PartUsingZoltanWithWeights(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes, bool outflag);

    Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map> roweles, Teuchos::RCP<Epetra_Map>& rownodes,
        Teuchos::RCP<Epetra_Comm> comm, bool outflag);

    void WeightedRepartitioning(Teuchos::RCP<DRT::Discretization> dis, bool assigndegreesoffreedom,
        bool initelements, bool doboundaryconditions);

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H
