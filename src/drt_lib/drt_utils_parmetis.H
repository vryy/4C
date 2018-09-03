/*!----------------------------------------------------------------------
\file drt_utils_parmetis.H

\brief A collection of helper methods for namespace DRT

\level 0

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef UTILS_PARMETIS_H
#define UTILS_PARMETIS_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"


namespace DRT
{
  /// forward declarations
  class Discretization;

  /// collection of useful and important but less central stuff
  namespace UTILS
  {
    /// redistribute a discretization using ParMetis
    /*!
      Actually, only the \p rownodes and \p colnodes maps are created.
     */

    void PartUsingParMetis(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map> roweles,
        Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<Epetra_Map>& colnodes,
        Teuchos::RCP<Epetra_Comm> comm, bool outflag, int parts = -1);

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
