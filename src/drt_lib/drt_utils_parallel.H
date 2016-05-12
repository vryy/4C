/*!----------------------------------------------------------------------
\file drt_utils_parallel.H

\brief A collection of helper methods for namespace DRT

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_PARALLEL_H
#define UTILS_PARALLEL_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"


namespace DRT
{
/// forward declarations
class Discretization;

/// collection of useful and important but less central stuff
namespace UTILS
{
  //! Given a nodal row and coloumn map, this routines creates the corresponding element maps and call FillComplete()
  void RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    );

  /// print information about parallel mesh distribution
  void PrintParallelDistribution(const DRT::Discretization& dis);

  /*!
  \brief Get a column vector made of a row vector.

  Using this method, a reference to a vector is returned.
  If the vector is supplied in dis.DofColMap() the vector itself will be returned.
  If the vector is NOT supplied in dis.DofColMap(), but in dis.DofRowMap(),
  a new vector with column map is allocated and the supplied vector is exported to it and returned.

  Note that the very same functionality is used in Discretization::SetState()!

  \param name (in): discretization
  \param state (in): vector of some data
  */
  Teuchos::RCP<const Epetra_Vector> GetColVersionOfRowVector(
      const Teuchos::RCP<const DRT::Discretization> dis,
      const Teuchos::RCP<const Epetra_Vector> state,
      const int nds = 0);


  ///recompute nodecolmap of standard discretization to include all nodes as of subdicretization
  Teuchos::RCP<Epetra_Map> ComputeNodeColMap
  (
      const Teuchos::RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
      const Teuchos::RCP<DRT::Discretization> subdis ///< subdiscretization prescribing ghosting
  );

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_PARALLEL_H
