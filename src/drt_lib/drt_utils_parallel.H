/*!----------------------------------------------------------------------
\file drt_utils_parallel.H
\brief A collection of helper methods for namespace DRT

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
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
      const Teuchos::RCP<const Epetra_Vector> state);


  ///recompute nodecolmap of standard discretization to include all nodes as of subdicretization
  Teuchos::RCP<Epetra_Map> ComputeNodeColMap
  (
      const Teuchos::RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
      const Teuchos::RCP<DRT::Discretization> subdis ///< subdiscretization prescribing ghosting
  );

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_PARALLEL_H




