/*!----------------------------------------------------------------------
\file drt_utils_parmetis.H
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
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

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

  void PartUsingParMetis(Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<Epetra_Map> roweles,
      Teuchos::RCP<Epetra_Map>& rownodes,
      Teuchos::RCP<Epetra_Map>& colnodes,
      Teuchos::RCP<Epetra_Comm> comm,
      bool outflag,
      int parts = -1);

  void PartUsingZoltanWithWeights(Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<Epetra_Map>& rownodes,
      Teuchos::RCP<Epetra_Map>& colnodes,
      bool outflag);

  Teuchos::RCP<const Epetra_CrsGraph> BuildGraph(Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<Epetra_Map> roweles,
      Teuchos::RCP<Epetra_Map>& rownodes,
      Teuchos::RCP<Epetra_Comm> comm,
      bool outflag);

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H




