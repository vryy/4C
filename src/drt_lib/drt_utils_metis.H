/*!----------------------------------------------------------------------
\file drt_utils_metis.H
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
#ifndef UTILS_METIS_H
#define UTILS_METIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_CrsGraph.h"
#include "Epetra_Vector.h"

namespace DRT
{
/// collection of useful and important but less central stuff
namespace UTILS
{

/*!
\brief Partition a graph with metis

Take an input graph with any row map and a vector of point weights with the same row
map and redistribute the graph using the partitioning library metis.
The output graph has a row map according to the new distribution.

\param graph : (in) any graph
\param weights: (int) a vector of weight factors for each row of graph.
                Map of weights must be pointwise equal to graph.RowMap().

\return graph with new distribution
*/
Teuchos::RCP<Epetra_CrsGraph> PartGraphUsingMetis(const Epetra_CrsGraph& graph,
                                                 const Epetra_Vector& weights);


} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H




