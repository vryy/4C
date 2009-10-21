/*!----------------------------------------------------------------------
\file fluid_rotsym_periodicbc_utils.H

\brief Methods to apply rotationally symmetric periodic boundary conditions
       for fluid problems on element level

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef FLUID_ROTSYM_PERIODICBC_UTILS_H
#define FLUID_ROTSYM_PERIODICBC_UTILS_H


#include "Epetra_Vector.h"
#include "../drt_lib/drt_node.H"

using namespace std;
using namespace Teuchos;

/*!
\brief  This class manages local transformations(rotation) of velocity fields

        o used for rotationally symmetric boundary conditions
        o used for rotationally symmetric boundary conditions

\author gjb

 */

namespace FLD
{

//! return specific component of rotated vector result
double GetComponentOfRotatedVectorField
(const int idf, ///< the component id 0,1 (and 2)
 const RCP<const Epetra_Vector> proc0data, ///< the data vector
 const int lid, ///< local dof id (w.r.t proc0data) of component idf
 const double rotangle ///< angle of rotation (RAD)
  );

/// Is given node a slave node of rotationally symmetric periodic boundary conditions?
/*!
  \param node (in) the node
  \param rotangle (out) the angle of rotation (RAD!)
*/
bool IsSlaveNodeOfRotSymPBC(const DRT::Node* node, double& rotangle);

} // namespace FLD

#endif  // #ifndef FLUID_ROTSYM_PERIODICBC_UTILS_H
#endif  // #ifdef CCADISCRET
