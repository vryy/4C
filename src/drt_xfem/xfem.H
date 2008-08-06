/*!
\file xfem.H

\brief definitions used in xfem context

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFEM_H
#define XFEM_H

#include "../drt_geometry/vector_definitions.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"


typedef Epetra_SerialDenseVector EpVec; ///< a short-cut for a one-dimensional Epetra vector
typedef Epetra_SerialDenseMatrix EpMat; ///< a short-cut for a two-dimensional Epetra matrix

//! namespace for all XFEM related stuff
namespace XFEM
{
  typedef BlitzVec3 PosX;             ///< position in physical coordinates (x,y,z)
  typedef BlitzVec3 PosXiDomain;      ///< position in element coordinates of domain element (r,s,t)
  typedef BlitzVec2 PosXiBoundary;    ///< position in element coordinates of boundary element (r,s)
  typedef BlitzVec3 PosEtaDomain;     ///< position in local coordinates of a domain integrationcell (r,s,t)
  typedef BlitzVec2 PosEtaBoundary;   ///< position in local coordinates of a boundary integrationcell (r,s)
}

#endif  // #ifndef XFEM_H
#endif  // #ifdef CCADISCRET
