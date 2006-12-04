/*!----------------------------------------------------------------------
\file linalg_utils.H
\brief A collection of helper methods for namespace LINALG

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef LINALG_UTILS_H
#define LINALG_UTILS_H


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{
  /*!
  \brief Create a new Epetra_CrsMatrix and return RefcountPtr to it
  
  \param rowmap (in): row map of matrix
  \param npr (in): estimated number of entries per row.
                   (need not be exact better should be too big rather then too small)
  */
  RefCountPtr<Epetra_CrsMatrix> CreateMatrix(const Epetra_Map& rowmap, const int npr);
  
  /*!
  \brief Create a new Epetra_Vector and return RefcountPtr to it
  
  \param rowmap (in): row map of vector
  \param npr (in): initializa vector to zero upon construction
  */
  RefCountPtr<Epetra_Vector> CreateVector(const Epetra_Map& rowmap, const bool init);

} // namespace LINALG


#endif  // #ifndef LINALG_UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
