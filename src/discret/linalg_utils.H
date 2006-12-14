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
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
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
  \param init (in): initializa vector to zero upon construction
  */
  RefCountPtr<Epetra_Vector> CreateVector(const Epetra_Map& rowmap, const bool init);

  /*!
  \brief Export a vector to a different map
  
  Values of source are copied to target where maps don't have to match.
  Prerequisite: Either the map of source OR the map of target has to be unique
                (will be tested)
  \warning When source is overlapping (and therefore target is unique), values
           in the overlapping region are inserted into the target on a first come 
           first serve basis, meaning they should be equal in the source to 
           be deterministic
  \param source (in) : source vector values are taken from
  \param target (out): target vector values will be inserted in
  */
  void Export(const Epetra_Vector& source, Epetra_Vector& target);

  /*!
  \brief Invert a symmetric dim*dim square matrix
  
  \param A (in/out): Matrix to be inverted
  \param dim (in) :  Dimension of matrix
  */
  void SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

} // namespace LINALG


#endif  // #ifndef LINALG_UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
