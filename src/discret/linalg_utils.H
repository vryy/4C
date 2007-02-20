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
#include "Epetra_SerialDenseVector.h"
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

  /*!
  \brief Assemble a Epetra_SerialDenseMatrix into a Epetra_CrsMatrix
  
  This is an individual call. Matrix A is not completed upon return.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator from matrix A to determine ownerships.
  Local matrix Aele has to be square.
  
  \param A (out)   : Sparse matrix to be assembled on
  \param Aele (in) : dense matrix to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owners of gids
  */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele, 
                const vector<int>& lm, const vector<int>& lmowner);

  /*!
  \brief Assemble a Epetra_SerialDenseVector into a Epetra_Vector
  
  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator from vector V to determine ownerships.
  
  \param V (out)   : Sparse matrix to be assembled on
  \param Vele (in) : dense matrix to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owners of gids
  */
  void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele, 
                const vector<int>& lm, const vector<int>& lmowner);


  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix
  */
  void Complete(Epetra_CrsMatrix& A);

  /*!
  \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA
  
  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Range, Row and Domain maps of A(^T) and B have to match.
  
  Note that this is a true parallel add, even in the transposed case!
  
  \param A          (in)     : Matrix to add to B (must have Filled()==true)
  \param transposeA (in)     : flag indicating whether transposed of A should be used
  \param scalarA    (in)     : scaling factor for A
  \param B          (in/out) : Matrix to be added to (must have Filled()==false)
  \param scalarB    (in)     : scaling factor for B
  */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
           Epetra_CrsMatrix& B, const double scalarB);

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations
    
  
  \param A         (in/out) : Matrix to of Ax=b
  \param x         (in/out) : initial guess vector of Ax=b
  \param b         (in/out) : rhs vector of Ax=b
  \param dbcval    (in)     : vector holding prescribed dirichlet values
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                              RefCountPtr<Epetra_Vector>&      x,
                              RefCountPtr<Epetra_Vector>&      b,
                              const RefCountPtr<Epetra_Vector> dbcval,
                              const RefCountPtr<Epetra_Vector> dbctoggle);

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations
    
  
  \param A         (in/out) : Matrix to of Ax=b
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                              const RefCountPtr<Epetra_Vector> dbctoggle);

} // namespace LINALG


#endif  // #ifndef LINALG_UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
