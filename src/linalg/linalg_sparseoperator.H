/*!----------------------------------------------------------------------
\file linalg_sparseoperator.H

\brief Base class for implementation of sparse operations (including sparse
       matrices)

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef LINALG_SPARSEOPERATOR_H_
#define LINALG_SPARSEOPERATOR_H_

#include <vector>

#include <Epetra_Operator.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_Vector;

namespace LINALG
{
  // forward declarations
  class BlockSparseMatrixBase;
  class SparseMatrixBase;
  class SparseMatrix;

  /*! \enum LINALG::DataAccess
   *  \brief Handling of data access (Copy or View)
   *
   *  If set to LINALG::Copy, user data will be copied at construction.
   *  If set to LINALG::View, user data will be encapsulated and used throughout
   *  the life of the object.
   *
   *  \note A separate LINALG::DataAccess is necessary in order to resolve
   *  possible ambiguity conflicts with the Epetra_DataAccess.
   *
   *  Use LINALG::DataAccess for construction of any LINALG matrix object.
   *  Use plain 'Copy' or 'View' for construction of any Epetra matrix object.
   *
   *  \author mayr.mt \date 10/2015
   */
  enum DataAccess
  {
    Copy,  ///< deep copy
    View   ///< reference to original data
  };

  /// Linear operator interface enhanced for use in FE simulations
  /*!

    The point in FE simulations is that you have to assemble (element)
    contributions to the global matrix, apply Dirichlet conditions in some way
    and finally solve the completed system of equations.

    Here we have an interface that has different implementations. The obvious
    one is the SparseMatrix, a single Epetra_CrsMatrix in a box, another one
    is BlockSparseMatrix, a block matrix build from a list of SparseMatrix.

    \author u.kue
    \date 02/08
   */
  class SparseOperator : public Epetra_Operator
  {
   public:
    /// return the internal Epetra_Operator
    /*!
      By default the SparseOperator is its own Epetra_Operator. However
      subclasses might have a better connection to Epetra.

      \warning Only low level solver routines are interested in the internal
      Epetra_Operator.
     */
    virtual Teuchos::RCP<Epetra_Operator> EpetraOperator() { return Teuchos::rcp(this, false); }

    /// set matrix to zero
    virtual void Zero() = 0;

    /// throw away the matrix and its graph and start anew
    virtual void Reset() = 0;

    /// Assemble a Epetra_SerialDenseMatrix into a matrix with striding
    /*!

    This is an individual call.  Will only assemble locally and will never
    do any commmunication.  All values that cannot be assembled locally will
    be ignored.  Will use the communicator and rowmap from matrix to
    determine ownerships.  Local matrix Aele has to be square.

    If matrix is Filled(), it stays so and you can only assemble to places
    already masked. An attempt to assemble into a non-existing place is a
    grave mistake.

    If matrix is not Filled(), the matrix is enlarged as required.

    \note Assembling to a non-Filled() matrix is much more expensive than to
    a Filled() matrix. If the sparse mask does not change it pays to keep
    the matrix around and assemble into the Filled() matrix.

    The first parameter \p eid is purely for performance enhancements. Plain
    sparse matrices do not know about finite elements and do not use the
    element id at all. However, BlockSparseMatrix might be created with
    specialized, problem specific assembling strategies. And these strategies
    might gain considerable performance advantages from knowing the element
    id.

    \param eid (in) : element gid
    \param Aele (in) : dense matrix to be assembled
    \param lm (in) : vector with gids
    \param lmowner (in) : vector with owner procs of gids
    */
    virtual void Assemble(int eid, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lm,
        const std::vector<int>& lmowner)
    {
      Assemble(eid, lmstride, Aele, lm, lmowner, lm);
    }

#if 0
    /// Assemble a Epetra_SerialDenseMatrix into a matrix
    /*!

    This is an individual call.  Will only assemble locally and will never
    do any commmunication.  All values that cannot be assembled locally will
    be ignored.  Will use the communicator and rowmap from matrix to
    determine ownerships.  Local matrix Aele has to be square.

    If matrix is Filled(), it stays so and you can only assemble to places
    already masked. An attempt to assemble into a non-existing place is a
    grave mistake.

    If matrix is not Filled(), the matrix is enlarged as required.

    \note Assembling to a non-Filled() matrix is much more expensive than to
    a Filled() matrix. If the sparse mask does not change it pays to keep
    the matrix around and assemble into the Filled() matrix.

    The first parameter \p eid is purely for performance enhancements. Plain
    sparse matrices do not know about finite elements and do not use the
    element id at all. However, BlockSparseMatrix might be created with
    specialized, problem specific assembling strategies. And these strategies
    might gain considerable performance advantages from knowing the element
    id.

    \param eid (in) : element gid
    \param Aele (in) : dense matrix to be assembled
    \param lm (in) : vector with gids
    \param lmowner (in) : vector with owner procs of gids
    */
    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(eid,Aele,lm,lmowner,lm);
    }
#endif

    /// Assemble a Epetra_SerialDenseMatrix into a matrix with striding
    /*!

      This is an individual call.
      Will only assemble locally and will never do any communication.
      All values that can not be assembled locally will be ignored.
      Will use the communicator and rowmap from matrix A to determine ownerships.
      Local matrix Aele may be \b square or \b rectangular.

      If matrix is Filled(), it stays so and you can only assemble to places
      already masked. An attempt to assemble into a non-existing place is a
      grave mistake.

      If matrix is not Filled(), the matrix is enlarged as required.

      \note Assembling to a non-Filled() matrix is much more expensive than to
      a Filled() matrix. If the sparse mask does not change it pays to keep
      the matrix around and assemble into the Filled() matrix.

      \note The user must provide an \b additional input vector 'lmcol'
      containing the column gids for assembly seperately!

      The first parameter \p eid is purely for performance enhancements. Plain
      sparse matrices do not know about finite elements and do not use the
      element id at all. However, BlockSparseMatrix might be created with
      specialized, problem specific assembling strategies. And these
      strategies might gain considerable performance advantages from knowing
      the element id.

      \param eid (in) : element gid
      \param Aele (in)       : dense matrix to be assembled
      \param lmrow (in)      : vector with row gids
      \param lmrowowner (in) : vector with owner procs of row gids
      \param lmcol (in)      : vector with column gids
    */
    virtual void Assemble(int eid, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol) = 0;

#if 0
    /// Assemble a Epetra_SerialDenseMatrix into a matrix
    /*!

      This is an individual call.
      Will only assemble locally and will never do any commmunication.
      All values that can not be assembled locally will be ignored.
      Will use the communicator and rowmap from matrix A to determine ownerships.
      Local matrix Aele may be \b square or \b rectangular.

      If matrix is Filled(), it stays so and you can only assemble to places
      already masked. An attempt to assemble into a non-existing place is a
      grave mistake.

      If matrix is not Filled(), the matrix is enlarged as required.

      \note Assembling to a non-Filled() matrix is much more expensive than to
      a Filled() matrix. If the sparse mask does not change it pays to keep
      the matrix around and assemble into the Filled() matrix.

      \note The user must provide an \b additional input vector 'lmcol'
      containing the column gids for assembly seperately!

      The first parameter \p eid is purely for performance enhancements. Plain
      sparse matrices do not know about finite elements and do not use the
      element id at all. However, BlockSparseMatrix might be created with
      specialized, problem specific assembling strategies. And these
      strategies might gain considerable performance advantages from knowing
      the element id.

      \param eid (in) : element gid
      \param Aele (in)       : dense matrix to be assembled
      \param lmrow (in)      : vector with row gids
      \param lmrowowner (in) : vector with owner procs of row gids
      \param lmcol (in)      : vector with column gids
    */
    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol) = 0;
#endif

    /// single value assemble using gids
    virtual void Assemble(double val, int rgid, int cgid) = 0;

    /// If Complete() has been called, this query returns true, otherwise it returns false.
    virtual bool Filled() const = 0;

    /// Call FillComplete on a matrix
    virtual void Complete() = 0;

    /// Call FillComplete on a matrix (for rectangular and square matrices)
    virtual void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap) = 0;

    /// Undo a previous Complete() call
    virtual void UnComplete() = 0;

    /// Apply dirichlet boundary condition to a matrix
    virtual void ApplyDirichlet(
        const Teuchos::RCP<const Epetra_Vector> dbctoggle, bool diagonalblock = true) = 0;

    /// Apply dirichlet boundary condition to a matrix
    ///
    ///  This method blanks the rows associated with Dirichlet DOFs
    ///  and puts a 1.0 at the diagonal entry if diagonalblock==true.
    ///  Only the rows are blanked, the columns are not touched.
    ///  We are left with a non-symmetric matrix, if the original
    ///  matrix was symmetric. However, the blanking of columns is computationally
    ///  quite expensive, because the matrix is stored in a sparse and distributed
    ///  manner.
    virtual void ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock = true) = 0;

    /** \brief Return TRUE if all Dirichlet boundary conditions have been applied
     *  to this matrix
     *
     *  \param (in) dbcmap: DBC map holding all dbc dofs
     *  \param (in) diagonalblock: Is this matrix a diagonalblock of a blocksparsematrix?
     *                             If it is only one block/matrix, this boolean should be TRUE.
     *  \param (in) trafo: pointer to an optional trafo matrix (see LocSys).
     *
     *  \author hiermeier \date 01/18 */
    virtual bool IsDbcApplied(const Epetra_Map& dbcmap, bool diagonalblock = true,
        const LINALG::SparseMatrix* trafo = NULL) const = 0;

    /// Returns the Epetra_Map object associated with the (full) domain of this operator.
    virtual const Epetra_Map& DomainMap() const = 0;

    /// Add one operator to another
    virtual void Add(const LINALG::SparseOperator& A, const bool transposeA, const double scalarA,
        const double scalarB) = 0;

    /// Add one SparseMatrixBase to another
    virtual void AddOther(LINALG::SparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB) const = 0;

    /// Add one BlockSparseMatrix to another
    virtual void AddOther(LINALG::BlockSparseMatrixBase& A, const bool transposeA,
        const double scalarA, const double scalarB) const = 0;

    /// Multiply all values by a constant value (in place: A <- ScalarConstant * A).
    virtual int Scale(double ScalarConstant) = 0;

    /// Matrix-vector product
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
  };


}  // namespace LINALG

#endif
/*LINALG_SPARSEOPERATOR_H_*/
