/*!----------------------------------------------------------------------
\file linalg_sparseoperator.H

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
#ifdef CCADISCRET

#ifndef LINALG_SPARSEOPERATOR_H_
#define LINALG_SPARSEOPERATOR_H_

#include <vector>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

namespace LINALG
{
  // forward declarations
  class BlockSparseMatrixBase;
  class SparseMatrix;


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
    virtual Teuchos::RCP<Epetra_Operator> EpetraOperator() { return Teuchos::rcp(this,false); }

    /// set matrix to zero
    virtual void Zero() = 0;

    /// throw away the matrix and its graph and start anew
    virtual void Reset() = 0;

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
    virtual void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle, bool diagonalblock=true) = 0;

    /// Apply dirichlet boundary condition to a matrix
    ///
    ///  This method blanks the rows associated with Dirichlet DOFs
    ///  and puts a 1.0 at the diagonal entry if diagonlblock==true.
    ///  Only the rows are blanked, the columns are not touched.
    ///  We are left with a non-symmetric matrix, if the original
    ///  matrix was symmtric. However, the blanking of columns is computionally
    ///  quite expensive, because the matrix is stored in a sparse and distributed
    ///  manner.
    virtual void ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock=true) = 0;

    /// Returns the Epetra_Map object associated with the (full) domain of this operator.
    virtual const Epetra_Map& DomainMap() const = 0;

    /// Add one operator to another
    virtual void Add(const LINALG::SparseOperator& A, const bool transposeA, const double scalarA,  const double scalarB) = 0;

    /// Add one SparseMatrix to another
    virtual const void AddOther(LINALG::SparseMatrix& A, const bool transposeA, const double scalarA,  const double scalarB) const = 0;

    /// Add one BlockSparseMatrix to another
    virtual const void AddOther(LINALG::BlockSparseMatrixBase& A, const bool transposeA, const double scalarA,  const double scalarB) const = 0;

    /// Matrix-vector product
    virtual int Multiply(bool TransA, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const = 0;

};


}

#endif

#endif /*LINALG_SPARSEOPERATOR_H_*/
