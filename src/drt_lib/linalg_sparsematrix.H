/*!----------------------------------------------------------------------
\file linalg_sparsematrix.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LINALG_SPARSEMATRIX_H_
#define LINALG_SPARSEMATRIX_H_


#include <Epetra_CrsMatrix.h>
#include <Epetra_Comm.h>


#include "linalg_sparseoperator.H"
#include "linalg_mapextractor.H"
#include "drt_dserror.H"



namespace LINALG
{

  // forward delcarations:
  class BlockSparseMatrixBase;
  template <class Strategy> class BlockSparseMatrix;
  
  /// A single sparse matrix enhanced with features for FE simulations
  /*!

    A single sparse matrix. Internally we have an Epetra_CrsMatrix. So we have
    all the glory of a fully parallel and fast sparse matrix. The added value
    is twofold. For one thing there are the FE specific operations. Assemble()
    adds an (element) matrix to the (global) sparse matrix and
    ApplyDirichlet() modifies the matrix to contain just ones on Dirichlet
    rows (the columns in other rows are not touched, so the matrix becomes
    unsymmetric).

    The second gain are the different states this matrix can be in. You can
    set explicitdirichlet==true in order to modify the matrix graph in each
    ApplyDirichlet() call to contain just the diagonal entries on those rows
    -- this essentially copies the matrix. (ML gains a lot from completely
    Dirichlet-constrained rows.) With explicitdirichlet==false the matrix
    graph is not touched, instead the Dirichlet rows are filled with zeros.

    With savegraph==true you specify that you want to keep the original matrix
    graph before you apply Dirichlet conditions. This way you can call Zero()
    and get an already Filled() matrix. You cannot alter its graph afterwards,
    but Assemble() is much faster if your matrix is already Filled(). Of course
    you can always Reset() you matrix, that is throw away the matrix graph and
    start anew with an empty matrix.

    \note A large part of the SparseMatrix interface consists of methods from
    the internal Epetra_CrsMatrix. If there are methods in Epetra_CrsMatrix
    and not in SparseMatrix that you would like to call (for legitimate
    reasons!) please add them to the SparseMatrix.

    \author u.kue
    \date 02/08
   */
class SparseMatrix : public SparseOperator
{
  public:

    /// construction of sparse matrix
    SparseMatrix(const Epetra_Map& rowmap, const int npr, bool explicitdirichlet=true, bool savegraph=false);

    /// construction of sparse matrix
    /*!
      Implicit construction from a Epetra_CrsMatrix

      Makes a copy of the Epetra_CrsMatrix.

      \pre matrix.Filled()==true
     */
    explicit SparseMatrix(const Epetra_CrsMatrix& matrix, bool explicitdirichlet=true, bool savegraph=false);

    /// construction of sparse matrix
    /*!
      Implicit construction from a Epetra_CrsMatrix

      \warning Does not make a copy of the Epetra_CrsMatrix. The RCP is copied
      internally. So this constructor creates a new view on the
      Epetra_CrsMatrix.
     */
    explicit SparseMatrix(Teuchos::RCP<Epetra_CrsMatrix> matrix, bool explicitdirichlet=true, bool savegraph=false);

    /// construction of a diagonal matrix from a vector
    /*!
      Creates diagonal matrix with range and domain map euqal to vector map.
      Sets diagonal values from vector and does NOT call Complete() on matrix

      Does not make a copy of the Epetra_CrsMatrix.
     */
    explicit SparseMatrix(const Epetra_Vector& diag, bool explicitdirichlet=true, bool savegraph=false);

    /// Copy constructor. Deep copy or view on matrix.
    /*!
      \warning A view assignment will have your matrix use the same internal
      data as the original matrix. Changes on one side will affect the
      other. However, some methods like Zero() or Reset() can, depending on
      the SparseMatrix flags, cut the connection. Do not rely on the view if
      you change one of these matrices! The view assignment is meant to
      provide a slim copy operation that transferes ownership from one matrix
      to the other before the original matrix is destroyed. Do a deep copy if
      both matrices are meant to live on.

      \param mat matrix to assign from
      \param access how to treat this assignment: Copy or View
     */
    SparseMatrix(const SparseMatrix& mat, Epetra_DataAccess access=Copy);

    /// Destructor
    virtual ~SparseMatrix();

    /// Assignment operator. Makes a deep copy.
    SparseMatrix& operator=(const SparseMatrix& mat);

    /// Assignment method. Deep copy or view on matrix.
    /*!
      Explicit method for the assignment operator. You can make an explicit
      copy of the internal Epetra_CrsMatrix or have a second view on it.

      \warning A view assignment will have your matrix use the same internal
      data as the original matrix. Changes on one side will affect the
      other. However, some methods like Zero() or Reset() can, depending on
      the SparseMatrix flags, cut the connection. Do not rely on the view if
      you change one of these matrices! The view assignment is meant to
      provide a slim copy operation that transferes ownership from one matrix
      to the other before the original matrix is destroyed. Do a deep copy if
      both matrices are meant to live on.

      \param access how to treat this assignment: Copy or View
      \param mat matrix to assign from
     */
    void Assign(Epetra_DataAccess access, const SparseMatrix& mat);

    /// return the internal Epetra_Operator
    /*!
      The internal Epetra_Operator here is the internal Epetra_CrsMatrix. This
      way the solver can down-cast to Epetra_CrsMatrix and access the matrix
      rows directly.

      \note This method is here for performance reasons.
     */
    virtual Teuchos::RCP<Epetra_Operator> EpetraOperator() { return sysmat_; }

    /// return the internal Epetra_CrsMatrix (you should not need this!)
    Teuchos::RCP<Epetra_CrsMatrix> EpetraMatrix() { return sysmat_; }

    /// return the internal Epetra_CrsMatrix (you should not need this!)
    Teuchos::RCP<Epetra_CrsMatrix> EpetraMatrix() const { return sysmat_; }

    /** \name FE methods */
    //@{

    /// set all matrix entries to zero
    void Zero();

    /// throw away the matrix and its graph and start anew
    void Reset();

    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(eid,Aele,lm,lmowner,lm);
    }

    void Assemble(int eid,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    /// single value assemble used by BlockSparseMatrix
    void Assemble(double val, int rgid, int cgid);

    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void UnComplete();

    void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle, bool diagonalblock=true);

    /// create matrix that contains all Dirichlet lines from my
    Teuchos::RCP<SparseMatrix> ExtractDirichletLines(const Teuchos::RCP<Epetra_Vector> dbctoggle);

    //@}

    /** \name Matrix Properties Query Methods */
    //@{

    /// If Complete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const { return sysmat_->Filled(); }

    /// Whether Dirichlet conditions should result in a trimmed graph row
    /*!
      ML requires rows of length 1 to recognize Dirichlet lines. However it is
      an expensive operation to apply Dirichlet conditions in this case.
     */
    bool ExplicitDirichlet() const { return explicitdirichlet_; }

    /// Whether the matrix graph should be saved when the matrix is zeroed
    /*!
      Saving the graph will result in constructing new matrices in Filled()
      state. This speeds up assembling but limits assembling to the current
      graph.
     */
    bool SaveGraph() const { return savegraph_; }

    //@}

    /** \name Attribute set methods */
    //@{

    /// If set true, transpose of this operator will be applied.
    virtual int SetUseTranspose(bool UseTranspose);

    //@}

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Returns the infinity norm of the global matrix.
    double NormInf() const;

    /// Returns the one norm of the global matrix.
    double NormOne() const;

    /// Returns the frobenius norm of the global matrix.
    double NormFrobenius() const;

    //@}

    /** \name Attribute access functions */
    //@{

    /// Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumEntries() const;

    /// Returns the Epetra_Map object associated with the rows of this matrix.
    const Epetra_Map& RowMap() const { return sysmat_->RowMap(); }

    /// Returns the Epetra_Map object that describes the set of column-indices that appear in each processor's locally owned matrix rows.
    const Epetra_Map& ColMap() const { return sysmat_->ColMap(); }

    /// Returns the Epetra_Map object associated with the domain of this matrix operator.
    const Epetra_Map& DomainMap() const { return sysmat_->DomainMap(); }

    /// Returns the Epetra_Map object associated with the range of this matrix operator.
    const Epetra_Map& RangeMap() const { return sysmat_->RangeMap(); }

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    /// Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    /// Returns true if the this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    /// Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm& Comm() const;

    /// Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map& OperatorDomainMap() const;

    /// Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map& OperatorRangeMap() const;

    //@}

    /** \name Computational methods */
    //@{

    /// Returns the result of a matrix multiplied by a Epetra_Vector x in y.
    int Multiply(bool TransA, const Epetra_Vector &x, Epetra_Vector &y) const;

    /// Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_MultiVector X in Y.
    int Multiply(bool TransA, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Scales the Epetra_CrsMatrix on the left with a Epetra_Vector x.
    int LeftScale(const Epetra_Vector &x);

    /// Scales the Epetra_CrsMatrix on the right with a Epetra_Vector x.
    int RightScale(const Epetra_Vector &x);

    //@}

    /** \name Insertion/Replace/SumInto methods */
    //@{

    /// Initialize all values in the matrix with constant value.
    int PutScalar(double ScalarConstant);

    /// Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
    int Scale(double ScalarConstant);

    /// Replaces diagonal values of the matrix with those in the user-provided vector.
    int ReplaceDiagonalValues(const Epetra_Vector &Diagonal);

    //@}

    /** \name Extraction methods */
    //@{

    /// Returns a copy of the main diagonal in a user-provided vector.
    int ExtractDiagonalCopy(Epetra_Vector &Diagonal) const;

    //@}

    /** \name Utility functions */
    //@{

    /// Compute transposed matrix explicitly
    /*!
      \note This is an expensive operation!
     */
    Teuchos::RCP<SparseMatrix> Transpose();

    /// Add a (transposed) Epetra_CrsMatrix to another: (*this) = (*this)*scalarB + A(^T)*scalarA
    /*!

    Add one matrix to another. the matrix (*this) to be added to must not be
    completed. Sparsity patterns of A and (*this) need not match and A and (*this) can be
    nonsymmetric in value and pattern.  Row map of A has to be a
    processor-local subset of the row map of (*this).

    \note This is a true parallel add, even in the transposed case!

    \param A          (in)     : Matrix to add to B (must have Filled()==true)
    \param transposeA (in)     : flag indicating whether transposed of A should be used
    \param scalarA    (in)     : scaling factor for A
    \param scalarB    (in)     : scaling factor for B
    */
    void Add(const SparseMatrix& A,
             const bool transposeA,
             const double scalarA,
             const double scalarB);

    /// Add a (transposed) Epetra_CrsMatrix to another: (*this) = (*this)*scalarB + A(^T)*scalarA
    void Add(const Epetra_CrsMatrix& A,
             const bool transposeA,
             const double scalarA,
             const double scalarB);

    /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
    /*!
      Multiply one matrix with another. Both matrices must be completed. Sparsity
      Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

      \note This is a true parallel multiplication, even in the transposed case.

      \note Does not call complete on C upon exit.

      \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
      \param transA         (in)     : flag indicating whether transposed of A should be used
      \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
      \param transB         (in)     : flag indicating whether transposed of B should be used
      \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
      \return Matrix product A(^T)*B(^T)
    */
    friend Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                               const SparseMatrix& B, bool transB,
                                               bool completeoutput = true);

    /*!
      \brief Split matrix in either 2x2 or 3x3 blocks (3x3 not yet impl.)

      Split given matrix 'this' into 2x2 block matrix and
      return result as templated BlockSparseMatrix.
      The MultiMapExtractor's provided have to be 2 and 2 maps,
      otherwise this method will throw an error.

      \note This is an expensive operation!

      \note This method will NOT call Complete() on the output
            BlockSparseMatrix.
     */
    template <class Strategy> Teuchos::RCP<LINALG::BlockSparseMatrix<Strategy> >
      Split(const MultiMapExtractor& domainmaps,const MultiMapExtractor& rangemaps);

    /// binary dump of matrix for debugging
    void Dump(std::string filename);

    //@}

  private:

    /// Split matrix in 2x2 blocks, where main diagonal blocks have to be square
    /*!
       Used by public Split, does not call Complete() on output matrix.
     */
    void Split2x2(BlockSparseMatrixBase& Abase);

    /// internal epetra matrix
    Teuchos::RCP<Epetra_CrsMatrix> sysmat_;

    /// whether to modify the matrix graph on apply Dirichlet
    bool explicitdirichlet_;

    /// whether to save the graph and assemble to a filled matrix next time
    bool savegraph_;

    /// (estimated) maximum number of entries per row
    int maxnumentries_;

    /// saved graph (if any)
    Teuchos::RCP<Epetra_CrsGraph> graph_;
  };

  /// output of SparseMatrix
  ostream& operator << (ostream& os, const LINALG::SparseMatrix& mat);

  
  Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                      const SparseMatrix& B, bool transB,
                                      bool completeoutput);
  
  
  /// merge 2x2 system
  Teuchos::RCP<SparseMatrix> Merge(const LINALG::SparseMatrix& Aii,
      const LINALG::SparseMatrix& Aig,
      const LINALG::SparseMatrix& Agi,
      const LINALG::SparseMatrix& Agg);


}  // end of LINALG namespace



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*!
  \brief Split matrix in either 2x2 or 3x3 blocks (3x3 not yet impl.)

  Split given matrix 'this' into 2x2 block matrix and
  return result as templated BlockSparseMatrix.
  The MultiMapExtractor's provided have to be 2 and 2 maps,
  otherwise this method will throw an error.

  \note This is an expensive operation!

  \note This method will NOT call Complete() on the output
        BlockSparseMatrix.
 */
template <class Strategy> Teuchos::RCP<LINALG::BlockSparseMatrix<Strategy> >
LINALG::SparseMatrix::Split(
    const MultiMapExtractor& domainmaps,
    const MultiMapExtractor& rangemaps)
{
  if (domainmaps.NumMaps()==2 && rangemaps.NumMaps()==2)
  {
    const int npr = EpetraMatrix()->MaxNumEntries();
    Teuchos::RCP<BlockSparseMatrix<Strategy> > blockA =
      rcp(new LINALG::BlockSparseMatrix<Strategy>(domainmaps,rangemaps,npr,
          explicitdirichlet_,savegraph_));
    this->Split2x2(*blockA);
    return blockA;
  }
  else
    dserror("Currently can only split matrix in 2x2 blocks");

  return Teuchos::null;
}


#endif
#endif
