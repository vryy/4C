#ifdef CCADISCRET

#ifndef LINALG_SYSTEMMATRIX_H
#define LINALG_SYSTEMMATRIX_H

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include "linalg_mapextractor.H"
#include "drt_dserror.H"



namespace LINALG
{
  // forward delcarations:
  class BlockSparseMatrixBase;
  template <class Strategy> class BlockSparseMatrix;


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

    \param Aele (in) : dense matrix to be assembled
    \param lm (in) : vector with gids
    \param lmowner (in) : vector with owner procs of gids
    */
    virtual void Assemble(const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(Aele,lm,lmowner,lm);
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

      \param Aele (in)       : dense matrix to be assembled
      \param lmrow (in)      : vector with row gids
      \param lmrowowner (in) : vector with owner procs of row gids
      \param lmcol (in)      : vector with column gids
    */
    virtual void Assemble(const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol) = 0;

    /// Call FillComplete on a matrix (for square matrices only!)
    virtual void Complete() = 0;

    /// Call FillComplete on a matrix (for rectangular and square matrices)
    virtual void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap) = 0;

    /// Apply dirichlet boundary condition to a matrix
    virtual void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle) = 0;

  };


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
     */
    explicit SparseMatrix(const Epetra_CrsMatrix& matrix, bool explicitdirichlet=true, bool savegraph=false);

    /// construction of sparse matrix
    /*!
      Implicit construction from a Epetra_CrsMatrix

      Does not make a copy of the Epetra_CrsMatrix.
     */
    explicit SparseMatrix(Teuchos::RCP<Epetra_CrsMatrix> matrix, bool explicitdirichlet=true, bool savegraph=false);

    /// construction of a diagonal matrix from a vector
    /*!
      Creates diagonal matrix with range and domain map euqal to vector map.
      Sets diagonal values from vector and does NOT call Complete() on matrix

      Does not make a copy of the Epetra_CrsMatrix.
     */
    explicit SparseMatrix(const Epetra_Vector& diag, bool explicitdirichlet=true, bool savegraph=false);

    /// copy constructor
    SparseMatrix(const SparseMatrix& mat);

    /// Destructor
    virtual ~SparseMatrix();

    /// assignment operator
    SparseMatrix& operator=(const SparseMatrix& mat);

    /// return the internal Epetra_Operator
    /*!
      The internal Epetra_Operator here is the internal Epetra_CrsMatrix. This
      way the solver can down-cast to Epetra_CrsMatrix and access the matrix
      rows directly.

      \note This method is here for performance reasons.
     */
    virtual Teuchos::RCP<Epetra_Operator> EpetraOperator() { return sysmat_; }

    /// return the internal Epetra_CrsMatrix (you should not need this!)
    Teuchos::RCP<Epetra_CrsMatrix>& EpetraMatrix() { return sysmat_; }

    /** \name FE methods */
    //@{

    /// set all matrix entries to zero
    void Zero();

    /// throw away the matrix and its graph and start anew
    void Reset();

    virtual void Assemble(const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(Aele,lm,lmowner,lm);
    }

    void Assemble(const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    /// single value assemble used by BlockSparseMatrix
    void Assemble(double val, int rgid, int cgid);

    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle);

    //@}

    /** \name Matrix Properties Query Methods */
    //@{

    /// If Complete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const { return sysmat_->Filled(); }

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

    /// Multiply a (transposed) matrix with me (transposed): C = A(^T)*B(^T)
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
    friend Teuchos::RCP<SparseMatrix> LINALG::Multiply(const SparseMatrix& A, bool transA,
                                                       const SparseMatrix& B, bool transB,
                                                       bool completeoutput = true);

    /// Split matrix in either 2x2 or 3x3 blocks
    /*!
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

    /// saved graph (if any)
    Teuchos::RCP<Epetra_CrsGraph> graph_;
  };


  Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                      const SparseMatrix& B, bool transB,
                                      bool completeoutput);

  /// Internal base class of BlockSparseMatrix that contains the non-template stuff
  /*!

    This is where the bookkeeping of the BlockSparseMatrix happens. We use two
    MultiMapExtractor objects to store the FullRangeMap() and the
    FullDomainMap() along with their many partial RangeMap() and
    DomainMap(). Most of the required SparseOperator methods can simply be
    implemented in terms of the matrix blocks.

    \author u.kue
    \date 02/08
  */
  class BlockSparseMatrixBase : public SparseOperator
  {
  public:

    /// constructor
    /*!
      \param domainmaps domain maps for all blocks
      \param rangemaps range maps for all blocks
      \param npr estimated number of entries per row in each block
      \param explicitdirichlet whether to remove Dirichlet zeros from the
      matrix graphs in each block
      \param savegraph whether to save the matrix graphs of each block and
      recreate filled matrices the next time
     */
    BlockSparseMatrixBase(const MultiMapExtractor& domainmaps,
                          const MultiMapExtractor& rangemaps,
                          int npr,
                          bool explicitdirichlet=true,
                          bool savegraph=false);

    virtual ~BlockSparseMatrixBase() {}

    /** \name Block matrix access */
    //@{

    const SparseMatrix& Matrix(int r, int c) const { return blocks_[r*Cols()+c]; }

    SparseMatrix& Matrix(int r, int c) { return blocks_[r*Cols()+c]; }

    inline const SparseMatrix& operator()(int r, int c) const { return Matrix(r,c); }

    inline SparseMatrix& operator()(int r, int c) { return Matrix(r,c); }

    //@}

    /** \name FE methods */
    //@{

    void Zero();
    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle);

    //@}

    /** \name Matrix Properties Query Methods */
    //@{

    /// If Complete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const;

    //@}

    /** \name Block maps */
    //@{

    /// number of row blocks
    int Rows() const { return rangemaps_.NumMaps(); }

    /// number of column blocks
    int Cols() const { return domainmaps_.NumMaps(); }

    /// range map for given row block
    Epetra_Map& RangeMap(int r) const { return *rangemaps_.Map(r); }

    /// domain map for given column block
    Epetra_Map& DomainMap(int r) const { return *domainmaps_.Map(r); }

    /// total matrix range map with all blocks
    Epetra_Map& FullRangeMap() const { return *rangemaps_.FullMap(); }

    /// total matrix domain map with all blocks
    Epetra_Map& FullDomainMap() const { return *domainmaps_.FullMap(); }

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
    virtual double NormInf() const;

    //@}

    /** \name Attribute access functions */
    //@{

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

  private:

    /// the full domain map together with all partial domain maps
    MultiMapExtractor domainmaps_;

    /// the full range map together with all partial range maps
    MultiMapExtractor rangemaps_;

    /// row major matrix block storage
    std::vector<SparseMatrix> blocks_;
  };


  /// Block matrix consisting of SparseMatrix blocks
  /*!

    There are strange algorithms that need to split a large sparse matrix into
    blocks. Such things happen, e.g., in FSI calculations with internal and
    interface splits, in fluid projection preconditioners or in contact
    simulations with slave and master sides. Unfortunatelly spliting a huge
    sparse matrix in (possibly) many blocks is nontrivial and expensive. So
    the idea here is to assemble into a block matrix in the first place.

    The difficulty with this approach is the handling of ghost entries in a
    parallel matrix. It is hard (expensive) to figure out to which column
    block each particular ghost entry belongs. That is why this class is
    templated with a Strategy. There is a default implementation for this
    template parameter DefaultBlockMatrixStrategy, that handles the most
    general case. That is DefaultBlockMatrixStrategy finds the right column
    block be heavy communication. But if there is some knowledge available in
    a particular case, it is easy to implement a specify Strategy that does
    not need to communicate that much.

    \author u.kue
    \date 02/08
  */
  template <class Strategy>
  class BlockSparseMatrix : public BlockSparseMatrixBase,
                            public Strategy
  {
  public:
    BlockSparseMatrix(const MultiMapExtractor& domainmaps,
                      const MultiMapExtractor& rangemaps,
                      int npr=81);

    void Assemble(const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    void Complete();
  };


  /// default strategy implementation for block matrix
  /*!

    This default implementation solves the ghost entry problem by remembering
    all ghost entries during the assembly in a private map. Afterwards
    Complete() needs to be called that finds the appropiate block for each
    ghost entry by communication an finally assembles these entries.

    This is the most general, most expensive implementation. You are
    encouraged to provide your own Strategy implementation if you know your
    specific block structure.

    \author u.kue
    \date 02/08
   */
  class DefaultBlockMatrixStrategy
  {
  public:

    /// construct with a block matrix base
    explicit DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat);

    /// find row block to a given row gid
    int RowBlock(int lrow, int rgid);

    /// find column block to a given column gid
    int ColBlock(int rblock, int lcol, int cgid);

    /// assemble into the given block
    void Assemble(double val,
                  int lrow, int rgid, int rblock,
                  int lcol, int cgid, int cblock);

    /// assemble the remaining ghost entries
    void Complete();

  private:

    /// my block matrix base
    BlockSparseMatrixBase& mat_;

    /// all ghost entries stored by row,column
    std::map<int,std::map<int,double> > ghost_;
  };

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const LINALG::SparseMatrix& mat);
ostream& operator << (ostream& os, const LINALG::BlockSparseMatrixBase& mat);

} // end of namespace LINALG


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy> Teuchos::RCP<LINALG::BlockSparseMatrix<Strategy> >
LINALG::SparseMatrix::Split(const MultiMapExtractor& domainmaps,
                            const MultiMapExtractor& rangemaps)
{
  if (domainmaps.NumMaps()==2 && rangemaps.NumMaps()==2)
  {
    const int npr = EpetraMatrix()->MaxNumEntries();
    Teuchos::RCP<BlockSparseMatrix<Strategy> > blockA =
      rcp(new LINALG::BlockSparseMatrix<Strategy>(domainmaps,rangemaps,npr));
    this->Split2x2(*blockA);
    return blockA;
  }
  else
    dserror("Currently can only split in 2x2");

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
LINALG::BlockSparseMatrix<Strategy>::BlockSparseMatrix(const MultiMapExtractor& domainmaps,
                                                       const MultiMapExtractor& rangemaps,
                                                       int npr)
  : BlockSparseMatrixBase(domainmaps,rangemaps,npr),
    // this was necessary, otherwise ambigous with copy constructor of Strategy
    Strategy((LINALG::BlockSparseMatrixBase&)(*this))
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
void LINALG::BlockSparseMatrix<Strategy>::Assemble(const Epetra_SerialDenseMatrix& Aele,
                                                    const std::vector<int>& lmrow,
                                                    const std::vector<int>& lmrowowner,
                                                    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();

  const int myrank = Comm().MyPID();
  // neccessary for compilation with suse gcc  4.1.0
  Strategy* strategy = (Strategy*)this;

  // loop rows of local matrix
  for (int lrow=0; lrow<lrowdim; ++lrow)
  {

    // check ownership of row
    if (lmrowowner[lrow]!=myrank)
      continue;

    int rgid = lmrow[lrow];
    int rblock = strategy->RowBlock(lrow,rgid);

    for (int lcol=0; lcol<lcoldim; ++lcol)
    {
      double val = Aele(lrow,lcol);
      int cgid = lmcol[lcol];
      int cblock = strategy->ColBlock(rblock,lcol,cgid);

      strategy->Assemble(val,lrow,rgid,rblock,lcol,cgid,cblock);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
void LINALG::BlockSparseMatrix<Strategy>::Complete()
{
  Strategy::Complete();
  BlockSparseMatrixBase::Complete();
}


#endif
#endif
