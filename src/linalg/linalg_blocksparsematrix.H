/*----------------------------------------------------------------------*/
/*! \file

\brief block sparse matrix

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef LINALG_BLOCKSPARSEMATRIX_H_
#define LINALG_BLOCKSPARSEMATRIX_H_

#include "linalg_sparsematrix.H"


namespace LINALG
{
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
    BlockSparseMatrixBase(const MultiMapExtractor& domainmaps, const MultiMapExtractor& rangemaps,
        int npr, bool explicitdirichlet = true, bool savegraph = false);

    virtual ~BlockSparseMatrixBase() {}

    /// make a copy of me
    virtual Teuchos::RCP<BlockSparseMatrixBase> Clone(DataAccess access) = 0;

    /// destroy the underlying Epetra objects
    virtual bool Destroy(bool throw_exception_for_blocks = true);

    /// setup of block preconditioners
    /*!
      This method can be implemented by subclasses that implement
      ApplyInverse() to execute a block preconditioner on the matrix.
     */
    virtual void SetupPreconditioner() {}

    /// Merge block matrix into a SparseMatrix
    Teuchos::RCP<SparseMatrix> Merge(bool explicitdirichlet = true) const;

    /** \name Block matrix access */
    //@{

    /// return block (r,c)
    const SparseMatrix& Matrix(int r, int c) const { return blocks_[r * Cols() + c]; }

    /// return block (r,c)
    SparseMatrix& Matrix(int r, int c) { return blocks_[r * Cols() + c]; }

    /// return block (r,c)
    inline const SparseMatrix& operator()(int r, int c) const { return Matrix(r, c); }

    /// return block (r,c)
    inline SparseMatrix& operator()(int r, int c) { return Matrix(r, c); }

    /// assign SparseMatrix to block (r,c)
    /*!
      \note The maps of the block have to match the maps of the given matrix.
     */
    void Assign(int r, int c, DataAccess access, const SparseMatrix& mat);

    //@}

    /** \name FE methods */
    //@{

    void Zero();
    void Reset();
    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void UnComplete();

    void ApplyDirichlet(
        const Teuchos::RCP<const Epetra_Vector> dbctoggle, bool diagonalblock = true);

    void ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock = true);

    /// derived
    virtual bool IsDbcApplied(const Epetra_Map& dbcmap, bool diagonalblock = true,
        const LINALG::SparseMatrix* trafo = NULL) const;

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
    const Epetra_Map& RangeMap(int r) const { return *rangemaps_.Map(r); }

    /// domain map for given column block
    const Epetra_Map& DomainMap(int r) const { return *domainmaps_.Map(r); }

    /// total matrix range map with all blocks
    const Epetra_Map& FullRangeMap() const { return *rangemaps_.FullMap(); }

    /// total matrix domain map with all blocks
    const Epetra_Map& FullDomainMap() const { return *domainmaps_.FullMap(); }

    /// total matrix domain map with all blocks (this is needed for
    /// consistency with LINALG::SparseMatrix)
    const Epetra_Map& DomainMap() const { return *domainmaps_.FullMap(); }

    /// total matrix row map with all blocks
    /*!
      \pre Filled()==true
     */
    Epetra_Map& FullRowMap() const { return *fullrowmap_; }

    /// total matrix column map with all blocks
    /*!
      \pre Filled()==true
     */
    Epetra_Map& FullColMap() const { return *fullcolmap_; }

    //@}

    /** \name Attribute set methods */
    //@{

    /// If set true, transpose of this operator will be applied.
    virtual int SetUseTranspose(bool UseTranspose);

    //@}

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Resolve virtual function of parent class
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Add a (transposed) BlockSparseMatrix: (*this) = (*this)*scalarB + A(^T)*scalarA
    virtual void Add(const BlockSparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB);

    /// Resolve virtual function of parent class
    virtual void Add(
        const SparseOperator& A, const bool transposeA, const double scalarA, const double scalarB);

    /// Resolve virtual function of parent class
    virtual void AddOther(SparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB) const;

    /// Resolve virtual function of parent class
    virtual void AddOther(BlockSparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB) const;

    /// Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
    virtual int Scale(double ScalarConstant);

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

    /// access to domain map extractor in derived classes
    const MultiMapExtractor& DomainExtractor() const { return domainmaps_; }

    /// access to range map extractor in derived classes
    const MultiMapExtractor& RangeExtractor() const { return rangemaps_; }

    /// friend functions
    friend Teuchos::RCP<BlockSparseMatrixBase> Multiply(const BlockSparseMatrixBase& A, bool transA,
        const BlockSparseMatrixBase& B, bool transB, bool explicitdirichlet, bool savegraph,
        bool completeoutput);

   protected:
    /// extract a partial map extractor from the full map extractor
    void GetPartialExtractor(const MultiMapExtractor& full_extractor,
        const std::vector<unsigned>& block_ids, MultiMapExtractor& partial_extractor) const;

   private:
    /// the full domain map together with all partial domain maps
    MultiMapExtractor domainmaps_;

    /// the full range map together with all partial range maps
    MultiMapExtractor rangemaps_;

    /// row major matrix block storage
    std::vector<SparseMatrix> blocks_;

    /// full matrix row map
    Teuchos::RCP<Epetra_Map> fullrowmap_;

    /// full matrix column map
    Teuchos::RCP<Epetra_Map> fullcolmap_;

    /// see matrix as transposed
    bool usetranspose_;
  };



  /// Block matrix consisting of SparseMatrix blocks
  /*!
      There are strange algorithms that need to split a large sparse matrix into
      blocks. Such things happen, e.g., in FSI calculations with internal and
      interface splits, in fluid projection preconditioners or in contact
      simulations with slave and master sides. Unfortunately splitting a huge
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
  class BlockSparseMatrix : public BlockSparseMatrixBase, public Strategy
  {
   public:
    BlockSparseMatrix(const MultiMapExtractor& domainmaps, const MultiMapExtractor& rangemaps,
        int npr = 81, bool explicitdirichlet = true, bool savegraph = false);

    /// clone the full block sparse matrix

    /** Do not forget to call Complete() after cloning, even if you
     *  use LINALG::View! */
    Teuchos::RCP<BlockSparseMatrixBase> Clone(DataAccess access);

    /// clone only a part of the block sparse matrix
    /** Do not forget to call Complete() after cloning, even if you
     *  use LINALG::View!
     *
     *  \param[in] access : consider copy or view of block matrices
     *  \param[in] row_block_ids : ID's of the row blocks to clone
     *  \param[in] col_block_ids : ID's of the column blocks to clone
     *
     *  \author hiermeier \date 04/17 */
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> Clone(DataAccess access,
        const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids);

    /// just a dummy that switches from strided assembly to standard assembly
    void Assemble(int eid, const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& Aele,
        const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
        const std::vector<int>& lmcol)
    {
      const int myrank = Comm().MyPID();
      Strategy::Assemble(eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
    }

#if 0
  void Assemble(int eid,
      const Epetra_SerialDenseMatrix& Aele,
      const std::vector<int>& lmrow,
      const std::vector<int>& lmrowowner,
      const std::vector<int>& lmcol)
  {
    const int myrank = Comm().MyPID();
    Strategy::Assemble(eid,myrank,Aele,lmrow,lmrowowner,lmcol);
  }
#endif
    /// single value assemble
    /*!

       \warning This method is less useful here. We just need to make the
       compiler happy. We assume "element matrices" of size 1x1 here. Do not use
       this method if your strategy depends on the real position of the dof in
       the element matrix.

     */
    void Assemble(double val, int rgid, int cgid) { Strategy::Assemble(val, rgid, cgid); }

    void Complete();

   private:
    /** \brief internal clone method which provides the possibility to clone only
     *         a sub-set of all blocks
     *
     *  This method is not supposed to be called from outside! See public variant.
     *
     *  \param[in] access : consider copy or view of block matrices
     *  \param[in] row_block_ids : ID's of the row blocks to clone
     *  \param[in] col_block_ids : ID's of the column blocks to clone
     *  \param[in] domain_extractor : necessary domain extractor
     *  \param[in] range_extractor : necessary range extractor
     *
     *  \author hiermeier \date 04/17 */
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> Clone(DataAccess access,
        const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids,
        const MultiMapExtractor& domain_extractor, const MultiMapExtractor& range_extractor);
  };


  /// default strategy implementation for block matrix
  /*!

      This default implementation solves the ghost entry problem by remembering
      all ghost entries during the assembly in a private map. Afterwards
      Complete() needs to be called that finds the appropriate block for each
      ghost entry by communication an finally assembles these entries.

      This is the most general, most expensive implementation. You are
      encouraged to provide your own Strategy implementation if you know your
      specific block structure.

      \sa BlockSparseMatrix

      \author u.kue
      \date 02/08
   */
  class DefaultBlockMatrixStrategy
  {
   public:
    /// construct with a block matrix base
    explicit DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat);

    /// assemble into the given block using nodal strides
    void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol);
#if 0
  /// assemble into the given block
  void Assemble(int eid,
      int myrank,
      const Epetra_SerialDenseMatrix& Aele,
      const std::vector<int>& lmrow,
      const std::vector<int>& lmrowowner,
      const std::vector<int>& lmcol);
#endif
    /// assemble into the given block
    void Assemble(double val, int rgid, int cgid);

    /// assemble the remaining ghost entries
    void Complete();

   protected:
    /// assemble into the given block
    void Assemble(double val, int lrow, int rgid, int rblock, int lcol, int cgid, int cblock);

    /// find row block to a given row gid
    int RowBlock(int rgid);

    /// find column block to a given column gid
    int ColBlock(int rblock, int cgid);

    /// access to the block sparse matrix for subclasses
    BlockSparseMatrixBase& Mat() { return mat_; }

   private:
    /// my block matrix base
    BlockSparseMatrixBase& mat_;

    /// all ghost entries stored by row,column
    std::map<int, std::map<int, double>> ghost_;

    /// scratch array for identifying column information
    std::vector<std::vector<int>> scratch_lcols_;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/

  /// output of BlockSparseMatrixBase
  std::ostream& operator<<(std::ostream& os, const LINALG::BlockSparseMatrixBase& mat);


  //////////////////////////////////
  /// helper functions

  /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
  /*!
    Multiply one matrix with another. Both matrices must be completed.
    Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

    \note This is a true parallel multiplication, even in the transposed case.

    \note Does call complete on C upon exit by default.

    \note In this version the flags explicitdirichlet and savegraph must be handed in.
          Thus, they can be defined explicitly, while in the standard version of Multipliy()
          above, result matrix C automatically inherits these flags from input matrix A

    \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
    \param transA         (in)     : flag indicating whether transposed of A should be used
    \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
    \param transB         (in)     : flag indicating whether transposed of B should be used
    \param explicitdirichlet (in)  : flag deciding on explicitdirichlet flag of C
    \param savegraph      (in)     : flag deciding on savegraph flag of C
    \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon
    output \return Matrix product A(^T)*B(^T)
  */
  Teuchos::RCP<BlockSparseMatrixBase> Multiply(const BlockSparseMatrixBase& A, bool transA,
      const BlockSparseMatrixBase& B, bool transB, bool explicitdirichlet, bool savegraph,
      bool completeoutput = true);

  Teuchos::RCP<BlockSparseMatrix<DefaultBlockMatrixStrategy>> BlockMatrix2x2(
      SparseMatrix& A00, SparseMatrix& A01, SparseMatrix& A10, SparseMatrix& A11);

  //////////////////////////////////


}  // end of namespace LINALG


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
LINALG::BlockSparseMatrix<Strategy>::BlockSparseMatrix(const MultiMapExtractor& domainmaps,
    const MultiMapExtractor& rangemaps, int npr, bool explicitdirichlet, bool savegraph)
    : BlockSparseMatrixBase(domainmaps, rangemaps, npr, explicitdirichlet, savegraph),
      // this was necessary, otherwise ambiguous with copy constructor of Strategy
      Strategy((LINALG::BlockSparseMatrixBase&)(*this))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
Teuchos::RCP<LINALG::BlockSparseMatrixBase> LINALG::BlockSparseMatrix<Strategy>::Clone(
    DataAccess access)
{
  std::vector<unsigned> row_block_ids(Rows());
  for (unsigned i = 0; i < static_cast<unsigned>(Rows()); ++i) row_block_ids[i] = i;

  std::vector<unsigned> col_block_ids(Cols());
  for (unsigned i = 0; i < static_cast<unsigned>(Cols()); ++i) col_block_ids[i] = i;

  return Clone(access, row_block_ids, col_block_ids, DomainExtractor(), RangeExtractor());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
Teuchos::RCP<LINALG::BlockSparseMatrixBase> LINALG::BlockSparseMatrix<Strategy>::Clone(
    DataAccess access, const std::vector<unsigned>& row_block_ids,
    const std::vector<unsigned>& col_block_ids, const MultiMapExtractor& domain_extractor,
    const MultiMapExtractor& range_extractor)
{
  int npr = Matrix(0, 0).MaxNumEntries();
  bool explicitdirichlet = Matrix(0, 0).ExplicitDirichlet();
  bool savegraph = Matrix(0, 0).SaveGraph();
  Teuchos::RCP<BlockSparseMatrixBase> bsm = Teuchos::rcp(new BlockSparseMatrix<Strategy>(
      domain_extractor, range_extractor, npr, explicitdirichlet, savegraph));

  for (std::vector<unsigned>::const_iterator r = row_block_ids.begin(); r != row_block_ids.end();
       ++r)
  {
    for (std::vector<unsigned>::const_iterator c = col_block_ids.begin(); c != col_block_ids.end();
         ++c)
    {
      bsm->Matrix(*r, *c).Assign(access, Matrix(*r, *c));
    }
  }
  return bsm;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
Teuchos::RCP<LINALG::BlockSparseMatrixBase> LINALG::BlockSparseMatrix<Strategy>::Clone(
    DataAccess access, const std::vector<unsigned>& row_block_ids,
    const std::vector<unsigned>& col_block_ids)
{
  if (std::lower_bound(row_block_ids.begin(), row_block_ids.end(), static_cast<unsigned>(Rows())) !=
      row_block_ids.end())
    dserror("The partial row block ids exceed the maximal possible id!");

  if (std::lower_bound(col_block_ids.begin(), col_block_ids.end(), static_cast<unsigned>(Cols())) !=
      col_block_ids.end())
    dserror("The partial column block ids exceed the maximal possible id!");

  if (row_block_ids.size() == 0 or col_block_ids.size() == 0)
    dserror("The provided row/col block id vector has a length of zero!");

  // extract domain extractors
  MultiMapExtractor p_domain_extractor;
  this->GetPartialExtractor(DomainExtractor(), col_block_ids, p_domain_extractor);

  // extract range extractors
  MultiMapExtractor p_range_extractor;
  this->GetPartialExtractor(RangeExtractor(), row_block_ids, p_range_extractor);

  return Clone(access, row_block_ids, col_block_ids, p_domain_extractor, p_range_extractor);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
void LINALG::BlockSparseMatrix<Strategy>::Complete()
{
  Strategy::Complete();
  BlockSparseMatrixBase::Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline int LINALG::DefaultBlockMatrixStrategy::RowBlock(int rgid)
{
  int rows = mat_.Rows();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    if (mat_.RangeMap(rblock).MyGID(rgid))
    {
      return rblock;
    }
  }
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline int LINALG::DefaultBlockMatrixStrategy::ColBlock(int rblock, int cgid)
{
  int cols = mat_.Cols();
  for (int cblock = 0; cblock < cols; ++cblock)
  {
    SparseMatrix& matrix = mat_.Matrix(rblock, cblock);

    // If we have a filled matrix we know the column map already.
    if (matrix.Filled())
    {
      if (matrix.ColMap().MyGID(cgid))
      {
        return cblock;
      }
    }

    // otherwise we can get just the non-ghost entries right now
    else if (mat_.DomainMap(cblock).MyGID(cgid))
    {
      return cblock;
    }
  }

  // ghost entries in a non-filled matrix will have to be done later

  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void LINALG::DefaultBlockMatrixStrategy::Assemble(int eid, int myrank,
    const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();

  dsassert(
      static_cast<int>(scratch_lcols_.size()) == mat_.Rows(), "Unexpected number of block rows");

  for (int rblock = 0; rblock < mat_.Rows(); ++rblock)
  {
    scratch_lcols_[rblock].resize(lcoldim);
    std::fill(scratch_lcols_[rblock].begin(), scratch_lcols_[rblock].end(), -1);
  }

  // loop rows of local matrix
  for (int lrow = 0; lrow < lrowdim; ++lrow)
  {
    // check ownership of row
    if (lmrowowner[lrow] != myrank) continue;

    int rgid = lmrow[lrow];
    int rblock = RowBlock(rgid);

    if (scratch_lcols_[rblock][0] == -1)
      for (int lcol = 0; lcol < lcoldim; ++lcol)
        scratch_lcols_[rblock][lcol] = ColBlock(rblock, lmcol[lcol]);

    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      double val = Aele(lrow, lcol);
      int cgid = lmcol[lcol];

      Assemble(val, lrow, rgid, rblock, lcol, cgid, scratch_lcols_[rblock][lcol]);
    }
  }
}

#if 0
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline
void LINALG::DefaultBlockMatrixStrategy::Assemble(int eid,
    int myrank,
    const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();

  // loop rows of local matrix
  for (int lrow=0; lrow<lrowdim; ++lrow)
  {

    // check ownership of row
    if (lmrowowner[lrow]!=myrank)
      continue;

    int rgid = lmrow[lrow];
    int rblock = RowBlock(rgid);

    for (int lcol=0; lcol<lcoldim; ++lcol)
    {
      double val = Aele(lrow,lcol);
      int cgid = lmcol[lcol];
      int cblock = ColBlock(rblock,cgid);

      Assemble(val,lrow,rgid,rblock,lcol,cgid,cblock);
    }
  }
}
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void LINALG::DefaultBlockMatrixStrategy::Assemble(double val, int rgid, int cgid)
{
  int rblock = RowBlock(rgid);
  int cblock = ColBlock(rblock, cgid);

  Assemble(val, 0, rgid, rblock, 0, cgid, cblock);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void LINALG::DefaultBlockMatrixStrategy::Assemble(
    double val, int lrow, int rgid, int rblock, int lcol, int cgid, int cblock)
{
#ifdef DEBUG
  if (rblock == -1) dserror("no block entry found for row gid=%d", rgid);
#endif

  if (cblock > -1)
  {
    SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
    matrix.Assemble(val, rgid, cgid);
  }
  else
  {
    // ghost entry in non-filled matrix. Save for later insertion.
    ghost_[rgid][cgid] += val;
  }
}



#endif
