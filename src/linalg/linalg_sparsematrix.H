/*!----------------------------------------------------------------------
\file linalg_sparsematrix.H

\brief Declaration of general BACI sparse matrix class

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef LINALG_SPARSEMATRIX_H_
#define LINALG_SPARSEMATRIX_H_

#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Comm.h>


#include "linalg_sparsematrixbase.H"
#include "linalg_mapextractor.H"
#include "../drt_lib/drt_dserror.H"

#include <boost/static_assert.hpp>


namespace LINALG
{
  // forward declarations:
  class BlockSparseMatrixBase;
  template <class Strategy> class BlockSparseMatrix;

  /// A single sparse matrix enhanced with features for FE simulations
  /*!

    A single sparse matrix. Internally we have an Epetra_CrsMatrix
    (or the sub-class Epetra_FECrsMatrix). So we have
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

    If FE_MATRIX is set as a flag, the implementation is based on an Epetra_FECrsMatrix.
    Nonlocal matrix values can be assembled by invoking FEAssemble()-methods
    instead of Assemble()-methods. Internally this will cause
    the GlobalAssemble()-method to distribute the nonlocal values to the owning
    procs before FillComplete is called on the matrix.
    Since Epetra_FECrsMatrix is a sub-class of Epetra_CrsMatrix,
    all other functionality can be used as described above.
    All SparseMatrix constructors create Epetra_CrsMatrices by default.

    \note A large part of the SparseMatrix interface consists of methods from
    the internal Epetra_CrsMatrix. If there are methods in Epetra_CrsMatrix
    and not in SparseMatrix that you would like to call (for legitimate
    reasons!) please add them to the SparseMatrix.

    \author u.kue
    \date 02/08
   */
class SparseMatrix : public SparseMatrixBase
{
  public:
    /*!
      flag for the underlying Epetra matrix type; CRS_MATRIX means that the implementation of the
      linalg sparse matrix is based on an Epetra_CrsMatrix.
      If FE_MATRIX is chosen, the implementation is based on an Epetra_FECrsMatrix.
      This allows the assembly and handling of nonlocal values in addition to Epetra_CrsMatrix-functionality.
     */
    enum MatrixType {CRS_MATRIX, FE_MATRIX};

    /// The following dummy templated constructors will catch all attempts to call any constructor with arguments
    /// that require implicit casting. Therefore you have to instantiate this class with exactly the types
    /// matching one of the constructors
#if __cplusplus >= 201103L
  template <typename... T> SparseMatrix(T...) = delete;
#else
  template <typename T0>                                                     SparseMatrix(T0)                 {BOOST_STATIC_ASSERT_MSG(sizeof(T0) == 0, "dummy templated constructor") BACI_ATTRIBUTE_UNUSED;}
  template <typename T0, typename T1>                                        SparseMatrix(T0, T1)             {BOOST_STATIC_ASSERT_MSG(sizeof(T0) == 0, "dummy templated constructor") BACI_ATTRIBUTE_UNUSED;}
  template <typename T0, typename T1, typename T2>                           SparseMatrix(T0, T1, T2)         {BOOST_STATIC_ASSERT_MSG(sizeof(T0) == 0, "dummy templated constructor") BACI_ATTRIBUTE_UNUSED;}
  template <typename T0, typename T1, typename T2, typename T3>              SparseMatrix(T0, T1, T2, T3)     {BOOST_STATIC_ASSERT_MSG(sizeof(T0) == 0, "dummy templated constructor") BACI_ATTRIBUTE_UNUSED;}
  template <typename T0, typename T1, typename T2, typename T3, typename T4> SparseMatrix(T0, T1, T2, T3, T4) {BOOST_STATIC_ASSERT_MSG(sizeof(T0) == 0, "dummy templated constructor") BACI_ATTRIBUTE_UNUSED;}
#endif

    /// construction of sparse matrix
    SparseMatrix(
      Teuchos::RCP<Epetra_CrsGraph> crsgraph,
      Teuchos::RCP<LINALG::MultiMapExtractor> dbcmaps);

    /// construction of sparse matrix
    SparseMatrix(
        const Epetra_Map&   rowmap,
        const int           npr,
        bool                explicitdirichlet=true,
        bool                savegraph=false,
        MatrixType          matrixtype = CRS_MATRIX);

    /// construction of sparse matrix using an individual estimate for number of non-zeros per row
    SparseMatrix(
        const Epetra_Map&   rowmap,
        std::vector<int> &  numentries,
        bool                explicitdirichlet=true,
        bool                savegraph=false,
        MatrixType          matrixtype = CRS_MATRIX);

    /// construction of sparse matrix
    /*!
       Makes either a deep copy of the Epetra_CrsMatrix or Epetra_FECrsMatrix.
       (Note: \pre matrix.Filled()==true)
       or an implicit construction from a Epetra_CrsMatrix or Epetra_FECrsMatrix
       where the Teuchos::RCP is copied internally leading to a new view on the
       Epetra_CrsMatrix or Epetra_FECrsMatrix.
     */
    SparseMatrix(
        Teuchos::RCP<Epetra_CrsMatrix>  matrix,
        DataAccess               access,
        bool                            explicitdirichlet=true,
        bool                            savegraph=false,
        MatrixType                      matrixtype = CRS_MATRIX);

    /// construction of a diagonal matrix from a vector
    /*!
      Creates diagonal matrix with range and domain map equal to vector map.
      Sets diagonal values from vector and does NOT call Complete() on matrix

      Allocates new memory for the Epetra_CrsMatrix or Epetra_FECrsMatrix
     */
    SparseMatrix(
        const Epetra_Vector&  diag,
        bool                  explicitdirichlet=true,
        bool                  savegraph=false,
        MatrixType            matrixtype = CRS_MATRIX);

    /// Copy constructor. Deep copy or view on matrix.
    /*!
      \warning A view assignment will have your matrix use the same internal
      data as the original matrix. Changes on one side will affect the
      other. However, some methods like Zero() or Reset() can, depending on
      the SparseMatrix flags, cut the connection. Do not rely on the view if
      you change one of these matrices! The view assignment is meant to
      provide a slim copy operation that transfers ownership from one matrix
      to the other before the original matrix is destroyed. Do a deep copy if
      both matrices are meant to live on.

      \param mat matrix to assign from
      \param access how to treat this assignment: Copy or View
     */
    SparseMatrix(const SparseMatrix& mat, DataAccess access = Copy);

    /// Destructor
    virtual ~SparseMatrix();

    /// Assignment operator. Makes a deep copy.
    SparseMatrix& operator=(const SparseMatrix& mat);

    /// Assignment method. Deep copy or view on matrix.
    /*!
      Explicit method for the assignment operator. You can make an explicit
      copy of the internal Epetra_CrsMatrix (or Epetra_FECrsMatrix)
      or have a second view on it.

      \warning A view assignment will have your matrix use the same internal
      data as the original matrix. Changes on one side will affect the
      other. However, some methods like Zero() or Reset() can, depending on
      the SparseMatrix flags, cut the connection. Do not rely on the view if
      you change one of these matrices! The view assignment is meant to
      provide a slim copy operation that transfers ownership from one matrix
      to the other before the original matrix is destroyed. Do a deep copy if
      both matrices are meant to live on.

      \param access how to treat this assignment: Copy or View
      \param mat matrix to assign from
     */
    void Assign(DataAccess access, const SparseMatrix& mat);

    /** \name FE methods */
    //@{

    /// set all matrix entries to zero
    void Zero();

    /// throw away the matrix and its graph and start anew
    void Reset();

    /// destroy the underlying Epetra objects
    virtual bool Destroy( bool throw_exception = true );

    /// assemble method for Epetra_CrsMatrices, if ONLY local values are assembled
    virtual void Assemble(int eid,
                          const std::vector<int>& lmstride,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(eid,lmstride,Aele,lm,lmowner,lm);
    }

    /// assemble method for Epetra_CrsMatrices, if ONLY local values are assembled
    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(eid,Aele,lm,lmowner,lm);
    }

    /// assemble method for Epetra_CrsMatrices, if ONLY local values are assembled
    void Assemble(int eid,
                  const std::vector<int>& lmstride,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    /// assemble method for Epetra_CrsMatrices, if ONLY local values are assembled
    void Assemble(int eid,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    /// single value assemble used by BlockSparseMatrix
    void Assemble(double val, int rgid, int cgid);

    /*!
      Assemble method for an Epetra_FECrsMatrix.
      This method is also able to handle the assembly of nonlocal values.
      It sets the doGlobalAssemble-flag to true and causes the
      GlobalAssemble() method to redistribute the non-local
      values to their owning procs, such that FillComplete can be safely
      called on this matrix.

      NOTE: This methods checks if rowowner == myrank. Only in this case
      values are set. This is needed if the method is called in a loop over
      column elements (which is the standard in baci) to avoid multiple same entries.
     */
    void FEAssemble(
            const Epetra_SerialDenseMatrix& Aele,
            const std::vector<int>& lmrow,
            const std::vector<int>& lmrowowner,
            const std::vector<int>& lmcol);

    /*!
      Assemble method for an Epetra_FECrsMatrices.
      This method is also able to handle the assembly of nonlocal values.
      It sets the doGlobalAssemble-flag to true and causes the
      GlobalAssemble() method to redistribute the non-local
      values to their owning procs, such that FillComplete can be safely
      called on this matrix.
     */
    void FEAssemble(
            const Epetra_SerialDenseMatrix& Aele,
            const std::vector<int>& lmrow,
            const std::vector<int>& lmcol);

    /*!
      Assemble method for an Epetra_FECrsMatrices.
      This method is also able
      to handle the assembly of nonlocal values.
      It sets the doGlobalAssemble-flag to true and causes the
      GlobalAssemble -method() to redistribute the non-local
      values to their owning procs, such that FillComplete can be savely
      called on this matrix.
     */
    void FEAssemble(double val, int rgid, int cgid);


    /*!
      The GlobalAssembleMethod() distributes nonlocal values to their owning procs
      for Epetra_FECrsMatrices.
      Afterwards Fillcomplete is called such as for Epetra_CrsMatrices.
     */
    void Complete();

    /*!
      The GlobalAssembleMethod() distributes nonlocal values to their owning procs
      for Epetra_FECrsMatrices.
      Afterwards Fillcomplete is called such as for Epetra_CrsMatrices.
     */
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void UnComplete();

    void ApplyDirichlet(const Teuchos::RCP<const Epetra_Vector> dbctoggle, bool diagonalblock=true);

    /// Apply dirichlet boundary condition to a matrix.
    ///
    ///  This method blanks the rows associated with Dirichlet DOFs
    ///  and puts a 1.0 at the diagonal entry if diagonlblock==true.
    ///  Only the rows are blanked, the columns are not touched.
    ///  We are left with a non-symmetric matrix, if the original
    ///  matrix was symmetric. However, the blanking of columns is computationally
    ///  quite expensive, because the matrix is stored in a sparse and distributed
    ///  manner.
    void ApplyDirichlet(const Epetra_Map& dbctoggle, bool diagonalblock=true);

    /// Apply dirichlet boundary condition to a matrix using a #trafo matrix
    ///
    /// This method the same as the method #ApplyDirichlet, but instead of
    /// 1.0 on diagonal, the corresponding row of #trafo is inserted. This
    /// is needed to treat efficiently Dirichlet BCs with local co-ordinate systems.
    /// The transformation matrix #trafo basically holds rotation matrices
    /// for the DOFs of the nodes.
    void ApplyDirichletWithTrafo(Teuchos::RCP<const LINALG::SparseMatrix> trafo,
                                 const Epetra_Map& dbctoggle,
                                 bool diagonalblock=true,
                                 bool complete=true);

    /// create matrix that contains all Dirichlet lines from my
    Teuchos::RCP<SparseMatrix> ExtractDirichletRows(const Teuchos::RCP<Epetra_Vector> dbctoggle);

    /// create matrix that contains all Dirichlet lines from my
    Teuchos::RCP<SparseMatrix> ExtractDirichletRows(const Epetra_Map& dbctoggle);

    //@}

    /** \name Matrix Properties Query Methods */
    //@{

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

    /// Return matrix type
    MatrixType GetMatrixtype() const {return matrixtype_; }

    //@}

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /** \name Utility functions */
    //@{

    /// Compute transposed matrix explicitly
    /*!
      \note This is an expensive operation!
     */
    Teuchos::RCP<SparseMatrix> Transpose();

    /// derived
    virtual int ReplaceRowMap(const Epetra_BlockMap& newmap);

  using SparseMatrixBase::Add;

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
    void Add(const SparseMatrixBase& A,
             const bool transposeA,
             const double scalarA,
             const double scalarB);

    /// Add a (transposed) Epetra_CrsMatrix to another: (*this) = (*this)*scalarB + A(^T)*scalarA
    void Add(const Epetra_CrsMatrix& A,
             const bool transposeA,
             const double scalarA,
             const double scalarB);

    /// Put a Epetra_CrsMatrix (partially) onto another: (*this) = B*scalarB
    /*!

    Put one matrix onto another. The matrix (*this) to be added to must not be
    completed. Sparsity patterns of A and (*this) need not match and A and (*this) can be
    nonsymmetric in value and pattern.  Row map of A has to be a
    processor-local subset of the row map of (*this).

    \param A          (in)     : Matrix to add to this (must have Filled()==true)
    \param scalarA    (in)     : scaling factor for #A
    \param rowmap     (in)     : to put selectively on rows in #rowmap (inactive if ==Teuchos::null)
    */
    void Put(const LINALG::SparseMatrix& A,
             const double scalarA,
             Teuchos::RCP<const Epetra_Map> rowmap);

    /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
    /*!
      Multiply one matrix with another. Both matrices must be completed.
      Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

      \note This is a true parallel multiplication, even in the transposed case.

      \note Does call complete on C upon exit by default.

      \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
      \param transA         (in)     : flag indicating whether transposed of A should be used
      \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
      \param transB         (in)     : flag indicating whether transposed of B should be used
      \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
      \return Matrix product A(^T)*B(^T)
    */
    friend Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                               const SparseMatrix& B, bool transB,
                                               bool completeoutput);


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
      \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
      \return Matrix product A(^T)*B(^T)
    */
    friend Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                               const SparseMatrix& B, bool transB,
                                               bool explicitdirichlet,
                                               bool savegraph,
                                               bool completeoutput);

    /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
    /*!
      Multiply one matrix with another. Both matrices must be completed.
      Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

      \note This is a true parallel multiplication, even in the transposed case.

      \note Does call complete on C upon exit by default.

      \note Uses ML as multiplication kernel, not EpetraExt.

      \note In this version the flags explicitdirichlet and savegraph must be handed in.
            Thus, they can be defined explicitly, while in the standard version of Multipliy()
            above, result matrix C automatically inherits these flags from input matrix A

      \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
      \param transA         (in)     : flag indicating whether transposed of A should be used
      \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
      \param transB         (in)     : flag indicating whether transposed of B should be used
      \param explicitdirichlet (in)  : flag deciding on explicitdirichlet flag of C
      \param savegraph      (in)     : flag deciding on savegraph flag of C
      \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
      \return Matrix product A(^T)*B(^T)
    */
    friend Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A, bool transA,
                                                 const SparseMatrix& B, bool transB,
                                                 bool explicitdirichlet,
                                                 bool savegraph,
                                                 bool completeoutput);

    /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
    /*!
      Multiply one matrix with another. Both matrices must be completed.
      Respective Range, Row and Domain maps of A and B have to match.

      \note This is a true parallel multiplication, even in the transposed case.

      \note Does call complete on C upon exit by default.

      \note Uses ML as multiplication kernel, not EpetraExt.

      \note In this version the flags explicitdirichlet and savegraph must be handed in.
            Thus, they can be defined explicitly, while in the standard version of Multipliy()
            above, result matrix C automatically inherits these flags from input matrix A

      \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
      \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
      \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
      \return Matrix product A*B
    */
    friend Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
                                                 const SparseMatrix& B,
                                                 bool complete);


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
      Split(const MultiMapExtractor& domainmaps,const MultiMapExtractor& rangemaps) const;

    /// binary dump of matrix for debugging
    void Dump(std::string filename);

    /// load a dump
    void Load(const Epetra_Comm& comm, std::string& filename);

    //@}

  private:

    /// Split matrix in 2x2 blocks, where main diagonal blocks have to be square
    /*!
       Used by public Split, does not call Complete() on output matrix.
     */
    void Split2x2(BlockSparseMatrixBase& Abase) const;

    /// Split matrix in MxN blocks
    /*!
       Used by public Split, does not call Complete() on output matrix.
     */
    void SplitMxN(BlockSparseMatrixBase& Abase) const;

    /// saved graph (if any)
    Teuchos::RCP<Epetra_CrsGraph> graph_;

    /// Dirichlet row map (if known)
    Teuchos::RCP<LINALG::MultiMapExtractor> dbcmaps_;

    /// whether to modify the matrix graph on apply Dirichlet
    bool explicitdirichlet_;

    /// whether to save the graph and assemble to a filled matrix next time
    bool savegraph_;

    /// matrix type (Epetra_CrsMatrix or Epetra_FECrsMatrix)
    MatrixType matrixtype_;

  };

  /// output of SparseMatrix
  std::ostream& operator << (std::ostream& os, const LINALG::SparseMatrix& mat);

  /// matrix matrix multiplication
  Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                      const SparseMatrix& B, bool transB,
                                      bool completeoutput = true);

  /// matrix matrix multiplication with explicitly defined flags
  Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA,
                                      const SparseMatrix& B, bool transB,
                                      bool explicitdirichlet,
                                      bool savegraph,
                                      bool completeoutput = true);

  /// merge 2x2 system
  Teuchos::RCP<SparseMatrix> Merge(const LINALG::SparseMatrix& Aii,
      const LINALG::SparseMatrix& Aig,
      const LINALG::SparseMatrix& Agi,
      const LINALG::SparseMatrix& Agg);


  Teuchos::RCP<LINALG::SparseMatrix> Eye(const Epetra_Map& map);

}  // end of LINALG namespace



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*!
  \brief Split matrix in either 2x2 or NxN blocks (with N>2)

  Split given matrix 'this' into 2x2 or NxN block matrix and
  return result as templated BlockSparseMatrix.
  The MultiMapExtractor's provided have to be 2x2 or NxN maps,
  otherwise this method will throw an error.

  \note This is an expensive operation!

  \note This method will NOT call Complete() on the output
        BlockSparseMatrix.
 */
template <class Strategy> Teuchos::RCP<LINALG::BlockSparseMatrix<Strategy> >
LINALG::SparseMatrix::Split(
    const MultiMapExtractor& domainmaps,
    const MultiMapExtractor& rangemaps) const
{
  // initialize resulting BlockSparseMatrix. no need to provide estimates of nonzeros because
  // all entries will be inserted at once anyway
  Teuchos::RCP<BlockSparseMatrix<Strategy> > blockA =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<Strategy>(domainmaps,rangemaps,0,explicitdirichlet_,savegraph_));

  // perform matrix splitting
  if(domainmaps.NumMaps() == 2 && rangemaps.NumMaps() == 2)
    this->Split2x2(*blockA);
  else if(domainmaps.NumMaps() > 0 && rangemaps.NumMaps() > 0)
    this->SplitMxN(*blockA);
  else
    dserror("Invalid number %d of row blocks or %d of column blocks for splitting operation!",rangemaps.NumMaps(),domainmaps.NumMaps());

  // return resulting BlockSparseMatrix
  return blockA;
}

#endif
