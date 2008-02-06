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
#include "Epetra_LAPACK.h"
#include "Teuchos_RefCountPtr.hpp"
#include "drt_exporter.H"

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
                   (need not be exact, better should be too big rather then too small)
  */
  RefCountPtr<Epetra_CrsMatrix> CreateMatrix(const Epetra_Map& rowmap, const int npr);

  /*!
  \brief Create a new Epetra_Vector and return RefcountPtr to it

  \param rowmap (in): row map of vector
  \param init (in): initializa vector to zero upon construction
  */
  RefCountPtr<Epetra_Vector> CreateVector(const Epetra_Map& rowmap, const bool init = true);

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
  void Export(const Epetra_MultiVector& source, Epetra_MultiVector& target);

  /*!
  \brief Make Epetra_SerialDenseMatrix symmetric by averaging upper and lower traingular indices

  \param A (in/out): Matrix to be symmetrised
  */
  void SymmetriseMatrix(Epetra_SerialDenseMatrix& A);

  /*!
  \brief Compute all eigenvalues and, optionally, eigenvectors
   of a real symmetric matrix A

  \param A (in/out): Matrix to be analysed, if eigv=true stores eigenvectors
  \param L (in/out): Vector of eigenvalues in ascending order
  \param dim (in) :  Dimension of matrix
  \param eigv (in):  flag to evaluate also eigenvectors
  */
  void SymmetricEigen(Epetra_SerialDenseMatrix& A,
                      Epetra_SerialDenseVector& L,const int dim, const char eigv='N');

  /*!
  \brief Multiply two 2nd order tensors A x B and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with A x B
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding A x B
  */
  void ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);

  /*!
  \brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with (A x B + B x A)
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding (A x B + B x A)
  */
  void ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);

  /*!
  \brief Multiply two 2nd order tensors A o B and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with A o B
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding A o B
  */
  void ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);
  /*!
  \brief Invert a symmetric dim*dim square matrix

  \param A (in/out): Matrix to be inverted
  \param dim (in) :  Dimension of matrix
  */
  void SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
  \brief Invert a nonsymmetric dim*dim square matrix

  \param A (in/out): Matrix to be inverted
  \param dim (in) :  Dimension of matrix
  */
  void NonSymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
  \brief Assemble a Epetra_SerialDenseMatrix into a Epetra_CrsMatrix

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator and rowmap from matrix A to determine ownerships.
  Local matrix Aele has to be square.

  If matrix A is Filled(), it stays so and you can only assemble to
  places already masked. An attempt to assemble into a non-existing
  place is a grave mistake.

  If matrix A is not Filled(), the matrix is enlarged as required.

  \note Assembling to a non-Filled() matrix is much more expensive
        than to a Filled() matrix. If the sparse mask
        does not change it pays to keep the matrix around and assemble
        into the Filled() matrix.

  \param A (out)   : Sparse matrix to be assembled on
  \param Aele (in) : dense matrix to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owner procs of gids
  */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
                const vector<int>& lm, const vector<int>& lmowner);

  /*!
  \brief Assemble a Epetra_SerialDenseMatrix into a Epetra_CrsMatrix

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator and rowmap from matrix A to determine ownerships.
  Local matrix Aele may be \b square or \b rectangular.

	This version of 'Assemble' does not work for a matrix A that is already
	Filled()! If matrix A is not Filled(), it will be enlarged as required.

	\note The user must provide an \b additional input vector 'lmcol'
				containing the column gids for assembly seperately!

  \param A (out)         : Sparse matrix to be assembled on
  \param Aele (in)       : dense matrix to be assembled
  \param lmrow (in)      : vector with row gids
  \param lmrowowner (in) : vector with owner procs of row gids
  \param lmcol (in)      : vector with column gids
  */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
                const vector<int>& lmrow, const vector<int>& lmrowowner,
                const vector<int>& lmcol);

  /*!
  \brief Assemble a Epetra_SerialDenseVector into a Epetra_Vector

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator from vector V to determine ownerships.

  \param V (out)   : Vector to be assembled on
  \param Vele (in) : dense vector to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owner procs of gids
  */
  void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                const vector<int>& lm, const vector<int>& lmowner);

  /*!
  \brief Assemble a Epetra_SerialDenseVector into a Epetra_MultiVector

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator from vector V to determine ownerships.

  \param V (out)   : Vector to be assembled on
  \param n (in)   : column index of MultiVector to be assembled on
  \param Vele (in) : dense vector to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owner procs of gids
  */
  void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                const vector<int>& lm, const vector<int>& lmowner);


  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)
  */
  void Complete(Epetra_CrsMatrix& A);

  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)

  This is the RCP wrapper of the above method.
  */
  inline void Complete(RCP<Epetra_CrsMatrix> A)
  { LINALG::Complete(*A); return; }

  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)
  */
  void Complete(Epetra_CrsMatrix& A, const Epetra_Map& domainmap, const Epetra_Map& rangemap);

  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)

  This is the RCP wrapper of the above method.
  */
  inline void Complete(RCP<Epetra_CrsMatrix> A, const Epetra_Map& domainmap, const Epetra_Map& rangemap)
  { LINALG::Complete(*A,domainmap,rangemap); return; }

  /*!
  \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Row map of A has to be a processor-local subset of the row map of B.


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
  \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Row map of A has to be a processor-local subset of the row map of B.


  Note that this is a true parallel add, even in the transposed case!
  This is the RCP wrapper of the above method.

  \param A          (in)     : Matrix to add to B (must have Filled()==true)
  \param transposeA (in)     : flag indicating whether transposed of A should be used
  \param scalarA    (in)     : scaling factor for A
  \param B          (in/out) : Matrix to be added to (must have Filled()==false)
  \param scalarB    (in)     : scaling factor for B
  */
  inline void Add(const RCP<Epetra_CrsMatrix> A, const bool transposeA, const double scalarA,
           RCP<Epetra_CrsMatrix> B, const double scalarB)
  { LINALG::Add(*A,transposeA,scalarA,*B,scalarB); return; }

  /*!
  \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

  Returns RCP to the transposed matrix of the input matrix A.

  \param A          (in)     : Matrix to transpose (must have Filled()==true)

  */
  RCP<Epetra_CrsMatrix> Transpose(const Epetra_CrsMatrix& A);

  /*!
  \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

  Returns RCP to the transposed matrix of the input matrix A.
  This is the RCP wrapper of the above method.

  \param A          (in)     : Matrix to transpose (must have Filled()==true)

  */
  inline RCP<Epetra_CrsMatrix> Transpose(const RCP<Epetra_CrsMatrix> A)
  { return LINALG::Transpose(*A); }

  /*!
  \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

  Multiply one matrix with another. Both matrices must be completed. Sparsity
  Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

  Note that this is a true parallel multiplication, even in the transposed case!

  \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
  \param transA     (in)     : flag indicating whether transposed of A should be used
  \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
  \param transB     (in)     : flag indicating whether transposed of B should be used
  \return Matrix product A(^T)*B(^T)
  */
  RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
                                 const Epetra_CrsMatrix& B, bool transB);

  /*!
  \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

  Multiply one matrix with another. Both matrices must be completed. Sparsity
  Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

  Note that this is a true parallel multiplication, even in the transposed case!
  This is the RCP wrapper of the above method.

  \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
  \param transA     (in)     : flag indicating whether transposed of A should be used
  \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
  \param transB     (in)     : flag indicating whether transposed of B should be used
  \return Matrix product A(^T)*B(^T)
  */
  inline RCP<Epetra_CrsMatrix> Multiply(const RCP<Epetra_CrsMatrix>& A, bool transA,
                                        const RCP<Epetra_CrsMatrix>& B, bool transB)
  { return Multiply(*A,transA,*B,transB); }

  /*!
  \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

  Multiply one matrix with another. All input matrices must be completed. Sparsity
  Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

  Note that this is a true parallel multiplication, even in the transposed case!

  \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
  \param transA     (in)     : flag indicating whether transposed of A should be used
  \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
  \param transB     (in)     : flag indicating whether transposed of B should be used
  \param C          (in)     : Matrix C (must have Filled()==true)
  \param transC     (in)     : flag indicating whether transposed of C should be used
  \return Matrix product A(^T)*B(^T)*C(^T)
  */
  RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
                                 const Epetra_CrsMatrix& B, bool transB,
                                 const Epetra_CrsMatrix& C, bool transC);

  /*!
  \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

  Multiply one matrix with another. All input matrices must be completed. Sparsity
  Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

  Note that this is a true parallel multiplication, even in the transposed case!
  This is the RCP wrapper of the above method.

  \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
  \param transA     (in)     : flag indicating whether transposed of A should be used
  \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
  \param transB     (in)     : flag indicating whether transposed of B should be used
  \param C          (in)     : Matrix C (must have Filled()==true)
  \param transC     (in)     : flag indicating whether transposed of C should be used
  \return Matrix product A(^T)*B(^T)*C(^T)
  */
  inline RCP<Epetra_CrsMatrix> Multiply(const RCP<Epetra_CrsMatrix>& A, bool transA,
                                        const RCP<Epetra_CrsMatrix>& B, bool transB,
                                        const RCP<Epetra_CrsMatrix>& C, bool transC)
  { return Multiply(*A,transA,*B,transB,*C,transC); }

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations

  Modifies a system of equations such that dirichlet boundary conditions are enforced.
  Prescribed dirichlet BC values are supplied in dbcval and dbctoggle, where
  a prescribed value is dbcval[i] and dbctoggle[i] = 1.0. No BC is enforced in
  all places where dbctoggle[i] != 1.0.<br>
  Let us denote the \f$ A_{2 \times 2} \f$ blocks of \f$A\f$ by
  \f$A_{ff}, A_{fD}, A_{Df}, A_{DD}\f$, where \f$f\f$ stands for 'free' and
  \f$D\f$ stands for 'Dirichlet BC'. Then, after a call to this method<br>

  \f$ A_{ff} = A_{ff}, \f$<br>
  \f$ A_{fD} = A_{fD}, \f$<br>
  \f$ A_{Df} = 0_{Df}, \f$<br>
  \f$ A_{DD} = I_{DD}, \f$<br>
  \f$ x_{D} = dbcval_{D}, \f$<br>
  \f$ b_{D} = dbcval_{D} \f$<br>

  and

  \f$ A_{ff} x_f + A_{fD} x_D = b_f \f$<br>
  \f$ 0 x_f + I_{DD} x_D = x_D \f$.<br>

  \note The matrix is then nonsymmetric. When using iterative methods on this
        linear system of equations that depend on the symmetry of the matrix (such as e.g. CG),
        the initial guess supplied to the solver has to be exact at the
        Dirichlet BCs. This should be easy, as the values at the Dirichlet BCs
        are known.

  \note The mask of matrix \f$A\f$ is not modified. That is the
        entries in \f$A_{Df}\f$ and \f$A_{DD}\f$ are set to zero, not
        removed. This way the matrix can be reused in the next step.

  \param A         (in/out) : Matrix of Ax=b
  \param x         (in/out) : initial guess vector x of Ax=b
  \param b         (in/out) : rhs vector b of Ax=b
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

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations


  \param x (in/out)         : vector x of Ax=b
  \param b (in/out)         : vector b of Ax=b
  \param dbcval (in)        : vector holding values that are supposed to be prescribed
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_Vector>&      x,
                              RefCountPtr<Epetra_Vector>&      b,
                              const RefCountPtr<Epetra_Vector> dbcval,
                              const RefCountPtr<Epetra_Vector> dbctoggle);


  /*!
  \brief split a matrix into a 2x2 block system where the rowmap of one of the blocks is given

  Splits a given matrix into a 2x2 block system where the rowmap of one of the blocks is given
  on input. Blocks A11 and A22 are assumed to be square.
  All values on entry have to be Teuchos::null except the given rowmap and matrix A.
  Note that either A11rowmap or A22rowmap or both have to be nonzero. In case
  both rowmaps are supplied they have to be an exact and nonoverlapping split of A->RowMap().
  Matrix blocks are FillComplete() on exit.

  \param A         : Matrix A on input
  \param A11rowmap : rowmap of A11 or null
  \param A22rowmap : rowmap of A22 or null
  \param A11       : on exit matrix block A11
  \param A12       : on exit matrix block A12
  \param A21       : on exit matrix block A21
  \param A22       : on exit matrix block A22
  */
  bool SplitMatrix2x2(RefCountPtr<Epetra_CrsMatrix> A,
                      RefCountPtr<Epetra_Map>& A11rowmap,
                      RefCountPtr<Epetra_Map>& A22rowmap,
                      RefCountPtr<Epetra_CrsMatrix>& A11,
                      RefCountPtr<Epetra_CrsMatrix>& A12,
                      RefCountPtr<Epetra_CrsMatrix>& A21,
                      RefCountPtr<Epetra_CrsMatrix>& A22);

  /*!
  \brief split a matrix into a 2x2 block system

  Splits a given matrix into a 2x2 block system. All values on entry have to be
  Teuchos::null except the given rowmap(s) / domainmap(s) and matrix A.
  Note that either A11rowmap or A22rowmap or both have to be nonzero!
  Note that either A11domainmap or A22domainmap or both have to be nonzero!
  In case both rowmaps / domainmaps are supplied they have to be an exact and
  nonoverlapping split of A->RowMap() / A->DomainMap().
  Matrix blocks are FillComplete() on exit.

  \param A            : Matrix A on input
  \param A11rowmap    : rowmap of A11 or null
  \param A22rowmap    : rowmap of A22 or null
  \param A11domainmap : domainmap of A11 or null
  \param A22domainmap : domainmap of A22 or null
  \param A11          : on exit matrix block A11
  \param A12          : on exit matrix block A12
  \param A21          : on exit matrix block A21
  \param A22          : on exit matrix block A22
  */
  bool SplitMatrix2x2(RCP<Epetra_CrsMatrix> A,
  										RCP<Epetra_Map>& A11rowmap,
                      RCP<Epetra_Map>& A22rowmap,
                      RCP<Epetra_Map>& A11domainmap,
                      RCP<Epetra_Map>& A22domainmap,
                      RCP<Epetra_CrsMatrix>& A11,
                      RCP<Epetra_CrsMatrix>& A12,
                      RCP<Epetra_CrsMatrix>& A21,
                      RCP<Epetra_CrsMatrix>& A22);

  /*!
  \brief split a rowmap of matrix A

  splits A->RowMap() into 2 maps and returns them, where one of the rowmaps has
  to be given on input

  \param Amap      : Map to split on input
  \param Agiven    : on entry submap that is given and part of Amap
  \return the remainder map of Amap that is not overlapping with Agiven
  */
  Epetra_Map* SplitMap(const Epetra_Map& Amap,
                       const Epetra_Map& Agiven);

  /*!
  \brief merges two given Epetra_Maps

  merges input map1 and input map2, both of which have to be unique,
  but may be overlapping, to a new map and returns RCP to it.

  \param map1      : one map to be merged
  \param map2      : the other map to be merged
  \return the (sorted) merged map of input maps map1 and map2
  */
  RCP<Epetra_Map> MergeMap(const Epetra_Map& map1,
                           const Epetra_Map& map2);

  /*!
  \brief merges two given Epetra_Maps

  merges input map1 and input map2 (given as RCP), both of which
  have to be unique, but may be overlapping, to a new map and returns
  RCP to it. The case that one or both input RCPs are null is
  detected and handled appropriately.

  \param map1      : one map to be merged
  \param map2      : the other map to be merged
  \return the (sorted) merged map of input maps map1 and map2
  */
  RCP<Epetra_Map> MergeMap(const RCP<Epetra_Map>& map1,
                           const RCP<Epetra_Map>& map2);

  /*!
  \brief split a vector into 2 non-overlapping pieces
  */
  bool SplitVector(const Epetra_Vector& x,
                   const Epetra_Map& x1map,
                   Epetra_Vector*&   x1,
                   const Epetra_Map& x2map,
                   Epetra_Vector*&   x2);

  /*!
  \brief Gather information of type vector<T> on a subset of processors

  This template gathers information provided in sdata on a subset of
  processors tprocs, where the length of the array tprocs is ntargetprocs.
  The redistributed data is returned in rdata which has appropiate size
  on output (size of rdata on input is arbitrary). ntargetprocs can be
  one to reduce data to one proc, it also can be equal to the total number
  of processors to make sdata redundant on all procs.

  \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
        except for that the Epetra version demands the data to be of constant
        size over all procs which this method does not require!

  \param sdata (in) : Information to be gathered on tprocs.
                      Length of sdata can be different on every proc.
  \param rdata (out): Information from sdata gathered on a subset of procs.
                      size of rdata on input is arbitrary, it is exact on output.
  \param ntargetprocs (in): length of tprocs
  \param tprocs (in): vector of procs ids the information in sdata shall be
                      gathered on.
  \param comm (in):   communicator to be used.


  */
  template<typename T> void Gather(vector<T>&         sdata,
                                   vector<T>&         rdata,
                                   const int          ntargetprocs,
                                   const int*         tprocs,
                                   const Epetra_Comm& comm)
  {
    const int myrank  = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc==1)
    {
      rdata = sdata;
      return; // nothing to do in serial
    }
    // build a map of data
    map<int,vector<T> > datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc,1,&myrank,0,comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i=0; i<ntargetprocs; ++i)
      if (tprocs[i]==myrank)
      {
        iamtarget = true;
        break;
      }
    vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i=0; i<numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1,tnummyelements,&targetvec[0],0,comm);
    // build an exporter and export data
    DRT::Exporter exporter(source,target,comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    int count=0;
    map<int,vector<int> >::iterator curr;
    for (curr=datamap.begin(); curr != datamap.end(); ++curr)
    {
      vector<T>& current = curr->second;
      const int size = (int)current.size();
      rdata.resize((int)rdata.size()+size);
      for (int i=0; i<size; ++i)
      rdata[count+i] = current[i];
      count += size;
    }
    return;
  }

  /// Create an allreduced vector of gids from the given Epetra_Map
  /*!
    We have nodes and elements with unique but otherwise arbitrary
    global ids. On rare occations, however, we need to allreduce a
    particular map to one or more processors. This is a building block
    for such occations. We allreduce the gids of the given Epetra_Map
    into a vector ordered by processor number.

    It is assumed that the given map does not overlap!

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param rredundant (o) redundant vector of global ids
    \param emap (i) distributed Epetra_Map

    \author u.kue
    \date 05/07
   */
  void AllreduceEMap(vector<int>& rredundant, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /*!
    We have nodes and elements with unique but otherwise arbitrary
    global ids. But unfortunately we need an allreduced vector of dof
    numbers during the dof assignment phase. In order to use such a
    vector we need to map from global ids to vector indexes. Here we
    provide that map.

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param idxmap (o) map from global ids to (redundant) vector indexes
    \param emap (i) distributed Epetra_Map

    \author u.kue
    \date 05/07
   */
  void AllreduceEMap(map<int,int>& idxmap, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /// on a distinct processor, all other procs create empty maps instead.
  /*!
    This method is currently used within the parallel post_drt_ensight
    filter in order to import all values stored in a distributed Epetra_Vector
    to processor 0 for writing them into file.

    \note see also documentation for the usual AllreduceEMap methods

    \param emap (i) any distributed Epetra_Map
    \param pid (i)  processor id where you want to have the allreduced map
    				exclusively
    \author gjb
    \date 11/07
   */
  RCP<Epetra_Map> AllreduceEMap(const Epetra_Map& emap, const int pid);

  /// find position of my map elements in a consecutive vector
  /*!
    The idea is to put the entries of a given map into a redundant
    vector, ordered by processor number. The map is assumed to be
    nonoverlapping. Here we figure out the index of our first entry in
    that vector.

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param nummyelements (i) number of elements on this proc
    \param comm (i) communicator

    \return vector position of first entry on each processor

    \author u.kue
    \date 05/07
   */
  int FindMyPos(int nummyelements, const Epetra_Comm& comm);

  /// Communication between all pairs of processes, with distinct data for each.
  /*!
    Sends a different vector<int> to each processes. The size of each vector may
    be different, zero-length vectors are allowed.
    Communication is implemented with the MPI function MPI_Alltoallv.

    \param comm (i) communicator
    \param send (i) vector of length comm.NumProc(), j-th element to be send to
    j-th processor.
    \param recv (o) vector of length comm.NumProc(), j-th element received from
    j-th processor.

    \author h.kue
    \date 09/07
   */
  void AllToAllCommunication( const Epetra_Comm& comm,
                              const vector< vector<int> >& send, vector< vector<int> >& recv );

} // namespace LINALG


#endif  // #ifndef LINALG_UTILS_H
#endif  // #ifdef CCADISCRET
