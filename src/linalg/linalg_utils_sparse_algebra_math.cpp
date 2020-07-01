/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic mathematical methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_sparse_algebra_math.H"
#include "../drt_lib/drt_dserror.H"
#include <EpetraExt_Transpose_RowMatrix.h>
#include <EpetraExt_MatrixMatrix.h>

namespace LINALG
{
  namespace
  {
    /*----------------------------------------------------------------------*
     |  Internal optimized matrix addition                 kronbichler 11/15|
     |  B += A(transposed)*scalarA                                          |
     |  Return value: number of local rows in A added successfully          |
     |             (in case B must be uncompleted, this must be remembered) |
     *----------------------------------------------------------------------*/
    int DoAdd(const Epetra_CrsMatrix& A, const double scalarA, Epetra_CrsMatrix& B,
        const double scalarB, const int startRow = 0)
    {
      if (!A.Filled()) dserror("Internal error, matrix A must have called FillComplete()");

      const int NumMyRows = A.NumMyRows();

      // Case 1 where matrix B is filled. In that case, we can attempt to add in local indices,
      // much faster than the global indices... :-)
      if (B.Filled())
      {
        if (startRow != 0) dserror("Internal error. Not implemented.");

        // step 1: get the indexing from A to B in a random-access array
        std::vector<int> AcolToBcol(A.ColMap().NumMyElements());
        for (int i = 0; i < A.ColMap().NumMyElements(); ++i)
          AcolToBcol[i] = B.ColMap().LID(A.ColMap().GID(i));

        std::vector<int> indicesInB(A.MaxNumEntries());

        // step 2: loop over all local rows in matrix A and attempt the addition in local index
        // space
        for (int i = 0; i < NumMyRows; ++i)
        {
          const int myRowB = B.RowMap().LID(A.RowMap().GID(i));
          if (myRowB == -1)
            dserror(
                "LINALG::Add: The row map of matrix B must be a superset of the row map of Matrix "
                "A.");

          // extract views of both the row in A and in B
          double *valuesA = 0, *valuesB = 0;
          int *indicesA = 0, *indicesB = 0;
          int NumEntriesA = -1, NumEntriesB = -1;
          A.ExtractMyRowView(i, NumEntriesA, valuesA, indicesA);
          B.ExtractMyRowView(myRowB, NumEntriesB, valuesB, indicesB);

          // check if we can identify all indices from matrix A in matrix B. quite a lot of
          // indirect addressing in here, but these are all local operations and thus pretty fast
          bool failed = false;
          indicesInB.clear();
          for (int jA = 0, jB = 0; jA < NumEntriesA; ++jA)
          {
            const int col = AcolToBcol[indicesA[jA]];
            if (col == -1)
            {
              failed = true;
              break;
            }
            while (jB < NumEntriesB && indicesB[jB] < col) ++jB;

            // did not find index in linear search (re-indexing from A.ColMap() to B.ColMap()
            // might pass through the indices differently), try binary search
            if (indicesB[jB] != col)
              jB = std::lower_bound(&indicesB[0], &indicesB[0] + NumEntriesB, col) - &indicesB[0];

            // not found, sparsity pattern of B does not contain the index from A -> terminate
            if (indicesB[jB] != col)
            {
              failed = true;
              break;
            }

            indicesInB.push_back(jB);
          }

          if (failed)
            return i;
          else
          {
            for (int j = 0; j < NumEntriesA; ++j) valuesB[indicesInB[j]] += scalarA * valuesA[j];
          }
        }

        return NumMyRows;
      }

      // case 2 where B is not filled -> good old slow addition in global indices

      // Loop over Aprime's rows and sum into
      std::vector<int> Indices(A.MaxNumEntries());
      std::vector<double> Values(A.MaxNumEntries());

      // Continue with i from the previous attempt
      for (int i = startRow; i < NumMyRows; ++i)
      {
        const int Row = A.GRID(i);
        int NumEntries = 0;
        int ierr = A.ExtractGlobalRowCopy(Row, Values.size(), NumEntries, &Values[0], &Indices[0]);
        if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", ierr);
        if (scalarA != 1.0)
          for (int j = 0; j < NumEntries; ++j) Values[j] *= scalarA;
        for (int j = 0; j < NumEntries; ++j)
        {
          int err = B.SumIntoGlobalValues(Row, 1, &Values[j], &Indices[j]);
          if (err < 0 || err == 2) err = B.InsertGlobalValues(Row, 1, &Values[j], &Indices[j]);
          if (err < 0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
        }
      }

      return NumMyRows;
    }
  }  // namespace
}  // namespace LINALG



/*----------------------------------------------------------------------*
 |  Add a sparse matrix to another                     kronbichler 11/15|
 |  B = B*scalarB + A(transposed)*scalarA                               |
 *----------------------------------------------------------------------*/
void LINALG::Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
    LINALG::SparseMatrixBase& B, const double scalarB)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");

  Epetra_CrsMatrix* Aprime = NULL;
  Teuchos::RCP<EpetraExt::RowMatrix_Transpose> Atrans = Teuchos::null;
  if (transposeA)
  {
    // Atrans = Teuchos::rcp(new EpetraExt::RowMatrix_Transpose(false,NULL,false));
    Atrans = Teuchos::rcp(new EpetraExt::RowMatrix_Transpose());
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
  {
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
  }

  if (scalarB == 0.)
    B.PutScalar(0.0);
  else if (scalarB != 1.0)
    B.Scale(scalarB);

  int rowsAdded = DoAdd(*Aprime, scalarA, *B.EpetraMatrix(), scalarB);
  int localSuccess = rowsAdded == Aprime->RowMap().NumMyElements();
  int globalSuccess = 0;
  B.Comm().MinAll(&localSuccess, &globalSuccess, 1);
  if (!globalSuccess)
  {
    if (!B.Filled()) dserror("Unexpected state of B (expected: B not filled, got: B filled)");

    // not successful -> matrix structure must be un-completed to be able to add new
    // indices.
    B.UnComplete();
    DoAdd(*Aprime, scalarA, *B.EpetraMatrix(), scalarB, rowsAdded);
    B.Complete();
  }
}



/*----------------------------------------------------------------------*
 |  Add a sparse matrix to another                           mwgee 12/06|
 |  B = B*scalarB + A(transposed)*scalarA                               |
 *----------------------------------------------------------------------*/
void LINALG::Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
    Epetra_CrsMatrix& B, const double scalarB)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");

  Epetra_CrsMatrix* Aprime = NULL;
  Teuchos::RCP<EpetraExt::RowMatrix_Transpose> Atrans = Teuchos::null;
  if (transposeA)
  {
    // Atrans = Teuchos::rcp(new EpetraExt::RowMatrix_Transpose(false,NULL,false));
    Atrans = Teuchos::rcp(new EpetraExt::RowMatrix_Transpose());
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
  {
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
  }

  if (scalarB == 0.)
    B.PutScalar(0.0);
  else if (scalarB != 1.0)
    B.Scale(scalarB);

  int rowsAdded = DoAdd(*Aprime, scalarA, B, scalarB);
  if (rowsAdded != Aprime->RowMap().NumMyElements())
    dserror("LINALG::Add: Could not add all entries from A into B in row %d",
        Aprime->RowMap().GID(rowsAdded));
}

/*----------------------------------------------------------------------*
 | Transpose matrix A                                         popp 02/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> LINALG::Transpose(const Epetra_CrsMatrix& A)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");

  Teuchos::RCP<EpetraExt::RowMatrix_Transpose> Atrans =
      Teuchos::rcp(new EpetraExt::RowMatrix_Transpose(/*false,NULL,false*/));
  Epetra_CrsMatrix* Aprime =
      &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));

  return Teuchos::rcp(new Epetra_CrsMatrix(*Aprime));
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 01/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> LINALG::Multiply(
    const Epetra_CrsMatrix& A, bool transA, const Epetra_CrsMatrix& B, bool transB, bool complete)
{
  /* ATTENTION (Q1/2013 and later)
   *
   * Be careful with LINALG::Multiply when using Trilinos Q1/2013.
   * The new EpetraExt MMM routines are very fast, but not well tested and
   * sometimes crash. To be on the safe side consider to use LINALG::MLMultiply
   * which is based on ML and well tested over several years.
   *
   * For improving speed of MM operations (even beyond MLMultiply) we probably
   * should rethink the design of our Multiply routines. As a starting point we
   * should have a look at the new variants of the EpetraExt routines. See also
   * the following personal communication with Chris Siefert
   *
   * "There are basically 3 different matrix-matrix-multiplies in the new EpetraExt MMM:
   *
   * 1) Re-use existing C.  I think this works, but it isn't well tested.
   * 2) Start with a clean C and FillComplete.  This is the one we're using in MueLu that works
   * fine. 3) Any other case.  This is the one you found a bug in.
   *
   * I'll try to track the bug in #3 down, but you really should be calling #2 (or #1) if at all
   * possible.
   *
   * -Chris"
   *
   */

  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  // do a very coarse guess of nonzeros per row (horrible memory consumption!)
  // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
  // a first guess for the bandwidth of C leading to much less memory allocation:
  const int guessnpr = std::max(A.MaxNumEntries(), B.MaxNumEntries());

  // create resultmatrix with correct rowmap
  Epetra_CrsMatrix* C = NULL;
  if (!transA)
    C = new Epetra_CrsMatrix(::Copy, A.OperatorRangeMap(), guessnpr, false);
  else
    C = new Epetra_CrsMatrix(::Copy, A.OperatorDomainMap(), guessnpr, false);

  int err = EpetraExt::MatrixMatrix::Multiply(A, transA, B, transB, *C, complete);
  if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = &d", err);

  return Teuchos::rcp(C);
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B*C                                   mwgee 02/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> LINALG::Multiply(const Epetra_CrsMatrix& A, bool transA,
    const Epetra_CrsMatrix& B, bool transB, const Epetra_CrsMatrix& C, bool transC, bool complete)
{
  Teuchos::RCP<Epetra_CrsMatrix> tmp = LINALG::Multiply(B, transB, C, transC, true);
  return LINALG::Multiply(A, transA, *tmp, false, complete);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SymmetriseMatrix(Epetra_SerialDenseMatrix& A)
{
  const int n = A.N();
  if (n != A.M()) dserror("Cannot symmetrize non-square matrix");
  // do not make deep copy of A, matrix addition and full scaling just to sym it
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
    {
      const double aver = 0.5 * (A(i, j) + A(j, i));
      A(i, j) = A(j, i) = aver;
    }
  return;
}
