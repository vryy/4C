/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of sparse matrix printing methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/
#include <fstream>

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_sparse_algebra_print.H"
#include "../drt_lib/drt_dserror.H"
#include <Ifpack_AdditiveSchwarz.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintSparsityToPostscript(const Epetra_RowMatrix& A)
{
  Ifpack_PrintSparsity(A);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintMatrixInMatlabFormat(
    std::string fname, const Epetra_CrsMatrix& A, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  int MyPID = A.RowMap().Comm().MyPID();
  int NumProc = A.RowMap().Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < NumProc; iproc++)
  {
    if (MyPID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      int NumMyRows1 = A.NumMyRows();
      int MaxNumIndices = A.MaxNumEntries();
      int* Indices = new int[MaxNumIndices];
      double* Values = new double[MaxNumIndices];
      int NumIndices;
      int i, j;

      for (i = 0; i < NumMyRows1; i++)
      {
        int Row = A.GRID(i);  // Get global row number
        A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);

        for (j = 0; j < NumIndices; j++)
        {
          os << std::setw(10) << Row + 1;         // increase index by one for matlab
          os << std::setw(10) << Indices[j] + 1;  // increase index by one for matlab
          os << std::setw(30) << std::setprecision(16) << std::scientific << Values[j];
          os << std::endl;
        }
      }

      delete[] Indices;
      delete[] Values;

      os << std::flush;

      // close file
      os.close();
    }
    // Do a few global ops to give I/O a chance to complete
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
  }

  // just to be sure
  if (os.is_open()) os.close();

  // have fun with your Matlab matrix
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintSerialDenseMatrixInMatlabFormat(
    std::string fname, const Epetra_SerialDenseMatrix& A, const bool newfile)
{
  // The following source code has been adapted from the PrintMatrixInMatlabFormat
  // method in order to also print a Epetra_SerialDenseMatrix.

  std::ofstream os;

  // open file for writing
  if (newfile)
    os.open(fname.c_str(), std::fstream::trunc);
  else
    os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

  int NumMyRows = A.RowDim();
  int NumMyColumns = A.ColDim();

  for (int i = 0; i < NumMyRows; i++)
  {
    for (int j = 0; j < NumMyColumns; j++)
    {
      os << std::setw(10) << i + 1;  // increase index by one for matlab
      os << std::setw(10) << j + 1;  // increase index by one for matlab
      os << std::setw(30) << std::setprecision(16) << std::scientific << A(i, j);
      os << std::endl;
    }
  }

  os << std::flush;

  // close file
  os.close();

  // just to be sure
  if (os.is_open()) os.close();

  // have fun with your Matlab matrix
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintBlockMatrixInMatlabFormat(std::string fname, const BlockSparseMatrixBase& A)
{
  // For each sub-matrix of A use the existing printing method
  for (int r = 0; r < A.Rows(); r++)
  {
    for (int c = 0; c < A.Cols(); c++)
    {
      const LINALG::SparseMatrix& M = A.Matrix(r, c);
      LINALG::PrintMatrixInMatlabFormat(fname, *(M.EpetraMatrix()), ((r == 0) && (c == 0)));
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintVectorInMatlabFormat(
    std::string fname, const Epetra_Vector& V, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  int MyPID = V.Map().Comm().MyPID();
  int NumProc = V.Map().Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < NumProc; iproc++)  // loop over all processors
  {
    if (MyPID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      int NumMyElements1 = V.Map().NumMyElements();
      int MaxElementSize1 = V.Map().MaxElementSize();
      int* MyGlobalElements1 = V.Map().MyGlobalElements();
      int* FirstPointInElementList1(NULL);
      if (MaxElementSize1 != 1) FirstPointInElementList1 = V.Map().FirstPointInElementList();
      double** A_Pointers = V.Pointers();

      for (int i = 0; i < NumMyElements1; i++)
      {
        for (int ii = 0; ii < V.Map().ElementSize(i); ii++)
        {
          int iii;
          if (MaxElementSize1 == 1)
          {
            os << std::setw(10) << MyGlobalElements1[i] + 1;  // add +1 for Matlab convention
            iii = i;
          }
          else
          {
            os << std::setw(10) << MyGlobalElements1[i] << "/" << std::setw(10) << ii;
            iii = FirstPointInElementList1[i] + ii;
          }

          os << std::setw(30) << std::setprecision(16)
             << A_Pointers[0][iii];  // print out values of 1. vector (only Epetra_Vector supported,
                                     // no Multi_Vector)
          os << std::endl;
        }
      }
      os << std::flush;
    }
    // close file
    os.close();

    // Do a few global ops to give I/O a chance to complete
    V.Map().Comm().Barrier();
    V.Map().Comm().Barrier();
    V.Map().Comm().Barrier();
  }

  // just to be sure
  if (os.is_open()) os.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintMapInMatlabFormat(std::string fname, const Epetra_Map& map, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  int MyPID = map.Comm().MyPID();
  int NumProc = map.Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < NumProc; iproc++)  // loop over all processors
  {
    if (MyPID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      int NumMyElements1 = map.NumMyElements();
      int MaxElementSize1 = map.MaxElementSize();
      int* MyGlobalElements1 = map.MyGlobalElements();

      for (int i = 0; i < NumMyElements1; i++)
      {
        for (int ii = 0; ii < map.ElementSize(i); ii++)
        {
          if (MaxElementSize1 == 1)
          {
            os << std::setw(10) << MyGlobalElements1[i] + 1;
          }
          else
          {
            os << std::setw(10) << MyGlobalElements1[i] + 1 << "/" << std::setw(10) << ii;
          }
          os << std::endl;
        }
      }
      os << std::flush;
    }
    // close file
    os.close();

    // Do a few global ops to give I/O a chance to complete
    map.Comm().Barrier();
    map.Comm().Barrier();
    map.Comm().Barrier();
  }

  // just to be sure
  if (os.is_open()) os.close();

  return;
}
