/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of sparse matrix printing methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/
#include <fstream>

#include "linalg_blocksparsematrix.H"
#include "linalg_utils_sparse_algebra_print.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <MueLu_UseDefaultTypes.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintSparsityToPostscript(const Epetra_RowMatrix& A) { Ifpack_PrintSparsity(A); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintMatrixInMatlabFormat(
    const std::string& fname, const Epetra_CrsMatrix& A, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  const int my_PID = A.RowMap().Comm().MyPID();
  const int num_proc = A.RowMap().Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    if (my_PID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      int num_my_rows = A.NumMyRows();
      int max_num_inidces = A.MaxNumEntries();
      std::vector<int> indices(max_num_inidces);
      std::vector<double> values(max_num_inidces);
      int num_indices;

      for (int row_lid = 0; row_lid < num_my_rows; row_lid++)
      {
        int row_gid = A.GRID(row_lid);  // Get global row number
        A.ExtractGlobalRowCopy(row_gid, max_num_inidces, num_indices, &values[0], &indices[0]);

        for (int col_lid = 0; col_lid < num_indices; col_lid++)
        {
          os << std::setw(10) << row_gid + 1;           // increase index by one for matlab
          os << std::setw(10) << indices[col_lid] + 1;  // increase index by one for matlab
          os << std::setw(30) << std::setprecision(16) << std::scientific << values[col_lid];
          os << std::endl;
        }
      }

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
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintBlockMatrixInMatlabFormat(
    const std::string& fname, const BlockSparseMatrixBase& A)
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
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintVectorInMatlabFormat(
    const std::string& fname, const Epetra_Vector& V, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  const int my_PID = V.Map().Comm().MyPID();
  const int num_proc = V.Map().Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < num_proc; iproc++)  // loop over all processors
  {
    if (my_PID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      const int num_my_elements = V.Map().NumMyElements();
      const int max_element_size = V.Map().MaxElementSize();
      int* my_global_elements = V.Map().MyGlobalElements();
      double** A_Pointers = V.Pointers();

      for (int lid = 0; lid < num_my_elements; lid++)
      {
        if (max_element_size == 1)
        {
          os << std::setw(10) << my_global_elements[lid] + 1;  // add +1 for Matlab convention

          os << std::setw(30) << std::setprecision(16) << A_Pointers[0][lid]
             << std::endl;  // print out values of 1. vector (no Multi_Vector)
        }
        else
        {
          int* first_point_in_element_list = V.Map().FirstPointInElementList();
          for (int ele_lid = 0; ele_lid < V.Map().ElementSize(lid); ele_lid++)
          {
            os << std::setw(10) << my_global_elements[lid] << "/" << std::setw(10) << ele_lid;

            os << std::setw(30) << std::setprecision(16)
               << A_Pointers[0][first_point_in_element_list[lid] + ele_lid]
               << std::endl;  // print out values of 1. vector (no Multi_Vector)
          }
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
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintMapInMatlabFormat(
    const std::string& fname, const Epetra_Map& map, const bool newfile)
{
  // The following source code has been adapted from the Print() method
  // of the Epetra_CrsMatrix class (see "Epetra_CrsMatrix.cpp").

  int my_PID = map.Comm().MyPID();
  int num_proc = map.Comm().NumProc();

  std::ofstream os;

  for (int iproc = 0; iproc < num_proc; iproc++)  // loop over all processors
  {
    if (my_PID == iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(), std::fstream::trunc);
      else
        os.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      const int num_my_elements = map.NumMyElements();
      const int max_element_size = map.MaxElementSize();
      int* my_global_elements = map.MyGlobalElements();

      for (int lid = 0; lid < num_my_elements; lid++)
      {
        for (int ele_lid = 0; ele_lid < map.ElementSize(lid); ele_lid++)
        {
          if (max_element_size == 1)
          {
            os << std::setw(10) << my_global_elements[lid] + 1;
          }
          else
          {
            os << std::setw(10) << my_global_elements[lid] + 1 << "/" << std::setw(10) << ele_lid;
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::WriteEpetraCrsMatrixAsXpetra(
    const std::string& filename, Teuchos::RCP<Epetra_CrsMatrix> matrix)
{
#include <Xpetra_UseShortNames.hpp>  // Include in scope to avoid clash with namespace IO
  using Teuchos::rcp;
  using Teuchos::RCP;

  RCP<CrsMatrix> ACrs = rcp(new EpetraCrsMatrix(matrix));
  RCP<CrsMatrixWrap> ACrsWrap = rcp(new CrsMatrixWrap(ACrs));
  RCP<Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(ACrsWrap);

  Xpetra::IO<double, int, int, Node>::Write(filename, *A);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::WriteEpetraMultiVectorAsXpetra(
    const std::string& filename, Teuchos::RCP<Epetra_MultiVector> vec)
{
  Xpetra::IO<double, int, int, Node>::Write(filename, *Xpetra::toXpetra<int, Node>(vec));
}
