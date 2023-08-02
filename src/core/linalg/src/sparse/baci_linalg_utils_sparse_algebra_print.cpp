/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of sparse matrix printing methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/
#include "baci_linalg_utils_sparse_algebra_print.H"

#include "baci_linalg_blocksparsematrix.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <MueLu_UseDefaultTypes.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include <fstream>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::PrintSparsityToPostscript(const Epetra_RowMatrix& A) { Ifpack_PrintSparsity(A); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::PrintMatrixInMatlabFormat(
    const std::string& filename, const Epetra_CrsMatrix& sparsematrix, const bool newfile)
{
  const auto& comm = sparsematrix.Comm();

  const int my_PID = comm.MyPID();
  const int num_proc = comm.NumProc();

  // loop over all procs and send row data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_rows_iproc = sparsematrix.NumMyRows();
    comm.Broadcast(&num_rows_iproc, 1, iproc);

    for (int row_lid_iproc = 0; row_lid_iproc < num_rows_iproc; ++row_lid_iproc)
    {
      // get gid of this row and communicate to all procs
      int row_gid_iproc = iproc == my_PID ? sparsematrix.GRID(row_lid_iproc) : 0;
      comm.Broadcast(&row_gid_iproc, 1, iproc);

      // get indices and values of this row and communicate to all procs
      int num_indices_iproc;
      std::vector<int> indices_iproc;
      std::vector<double> values_iproc;
      if (iproc == my_PID)
      {
        const int max_num_inidces = sparsematrix.MaxNumEntries();
        indices_iproc.resize(max_num_inidces);
        values_iproc.resize(max_num_inidces);

        sparsematrix.ExtractGlobalRowCopy(row_gid_iproc, max_num_inidces, num_indices_iproc,
            values_iproc.data(), indices_iproc.data());
      }
      comm.Broadcast(&num_indices_iproc, 1, iproc);
      values_iproc.resize(num_indices_iproc);
      indices_iproc.resize(num_indices_iproc);

      comm.Broadcast(values_iproc.data(), num_indices_iproc, iproc);
      comm.Broadcast(indices_iproc.data(), num_indices_iproc, iproc);

      if (my_PID == 0)
      {
        std::ofstream os;
        // create new file
        if (newfile and iproc == 0 and row_lid_iproc == 0)
          os.open(filename.c_str(), std::fstream::trunc);
        else
          os.open(filename.c_str(), std::fstream::ate | std::fstream::app);

        for (int col_idx = 0; col_idx < num_indices_iproc; col_idx++)
        {
          os << std::setw(10) << row_gid_iproc + 1;           // increase index by one for matlab
          os << std::setw(10) << indices_iproc[col_idx] + 1;  // increase index by one for matlab
          os << std::setw(30) << std::setprecision(16) << std::scientific << values_iproc[col_idx];
          os << std::endl;
        }
        os << std::flush;
      }
    }
    // wait, until proc 0 has written
    comm.Barrier();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::PrintBlockMatrixInMatlabFormat(
    const std::string& filename, const BlockSparseMatrixBase& blockmatrix)
{
  for (int row = 0; row < blockmatrix.Rows(); row++)
  {
    for (int col = 0; col < blockmatrix.Cols(); col++)
    {
      const auto& sparsematrix = blockmatrix.Matrix(row, col);
      CORE::LINALG::PrintMatrixInMatlabFormat(
          filename, *(sparsematrix.EpetraMatrix()), ((row == 0) && (col == 0)));
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::PrintVectorInMatlabFormat(
    const std::string& filename, const Epetra_Vector& vector, const bool newfile)
{
  const auto& comm = vector.Comm();

  const int my_PID = comm.MyPID();
  const int num_proc = comm.NumProc();

  // loop over all procs and send data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_elements_iproc = vector.Map().NumMyElements();
    int max_element_size_iproc = vector.Map().MaxElementSize();

    comm.Broadcast(&num_elements_iproc, 1, iproc);
    comm.Broadcast(&max_element_size_iproc, 1, iproc);

    std::vector<int> global_elements_iproc(num_elements_iproc);
    std::vector<double> values_iproc(num_elements_iproc);
    std::vector<int> first_point_in_element_list_iproc(num_elements_iproc);

    if (iproc == my_PID)
    {
      for (int i = 0; i < num_elements_iproc; ++i)
      {
        global_elements_iproc[i] = vector.Map().MyGlobalElements()[i];
        values_iproc[i] = vector.Values()[i];
        first_point_in_element_list_iproc[i] = vector.Map().FirstPointInElementList()[i];
      }
    }

    comm.Broadcast(global_elements_iproc.data(), num_elements_iproc, iproc);
    comm.Broadcast(values_iproc.data(), num_elements_iproc, iproc);
    comm.Broadcast(first_point_in_element_list_iproc.data(), num_elements_iproc, iproc);

    if (my_PID == 0)
    {
      std::ofstream os;
      if (newfile and iproc == 0)
        os.open(filename.c_str(), std::fstream::trunc);
      else
        os.open(filename.c_str(), std::fstream::ate | std::fstream::app);

      for (int lid = 0; lid < num_elements_iproc; lid++)
      {
        if (max_element_size_iproc == 1)
        {
          os << std::setw(10) << global_elements_iproc[lid] + 1;  // add +1 for Matlab convention

          os << std::setw(30) << std::setprecision(16) << values_iproc[lid]
             << std::endl;  // print out values of 1. vector (no Multi_Vector)
        }
        else
        {
          for (int ele_lid = 0; ele_lid < vector.Map().ElementSize(lid); ele_lid++)
          {
            os << std::setw(10) << global_elements_iproc[lid] << "/" << std::setw(10) << ele_lid;

            os << std::setw(30) << std::setprecision(16)
               << values_iproc[first_point_in_element_list_iproc[lid] + ele_lid]
               << std::endl;  // print out values of 1. vector (no Multi_Vector)
          }
        }
        os << std::flush;
      }
    }
    // wait, until proc 0 has written
    comm.Barrier();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::PrintMapInMatlabFormat(
    const std::string& filename, const Epetra_Map& map, const bool newfile)
{
  const auto& comm = map.Comm();

  const int my_PID = comm.MyPID();
  const int num_proc = comm.NumProc();

  // loop over all procs and send data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_elements_iproc = map.NumMyElements();
    int max_element_size_iproc = map.MaxElementSize();

    comm.Broadcast(&num_elements_iproc, 1, iproc);
    comm.Broadcast(&max_element_size_iproc, 1, iproc);

    std::vector<int> global_elements_iproc(num_elements_iproc);

    if (iproc == my_PID)
    {
      for (int i = 0; i < num_elements_iproc; ++i)
        global_elements_iproc[i] = map.MyGlobalElements()[i];
    }
    comm.Broadcast(global_elements_iproc.data(), num_elements_iproc, iproc);

    if (my_PID == 0)
    {
      std::ofstream os;
      if (newfile and iproc == 0)
        os.open(filename.c_str(), std::fstream::trunc);
      else
        os.open(filename.c_str(), std::fstream::ate | std::fstream::app);

      for (int lid = 0; lid < num_elements_iproc; lid++)
      {
        for (int ele_lid = 0; ele_lid < map.ElementSize(lid); ele_lid++)
        {
          if (max_element_size_iproc == 1)
          {
            os << std::setw(10) << global_elements_iproc[lid] + 1;
          }
          else
          {
            os << std::setw(10) << global_elements_iproc[lid] + 1 << "/" << std::setw(10)
               << ele_lid;
          }
          os << std::endl;
        }
      }
      os << std::flush;
    }
    // wait, until proc 0 has written
    comm.Barrier();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::WriteEpetraCrsMatrixAsXpetra(
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
void CORE::LINALG::WriteEpetraMultiVectorAsXpetra(
    const std::string& filename, Teuchos::RCP<Epetra_MultiVector> vec)
{
  Xpetra::IO<double, int, int, Node>::Write(filename, *Xpetra::toXpetra<int, Node>(vec));
}
