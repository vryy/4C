// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_utils_sparse_algebra_io.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_linalg_utils_exceptions.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::print_matrix_in_matlab_format(
    const std::string& filename, const Core::LinAlg::SparseMatrix& sparsematrix, const bool newfile)
{
  const auto& comm = sparsematrix.get_comm();

  const int my_PID = Core::Communication::my_mpi_rank(comm);
  const int num_proc = Core::Communication::num_mpi_ranks(comm);

  // loop over all procs and send row data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_rows_iproc = sparsematrix.num_my_rows();
    Core::Communication::broadcast(&num_rows_iproc, 1, iproc, comm);

    for (int row_lid_iproc = 0; row_lid_iproc < num_rows_iproc; ++row_lid_iproc)
    {
      // get gid of this row and communicate to all procs
      int row_gid_iproc = iproc == my_PID ? sparsematrix.global_row_index(row_lid_iproc) : 0;
      Core::Communication::broadcast(&row_gid_iproc, 1, iproc, comm);

      // get indices and values of this row and communicate to all procs
      int num_indices_iproc;
      std::vector<int> indices_iproc;
      std::vector<double> values_iproc;
      if (iproc == my_PID)
      {
        const int max_num_inidces = sparsematrix.max_num_entries();
        indices_iproc.resize(max_num_inidces);
        values_iproc.resize(max_num_inidces);

        sparsematrix.extract_global_row_copy(row_gid_iproc, max_num_inidces, num_indices_iproc,
            values_iproc.data(), indices_iproc.data());
      }
      Core::Communication::broadcast(&num_indices_iproc, 1, iproc, comm);
      values_iproc.resize(num_indices_iproc);
      indices_iproc.resize(num_indices_iproc);

      Core::Communication::broadcast(values_iproc.data(), num_indices_iproc, iproc, comm);
      Core::Communication::broadcast(indices_iproc.data(), num_indices_iproc, iproc, comm);

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
    Core::Communication::barrier(comm);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::print_block_matrix_in_matlab_format(
    const std::string& filename, const BlockSparseMatrixBase& blockmatrix)
{
  for (int row = 0; row < blockmatrix.rows(); row++)
  {
    for (int col = 0; col < blockmatrix.cols(); col++)
    {
      const auto& sparsematrix = blockmatrix.matrix(row, col);
      print_matrix_in_matlab_format(filename, sparsematrix, ((row == 0) && (col == 0)));
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::print_vector_in_matlab_format(
    const std::string& filename, const Core::LinAlg::Vector<double>& vector, const bool newfile)
{
  const auto& comm = vector.get_comm();

  const int my_PID = Core::Communication::my_mpi_rank(comm);
  const int num_proc = Core::Communication::num_mpi_ranks(comm);

  // loop over all procs and send data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_elements_iproc = vector.get_map().num_my_elements();
    int max_element_size_iproc = vector.get_map().max_element_size();

    Core::Communication::broadcast(&num_elements_iproc, 1, iproc, comm);
    Core::Communication::broadcast(&max_element_size_iproc, 1, iproc, comm);

    std::vector<int> global_elements_iproc(num_elements_iproc);
    std::vector<double> values_iproc(num_elements_iproc);
    std::vector<int> first_point_in_element_list_iproc(num_elements_iproc);

    if (iproc == my_PID)
    {
      for (int i = 0; i < num_elements_iproc; ++i)
      {
        global_elements_iproc[i] = vector.get_map().my_global_elements()[i];
        values_iproc[i] = vector.get_values()[i];
        first_point_in_element_list_iproc[i] =
            vector.get_map().get_epetra_block_map().FirstPointInElementList()[i];
      }
    }

    Core::Communication::broadcast(global_elements_iproc.data(), num_elements_iproc, iproc, comm);
    Core::Communication::broadcast(values_iproc.data(), num_elements_iproc, iproc, comm);
    Core::Communication::broadcast(
        first_point_in_element_list_iproc.data(), num_elements_iproc, iproc, comm);

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
          for (int ele_lid = 0; ele_lid < vector.get_map().element_size(lid); ele_lid++)
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
    Core::Communication::barrier(comm);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::print_map_in_matlab_format(
    const std::string& filename, const Core::LinAlg::Map& map, const bool newfile)
{
  const auto& comm = map.get_comm();

  const int my_PID = Core::Communication::my_mpi_rank(comm);
  const int num_proc = Core::Communication::num_mpi_ranks(comm);

  // loop over all procs and send data to proc 0
  for (int iproc = 0; iproc < num_proc; iproc++)
  {
    int num_elements_iproc = map.num_my_elements();
    int max_element_size_iproc = map.max_element_size();

    Core::Communication::broadcast(&num_elements_iproc, 1, iproc, comm);
    Core::Communication::broadcast(&max_element_size_iproc, 1, iproc, comm);

    std::vector<int> global_elements_iproc(num_elements_iproc);

    if (iproc == my_PID)
    {
      for (int i = 0; i < num_elements_iproc; ++i)
        global_elements_iproc[i] = map.my_global_elements()[i];
    }
    Core::Communication::broadcast(global_elements_iproc.data(), num_elements_iproc, iproc, comm);

    if (my_PID == 0)
    {
      std::ofstream os;
      if (newfile and iproc == 0)
        os.open(filename.c_str(), std::fstream::trunc);
      else
        os.open(filename.c_str(), std::fstream::ate | std::fstream::app);

      for (int lid = 0; lid < num_elements_iproc; lid++)
      {
        for (int ele_lid = 0; ele_lid < map.element_size(lid); ele_lid++)
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
    Core::Communication::barrier(comm);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::SparseMatrix Core::LinAlg::read_matrix_market_file_as_sparse_matrix(
    const std::string& filename, MPI_Comm comm)
{
  Epetra_CrsMatrix* crs_marix;

  ASSERT_EPETRA_CALL(EpetraExt::MatrixMarketFileToCrsMatrix(
      filename.c_str(), Core::Communication::as_epetra_comm(comm), crs_marix));

  return SparseMatrix(std::shared_ptr<Epetra_CrsMatrix>(crs_marix), DataAccess::Share);
}

Core::LinAlg::MultiVector<double> Core::LinAlg::read_matrix_market_file_as_multi_vector(
    const std::filesystem::path& filename, const Core::LinAlg::Map& map)
{
  // Epetra Matrix market file reader does not support complex maps.
  // Hence, read on rank 0 and distribute to the target ranks
  int rank = Core::Communication::my_mpi_rank(map.get_comm());
  Map rank0_map{
      map.num_global_elements(), rank == 0 ? map.num_global_elements() : 0, 0, map.get_comm()};


  Epetra_MultiVector* multi_vector_raw = nullptr;

  ASSERT_EPETRA_CALL(EpetraExt::MatrixMarketFileToMultiVector(
      filename.string().c_str(), rank0_map.get_epetra_block_map(), multi_vector_raw));
  FOUR_C_ASSERT_ALWAYS(multi_vector_raw != nullptr, "Failed to read matrix market file {}.",
      filename.relative_path().string());

  // make sure to delete the raw pointers from above
  std::unique_ptr<Epetra_MultiVector> auto_deleter(multi_vector_raw);

  MultiVector<double> mv_serial(*multi_vector_raw);

  FOUR_C_ASSERT_ALWAYS(mv_serial.global_length() == map.num_global_elements(),
      "The matrix in the matrix market file has {} rows, but the provided map has {} elements.",
      multi_vector_raw->GlobalLength(), map.num_global_elements());

  MultiVector<double> mv_distributed(map, mv_serial.num_vectors());

  export_to(mv_serial, mv_distributed);

  return mv_distributed;
}

FOUR_C_NAMESPACE_CLOSE
