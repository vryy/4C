// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_IO_HPP
#define FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_IO_HPP

#include "4C_config.hpp"

#include "4C_linalg_graph.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"

#include <filesystem>

FOUR_C_NAMESPACE_OPEN

// Forward declarations
namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Core::LinAlg
{
  //! Print content of @p sparsematrix in Matlab format to file @p filename. Create new file or
  //! overwrite existing one if @p newfile is true
  FOUR_C_API(FOUR_C_CORE)
  void print_matrix_in_matlab_format(const std::string& filename,
      const Core::LinAlg::SparseMatrix& sparsematrix, const bool newfile = true);

  //! Print content of @p blockmatrix in Matlab format to file @p filename
  FOUR_C_API(FOUR_C_CORE)
  void print_block_matrix_in_matlab_format(
      const std::string& filename, const BlockSparseMatrixBase& blockmatrix);

  //! Print content of @p vector in Matlab format to file @p filename. Create new file or overwrite
  //! existing one if @p newfile is true
  FOUR_C_API(FOUR_C_CORE)
  void print_vector_in_matlab_format(const std::string& filename,
      const Core::LinAlg::Vector<double>& vector, const bool newfile = true);

  //! Print content of @p map in Matlab format to file @p filename. Create new file or overwrite
  //! existing one if @p newfile is true
  FOUR_C_API(FOUR_C_CORE)
  void print_map_in_matlab_format(
      const std::string& filename, const Core::LinAlg::Map& map, const bool newfile = true);

  /**
   * @brief Read a Matrix Market file and construct a distributed sparse matrix.
   *
   * This function reads a matrix stored in the standard **Matrix Market (.mtx)**
   * format from disk and converts it into a Core::LinAlg::SparseMatrix.
   * The matrix is distributed across MPI processes according to the communicator provided.
   *
   * The function assumes the Matrix Market file represents a sparse matrix and that all MPI ranks
   * in the communicator collectively participate in the read and construction process.
   *
   * @param filename Path to the Matrix Market (.mtx) file to be read.
   * @param comm MPI communicator over which the sparse matrix will be distributed.
   *
   * @return A Core::LinAlg::SparseMatrix containing the data read from the file, distributed
   *         across the given MPI communicator.
   */
  FOUR_C_API(FOUR_C_CORE)
  Core::LinAlg::SparseMatrix read_matrix_market_file_as_sparse_matrix(
      const std::string& filename, MPI_Comm comm);


  /*! \brief Read a multi vector from a Matrix Market format (ASCII)
   *
   * ### Matrix Market format (ASCII)
   * The file can also be created from Python using `scipy.io.mmwrite`, for example:
   * @code{.py}
   * from scipy.io import mmwrite
   * mmwrite("matrix.mtx", projmatrix, precision=10)
   * @endcode
   *
   * @param filename  Path to the input file containing the matrix.
   * @param map Map describing the distribution of the multi vector.
   * @return MultiVector containing the data read from the Matrix Market file, distributed according
   * to \p map.
   */
  FOUR_C_API(FOUR_C_CORE)
  Core::LinAlg::MultiVector<double> read_matrix_market_file_as_multi_vector(
      const std::filesystem::path& filename, const Core::LinAlg::Map& map);
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
