// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_UTILS_HPP
#define FOUR_C_COMM_UTILS_HPP


#include "4C_config.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_DefaultMpiComm.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  // forward declaration
  class Communicators;

  /**
   * The known types for nested parallelism.
   */
  enum class NestedParallelismType
  {
    every_group_read_dat_file,
    separate_dat_files,
    no_nested_parallelism
  };

  //! create a local and a global communicator for the problem
  std::shared_ptr<Communicators> create_comm(std::vector<std::string> argv);

  /*! \brief debug routine to compare vectors from different parallel 4C runs
   *
   * You can add Core::Communication::AreDistributedVectorsIdentical in your code which will lead to
   * a comparison of the given vector for different executables and/or configurations. Command for
   * using this feature: \n mpirun -np 1 ./4C -nptype=diffgroup0 input.dat xxx_ser : -np 3
   * ./other-4C -nptype=diffgroup1 other-input.dat xxx_par \n Do not forget to include the header
   * (#include "4C_comm_utils.hpp"), otherwise it won't compile.
   *
   * A further nice option is to compare results from different executables used for
   * running the same simulation.
   *
   * \note You need to add the AreDistributedVectorsIdentical method in both executables at the same
   * position in the code
   *
   * \param communicators (in): communicators containing local and global comm
   * \param vec           (in): vector to compare
   * \param name          (in): user given name for the vector (needs to match within gcomm)
   * \param tol           (in): comparison tolerance for infinity norm
   * \return boolean to indicate if compared vectors are identical
   */
  bool are_distributed_vectors_identical(const Communicators& communicators,
      const Core::LinAlg::MultiVector<double>& vec, const char* name, double tol = 1.0e-14);

  /*! \brief debug routine to compare sparse matrices from different parallel 4C runs
   *
   * You can add Core::Communication::AreDistributedSparseMatricesIdentical in your code which will
   * lead to a comparison of the given sparse matrices for different executables and/or
   * configurations. Command for using this feature: \n mpirun -np 1 ./4C -nptype=diffgroup0
   * input.dat xxx_ser : -np 3 ./other-4C -nptype=diffgroup1 other-input.dat xxx_par \n Do not
   * forget to include the header (#include "4C_comm_utils.hpp"), otherwise it won't compile.
   *
   * A further nice option is to compare results from different executables used for
   * running the same simulation.
   *
   * \note You need to add the AreDistributedSparseMatricesIdentical method in both executables at
   * the same position in the code.
   *
   * \note From Core::LinAlg::SparseOperator to CrsMatrix, just do:
   * std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(yoursparseoperator)->EpetraMatrix()
   *
   * \param communicators (in): communicators containing local and global comm
   * \param matrix        (in): matrix to compare
   * \param name          (in): user given name for the matrix (needs to match within gcomm)
   * \param tol           (in): comparison tolerance for infinity norm
   * \return boolean to indicate if compared vectors are identical
   */
  bool are_distributed_sparse_matrices_identical(const Communicators& communicators,
      Epetra_CrsMatrix& matrix, const char* name, double tol = 1.0e-14);

  //! transform MPI_Comm to Teuchos::Comm, std::shared_ptr version
  template <class Datatype>
  std::shared_ptr<const Teuchos::Comm<Datatype>> to_teuchos_comm(MPI_Comm comm)
  {
    return std::make_shared<Teuchos::MpiComm<Datatype>>(comm);
  }


  class Communicators
  {
   public:
    Communicators(int groupId, int ngroup, std::map<int, int> lpidgpid, MPI_Comm lcomm,
        MPI_Comm gcomm, NestedParallelismType npType);

    /// return group id
    int group_id() const { return group_id_; }

    /// return number of groups
    int num_groups() const { return ngroup_; }

    /// return group size
    int group_size() const { return Core::Communication::num_mpi_ranks(lcomm_); }

    /// return global processor id of local processor id
    int gpid(int LPID) { return lpidgpid_[LPID]; }

    /// return local processor id of global processor id if GPID is in this group
    int lpid(int GPID);

    /// return local communicator
    MPI_Comm local_comm() const { return lcomm_; }

    /// return local communicator
    MPI_Comm global_comm() const { return gcomm_; }

    /// set a sub group communicator
    void set_sub_comm(MPI_Comm subcomm);

    /// return sub group communicator
    MPI_Comm sub_comm() const { return subcomm_; }

    /// return nested parallelism type
    NestedParallelismType np_type() const { return np_type_; }

   private:
    /// group id
    int group_id_;

    /// number of groups
    int ngroup_;

    /// map from local processor ids to global processor ids
    std::map<int, int> lpidgpid_;

    /// local communicator
    MPI_Comm lcomm_;

    /// global communicator
    MPI_Comm gcomm_;

    /// sub communicator
    MPI_Comm subcomm_;

    /// nested parallelism type
    NestedParallelismType np_type_;
  };


}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
