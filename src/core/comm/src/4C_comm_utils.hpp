/*----------------------------------------------------------------------*/
/*! \file
\brief Helper class for everything that deals with communication, e.g.
       MPI, Epetra_Comm and further communicators
\level 0
*/


#ifndef FOUR_C_COMM_UTILS_HPP
#define FOUR_C_COMM_UTILS_HPP


#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>
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
  Teuchos::RCP<Communicators> CreateComm(std::vector<std::string> argv);

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
  bool AreDistributedVectorsIdentical(const Communicators& communicators,
      Teuchos::RCP<const Epetra_MultiVector> vec, const char* name, double tol = 1.0e-14);

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
   * Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(yoursparseoperator)->EpetraMatrix()
   *
   * \param communicators (in): communicators containing local and global comm
   * \param matrix        (in): matrix to compare
   * \param name          (in): user given name for the matrix (needs to match within gcomm)
   * \param tol           (in): comparison tolerance for infinity norm
   * \return boolean to indicate if compared vectors are identical
   */
  bool AreDistributedSparseMatricesIdentical(const Communicators& communicators,
      Teuchos::RCP<Epetra_CrsMatrix> matrix, const char* name, double tol = 1.0e-14);

  //! transform Epetra_Comm to Teuchos::Comm, Teuchos::RCP version
  template <class Datatype>
  Teuchos::RCP<const Teuchos::Comm<Datatype>> toTeuchosComm(const Epetra_Comm& comm)
  {
    try
    {
      const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
      Teuchos::RCP<Teuchos::MpiComm<Datatype>> mpicomm =
          Teuchos::rcp(new Teuchos::MpiComm<Datatype>(Teuchos::opaqueWrapper(mpiComm.Comm())));
      return Teuchos::rcp_dynamic_cast<const Teuchos::Comm<Datatype>>(mpicomm);
    }
    catch (std::bad_cast& b)
    {
      FOUR_C_THROW(
          "Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm "
          "object is unknown");
    }
    FOUR_C_THROW(
        "Something went wrong with converting an Epetra_Comm to a Teuchos communicator! You should "
        "not be here!");
    return Teuchos::null;
  }


  class Communicators
  {
   public:
    Communicators(int groupId, int ngroup, std::map<int, int> lpidgpid,
        Teuchos::RCP<Epetra_Comm> lcomm, Teuchos::RCP<Epetra_Comm> gcomm,
        NestedParallelismType npType);

    /// return group id
    int group_id() const { return group_id_; }

    /// return number of groups
    int num_groups() const { return ngroup_; }

    /// return group size
    int group_size() const { return lcomm_->NumProc(); }

    /// return global processor id of local processor id
    int gpid(int LPID) { return lpidgpid_[LPID]; }

    /// return local processor id of global processor id if GPID is in this group
    int lpid(int GPID);

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> local_comm() const { return lcomm_; }

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> global_comm() const { return gcomm_; }

    /// set a sub group communicator
    void set_sub_comm(Teuchos::RCP<Epetra_Comm> subcomm);

    /// return sub group communicator
    Teuchos::RCP<Epetra_Comm> sub_comm() const { return subcomm_; }

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
    Teuchos::RCP<Epetra_Comm> lcomm_;

    /// global communicator
    Teuchos::RCP<Epetra_Comm> gcomm_;

    /// sub communicator
    Teuchos::RCP<Epetra_Comm> subcomm_;

    /// nested parallelism type
    NestedParallelismType np_type_;
  };


}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
