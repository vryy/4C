/*----------------------------------------------------------------------*/
/*! \file
\brief Helper class for everything that deals with communication, e.g.
       MPI, Epetra_Comm and further communicators
\level 0
*/


#ifndef BACI_COMM_UTILS_H
#define BACI_COMM_UTILS_H


#include "baci_config.H"

#include "baci_utils_exceptions.H"

#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_DefaultMpiComm.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::COMM
{
  // forward declaration
  class Communicators;

  /**
   * The known types for nested parallelism.
   */
  enum class NestedParallelismType
  {
    copy_dat_file,
    every_group_read_dat_file,
    separate_dat_files,
    no_nested_parallelism
  };

  //! create a local and a global communicator for the problem
  Teuchos::RCP<Communicators> CreateComm(std::vector<std::string> argv);

  /*! \brief debug routine to compare vectors from different parallel baci runs
   *
   * You can add CORE::COMM::AreDistributedVectorsIdentical in your code which will lead to a
   * comparison of the given vector for different executables and/or configurations.
   * Command for using this feature: \n
   * mpirun -np 1 ./baci-release -nptype=diffgroup0 input.dat xxx_ser : -np 3 ./other-baci-release
   * -nptype=diffgroup1 other-input.dat xxx_par \n
   * Do not forget to include the header (#include "baci_comm_utils.H"), otherwise it won't compile.
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

  /*! \brief debug routine to compare sparse matrices from different parallel baci runs
   *
   * You can add CORE::COMM::AreDistributedSparseMatricesIdentical in your code which will lead to a
   * comparison of the given sparse matrices for different executables and/or configurations.
   * Command for using this feature: \n
   * mpirun -np 1 ./baci-release -nptype=diffgroup0 input.dat xxx_ser : -np 3 ./other-baci-release
   * -nptype=diffgroup1 other-input.dat xxx_par \n
   * Do not forget to include the header (#include "baci_comm_utils.H"), otherwise it won't compile.
   *
   * A further nice option is to compare results from different executables used for
   * running the same simulation.
   *
   * \note You need to add the AreDistributedSparseMatricesIdentical method in both executables at
   * the same position in the code.
   *
   * \note From CORE::LINALG::SparseOperator to CrsMatrix, just do:
   * Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(yoursparseoperator)->EpetraMatrix()
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
  template <class datatype>
  Teuchos::RCP<const Teuchos::Comm<datatype>> toTeuchosComm(const Epetra_Comm& comm)
  {
    try
    {
      const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
      Teuchos::RCP<Teuchos::MpiComm<datatype>> mpicomm =
          Teuchos::rcp(new Teuchos::MpiComm<datatype>(Teuchos::opaqueWrapper(mpiComm.Comm())));
      return Teuchos::rcp_dynamic_cast<const Teuchos::Comm<datatype>>(mpicomm);
    }
    catch (std::bad_cast& b)
    {
      dserror(
          "Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm "
          "object is unknown");
    }
    dserror(
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
    int GroupId() const { return groupId_; }

    /// return number of groups
    int NumGroups() const { return ngroup_; }

    /// return group size
    int GroupSize() const { return lcomm_->NumProc(); }

    /// return global processor id of local processor id
    int GPID(int LPID) { return lpidgpid_[LPID]; }

    /// return local processor id of global processor id if GPID is in this group
    int LPID(int GPID);

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> LocalComm() const { return lcomm_; }

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> GlobalComm() const { return gcomm_; }

    /// set a sub group communicator
    void SetSubComm(Teuchos::RCP<Epetra_Comm> subcomm);

    /// return sub group communicator
    Teuchos::RCP<Epetra_Comm> SubComm() const { return subcomm_; }

    /// return nested parallelism type
    NestedParallelismType NpType() const { return npType_; }

   private:
    /// group id
    int groupId_;

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
    NestedParallelismType npType_;
  };


}  // namespace CORE::COMM

BACI_NAMESPACE_CLOSE

#endif
