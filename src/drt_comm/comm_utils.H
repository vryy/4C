/*----------------------------------------------------------------------*/
/*!
\file comm_utils.H

\brief Helper class for everything that deals with communication, e.g.
       MPI, Epetra_Comm and further communicators

\level 0

\maintainer Martin Kronbichler

*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 01/12 |
 *----------------------------------------------------------------------*/

#ifndef COMM_UTILS_H
#define COMM_UTILS_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include <Teuchos_DefaultMpiComm.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem_enums.H"

namespace DRT
{
  class Discretization;
  class Problem;
}  // namespace DRT

/*----------------------------------------------------------------------*
 |                                                          ghamm 01/12 |
 *----------------------------------------------------------------------*/
namespace COMM_UTILS
{
  // forward declaration
  class NestedParGroup;

  //! create a local and a global communicator for the problem
  void CreateComm(int argc, char** argv);

#if (0)
  //! distribute discretizations from one group to the others
  void BroadcastDiscretizations(const int bgroup);
#endif

  //! distribute discretizations from one group to the others with ponzi scheme
  void BroadcastDiscretizations(int instance = 0);

  //! distribute a discretization from one group to one other
  void NPDuplicateDiscretization(const int sgroup, const int rgroup,
      Teuchos::RCP<NestedParGroup> group, Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<Epetra_MpiComm> icomm);

  //! distribute a discretization from one group to one other
  //! if the groups have the same size. In this special case it is desirable
  //! that the parallel distribution is the same for all groups. Hence,
  //! a different braodcasting scheme is used
  void NPDuplicateDiscretizationEqualGroupSize(const int sgroup, const int rgroup,
      Teuchos::RCP<NestedParGroup> group, Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<Epetra_MpiComm> icomm);

  //! comparison of two vectors from parallel baci runs for debugging reason using the given
  //! tolerance for each entry
  void CompareVectors(
      Teuchos::RCP<const Epetra_MultiVector> vec, const char* name, double tol = 1.0e-14);

  //! comparison of two sparse matrices from parallel baci runs for debugging reason using the given
  //! tolerance for each entry
  void CompareSparseMatrices(
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


  class NestedParGroup
  {
   public:
    NestedParGroup(int groupId, int ngroup, std::map<int, int> lpidgpid,
        Teuchos::RCP<Epetra_Comm> lcomm, Teuchos::RCP<Epetra_Comm> gcomm, NP_TYPE npType);

    NestedParGroup(NestedParGroup& npgroup);

    virtual ~NestedParGroup(){};

    /// return group id
    int GroupId() { return groupId_; }

    /// return number of groups
    int NumGroups() { return ngroup_; }

    /// return group size
    int GroupSize() { return lcomm_->NumProc(); }

    /// return global processor id of local processor id
    int GPID(int LPID) { return lpidgpid_[LPID]; }

    /// return local processor id of global processor id if GPID is in this group
    int LPID(int GPID);

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> LocalComm() { return lcomm_; }

    /// return local communicator
    Teuchos::RCP<Epetra_Comm> GlobalComm() { return gcomm_; }

    /// set a sub group communicator
    void SetSubComm(Teuchos::RCP<Epetra_Comm> subcomm);

    /// return sub group communicator
    Teuchos::RCP<Epetra_Comm> SubComm() { return subcomm_; }

    /// return nested parallelism type
    NP_TYPE NpType() { return npType_; }

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
    NP_TYPE npType_;
  };


}  // namespace COMM_UTILS



/*----------------------------------------------------------------------*/
#endif  // COMM_UTILS_H
