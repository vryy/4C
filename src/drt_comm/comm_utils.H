/*----------------------------------------------------------------------*/
/*!
\file comm_utils.H

\brief Helper class for everything that deals with communication, e.g.
       MPI, Epetra_Comm and further communicators

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 01/12 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef COMM_UTILS_H
#define COMM_UTILS_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include <Teuchos_DefaultMpiComm.hpp>
#include <Epetra_MpiComm.h>
#include <mpi.h>

#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |                                                          ghamm 01/12 |
 *----------------------------------------------------------------------*/
/// COMM_UTILS
namespace COMM_UTILS
{
  //! create a local and a global communicator for the problem
  void CreateComm(int argc, char** argv);

  //! transform Epetra_Comm to Teuchos::Comm, RCP version
  template<class datatype>
  Teuchos::RCP<const Teuchos::Comm<datatype> > toTeuchosComm(const Epetra_Comm & comm)
  {
    try {
      const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
      Teuchos::RCP<Teuchos::MpiComm<datatype> > mpicomm =  Teuchos::rcp(new Teuchos::MpiComm<datatype>(Teuchos::opaqueWrapper(mpiComm.Comm())));
      return Teuchos::rcp_dynamic_cast<const Teuchos::Comm<datatype> >(mpicomm);
    }
    catch (std::bad_cast & b)
    {
      dserror("Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm object is unknown");
    }
    dserror("Something went wrong with converting an Epetra_Comm to a Teuchos communicator! You should not be here!");
    return Teuchos::null;
  }

} // namespace COMM_UTILS



/*----------------------------------------------------------------------*/
#endif  // COMM_UTILS_H
#endif  // CCADISCRET
