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

/*----------------------------------------------------------------------*
 |                                                          ghamm 01/12 |
 *----------------------------------------------------------------------*/
/// COMM_UTILS
namespace COMM_UTILS
{
  //! create a local and a global communicator for the problem
  void CreateComm(int argc, char** argv);

  //! transform Epetra_Comm to Teuchos::Comm, RCP version
  Teuchos::RCP<const Teuchos::Comm<Teuchos::Ordinal> > toTeuchosComm(const Epetra_Comm & comm);

} // namespace COMM_UTILS



/*----------------------------------------------------------------------*/
#endif  // COMM_UTILS_H
#endif  // CCADISCRET
