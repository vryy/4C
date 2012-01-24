/*----------------------------------------------------------------------*/
/*!
\file comm_utils.cpp

\brief Helper class for everything that deals with communication

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

#include <Epetra_MpiComm.h>

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include "comm_utils.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 | converts Epetra_Comm into Teuchos::Comm                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::Comm<int> > COMM_UTILS::toTeuchosComm(
  const Epetra_Comm& comm
  )
{
  try {
    const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
    Teuchos::RCP<Teuchos::MpiComm<int> > mpicomm =  Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm.Comm())));
    return Teuchos::rcp_dynamic_cast<const Teuchos::Comm<int> >(mpicomm);
  }
  catch (std::bad_cast & b) {}
  dserror("Something went wrong with converting the communicator! You should not be here!");

  return Teuchos::null;
}



/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
