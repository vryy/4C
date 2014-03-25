/*----------------------------------------------------------------------*/
/*!
\file crackUtils.H

\brief Utility functions for crack propagation problem

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef CRACKUTILS_H_
#define CRACKUTILS_H_

#include "../drt_lib/drt_discret.H"

namespace DRT {
namespace CRACK {
namespace UTILS{


  /*!
  \brief Update the given EpetraVector for the new map of discretization
   */
  void UpdateThisEpetraVectorCrack( const Teuchos::RCP<DRT::Discretization>& discret, Teuchos::RCP<Epetra_Vector>& vec,
                                  const std::map<int,int>& oldnewIds );

} /* namespace CRACK */
} /* namespace CRACK */
} /* namespace DRT */
#endif
