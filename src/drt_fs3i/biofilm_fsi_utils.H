/*----------------------------------------------------------------------*/
/*! \file

\brief utils for biofilm fs3i

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15249

 *----------------------------------------------------------------------*/

#ifndef BIOFILMFSIUTILS_H_
#define BIOFILMFSIUTILS_H_

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
}

namespace FS3I
{
  namespace BIOFILM
  {
    namespace UTILS
    {
      void ScatraChangeConfig(Teuchos::RCP<DRT::Discretization> scatradis,
          Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Vector> disp);

    } /* namespace UTILS */
  }   // namespace BIOFILM
} /* namespace FS3I */

#endif
