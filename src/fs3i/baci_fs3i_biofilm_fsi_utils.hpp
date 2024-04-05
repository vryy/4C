/*----------------------------------------------------------------------*/
/*! \file

\brief utils for biofilm fs3i

\level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP

#include "baci_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

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

BACI_NAMESPACE_CLOSE

#endif
