/*----------------------------------------------------------------------*/
/*! \file

\brief utils for biofilm fs3i

\level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}

namespace FS3I
{
  namespace BioFilm
  {
    namespace UTILS
    {
      void ScatraChangeConfig(Teuchos::RCP<Discret::Discretization> scatradis,
          Teuchos::RCP<Discret::Discretization> dis, Teuchos::RCP<Epetra_Vector> disp);

    } /* namespace UTILS */
  }   // namespace BioFilm
} /* namespace FS3I */

FOUR_C_NAMESPACE_CLOSE

#endif
