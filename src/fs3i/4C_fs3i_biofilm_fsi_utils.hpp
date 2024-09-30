/*----------------------------------------------------------------------*/
/*! \file

\brief utils for biofilm fs3i

\level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_UTILS_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FS3I
{
  namespace BioFilm
  {
    namespace UTILS
    {
      void scatra_change_config(Teuchos::RCP<Core::FE::Discretization> scatradis,
          Teuchos::RCP<Core::FE::Discretization> dis, Teuchos::RCP<Core::LinAlg::Vector> disp);

    } /* namespace UTILS */
  }   // namespace BioFilm
} /* namespace FS3I */

FOUR_C_NAMESPACE_CLOSE

#endif
