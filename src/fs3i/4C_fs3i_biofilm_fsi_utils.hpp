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
    namespace Utils
    {
      void scatra_change_config(Core::FE::Discretization& scatradis, Core::FE::Discretization& dis,
          Core::LinAlg::Vector<double>& disp);

    } /* namespace Utils */
  }   // namespace BioFilm
} /* namespace FS3I */

FOUR_C_NAMESPACE_CLOSE

#endif
