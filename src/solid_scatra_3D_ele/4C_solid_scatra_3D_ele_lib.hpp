#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_LIB_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_LIB_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_io_linedefinition.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  /*!
   * @brief Read the scatra implementation type from the a container of data to create the element
   *
   * @param container
   * @return Inpar::ScaTra::ImplType
   */
  Inpar::ScaTra::ImplType read_scatra_impl_type(const Core::IO::InputParameterContainer& container);
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
