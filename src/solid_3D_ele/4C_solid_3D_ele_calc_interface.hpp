/*! \file

\brief Interface of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_INTERFACE_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_INTERFACE_HPP


#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret::ELEMENTS
{
  struct StressIO
  {
    Inpar::Solid::StressType type;
    std::vector<char>& mutable_data;
  };

  struct StrainIO
  {
    Inpar::Solid::StrainType type;
    std::vector<char>& mutable_data;
  };

}  // namespace Discret::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif
