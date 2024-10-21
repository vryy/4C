#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_INTERFACE_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_INTERFACE_HPP


#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret::ELEMENTS
{
  struct CouplStressIO
  {
    Inpar::Solid::StressType type;
    std::vector<char>& mutable_data;
  };


}  // namespace Discret::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif