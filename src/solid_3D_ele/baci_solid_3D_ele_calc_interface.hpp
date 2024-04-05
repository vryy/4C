/*! \file

\brief Interface of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_INTERFACE_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_INTERFACE_HPP


#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"

BACI_NAMESPACE_OPEN


namespace DRT::ELEMENTS
{
  struct StressIO
  {
    INPAR::STR::StressType type;
    std::vector<char>& mutable_data;
  };

  struct StrainIO
  {
    INPAR::STR::StrainType type;
    std::vector<char>& mutable_data;
  };

}  // namespace DRT::ELEMENTS
BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_CALC_INTERFACE_H
