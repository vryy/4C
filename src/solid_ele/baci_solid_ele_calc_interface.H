/*! \file

\brief Interface of solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_INTERFACE_H
#define BACI_SOLID_ELE_CALC_INTERFACE_H


#include "baci_inpar_structure.H"


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
#endif  // SOLID_ELE_CALC_INTERFACE_H
