/*! \file
\brief Interface for all thermal material models

\level 3

*/

#ifndef FOUR_C_MAT_THERMO_HPP
#define FOUR_C_MAT_THERMO_HPP

#include "4C_config.hpp"

#include "4C_mat_material.hpp"
#include "4C_mat_trait_thermo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class ThermoMaterial : public Material, public TRAIT::Thermo
  {
  };
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
