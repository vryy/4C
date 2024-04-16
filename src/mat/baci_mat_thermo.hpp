/*! \file
\brief Interface for all thermal material models

\level 3

*/

#ifndef FOUR_C_MAT_THERMO_HPP
#define FOUR_C_MAT_THERMO_HPP

#include "baci_config.hpp"

#include "baci_mat_material.hpp"
#include "baci_mat_trait_thermo.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  class ThermoMaterial : public Material, public TRAIT::Thermo
  {
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
