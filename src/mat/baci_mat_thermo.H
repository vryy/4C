/*! \file
\brief Interface for all thermal material models

\level 3

*/

#ifndef BACI_MAT_THERMO_H
#define BACI_MAT_THERMO_H

#include "baci_config.H"

#include "baci_mat_material.H"
#include "baci_mat_trait_thermo.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  class ThermoMaterial : public Material, public TRAIT::Thermo
  {
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
