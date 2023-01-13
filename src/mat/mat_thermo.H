/*! \file
\brief Interface for all thermal material models

\level 3

*/

#ifndef MATERIAL_THERMO_H
#define MATERIAL_THERMO_H

#include "mat_material.H"
#include "mat_trait_thermo.H"

namespace MAT
{
  class ThermoMaterial : public Material, public TRAIT::Thermo
  {
  };
}  // namespace MAT

#endif
