/*! \file
\brief Interface for all thermo-mechanical material models

\level 3

*/

#ifndef FOUR_C_MAT_THERMOMECHANICAL_HPP
#define FOUR_C_MAT_THERMOMECHANICAL_HPP

#include "4C_config.hpp"

#include "4C_mat_so3_material.hpp"
#include "4C_mat_trait_thermo_solid.hpp"

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  /*!
   * Interface for all thermo-mechanical materials
   *
   * @note this interface  is inspired by the way materials are supposed to work in 4C, not all
   * materials are moved to this new interface
   */
  class ThermoMechanicalMaterial : public So3Material, public TRAIT::ThermoSolid
  {
  };

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif
