/*! \file
\brief Interface for all thermo-mechanical material models

\level 3

*/

#ifndef BACI_MAT_THERMOMECHANICAL_H
#define BACI_MAT_THERMOMECHANICAL_H

#include "baci_config.H"

#include "baci_mat_so3_material.H"
#include "baci_mat_trait_thermo_solid.H"

BACI_NAMESPACE_OPEN


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


BACI_NAMESPACE_CLOSE

#endif
