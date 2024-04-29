/*----------------------------------------------------------------------*/
/*! \file
\brief Factory for materials

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MATERIAL_FACTORY_HPP
#define FOUR_C_MAT_MATERIAL_FACTORY_HPP



#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_inpar_material.hpp"
#include "4C_material_base.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/// MAT: materials
namespace MAT
{
  const int NUM_STRESS_3D = 6;  ///< 6 stresses for 3D

  /// create element material object given the number of a material definition
  Teuchos::RCP<CORE::MAT::Material> Factory(int matnum  ///< material ID
  );

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
