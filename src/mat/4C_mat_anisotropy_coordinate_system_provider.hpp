// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_HPP
#define FOUR_C_MAT_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration

  class CylinderCoordinateSystemProvider;
  /*!
   * \brief Interface of a Coordinate System Provider
   */
  class CoordinateSystemProvider
  {
   public:
    virtual ~CoordinateSystemProvider() = default;

    /*!
     * \brief Returns the reference to the cylinder coordinate system (if present), otherwise
     * #Teuchos::null
     *
     * \return const Teuchos::RCP<CylinderCoordinateSystemProvider>
     */
    virtual Teuchos::RCP<const CylinderCoordinateSystemProvider> get_cylinder_coordinate_system()
        const = 0;
  };

  class CoordinateSystemHolder : public CoordinateSystemProvider
  {
   public:
    Teuchos::RCP<const CylinderCoordinateSystemProvider> get_cylinder_coordinate_system()
        const override
    {
      return cylinder_coordinate_system_;
    }

    void set_cylinder_coordinate_system_provider(
        Teuchos::RCP<const CylinderCoordinateSystemProvider> cylinderCoordinateSystem)
    {
      cylinder_coordinate_system_ = cylinderCoordinateSystem;
    }

   private:
    Teuchos::RCP<const CylinderCoordinateSystemProvider> cylinder_coordinate_system_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
