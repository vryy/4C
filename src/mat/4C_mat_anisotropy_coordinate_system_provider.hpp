/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an coordinate system provider

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_HPP
#define FOUR_C_MAT_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_HPP

#include "4C_config.hpp"

#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace MAT
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

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
