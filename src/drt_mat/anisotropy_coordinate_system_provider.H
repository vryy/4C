/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an coordinate system provider

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_H
#define BACI_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_H

#include <Teuchos_RCPDecl.hpp>

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
    virtual const Teuchos::RCP<const CylinderCoordinateSystemProvider> GetCylinderCoordinateSystem()
        const = 0;
  };

  class CoordinateSystemHolder : public CoordinateSystemProvider
  {
   public:
    const Teuchos::RCP<const CylinderCoordinateSystemProvider> GetCylinderCoordinateSystem()
        const override
    {
      return cylinderCoordinateSystem_;
    }

    void SetCylinderCoordinateSystemProvider(
        Teuchos::RCP<const CylinderCoordinateSystemProvider> cylinderCoordinateSystem)
    {
      cylinderCoordinateSystem_ = cylinderCoordinateSystem;
    }

   private:
    Teuchos::RCP<const CylinderCoordinateSystemProvider> cylinderCoordinateSystem_;
  };

}  // namespace MAT

#endif  // BACI_ANISOTROPY_COORDINATE_SYSTEM_PROVIDER_H
