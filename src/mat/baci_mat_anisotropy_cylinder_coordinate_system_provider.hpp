/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an cylinder coordinate system provider

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_HPP
#define BACI_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

BACI_NAMESPACE_OPEN


namespace MAT
{
  /*!
   * \brief Interface of a Cylinder Coordinate System Provider
   */
  class CylinderCoordinateSystemProvider
  {
   public:
    virtual ~CylinderCoordinateSystemProvider() = default;
    /*!
     * \brief Returns the unit vector pointing in radial direction
     *
     * \return unit vector pointing in radial direction
     */
    virtual const CORE::LINALG::Matrix<3, 1>& GetRad() const = 0;

    /*!
     * \brief Returns a reference to the unit vector pointing in axial direction
     *
     * \return const CORE::LINALG::Matrix<3, 1>& Reference to the unit vector pointing in axial
     * direction
     */
    virtual const CORE::LINALG::Matrix<3, 1>& GetAxi() const = 0;

    /*!
     * \brief Returns a reference to the unit vector pointing in circumferential direction
     *
     * \return const CORE::LINALG::Matrix<3, 1>& Reference to the unit vector pointing in
     * circumferential direction
     */
    virtual const CORE::LINALG::Matrix<3, 1>& GetCir() const = 0;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_H
