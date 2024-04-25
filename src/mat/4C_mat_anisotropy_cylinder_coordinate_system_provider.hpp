/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an cylinder coordinate system provider

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_HPP
#define FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


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

FOUR_C_NAMESPACE_CLOSE

#endif
