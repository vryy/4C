/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an cylinder coordinate system provider

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_H
#define BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_H

#include "../linalg/linalg_fixedsizematrix.H"


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
    virtual const LINALG::Matrix<3, 1>& GetRad() const = 0;

    /*!
     * \brief Returns a reference to the unit vector pointing in axial direction
     *
     * \return const LINALG::Matrix<3, 1>& Reference to the unit vector pointing in axial direction
     */
    virtual const LINALG::Matrix<3, 1>& GetAxi() const = 0;

    /*!
     * \brief Returns a reference to the unit vector pointing in circumferential direction
     *
     * \return const LINALG::Matrix<3, 1>& Reference to the unit vector pointing in circumferential
     * direction
     */
    virtual const LINALG::Matrix<3, 1>& GetCir() const = 0;
  };

}  // namespace MAT

#endif  // BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_PROVIDER_H
