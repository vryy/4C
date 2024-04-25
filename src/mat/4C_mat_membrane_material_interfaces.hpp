/*----------------------------------------------------------------------*/
/*! \file
\brief Interface classes for membrane materials

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MEMBRANE_MATERIAL_INTERFACES_HPP
#define FOUR_C_MAT_MEMBRANE_MATERIAL_INTERFACES_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  /*!
   * @brief Base interface for membrane materials in local membrane coordinate system
   *
   * This class provides an interface for membrane materials that are evaluated using a local
   * membrane coordinate system. The interface provides a method to compute the in-plane stresses of
   * membranes given the deformation in the local coordinate system.
   *
   */
  class MembraneMaterialLocalCoordinates
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~MembraneMaterialLocalCoordinates() = default;

    /*!
     * @brief Update material in local coordinates
     *
     * @param defgrd (in) : Deformation gradient in local coordinates
     * @param params (in) : Container for additional information
     * @param Q_trafo (in) : Transformation from local membrane orthonormal coordinates to global
     * coordinates
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void UpdateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd,
        Teuchos::ParameterList& params, const CORE::LINALG::Matrix<3, 3>& Q_trafo, int gp,
        int eleGID) = 0;

    /*!
     * @brief Evaluate stress response plus elasticity tensor for membranes assuming
     * incompressibility and plane stress
     *
     * @param defgrd (in) : Deformation gradient in local coordinates
     * @param cauchygreen (in) : right Cauchy-Green tensor in local coordinates
     * @param params (in) : Container for additional information
     * @param Q_trafo (in) : Transformation from local membrane orthonormal coordinates to global
     * coordinates
     * @param stress (out) : 2nd Piola Kirchhoff stress in stress-like Voigt notation in local
     * membrane coordinates
     * @param cmat (out) : elasticity tensor in local membrane coordinates
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void EvaluateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd,
        const CORE::LINALG::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
        const CORE::LINALG::Matrix<3, 3>& Q_trafo, CORE::LINALG::Matrix<3, 1>& stress,
        CORE::LINALG::Matrix<3, 3>& cmat, int gp, int eleGID) = 0;
  };

  /*!
   * @brief Base interface for membrane materials that are evaluated using the global coordinate
   * system.
   *
   */
  class MembraneMaterialGlobalCoordinates
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~MembraneMaterialGlobalCoordinates() = default;

    /*!
     * @brief Evaluate stress response plus elasticity tensor for membranes assuming
     * incompressibility and plane stress
     *
     * @param defgrd (in) : Deformation gradient in global coordinates
     * @param params (in) : Container for additional information
     * @param stress (out) : 2nd Piola Kirchhoff stress in stress-like Voigt notation in global
     * coordinates
     * @param cmat (out) : elasticity tensor in global coordinates
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void EvaluateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<3, 3>& stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) = 0;
  };

  /*!
   * @brief Base interface for membrane materials that have an inelastic deformation in thickness
   * direction.
   *
   * Usually, incompressibility is only satisfied by the elastic deformation. This interface let the
   * material interface decide the thickness of the membrane.
   *
   */
  class MembraneMaterialInelasticThickness
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~MembraneMaterialInelasticThickness() = default;

    /*!
     * @brief Evaluate membrane thickness stretch (mostly such that the elastic deformation is
     * incompressible)
     *
     * @param defgrd (in) : Deformation gradient in global coordinates
     * @param params (in) : Container for additional information
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     *
     * @return Adapted stretch in thickness direction
     */
    virtual double EvaluateMembraneThicknessStretch(const CORE::LINALG::Matrix<3, 3>& defgrd,
        Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
