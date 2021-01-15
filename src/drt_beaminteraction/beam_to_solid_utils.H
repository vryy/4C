/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for beam-to-solid interactions.

\level 3

*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_UTILS_H
#define BEAM_TO_SOLID_UTILS_H


// Forward declarations.
namespace INPAR
{
  namespace BEAMTOSOLID
  {
    enum class BeamToSolidRotationCoupling;
  }  // namespace BEAMTOSOLID
}  // namespace INPAR
namespace LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}  // namespace LINALG
namespace DRT
{
  class Element;
}  // namespace DRT


namespace BEAMINTERACTION
{
  /**
   * \brief Get the rotation vector of a triad constructed in the solid.
   * @param rot_coupling_type (in) Type of triad construction.
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVector(
      const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
      const LINALG::Matrix<3, 1, double>& xi,
      const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad.
   *
   * @param rot_coupling_type (in) Type of triad construction.
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3D(
      const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
      const LINALG::Matrix<3, 1, double>& xi,
      const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Perform a 2D polar decomposition of the deformation gradient and return the rotation
   * vector (2d) of R.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorPolarDecomposition2D(const LINALG::Matrix<3, 1, double>& xi,
      const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on a 2d deformation gradient and return the rotation
   * vector (2d) of said triad.
   *
   * @param rot_coupling_type (in) Type of triad construction.
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient2D(
      const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
      const LINALG::Matrix<3, 1, double>& xi,
      const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Check if the given solid deformation gradient as well as the given beam cross section
   * quaternion are plane with respect to the y-z plane.
   * @param deformation_gradient (in) Deformation gradient at a solid point solid.
   * @param quaternion_beam_ref (in) Quaternion of a beam cross section.
   */
  template <typename scalar_type>
  void CheckPlaneRotations(const LINALG::Matrix<3, 3, scalar_type> deformation_gradient,
      const LINALG::Matrix<4, 1, double>& quaternion_beam_ref);
}  // namespace BEAMINTERACTION

#endif
