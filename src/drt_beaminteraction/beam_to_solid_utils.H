/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for beam-to-solid interactions.

\level 3

*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_UTILS_H
#define BEAM_TO_SOLID_UTILS_H


#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidMortarManager;
  class BeamContactPair;
}  // namespace BEAMINTERACTION
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
  class SparseMatrix;
}  // namespace LINALG
namespace DRT
{
  class Element;
  class Discretization;
}  // namespace DRT
namespace LARGEROTATIONS
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}
class Epetra_Vector;
class Epetra_FEVector;


namespace BEAMINTERACTION
{
  /**
   * \brief Setup the triad interpolation scheme for the current triad and reference triad of the
   * given beam element.
   * @param discret (in) Discretization.
   * @param displacement_vector (in) Global displacement vector.
   * @param ele (in) Pointer to the beam element.
   * @param triad_interpolation_scheme (out) Interpolation of current triad field..
   * @param ref_triad_interpolation_scheme (out) Interpolation of reference triad field.
   */
  void GetBeamTriadInterpolationScheme(const DRT::Discretization& discret,
      const Teuchos::RCP<const Epetra_Vector>& displacement_vector, const DRT::Element* ele,
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>& triad_interpolation_scheme,
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>&
          ref_triad_interpolation_scheme);

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

  /**
   * \brief Assemble local mortar contributions from the classical mortar matrices D and M into the
   * global matrices.
   *
   * This function assumes that the mortar contributions are symmetric, i.e. global_GB =
   * global_FB^T and global_GS = global_FS^T.
   *
   * @param pair (in) The beam-to-solid pair.
   * @param discret (in)
   * @param mortar_manager (in)
   * @param global_GB (in/out)
   * @param global_GS (in/out)
   * @param global_FB (in/out)
   * @param global_FS (in/out)
   * @param global_constraint (in/out)
   * @param global_kappa (in/out)
   * @param global_lambda_active (in/out)
   * @param local_D (in) Local D matrix of the pair.
   * @param local_M (in) Local M matrix of the pair.
   * @param local_kappa (in) Local scaling vector of the pair.
   * @param local_constraint (in) Local constraint contributions of the pair.
   */
  template <typename beam, typename other, typename mortar>
  void AssembleLocalMortarContributions(const BEAMINTERACTION::BeamContactPair* pair,
      const DRT::Discretization& discret, const BeamToSolidMortarManager* mortar_manager,
      LINALG::SparseMatrix& global_GB, LINALG::SparseMatrix& global_GS,
      LINALG::SparseMatrix& global_FB, LINALG::SparseMatrix& global_FS,
      Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
      Epetra_FEVector& global_lambda_active,
      const LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
      const LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>& local_M,
      const LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
      const LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint);
}  // namespace BEAMINTERACTION

#endif
