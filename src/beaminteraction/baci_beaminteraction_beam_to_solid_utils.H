/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for beam-to-solid interactions.

\level 3

*/


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_UTILS_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_UTILS_H


#include <Teuchos_RCP.hpp>
#include <vector>


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
    enum class BeamToSolidMortarShapefunctions;
  }  // namespace BEAMTOSOLID
}  // namespace INPAR
namespace CORE::LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
  class SparseMatrix;
}  // namespace CORE::LINALG
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
namespace BEAMINTERACTION
{
  class BeamContactPair;
  class BeamToSolidMortarManager;
  class BeamToSolidSurfaceContactParams;
}  // namespace BEAMINTERACTION
class Epetra_Vector;
class Epetra_FEVector;


namespace BEAMINTERACTION
{
  /**
   * \brief Evaluate the penalty force depending on the gap function.
   * @param gap (in) Gap function value.
   * @return Penalty force.
   */
  template <typename scalar_type>
  scalar_type PenaltyForce(const scalar_type& gap,
      const Teuchos::RCP<const BeamToSolidSurfaceContactParams>& contact_params);

  /**
   * \brief Evaluate the penalty potential depending on the gap function.
   * @param gap (in) Gap function value.
   * @return Penalty potential.
   */
  template <typename scalar_type>
  scalar_type PenaltyPotential(const scalar_type& gap,
      const Teuchos::RCP<const BeamToSolidSurfaceContactParams>& contact_params);

  /**
   * \brief Get the number of Lagrange multiplicator values corresponding to the beam nodes and beam
   * element.
   * @param shape_function (in) Mortar shape function.
   * @param n_lambda_node (out) Number of Lagrange multiplicators per node.
   * @param n_lambda_element_ (out) Number of Lagrange multiplicators per element.
   */
  void MortarShapeFunctionsToNumberOfLagrangeValues(
      const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions shape_function,
      unsigned int& n_lambda_node, unsigned int& n_lambda_element);

  /**
   * \brief Get the GIDs of the Lagrange multiplicator unknowns for a beam-to-solid pair.
   * @param mortar_manager (in) Mortar manager for the beam-to-solid condition
   * @param contact_pair (in) Beam-to-solid contact pair
   * @param n_mortar_pos (in) Number of positional mortar DOFs associated with the pair
   * @param n_mortar_rot (in) Number of rotational mortar DOFs associated with the pair
   * @param lambda_gid_pos (out) GIDs of positional mortar DOFs associated with the pair
   * @param lambda_gid_rot (out) GIDs of rotational mortar DOFs associated with the pair
   */
  void GetMortarGID(const BeamToSolidMortarManager* mortar_manager,
      const BEAMINTERACTION::BeamContactPair* contact_pair, const unsigned int n_mortar_pos,
      const unsigned int n_mortar_rot, std::vector<int>* lambda_gid_pos,
      std::vector<int>* lambda_gid_rot);

  /**
   * \brief Setup the triad interpolation scheme for the current triad and reference triad of the
   * given beam element.
   * @param discret (in) Discretization.
   * @param displacement_vector (in) Global displacement vector.
   * @param ele (in) Pointer to the beam element.
   * @param triad_interpolation_scheme (out) Interpolation of current triad field..
   * @param ref_triad_interpolation_scheme (out) Interpolation of reference triad field.
   */
  void GetBeamTriadInterpolationScheme(const ::DRT::Discretization& discret,
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
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on the average vector of the deformed triad.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DGeneral(
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on cross section basis vectors.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on cross section basis vectors.
   *
   * @param F (in) Deformation gradient.
   * @param beam_ref_triad (in) Reference triad of the beam.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   */
  template <typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
      const CORE::LINALG::Matrix<3, 3, scalar_type>& F,
      const CORE::LINALG::Matrix<3, 3, double>& beam_ref_triad,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on the first basis vector of the deformed
   * triad.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   * @param element (in) Pointer to the solid element.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DBase1(
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction starts with a user-given base vector.
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
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

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
  void GetSolidRotationVectorPolarDecomposition2D(const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

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
      const CORE::LINALG::Matrix<3, 1, double>& xi,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
      const CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
      CORE::LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element = nullptr);

  /**
   * \brief Check if the given solid deformation gradient as well as the given beam cross section
   * quaternion are plane with respect to the y-z plane.
   * @param deformation_gradient (in) Deformation gradient at a solid point solid.
   * @param quaternion_beam_ref (in) Quaternion of a beam cross section.
   */
  template <typename scalar_type>
  void CheckPlaneRotations(const CORE::LINALG::Matrix<3, 3, scalar_type> deformation_gradient,
      const CORE::LINALG::Matrix<4, 1, double>& quaternion_beam_ref);

  /**
   * \brief Assemble local mortar contributions from the classical mortar matrices D and M into the
   * global matrices.
   *
   * This function assumes that the mortar contributions are symmetric, i.e. global_G_B =
   * global_FB_L^T and global_G_S = global_FS_L^T.
   *
   * @param pair (in) The beam-to-solid pair.
   * @param discret (in) Discretization
   * @param mortar_manager (in) Mortar manager for the beam-to-solid condition
   * @param global_G_B (in/out) Constraint equations derived w.r.t the beam DOFs
   * @param global_G_S (in/out) Constraint equations derived w.r.t the solid DOFs
   * @param global_FB_L (in/out) Beam force vector derived w.r.t the Lagrange multipliers
   * @param global_FS_L (in/out) Solid force vector derived w.r.t the Lagrange multipliers
   * @param global_constraint (in/out) Global constraint equations
   * @param global_kappa (in/out) Global penalty scaling vector equations
   * @param global_lambda_active (in/out) Global vector keeping track of active lagrange multipliers
   * @param local_D (in) Local D matrix of the pair.
   * @param local_M (in) Local M matrix of the pair.
   * @param local_kappa (in) Local scaling vector of the pair.
   * @param local_constraint (in) Local constraint contributions of the pair.
   * @param n_mortar_rot (int) Number of total rotational Lagrange multiplier DOFs per beam.
   */
  template <typename beam, typename other, typename mortar>
  void AssembleLocalMortarContributions(const BEAMINTERACTION::BeamContactPair* pair,
      const ::DRT::Discretization& discret, const BeamToSolidMortarManager* mortar_manager,
      CORE::LINALG::SparseMatrix& global_G_B, CORE::LINALG::SparseMatrix& global_G_S,
      CORE::LINALG::SparseMatrix& global_FB_L, CORE::LINALG::SparseMatrix& global_FS_L,
      Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
      Epetra_FEVector& global_lambda_active,
      const CORE::LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
      const CORE::LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>& local_M,
      const CORE::LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
      const CORE::LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint,
      const unsigned int n_mortar_rot = 0);
}  // namespace BEAMINTERACTION

#endif
