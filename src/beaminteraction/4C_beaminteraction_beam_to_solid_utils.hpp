/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for beam-to-solid interactions.

\level 3

*/


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_UTILS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_UTILS_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_element.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidMortarManager;
  class BeamContactPair;
}  // namespace BEAMINTERACTION
namespace Inpar
{
  namespace BeamToSolid
  {
    enum class BeamToSolidRotationCoupling;
    enum class BeamToSolidMortarShapefunctions;
  }  // namespace BeamToSolid
}  // namespace Inpar
namespace Core::LinAlg
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace LargeRotations
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
   * @param n_dim (in) Spatial dimension of Lagrange multiplicator field.
   * @return {n_lambda_node, n_lambda_element} Number of Lagrange multiplicators per node and per
   * element.
   */
  [[nodiscard]] std::pair<unsigned int, unsigned int> MortarShapeFunctionsToNumberOfLagrangeValues(
      const Inpar::BeamToSolid::BeamToSolidMortarShapefunctions shape_function,
      const unsigned int n_dim);

  /**
   * \brief Setup the triad interpolation scheme for the current triad and reference triad of the
   * given beam element.
   * @param discret (in) discretization.
   * @param displacement_vector (in) Global displacement vector.
   * @param ele (in) Pointer to the beam element.
   * @param triad_interpolation_scheme (out) Interpolation of current triad field..
   * @param ref_triad_interpolation_scheme (out) Interpolation of reference triad field.
   */
  void GetBeamTriadInterpolationScheme(const Core::FE::Discretization& discret,
      const Teuchos::RCP<const Epetra_Vector>& displacement_vector,
      const Core::Elements::Element* ele,
      LargeRotations::TriadInterpolationLocalRotationVectors<3, double>& triad_interpolation_scheme,
      LargeRotations::TriadInterpolationLocalRotationVectors<3, double>&
          ref_triad_interpolation_scheme);

  /**
   * \brief Get the rotation vector of a triad constructed in the solid.
   * @param rot_coupling_type (in) Type of triad construction.
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVector(
      const Inpar::BeamToSolid::BeamToSolidRotationCoupling& rot_coupling_type,
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on the average vector of the deformed triad.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DGeneral(
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

  /**
   * \brief Construct a solid triad depending on the deformation gradient and return the rotation
   * vector of said triad. The construction is based on cross section basis vectors.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

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
      const Core::LinAlg::Matrix<3, 3, scalar_type>& F,
      const Core::LinAlg::Matrix<3, 3, double>& beam_ref_triad,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

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
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3DBase1(
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

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
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient3D(
      const Inpar::BeamToSolid::BeamToSolidRotationCoupling& rot_coupling_type,
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

  /**
   * \brief Perform a 2D polar decomposition of the deformation gradient and return the rotation
   * vector (2d) of R.
   *
   * @param xi (in) Parameter coordinates in the solid.
   * @param q_solid_ref (in) Reference position of the solid.
   * @param q_solid (in) Displacement of the solid.
   * @param quaternion_beam_ref (in) Beam reference quaternion at the solid point.
   * @param psi_solid (out) Rotation vector of the constructed solid triad.
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorPolarDecomposition2D(const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

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
   */
  template <typename solid, typename scalar_type>
  void GetSolidRotationVectorDeformationGradient2D(
      const Inpar::BeamToSolid::BeamToSolidRotationCoupling& rot_coupling_type,
      const Core::LinAlg::Matrix<3, 1, double>& xi,
      const GEOMETRYPAIR::ElementData<solid, double>& q_solid_ref,
      const GEOMETRYPAIR::ElementData<solid, scalar_type>& q_solid,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref,
      Core::LinAlg::Matrix<3, 1, scalar_type>& psi_solid);

  /**
   * \brief Check if the given solid deformation gradient as well as the given beam cross section
   * quaternion are plane with respect to the y-z plane.
   * @param deformation_gradient (in) Deformation gradient at a solid point solid.
   * @param quaternion_beam_ref (in) Quaternion of a beam cross section.
   */
  template <typename scalar_type>
  void CheckPlaneRotations(const Core::LinAlg::Matrix<3, 3, scalar_type> deformation_gradient,
      const Core::LinAlg::Matrix<4, 1, double>& quaternion_beam_ref);

  /**
   * \brief Assemble local mortar contributions from the classical mortar matrices D and M into the
   * global matrices.
   *
   * This function assumes that the mortar contributions are symmetric, i.e. global_G_B =
   * global_FB_L^T and global_G_S = global_FS_L^T.
   *
   * @param pair (in) The beam-to-solid pair.
   * @param discret (in) discretization
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
      const Core::FE::Discretization& discret, const BeamToSolidMortarManager* mortar_manager,
      Core::LinAlg::SparseMatrix& global_G_B, Core::LinAlg::SparseMatrix& global_G_S,
      Core::LinAlg::SparseMatrix& global_FB_L, Core::LinAlg::SparseMatrix& global_FS_L,
      Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
      Epetra_FEVector& global_lambda_active,
      const Core::LinAlg::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
      const Core::LinAlg::Matrix<mortar::n_dof_, other::n_dof_, double>& local_M,
      const Core::LinAlg::Matrix<mortar::n_dof_, 1, double>& local_kappa,
      const Core::LinAlg::Matrix<mortar::n_dof_, 1, double>& local_constraint,
      const unsigned int n_mortar_rot = 0);
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
