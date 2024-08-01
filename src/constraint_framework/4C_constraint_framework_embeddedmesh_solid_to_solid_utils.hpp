/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for solid-to-solid interactions.

\level 3

*/

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_SOLID_TO_SOLID_UTILS_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_SOLID_TO_SOLID_UTILS_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_embeddedmesh_params.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

class Epetra_FEVector;
class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  class CutWizard;
  namespace CUT
  {
    class BoundaryCell;
    class Element;
    class Mesh;
    class Point;
    class Side;
    class VolumeCell;
  }  // namespace CUT
}  // namespace Core::Geo

namespace Core::LinAlg
{
  template <unsigned int rows, unsigned int cols, class ValueType>
  class Matrix;
}  // namespace Core::LinAlg

namespace DRT
{
  class Node;
  class Element;
  class Discretization;
}  // namespace DRT

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace CONSTRAINTS::EMBEDDEDMESH
{
  class SolidToSolidMortarManager;
  class SolidInteractionPair;

  /**
   * \brief Free function that prepares and performs the cut.
   * @param cutwizard (in) object of the cut library that performs the cut operation.
   * @param discret (in) Discretization
   */
  void prepare_and_perform_cut(Teuchos::RCP<Core::Geo::CutWizard> cutwizard,
      Teuchos::RCP<Core::FE::Discretization>& discret,
      CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams& embedded_mesh_coupling_params);

  /**
   * \brief Free function to get coupling pairs and background elements
   * @param cutwizard (in) object of the cut library that performs the cut operation.
   * @param embedded_mesh_params_ptr (in) pointer to the embeddedmesh parameters
   * @param discret (in) solid discretization
   * @param embedded_mesh_solid_interaction_pairs (out) embedded mesh coupling pairs
   * @param cut_elements_vector (out) vector of cut elements
   */
  void get_coupling_pairs_and_background_elements(Teuchos::RCP<Core::Geo::CutWizard>& cutwizard,
      CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams& params_ptr,
      Teuchos::RCP<Core::FE::Discretization>& discret,
      std::vector<Teuchos::RCP<CONSTRAINTS::EMBEDDEDMESH::SolidInteractionPair>>&
          embeddedmesh_coupling_pairs,
      std::vector<Core::Elements::Element*>& cut_elements_vector);

  /**
   * \brief Change integration rule of cut background elements
   * @param cut_elements_vector (in) vector of cut elements
   * @param cutwizard (in) object of the cut library that performs the cut operation.
   */
  void change_gauss_rule_of_cut_elements(std::vector<Core::Elements::Element*> cut_elements_vector,
      Teuchos::RCP<Core::Geo::CutWizard>& cutwizard);

  /**
   * \brief Get the number of Lagrange multiplier values corresponding to the solid nodes and
   * solid element.
   * @param shape_function (in) Mortar shape function.
   * @param n_lambda_node (out) Number of Lagrange multiplicators per node.
   */
  void mortar_shape_functions_to_number_of_lagrange_values(
      const Inpar::CONSTRAINTS::SolidToSolidMortarShapefunctions shape_function,
      unsigned int& n_lambda_node);

  /**
   * \brief Assemble local mortar contributions from the classical mortar matrices D and M into the
   * global matrices.
   *
   * This function assumes that the mortar contributions are symmetric, i.e. global_g_b =
   * global_fb_l^T and global_g_s = global_fs_l^T.
   *
   * @param pair (in) The beam-to-solid pair.
   * @param discret (in) Discretization
   * @param mortar_manager (in) Mortar manager for the solid-to-solid condition
   * @param global_g_bl (in/out) Constraint equations derived w.r.t the interface (from the boundary
   * layer) DOFs
   * @param global_g_bg (in/out) Constraint equations derived w.r.t the background solid DOFs
   * @param global_fbl_l (in/out) Interface force vector derived w.r.t the Lagrange multipliers
   * @param global_fbg_l (in/out) Background force vector derived w.r.t the Lagrange multipliers
   * @param global_constraint (in/out) Global constraint equations
   * @param global_kappa (in/out) Global penalty scaling vector equations
   * @param global_lambda_active (in/out) Global vector keeping track of active lagrange multipliers
   * @param local_D (in) Local D matrix of the pair.
   * @param local_M (in) Local M matrix of the pair.
   * @param local_kappa (in) Local scaling vector of the pair.
   * @param local_constraint (in) Local constraint contributions of the pair.
   */
  template <typename Interface, typename Background, typename Mortar>
  void assemble_local_mortar_contributions(
      const CONSTRAINTS::EMBEDDEDMESH::SolidInteractionPair* pair,
      const Core::FE::Discretization& discret,
      const CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager* mortar_manager,
      Core::LinAlg::SparseMatrix& global_g_bl, Core::LinAlg::SparseMatrix& global_g_bg,
      Core::LinAlg::SparseMatrix& global_fbl_l, Core::LinAlg::SparseMatrix& global_fbg_l,
      Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
      Epetra_FEVector& global_lambda_active,
      const Core::LinAlg::Matrix<Mortar::n_dof_, Interface::n_dof_, double>& local_D,
      const Core::LinAlg::Matrix<Mortar::n_dof_, Background::n_dof_, double>& local_M,
      const Core::LinAlg::Matrix<Mortar::n_dof_, 1, double>& local_kappa,
      const Core::LinAlg::Matrix<Mortar::n_dof_, 1, double>& local_constraint);

  /**
   * \brief Get the GIDs of the Lagrange multiplicator unknowns for a beam-to-solid pair.
   * @param mortar_manager (in) Mortar manager for the interface-to-background condition
   * @param interaction_pair (in) interface-to-background interaction pair
   * @param n_mortar_pos (in) Number of positional mortar DOFs associated with the pair
   * @param lambda_gid_pos (out) GIDs of positional mortar DOFs associated with the pair
   */
  void get_mortar_gid(const CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager* mortar_manager,
      const CONSTRAINTS::EMBEDDEDMESH::SolidInteractionPair* interaction_pair,
      const unsigned int n_mortar_pos, std::vector<int>* lambda_gid_pos);

  bool is_interface_node(Core::Nodes::Node const& node);

  bool is_interface_element_surface(Teuchos::RCP<Core::Elements::Element> ele);

  void get_current_element_displacement(Core::FE::Discretization const& discret,
      Core::Elements::Element const* ele, const Epetra_Vector& displacement_vector,
      std::vector<double>& eledisp);

  Core::FE::GaussIntegration create_gauss_integration_from_collection(
      std::vector<Core::FE::GaussIntegration>& intpoints_vector);

  /**
   * \brief Returns the shape function for the mortar Lagrange multipliers.
   */
  Inpar::CONSTRAINTS::SolidToSolidMortarShapefunctions define_shape_functions_lagrange_multipliers(
      Core::FE::CellType celltype);

}  // namespace CONSTRAINTS::EMBEDDEDMESH

FOUR_C_NAMESPACE_CLOSE

#endif