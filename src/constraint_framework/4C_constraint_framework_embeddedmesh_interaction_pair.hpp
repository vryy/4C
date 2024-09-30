/*----------------------------------------------------------------------*/
/*! \file

\brief one generic solid-to-solid interaction element pair

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_INTERACTION_PAIR_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_INTERACTION_PAIR_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_embeddedmesh_params.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_io_visualization_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

#include <unordered_set>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Cut
{
  class CutWizard;
  class BoundaryCell;
}  // namespace Cut

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT
namespace GEOMETRYPAIR
{
  class GeometryPair;
}  // namespace GEOMETRYPAIR

namespace CONSTRAINTS::EMBEDDEDMESH
{
  class SolidToSolidMortarManager;

  class SolidInteractionPair
  {
   public:
    /*!
    \brief Constructor
    */
    SolidInteractionPair(Teuchos::RCP<Core::Elements::Element> element1,
        Core::Elements::Element* element2,
        CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams& params_ptr,
        Teuchos::RCP<Cut::CutWizard> cutwizard_ptr,
        std::vector<Teuchos::RCP<Cut::BoundaryCell>>& boundary_cells);

    /*!
    \brief Destructor
    */
    virtual ~SolidInteractionPair() = default;

    //! @name Access methods
    inline CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams params() const { return params_; }

    /*!
    \brief Get first element
    */
    const Core::Elements::Element& element_1() const { return *element1_.get(); };

    /*!
    \brief Get second element
    */
    const Core::Elements::Element& element_2() const { return *element2_; };

    /*!
    \brief Get the number of segments
    */
    inline unsigned int get_num_segments() { return boundary_cells_.size(); }

    //! @name Visualization methods
    /*!
    \brief Get the Gauss points of element2_ after cut for visualization.
    */
    virtual void get_projected_gauss_rule_in_cut_element(
        const Core::IO::VisualizationData& cut_element_integration_points_visualization_data) = 0;

    /*!
    \brief Get the Gauss points of element1_ and element2_ for the evaluation of mortar matrices.
    */
    virtual void get_projected_gauss_rule_on_interface(int num_segment,
        const Core::IO::VisualizationData& background_integration_points_visualization_data,
        const Core::IO::VisualizationData& interface_integration_points_visualization_data) = 0;

    /*!
    \brief Get the Lagrange multiplier field evaluated on the interface nodes for visualization.
    */
    virtual void get_pair_visualization(
        const Core::IO::VisualizationData& lagrange_multipliers_visualization_data,
        Teuchos::RCP<Core::LinAlg::Vector<double>> lambda,
        const CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager* mortar_manager,
        Teuchos::RCP<std::unordered_set<int>> interface_tracker) = 0;


    //! @name Evaluation methods
    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling.
     * @param discret (in) Discretization, used to get the boundary layer GIDs.
     * @param mortar_manager (in) Mortar manager, used to get the Lagrange multiplier GIDs.
     * @param global_g_bl_ (in/out) Constraint equations derived w.r.t the boundary layer DOFs.
     * @param global_g_bg_ (in/out) Constraint equations derived w.r.t the background mesh DOFs.
     * @param global_fbl_l_ (in/out) Boundary layer force vector derived w.r.t the Lagrange
     * multipliers.
     * @param global_fbg_l_ (in/out) Background mesh force vector derived w.r.t the Lagrange
     * multipliers.
     * @param global_constraint (in/out) Global constraint vector.
     * @param global_kappa (in/out) Global scaling matrix.
     * @param global_lambda_active (in/out) Global vector with active Lagrange multipliers.
     * @param displacement_vector (in) Global displacement vector.
     */
    virtual void evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization& discret,
        const CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager* mortar_manager,
        Core::LinAlg::SparseMatrix& global_g_bl_, Core::LinAlg::SparseMatrix& global_g_bg_,
        Core::LinAlg::SparseMatrix& global_fbl_l_, Core::LinAlg::SparseMatrix& global_fbg_l_,
        Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
        Epetra_FEVector& global_lambda_active) = 0;

    /**
     * \brief Set the current element displacement.
     * @param discret (in) Discretization, used to get the boundary layer GIDs.
     * @param dispnp (in) Displacement at \f$t_{n}\f$
     */
    virtual void set_current_element_position(Core::FE::Discretization const& discret,
        const Core::LinAlg::Vector<double>& displacement_vector) = 0;

    //! Get the cutwizard
    Teuchos::RCP<Cut::CutWizard> get_cutwizard()
    {
      if (cutwizard_ptr_ == Teuchos::null) FOUR_C_THROW("The cut wizard hasn't been initialized!");
      return cutwizard_ptr_;
    }

    //! Get the boundary cells
    std::vector<Teuchos::RCP<Cut::BoundaryCell>> get_boundary_cells() { return boundary_cells_; }

   protected:
    //! embedded mesh parameter data container
    CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams params_;

   private:
    //! @name member variables
    //! first element of interacting pair
    const Teuchos::RCP<Core::Elements::Element> element1_;

    //! second element of interacting pair
    const Core::Elements::Element* element2_;

    //! pointer to the cutwizard
    Teuchos::RCP<Cut::CutWizard> cutwizard_ptr_ = Teuchos::null;

    //! boundary cells that are related to this coupling pair
    std::vector<Teuchos::RCP<Cut::BoundaryCell>> boundary_cells_;
  };
}  // namespace CONSTRAINTS::EMBEDDEDMESH

FOUR_C_NAMESPACE_CLOSE

#endif