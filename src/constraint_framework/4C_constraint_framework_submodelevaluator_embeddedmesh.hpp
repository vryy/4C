#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_EMBEDDEDMESH_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_EMBEDDEDMESH_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_submodelevaluator_base.hpp"
#include "4C_fem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONSTRAINTS::EMBEDDEDMESH
{
  class SolidToSolidMortarManager;
}

namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  class EmbeddedMeshConstraintManager : public ConstraintBase
  {
   public:
    /*!
    \brief Standard Constructor
    */
    EmbeddedMeshConstraintManager(Teuchos::RCP<Core::FE::Discretization> discret_ptr,
        const Core::LinAlg::Vector<double>& dispnp);

    //! @name Public evaluation methods

    /*!
     * \brief Reset the constraint stiffness matrix and delete node pairs
     */
    void reset() override
    {
      // Nothing implemented
    }

    /*! Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
     */
    bool evaluate_force_stiff(const Core::LinAlg::Vector<double>& displacement_vector,
        Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& global_state_ptr,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> me_stiff_ptr,
        Teuchos::RCP<Core::LinAlg::Vector<double>> me_force_ptr) override;

    //! Evaluate the matrices of the saddle-point system
    void evaluate_coupling_terms(Solid::TimeInt::BaseDataGlobalState& gstate) override
    {
      // Nothing implemented
    }

    //! derived
    void runtime_output_step_state(std::pair<double, int> output_time_and_step) override;

    //@}

   private:
    //! Pointer to the mortar manager. This object stores the relevant mortar matrices.
    Teuchos::RCP<CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager> mortar_manager_;
  };
}  // namespace CONSTRAINTS::SUBMODELEVALUATOR

FOUR_C_NAMESPACE_CLOSE

#endif  // FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_EMBEDDEDMESH_HPP
