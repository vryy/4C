/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all meshtying terms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MESHTYING_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MESHTYING_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CONTACT
{
  class Manager;
  class MtAbstractStrategy;
}  // namespace CONTACT

namespace Mortar
{
  class StrategyBase;
}  // namespace Mortar

namespace STR
{
  namespace MODELEVALUATOR
  {
    class MeshtyingData;

    /*! \brief Model evaluator for meshtying problems
     *
     */
    class Meshtying : public Generic
    {
     public:
      //! constructor
      Meshtying();


      /*! \brief Initialize class variables [derived]
       *
       * @param eval_data_ptr
       * @param gstate_ptr
       * @param gio_ptr
       * @param int_ptr
       * @param timint_ptr
       * @param dof_offset
       */
      void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TimeInt::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TimeInt::Base>& timint_ptr, const int& dof_offset) override;

      //! setup class variables [derived]
      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //!@{

      //! [derived]
      Inpar::STR::ModelType Type() const override { return Inpar::STR::model_meshtying; }

      //! [derived]
      void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) override;

      //! [derived]
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      //! [derived]
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! [derived]
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const Inpar::STR::PredEnum& pred_type) override{};

      //! [derived]
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! [derived]
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override{};

      //! [derived]
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override{};

      //! [derived]
      void run_post_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::Nln::Group& grp) override;

      //! [derived]
      void run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::Nln::Group& grp) override;

      //! [derived]
      void update_step_state(const double& timefac_n) override{};

      //! [derived]
      void update_step_element() override{};

      //! [derived]
      void determine_stress_strain() override{};

      //! [derived]
      void determine_energy() override{};

      //! [derived]
      void determine_optional_quantity() override{};

      //! [derived]
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override{};

      //! [derived]
      void reset_step_state() override{};

      //! [derived]
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! [derived]
      void post_output() override{};

      /*! \brief Reset model specific variables (without jacobian) [derived]
       *
       * Nothing to do in case of meshtying.
       *
       * \param[in] x Current full state vector
       */
      void Reset(const Epetra_Vector& x) override{};

      //! \brief Perform actions just before the Evaluate() call [derived]
      void pre_evaluate() override{};

      //! \brief Perform actions right after the Evaluate() call [derived]
      void post_evaluate() override{};

      //! @}

      //! @name Call-back routines
      //!@{

      Teuchos::RCP<const Core::LinAlg::SparseMatrix> GetJacobianBlock(
          const STR::MatBlockType bt) const;

      /** \brief Assemble the structural right-hand side vector
       *
       *  \param[in] without_these_models  Exclude all models defined in this vector
       *                                   during the assembly
       *  \param[in] apply_dbc             Apply Dirichlet boundary conditions
       *
       *  \author hiermeier \date 08/17 */
      Teuchos::RCP<Epetra_Vector> assemble_force_of_models(
          const std::vector<Inpar::STR::ModelType>* without_these_models = nullptr,
          const bool apply_dbc = false) const;

      virtual Teuchos::RCP<Core::LinAlg::SparseOperator> GetAuxDisplJacobian() const
      {
        return Teuchos::null;
      };

      void evaluate_weighted_gap_gradient_error();

      //! [derived]
      bool evaluate_force() override;

      //! [derived]
      bool evaluate_stiff() override;

      //! [derived]
      bool evaluate_force_stiff() override;

      /*!
      \brief Apply results of mesh initialization to the underlying problem discretization

      \note This is only necessary in case of a mortar method.

      \warning This routine modifies the reference coordinates of slave nodes at the meshtying
      interface.

      @param[in] Xslavemod Vector with modified nodal positions
      */
      void apply_mesh_initialization(Teuchos::RCP<const Epetra_Vector> Xslavemod);

      //!@}

      //! @name Accessors
      //!@{

      //! Returns a pointer to the underlying meshtying strategy object
      const Teuchos::RCP<CONTACT::MtAbstractStrategy>& StrategyPtr();

      //! Returns the underlying meshtying strategy object
      CONTACT::MtAbstractStrategy& Strategy();
      const CONTACT::MtAbstractStrategy& Strategy() const;

      //!@}

     protected:
     private:
      /// Set the correct time integration parameters within the meshtying strategy
      void set_time_integration_info(CONTACT::MtAbstractStrategy& strategy) const;

      //! meshtying strategy
      Teuchos::RCP<CONTACT::MtAbstractStrategy> strategy_ptr_;

      //! Mesh relocation for conservation of angular momentum
      Teuchos::RCP<Epetra_Vector> mesh_relocation_;
    };  // namespace MODELEVALUATOR

  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
