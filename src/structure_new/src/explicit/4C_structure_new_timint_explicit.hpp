/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP


#include "4C_config.hpp"

#include "4C_structure_new_expl_generic.hpp"
#include "4C_structure_new_nln_solver_generic.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace TimeInt
  {
    /** \brief Explicit time integration strategy
     *
     * \author Michael Hiermeier */
    class Explicit : public Base
    {
     public:
      //! constructor
      Explicit();

      void setup() override;

      int integrate() override;

      int integrate_step() override;

      void prepare_time_step() override;

      void update_state_incrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void determine_stress_strain() override;

      void evaluate() override;

      void evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void set_state(const Teuchos::RCP<Epetra_Vector>& x) override;

      void reset_step() override;

      Inpar::Solid::ConvergenceStatus solve() override;

      void prepare_partition_step() override;

      void update(double endtime) override;

      void print_step() override;

      Inpar::Solid::StcScale get_stc_algo() override;

      Teuchos::RCP<Core::LinAlg::SparseMatrix> get_stc_mat() override;

      Teuchos::RCP<const Epetra_Vector> initial_guess() override;

      Teuchos::RCP<const Epetra_Vector> get_f() const override;

      Teuchos::RCP<Epetra_Vector> freact() override;

      Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override;

      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

      void use_block_matrix(Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> rangemaps) override;

      ///@}

      //! @name Attribute access functions
      //@{

      enum Inpar::Solid::DynamicType method_name() const override;

      bool is_implicit() const override { return false; }

      bool is_explicit() const override { return true; }

      int method_steps() const override;

      int method_order_of_accuracy_dis() const override;

      int method_order_of_accuracy_vel() const override;

      double method_lin_err_coeff_dis() const override;

      double method_lin_err_coeff_vel() const override;

      //@}

     protected:
      Solid::EXPLICIT::Generic& expl_int()
      {
        check_init_setup();
        return *explint_ptr_;
      };

      Solid::Nln::SOLVER::Generic& nln_solver()
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

     private:
      //! ptr to the explicit time integrator object
      Teuchos::RCP<Solid::EXPLICIT::Generic> explint_ptr_;

      //! ptr to the non-linear solver object
      Teuchos::RCP<Solid::Nln::SOLVER::Generic> nlnsolver_ptr_;
    };
  }  // namespace TimeInt
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
