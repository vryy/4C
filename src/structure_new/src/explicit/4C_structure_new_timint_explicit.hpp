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

namespace STR
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


      void Setup() override;

      int Integrate() override;

      int IntegrateStep() override;

      void prepare_time_step() override;

      void update_state_incrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void determine_stress_strain() override;

      void Evaluate() override;

      void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void set_state(const Teuchos::RCP<Epetra_Vector>& x) override;

      void reset_step() override;

      Inpar::STR::ConvergenceStatus Solve() override;

      void prepare_partition_step() override;

      void Update(double endtime) override;

      void PrintStep() override;

      Inpar::STR::StcScale GetSTCAlgo() override;

      Teuchos::RCP<Core::LinAlg::SparseMatrix> GetSTCMat() override;

      Teuchos::RCP<const Epetra_Vector> initial_guess() override;

      Teuchos::RCP<const Epetra_Vector> GetF() const override;

      Teuchos::RCP<Epetra_Vector> Freact() override;

      Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix() override;

      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> BlockSystemMatrix() override;

      void use_block_matrix(Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> rangemaps) override;

      ///@}

      //! @name Attribute access functions
      //@{

      enum Inpar::STR::DynamicType MethodName() const override;

      bool IsImplicit() const override { return false; }

      bool IsExplicit() const override { return true; }

      int MethodSteps() const override;

      int method_order_of_accuracy_dis() const override;

      int method_order_of_accuracy_vel() const override;

      double method_lin_err_coeff_dis() const override;

      double method_lin_err_coeff_vel() const override;

      //@}

     protected:
      STR::EXPLICIT::Generic& expl_int()
      {
        check_init_setup();
        return *explint_ptr_;
      };

      STR::Nln::SOLVER::Generic& nln_solver()
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

     private:
      //! ptr to the explicit time integrator object
      Teuchos::RCP<STR::EXPLICIT::Generic> explint_ptr_;

      //! ptr to the non-linear solver object
      Teuchos::RCP<STR::Nln::SOLVER::Generic> nlnsolver_ptr_;
    };
  }  // namespace TimeInt
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
