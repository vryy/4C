/*----------------------------------------------------------------------*/
/*! \file
  \brief Lagrange multiplier function: solve a least squares problem to compute
  the Lagrange multiplier value dependent on the current displacement state

  \level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_HPP
#define FOUR_C_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_utils.hpp"
#include "4C_inpar_solver.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace CONTACT
{
  class ParamsInterface;

  namespace Aug
  {
    class DataContainer;
    class Strategy;

    class LagrangeMultiplierFunction
    {
     public:
      /// constructor
      LagrangeMultiplierFunction();


      void Init(const Strategy* const strategy, CONTACT::Aug::DataContainer& data);

      void Setup();

      void Redistribute();

      Teuchos::RCP<Epetra_Vector> Compute(const CONTACT::ParamsInterface& cparams);

      Teuchos::RCP<Epetra_Vector> FirstOrderDirection(
          const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dincr);

     private:
      Teuchos::RCP<Core::LinAlg::Solver> create_linear_solver(const int lin_sol_id,
          const Epetra_Comm& comm, enum Core::LinearSolver::SolverType& solver_type) const;

      void lin_solve(
          Core::LinAlg::SparseOperator& mat, Epetra_MultiVector& rhs, Epetra_MultiVector& sol);

      inline void check_init() const
      {
        if (not isinit_) FOUR_C_THROW("Call Init() first!");
      }

      inline void check_init_setup() const
      {
        check_init();
        if (not issetup_) FOUR_C_THROW("Call Setup() first!");
      }

      Teuchos::RCP<Epetra_Vector> get_structure_gradient(
          const CONTACT::ParamsInterface& cparams) const;

      void create_b_matrix();

      void assemble_gradient_b_matrix_contribution(
          const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const;

      void assemble_gradient_bb_matrix_contribution(
          const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const;

     private:
      bool isinit_;
      bool issetup_;

      const Strategy* strategy_;

      plain_interface_set interfaces_;

      Teuchos::RCP<DataContainer> data_;

      enum Core::LinearSolver::SolverType lin_solver_type_;

      Teuchos::RCP<Core::LinAlg::Solver> lin_solver_;

      // B-matrix
      Teuchos::RCP<Core::LinAlg::SparseMatrix> bmat_;

    };  // class LagrangeMultiplierFunction
  }     // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
