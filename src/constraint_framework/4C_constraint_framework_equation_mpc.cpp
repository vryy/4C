/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all multipoint constraint equation terms

\level 3
    */
/*-----------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_constraint_framework_equation_mpc.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::LinearCoupledEquation::evaluate_equation(
    Core::LinAlg::SparseMatrix& Q_dd, Core::LinAlg::SparseMatrix& Q_dL,
    Core::LinAlg::SparseMatrix& Q_Ld, Epetra_Vector& constraint_vector, const Epetra_Vector& D_np1)
{
  double constraintViolation = 0.;

  // Iterate over the elements (coefficient, rowId, dofId) in equationData.
  // Each element of equation data represents one term of the defined multipoint constraints
  // The rowId is equivalent to the Number of the equation
  for (const auto& [coefficient, rowId, dofId] : equation_data_)
  {
    // stiffness contribution
    // FixMe: Switch to FEAssemble
    Q_dL.assemble(coefficient, dofId, rowId);
    Q_Ld.assemble(coefficient, rowId, dofId);

    // force contribution
    constraintViolation = D_np1.Values()[dofId] * coefficient;
    constraint_vector.SumIntoGlobalValues(1, &constraintViolation, &rowId);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::get_number_of_mp_cs() const
{
  return n_dof_coupled_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::get_first_row_id() const
{
  return first_row_id_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::set_first_row_id(
    int global_row_id)
{
  first_row_id_ = global_row_id;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONSTRAINTS::SUBMODELEVALUATOR::LinearCoupledEquation::LinearCoupledEquation(
    int id, const std::vector<int>& dofs, std::vector<double> coefficients)
{
  Core::IO::cout(Core::IO::debug) << "\nLinear coupled equation saved (ID: " << id << ")\n ";
  Core::IO::cout(Core::IO::debug) << " 0 = ";  // #Todo

  set_first_row_id(id);

  for (std::vector<double>::size_type i = 0; i < coefficients.size(); ++i)
  {
    TermData term = {coefficients[i], id, dofs[i]};
    equation_data_.emplace_back(term);

    Core::IO::cout(Core::IO::debug) << " + " << coefficients[i] << " * d" << dofs[i];
  }
  Core::IO::cout(Core::IO::debug) << Core::IO::endl;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE