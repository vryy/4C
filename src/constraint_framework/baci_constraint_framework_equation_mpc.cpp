/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all multipoint constraint equation terms

\level 3
    */
/*-----------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_constraint_framework_equation_mpc.hpp"

#include "baci_io_pstream.hpp"
#include "baci_lib_discret.hpp"

BACI_NAMESPACE_OPEN
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::LinearCoupledEquation::EvaluateEquation(
    CORE::LINALG::SparseMatrix& Q_dd, CORE::LINALG::SparseMatrix& Q_dL,
    CORE::LINALG::SparseMatrix& Q_Ld, Epetra_Vector& constraint_vector, const Epetra_Vector& D_np1)
{
  double constraintViolation = 0.;

  // Iterate over the elements (coefficient, rowId, dofId) in equationData.
  // Each element of equation data represents one term of the defined multipoint constraints
  // The rowId is equivalent to the Number of the equation
  for (const auto& [coefficient, rowId, dofId] : equationData_)
  {
    // stiffness contribution
    // FixMe: Switch to FEAssemble
    Q_dL.Assemble(coefficient, dofId, rowId);
    Q_Ld.Assemble(coefficient, rowId, dofId);

    // force contribution
    constraintViolation = D_np1.Values()[dofId] * coefficient;
    constraint_vector.SumIntoGlobalValues(1, &constraintViolation, &rowId);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::GetNumberOfMPCs() const
{
  return nDofCoupled_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::GetFirstRowId() const
{
  return firstRowId_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::MultiPointConstraintEquationBase::SetFirstRowId(
    int global_row_id)
{
  firstRowId_ = global_row_id;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONSTRAINTS::SUBMODELEVALUATOR::LinearCoupledEquation::LinearCoupledEquation(
    int id, const std::vector<int>& dofs, std::vector<double> coefficients)
{
  IO::cout(IO::debug) << "\nLinear coupled equation saved (ID: " << id << ")\n ";
  IO::cout(IO::debug) << " 0 = ";  //#Todo

  SetFirstRowId(id);

  for (std::vector<double>::size_type i = 0; i < coefficients.size(); ++i)
  {
    TermData term = {coefficients[i], id, dofs[i]};
    equationData_.emplace_back(term);

    IO::cout(IO::debug) << " + " << coefficients[i] << " * d" << dofs[i];
  }
  IO::cout(IO::debug) << IO::endl;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE