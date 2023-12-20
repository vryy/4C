/*----------------------------------------------------------------------*/
/*! \file
  \brief Lagrange multiplier function: solve a least squares problem to compute
  the Lagrange multiplier value dependent on the current displacement state

  \level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H
#define BACI_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H

#include "baci_config.H"

#include "baci_contact_aug_utils.H"
#include "baci_inpar_solver.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Comm;

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  class ParamsInterface;

  namespace AUG
  {
    class DataContainer;
    class Strategy;

    class LagrangeMultiplierFunction
    {
     public:
      /// constructor
      LagrangeMultiplierFunction();


      void Init(const Strategy* const strategy, CONTACT::AUG::DataContainer& data);

      void Setup();

      void Redistribute();

      Teuchos::RCP<Epetra_Vector> Compute(const CONTACT::ParamsInterface& cparams);

      Teuchos::RCP<Epetra_Vector> FirstOrderDirection(
          const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dincr);

     private:
      Teuchos::RCP<CORE::LINALG::Solver> CreateLinearSolver(const int lin_sol_id,
          const Epetra_Comm& comm, enum INPAR::SOLVER::SolverType& solver_type) const;

      void LinSolve(
          CORE::LINALG::SparseOperator& mat, Epetra_MultiVector& rhs, Epetra_MultiVector& sol);

      inline void CheckInit() const
      {
        if (not isinit_) dserror("Call Init() first!");
      }

      inline void CheckInitSetup() const
      {
        CheckInit();
        if (not issetup_) dserror("Call Setup() first!");
      }

      Teuchos::RCP<Epetra_Vector> GetStructureGradient(
          const CONTACT::ParamsInterface& cparams) const;

      void CreateBMatrix();

      void AssembleGradientBMatrixContribution(
          const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const;

      void AssembleGradientBBMatrixContribution(
          const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const;

     private:
      bool isinit_;
      bool issetup_;

      const Strategy* strategy_;

      plain_interface_set interfaces_;

      Teuchos::RCP<DataContainer> data_;

      enum INPAR::SOLVER::SolverType lin_solver_type_;

      Teuchos::RCP<CORE::LINALG::Solver> lin_solver_;

      // B-matrix
      Teuchos::RCP<CORE::LINALG::SparseMatrix> bmat_;

    };  // class LagrangeMultiplierFunction
  }     // namespace AUG
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H
