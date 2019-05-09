/*!
  \brief Lagrange multiplier function: solve a least squares problem to compute
  the Lagrange multiplier value dependent on the current displacement state

  \level 3
  \maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H_

#include <Teuchos_RCP.hpp>
#include "contact_aug_utils.H"

#include "../drt_inpar/inpar_solver.H"
#include "../drt_lib/drt_dserror.H"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Comm;

namespace LINALG
{
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG

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

      /// destructor
      ~LagrangeMultiplierFunction(){};

      void Init(const Strategy* const strategy, CONTACT::AUG::DataContainer& data);

      void Setup();

      void Redistribute();

      Teuchos::RCP<Epetra_Vector> Compute(const CONTACT::ParamsInterface& cparams);

      Teuchos::RCP<Epetra_Vector> FirstOrderDirection(
          const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dincr);

     private:
      Teuchos::RCP<LINALG::Solver> CreateLinearSolver(const int lin_sol_id, const Epetra_Comm& comm,
          enum INPAR::SOLVER::SolverType& solver_type) const;

      void LinSolve(LINALG::SparseOperator& mat, Epetra_MultiVector& rhs, Epetra_MultiVector& sol);

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

      Teuchos::RCP<LINALG::Solver> lin_solver_;

      // B-matrix
      Teuchos::RCP<LINALG::SparseMatrix> bmat_;

    };  // class LagrangeMultiplierFunction
  }     // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_MULTIPLIER_FUNCTION_H_ */
