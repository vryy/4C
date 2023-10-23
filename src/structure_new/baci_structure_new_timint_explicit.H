/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_EXPLICIT_H
#define BACI_STRUCTURE_NEW_TIMINT_EXPLICIT_H


#include "baci_structure_new_expl_generic.H"
#include "baci_structure_new_nln_solver_generic.H"
#include "baci_structure_new_timint_base.H"

namespace STR
{
  namespace TIMINT
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

      void PrepareTimeStep() override;

      void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void DetermineStressStrain() override;

      void Evaluate() override;

      void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void SetState(const Teuchos::RCP<Epetra_Vector>& x) override;

      void ResetStep() override;

      INPAR::STR::ConvergenceStatus Solve() override;

      void PreparePartitionStep() override;

      void Update(double endtime) override;

      void PrintStep() override;

      INPAR::STR::STC_Scale GetSTCAlgo() override;

      Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override;

      Teuchos::RCP<const Epetra_Vector> InitialGuess() override;

      Teuchos::RCP<const Epetra_Vector> GetF() const override;

      Teuchos::RCP<Epetra_Vector> Freact() override;

      Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

      Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

      void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override;

      ///@}

     protected:
      STR::EXPLICIT::Generic& ExplInt()
      {
        CheckInitSetup();
        return *explint_ptr_;
      };

      STR::NLN::SOLVER::Generic& NlnSolver()
      {
        CheckInitSetup();
        return *nlnsolver_ptr_;
      };

     private:
      //! ptr to the explicit time integrator object
      Teuchos::RCP<STR::EXPLICIT::Generic> explint_ptr_;

      //! ptr to the non-linear solver object
      Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_ptr_;
    };
  }  // namespace TIMINT
}  // namespace STR

#endif  // STRUCTURE_NEW_TIMINT_EXPLICIT_H
