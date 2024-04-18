/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_EXPLICIT_HPP


#include "baci_config.hpp"

#include "baci_structure_new_expl_generic.hpp"
#include "baci_structure_new_nln_solver_generic.hpp"
#include "baci_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

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

      INPAR::STR::StcScale GetSTCAlgo() override;

      Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override;

      Teuchos::RCP<const Epetra_Vector> InitialGuess() override;

      Teuchos::RCP<const Epetra_Vector> GetF() const override;

      Teuchos::RCP<Epetra_Vector> Freact() override;

      Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

      Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

      void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override;

      ///@}

      //! @name Attribute access functions
      //@{

      enum INPAR::STR::DynamicType MethodName() const override;

      bool IsImplicit() const override { return false; }

      bool IsExplicit() const override { return true; }

      int MethodSteps() const override;

      int MethodOrderOfAccuracyDis() const override;

      int MethodOrderOfAccuracyVel() const override;

      double MethodLinErrCoeffDis() const override;

      double MethodLinErrCoeffVel() const override;

      //@}

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

FOUR_C_NAMESPACE_CLOSE

#endif
