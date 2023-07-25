/*-----------------------------------------------------------*/
/*! \file

\brief Implicit structural time integration strategy.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_IMPLICIT_H
#define BACI_STRUCTURE_NEW_TIMINT_IMPLICIT_H

#include "baci_structure_new_timint_implicitbase.H"  // base class

// forward declarations ...
namespace NOX::NLN::INNER::StatusTest
{
  class Generic;
}  // namespace NOX::NLN::INNER::StatusTest
namespace STR
{
  namespace IMPLICIT
  {
    class Generic;
  }  // namespace IMPLICIT
  namespace PREDICT
  {
    class Generic;
  }  // namespace PREDICT
  namespace NLN::SOLVER
  {
    class Generic;
    namespace INTERFACE
    {
      class Required;
    }  // namespace INTERFACE
  }    // namespace NLN::SOLVER

  namespace TIMINT
  {
    /** \brief Implicit time integration strategy
     *
     * \author Michael Hiermeier */
    class Implicit : public ImplicitBase
    {
     public:
      void Setup() override;

      int Integrate() override;

      int IntegrateStep() override;

      /// set the state of the nox group and the global state data container
      /// see class \ref ADAPTER::StructureNew for a detailed documentation.
      void SetState(const Teuchos::RCP<Epetra_Vector>& x) override;

      /*! \brief nonlinear solve
       *
       *  Do the nonlinear solve, i.e. (multiple) corrector,
       *  for the time step. All boundary conditions have
       *  been set. */
      INPAR::STR::ConvergenceStatus Solve() override;

      /** \brief Identify residual
       *
       *  This method does not predict the target solution but
       *  evaluates the residual and the stiffness matrix.
       *  In partitioned solution schemes, it is better to keep the current
       *  solution instead of evaluating the initial guess (as the predictor)
       *  does. */
      void PreparePartitionStep() override;

      //! Prepare time step
      void PrepareTimeStep() override;

      //! @name Accessors
      //! @{
      //! return the predictor
      [[nodiscard]] const STR::PREDICT::Generic& Predictor() const
      {
        CheckInitSetup();
        return *predictor_ptr_;
      }

      //! @}
      [[nodiscard]] Teuchos::RCP<const STR::NLN::SOLVER::Generic> GetNlnSolverPtr() const
      {
        return nlnsolver_ptr_;
      };

      //! do something in case nonlinear solution does not converge for some reason
      INPAR::STR::ConvergenceStatus PerformErrorAction(
          INPAR::STR::ConvergenceStatus nonlinsoldiv) override;


      //! check, if according to divercont flag time step size can be increased
      void CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status);

      //! returns pointer to generic implicit object
      const Teuchos::RCP<STR::IMPLICIT::Generic> ImplIntPtr()
      {
        CheckInitSetup();
        return implint_ptr_;
      };

      /// Update State Incrementally for coupled problems with monolithic approach
      void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void DetermineStressStrain() override;

      ///  Evaluate routine for coupled problems with monolithic approach
      void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;
      void Evaluate() override;

      /** \brief Print structural jacobian matrix into a text file for later use
       *  in MATLAB
       *
       *  This routine can be activated via the input parameter
       *  %STRUCT_JACOBIAN_MATLAB. See the corresponding inpar section for more
       *  details.
       *
       *  The text file can be found in the user-provided output folder using the
       *  following file name extension:
       *
       *  [OUTPUT-FOLDER]/[OUTPUT FILE NAME]_str_jacobian_step-[STEP]_nlniter-[NEWTON-ITERATION].mtl
       *
       *  \author hiermeier \date 06/17 */
      void PrintJacobianInMatlabFormat(const NOX::NLN::Group& curr_grp) const;

      /// compute the condition number of the tangential stiffness matrix
      void ComputeConditionNumber(const NOX::NLN::Group& grp) const;

     protected:
      //! returns the current solution group
      [[nodiscard]] const NOX::Abstract::Group& GetSolutionGroup() const override;

      //! returns the current solution group ptr
      Teuchos::RCP<NOX::Abstract::Group> SolutionGroupPtr() override;

      STR::IMPLICIT::Generic& ImplInt()
      {
        CheckInitSetup();
        return *implint_ptr_;
      };

      STR::PREDICT::Generic& Predictor()
      {
        CheckInitSetup();
        return *predictor_ptr_;
      };

      const Teuchos::RCP<STR::PREDICT::Generic> PredictorPtr()
      {
        CheckInitSetup();
        return predictor_ptr_;
      };

      [[nodiscard]] const STR::NLN::SOLVER::Generic& NlnSolver() const
      {
        CheckInitSetup();
        return *nlnsolver_ptr_;
      };

      STR::NLN::SOLVER::Generic& NlnSolver()
      {
        CheckInitSetup();
        return *nlnsolver_ptr_;
      };

      const Teuchos::RCP<STR::NLN::SOLVER::Generic> NlnSolverPtr()
      {
        CheckInitSetup();
        return nlnsolver_ptr_;
      };



     private:
      //! ptr to the implicit time integrator object
      Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

      //! ptr to the predictor object
      Teuchos::RCP<STR::PREDICT::Generic> predictor_ptr_;

      //! ptr to the non-linear solver object
      Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_ptr_;

      //! ptr to the nox group object
      Teuchos::RCP<NOX::Abstract::Group> grp_ptr_;
    };
  }  // namespace TIMINT
}  // namespace STR


#endif  // STRUCTURE_NEW_TIMINT_IMPLICIT_H
