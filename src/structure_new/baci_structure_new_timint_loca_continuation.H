/*-----------------------------------------------------------*/
/*! \file

\brief Library of Continuation Algorithms (LOCA) implementation.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_LOCA_CONTINUATION_H
#define BACI_STRUCTURE_NEW_TIMINT_LOCA_CONTINUATION_H


#include "baci_loca_continuation_nln_stepper.H"      // enum
#include "baci_solver_nonlin_nox_enum_lists.H"       // enum
#include "baci_structure_new_timint_implicitbase.H"  // base class

// forward declaration
namespace LOCA
{
  class ParameterVector;
  namespace MultiContinuation
  {
    class ConstraintInterface;
  }  // namespace MultiContinuation
  namespace Epetra
  {
    class Factory;
  }  // namespace Epetra
}  // namespace LOCA
namespace STR
{
  namespace TIMINT
  {
    /** \brief LOCA continuation strategy
     *
     * \author Michael Hiermeier */
    class LOCAContinuation : public ImplicitBase
    {
     public:
      //! constructor
      LOCAContinuation();

      //! destructor
      ~LOCAContinuation() override = default;

      //! derived setup routine
      void Setup() override;

      //! call of the stepper run() routine (derived)
      int Integrate() override;

      /// integrate the current step (derived)
      int IntegrateStep() override;

      /// Update State Incrementally for coupled problems with monolithic approach
      void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      ///  Evaluate routine for coupled problems with monolithic approach
      void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      /// don't update displacement but evaluate elements (implicit only)
      void Evaluate() override;

      void DetermineStressStrain() override { dserror("Currently unsupported!"); }

      /// set the state of the nox group and the global state data container
      void SetState(const Teuchos::RCP<Epetra_Vector>& x) override;

      /*! \brief nonlinear solve
       *
       *  Do the nonlinear solve, i.e. (multiple) corrector,
       *  for the time step. All boundary conditions have
       *  been set. */
      INPAR::STR::ConvergenceStatus Solve() override;

      //! start new time step (derived)
      void PrepareTimeStep() override;

      //! start new partitioned step
      void PreparePartitionStep() override;


     protected:
      //! Returns the current solution group
      const NOX::Abstract::Group& GetSolutionGroup() const override;

      //! Returns the current solution group ptr
      Teuchos::RCP<NOX::Abstract::Group> SolutionGroupPtr() override;

      //! Check for necessary parameters which have to be supplied by the user.
      void CheckForNecessaryParameters() const;

      //! combine the LOCA and NOX parameter lists
      void CombineParameterLists();

      //! create the LOCA parameter vector (continuation informations)
      Teuchos::RCP<LOCA::ParameterVector> CreateParameterVector(
          const Teuchos::ParameterList& p_loca);

     protected:
      //! parameter list which holds the loca and the nox parameters
      Teuchos::RCP<Teuchos::ParameterList> p_loca_nox_ptr_;

      //! LOCA constraints pointer (e.g. arc-length constraints object)
      Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> loca_constraints_ptr_;

      Teuchos::RCP<std::vector<std::string>> loca_constraint_param_names_ptr_;

      //! LOCA stepper (doing the main work)
      Teuchos::RCP<LOCA::NLN::Stepper> loca_stepper_ptr_;

      //! LOCA stepper status
      LOCA::Abstract::Iterator::IteratorStatus loca_stepper_status_;
    };  // class LOCAContinuation
  }     // namespace TIMINT
}  // namespace STR


#endif  // STRUCTURE_NEW_TIMINT_LOCA_CONTINUATION_H
