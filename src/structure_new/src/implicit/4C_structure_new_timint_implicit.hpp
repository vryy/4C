/*-----------------------------------------------------------*/
/*! \file

\brief Implicit structural time integration strategy.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICIT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICIT_HPP

#include "4C_config.hpp"

#include "4C_structure_new_timint_implicitbase.hpp"  // base class

FOUR_C_NAMESPACE_OPEN

// forward declarations ...
namespace STR
{
  namespace IMPLICIT
  {
    class Generic;
  }  // namespace IMPLICIT
  namespace Predict
  {
    class Generic;
  }  // namespace Predict
  namespace Nln::SOLVER
  {
    class Generic;
    namespace INTERFACE
    {
      class Required;
    }  // namespace INTERFACE
  }    // namespace Nln::SOLVER

  namespace TimeInt
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
      /// see class \ref Adapter::StructureNew for a detailed documentation.
      void set_state(const Teuchos::RCP<Epetra_Vector>& x) override;

      /*! \brief nonlinear solve
       *
       *  Do the nonlinear solve, i.e. (multiple) corrector,
       *  for the time step. All boundary conditions have
       *  been set. */
      Inpar::STR::ConvergenceStatus Solve() override;

      /** \brief Identify residual
       *
       *  This method does not predict the target solution but
       *  evaluates the residual and the stiffness matrix.
       *  In partitioned solution schemes, it is better to keep the current
       *  solution instead of evaluating the initial guess (as the predictor)
       *  does. */
      void prepare_partition_step() override;

      //! Prepare time step
      void prepare_time_step() override;

      //! @name Accessors
      //! @{
      //! return the predictor
      [[nodiscard]] const STR::Predict::Generic& Predictor() const
      {
        check_init_setup();
        return *predictor_ptr_;
      }

      //! @}
      [[nodiscard]] Teuchos::RCP<const STR::Nln::SOLVER::Generic> get_nln_solver_ptr() const
      {
        return nlnsolver_ptr_;
      };

      //! do something in case nonlinear solution does not converge for some reason
      Inpar::STR::ConvergenceStatus PerformErrorAction(
          Inpar::STR::ConvergenceStatus nonlinsoldiv) override;

      //! check, if according to divercont flag time step size can be increased
      void check_for_time_step_increase(Inpar::STR::ConvergenceStatus& status);

      //! returns pointer to generic implicit object
      Teuchos::RCP<STR::IMPLICIT::Generic> impl_int_ptr()
      {
        check_init_setup();
        return implint_ptr_;
      };

      /// Update State Incrementally for coupled problems with monolithic approach
      void update_state_incrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      void determine_stress_strain() override;

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
      void print_jacobian_in_matlab_format(const NOX::Nln::Group& curr_grp) const;

      /// compute the condition number of the tangential stiffness matrix
      void compute_condition_number(const NOX::Nln::Group& grp) const;

     protected:
      //! returns the current solution group
      [[nodiscard]] const ::NOX::Abstract::Group& get_solution_group() const override;

      //! returns the current solution group ptr
      Teuchos::RCP<::NOX::Abstract::Group> solution_group_ptr() override;

      STR::IMPLICIT::Generic& impl_int()
      {
        check_init_setup();
        return *implint_ptr_;
      };

      STR::Predict::Generic& predictor()
      {
        check_init_setup();
        return *predictor_ptr_;
      };

      Teuchos::RCP<STR::Predict::Generic> predictor_ptr()
      {
        check_init_setup();
        return predictor_ptr_;
      };

      [[nodiscard]] const STR::Nln::SOLVER::Generic& nln_solver() const
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

      STR::Nln::SOLVER::Generic& nln_solver()
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

      Teuchos::RCP<STR::Nln::SOLVER::Generic> nln_solver_ptr()
      {
        check_init_setup();
        return nlnsolver_ptr_;
      };

      //! @name Attribute access functions
      //@{

      //! Provide Name
      enum Inpar::STR::DynamicType method_name() const override;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      int method_steps() const override;

      //! Give local order of accuracy of displacement part
      int method_order_of_accuracy_dis() const override;

      //! Give local order of accuracy of velocity part
      int method_order_of_accuracy_vel() const override;

      //! Return linear error coefficient of displacements
      double method_lin_err_coeff_dis() const override;

      //! Return linear error coefficient of velocities
      double method_lin_err_coeff_vel() const override;

      ///@}

     private:
      //! ptr to the implicit time integrator object
      Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

      //! ptr to the predictor object
      Teuchos::RCP<STR::Predict::Generic> predictor_ptr_;

      //! ptr to the non-linear solver object
      Teuchos::RCP<STR::Nln::SOLVER::Generic> nlnsolver_ptr_;

      //! ptr to the nox group object
      Teuchos::RCP<::NOX::Abstract::Group> grp_ptr_;
    };
  }  // namespace TimeInt
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
