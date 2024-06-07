/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all explicit time integrators


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_EXPL_GENERIC_HPP
#define FOUR_C_STRUCTURE_NEW_EXPL_GENERIC_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_abstract_prepostoperator.hpp"
#include "4C_structure_new_integrator.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief A generic fully explicit time integrator for solid dynamics
     *
     * This serves as a base class for all fully explicit time integration schemes.
     */
    class Generic : public STR::Integrator
    {
     public:
      //! constructor
      Generic() = default;

      //! Setup (has to be implemented by the derived classes)
      void Setup() override;

      //! Apply the right hand side only (derived)
      bool apply_force(const Epetra_Vector& x, Epetra_Vector& f) override;

      //! \brief Apply the stiffness only (derived)
      bool apply_stiff(const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac) override;

      //! \brief Apply force and stiff at once (derived)
      bool apply_force_stiff(
          const Epetra_Vector& x, Epetra_Vector& f, Core::LinAlg::SparseOperator& jac) override;

      /*! \brief (derived)
       *
       */
      bool apply_correction_system(const enum NOX::Nln::CorrectionType type,
          const std::vector<Inpar::STR::ModelType>& constraint_models, const Epetra_Vector& x,
          Epetra_Vector& f, Core::LinAlg::SparseOperator& jac) override
      {
        return false;
      }

      /*! \brief Calculate characteristic/reference norms for forces (derived)
       */
      [[nodiscard]] double calc_ref_norm_force(
          const enum ::NOX::Abstract::Vector::NormType& type) const override;

      //! return the default step length
      [[nodiscard]] double get_default_step_length() const;

      //! compute the scaling operator for element based scaling using PTC (derived)
      void compute_jacobian_contributions_from_element_level_for_ptc(
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& scalingMatrixOpPtr) override;

      //! Assemble the right hand side
      bool assemble_force(Epetra_Vector& f,
          const std::vector<Inpar::STR::ModelType>* without_these_models = nullptr) const override;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      void pre_update() override{/* do nothing for now */};

      //! (derived)
      void update_constant_state_contributions() override;

      virtual void OutputStepstate() { ; };

      //! Update everything on element level after time step and after output (derived)
      void update_step_element() override;

      //! things that should be done after updating (derived)
      void post_update() override;

      /*! \brief Remove contributions from the structural right-hand side stemming
       *  from any condensation operations (typical example is contact) (derived)
       */
      void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const override;
      //! @}

      //! @name Attribute access functions
      //@{

      //! Provide Name
      [[nodiscard]] virtual enum Inpar::STR::DynamicType method_name() const = 0;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      [[nodiscard]] virtual int method_steps() const = 0;

      //! Give local order of accuracy of displacement part
      [[nodiscard]] virtual int method_order_of_accuracy_dis() const = 0;

      //! Give local order of accuracy of velocity part
      [[nodiscard]] virtual int method_order_of_accuracy_vel() const = 0;

      //! Return linear error coefficient of displacements
      [[nodiscard]] virtual double method_lin_err_coeff_dis() const = 0;

      //! Return linear error coefficient of velocities
      [[nodiscard]] virtual double method_lin_err_coeff_vel() const = 0;

      //! @}

     protected:
      //! reset the time step dependent parameters for the element evaluation [derived]
      void reset_eval_params() override;
    };
  }  // namespace EXPLICIT
}  // namespace STR

namespace NOX
{
  namespace Nln
  {
    namespace PrePostOp
    {
      namespace EXPLICIT
      {
        /*! \brief Explicit time integration helper class
         */
        class Generic : public NOX::Nln::Abstract::PrePostOperator
        {
         public:
          //! constructor
          Generic(const STR::EXPLICIT::Generic& expl)
              : default_step_(expl.get_default_step_length())
          {
          }

          //! derived
          void runPreIterate(const ::NOX::Solver::Generic& solver) override;

          //! derived
          void runPreSolve(const ::NOX::Solver::Generic& nlnSolver) override;

          //! derived
          void runPreComputeX(const NOX::Nln::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::Nln::Group& curr_grp) override;

          //! derived
          void runPostComputeX(const NOX::Nln::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::Nln::Group& curr_grp) override;

         private:
          //! default step length
          const double default_step_;
        };  // class Generic
      }     // namespace EXPLICIT
    }       // namespace PrePostOp
  }         // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
