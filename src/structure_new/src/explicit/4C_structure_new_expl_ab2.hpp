/*-----------------------------------------------------------*/
/*! \file

\brief Adams-Bashforth-2 time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_EXPL_AB2_HPP
#define FOUR_C_STRUCTURE_NEW_EXPL_AB2_HPP

#include "4C_config.hpp"

#include "4C_structure_new_expl_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief Adams-Bashforth-2 time integration for solid dynamics
     *
     */
    class AdamsBashforth2 : public Generic
    {
     public:
      //! constructor
      AdamsBashforth2();


      //! Setup class variables (derived)
      void Setup() override;

      //! Post setup operation (compute initial equilibrium state) (derived)
      void PostSetup() override;

      //! Set state variables (derived)
      void SetState(const Epetra_Vector& x) override;

      //! return integration factor (derived)
      [[nodiscard]] double GetIntParam() const override { return -1.0; }

      /*! \brief Add the viscous and mass contributions to the right hand side
       */
      void add_visco_mass_contributions(Epetra_Vector& f) const override;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       */
      void add_visco_mass_contributions(CORE::LINALG::SparseOperator& jac) const override;

      //! Update configuration after time step (derived)
      void UpdateStepState() override;

      //! (derived)
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      /*! read restart information of the different time integration schemes
       *  and model evaluators (derived) */
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! @name Attribute access functions
      //@{

      //! Return time integrator name
      [[nodiscard]] enum INPAR::STR::DynamicType MethodName() const override
      {
        return INPAR::STR::dyna_ab2;
      }

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a m-multistep method returns m
      [[nodiscard]] int MethodSteps() const override { return 2; }

      //! Give local order of accuracy of displacement part
      [[nodiscard]] int method_order_of_accuracy_dis() const override { return 2; }

      //! Give local order of accuracy of velocity part
      [[nodiscard]] int method_order_of_accuracy_vel() const override { return 2; }

      //! Return linear error coefficient of displacements
      [[nodiscard]] double method_lin_err_coeff_dis() const override;

      //! Return linear error coefficient of velocities
      [[nodiscard]] double method_lin_err_coeff_vel() const override
      {
        return method_lin_err_coeff_dis();
      }

      //@}

     private:
      //! viscous force vector F_viscous F_{viscous;n+1}
      Teuchos::RCP<Epetra_Vector> fvisconp_ptr_;

      //! viscous force vector F_viscous F_{viscous;n}
      Teuchos::RCP<Epetra_Vector> fviscon_ptr_;

      //! pointer to inertial force vector F_{inertial,n+1} at new time
      Teuchos::RCP<Epetra_Vector> finertianp_ptr_;

      //! pointer to inertial force vector F_{inertial,n} at last time
      Teuchos::RCP<Epetra_Vector> finertian_ptr_;
    };
  }  // namespace EXPLICIT
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
