/*-----------------------------------------------------------*/
/*! \file

\brief Central differences time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_EXPL_CENTRDIFF_HPP
#define FOUR_C_STRUCTURE_NEW_EXPL_CENTRDIFF_HPP

#include "4C_config.hpp"

#include "4C_structure_new_expl_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief Central differences time integration for solid dynamics
     *
     */
    class CentrDiff : public Generic
    {
     public:
      //! constructor
      CentrDiff();

      //! Setup class variables (derived)
      void setup() override;

      //! Post setup operation (compute initial equilibrium state) (derived)
      void post_setup() override;

      //! Set state variables (derived)
      void set_state(const Epetra_Vector& x) override;

      //! return integration factor (derived)
      [[nodiscard]] double get_int_param() const override { return -1.0; }

      /*! \brief Add the viscous and mass contributions to the right hand side
       */
      void add_visco_mass_contributions(Epetra_Vector& f) const override;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       */
      void add_visco_mass_contributions(Core::LinAlg::SparseOperator& jac) const override;

      //! Update configuration after time step (derived)
      void update_step_state() override;

      //! (derived)
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      /*! read restart information of the different time integration schemes
       *  and model evaluators (derived) */
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! @name Attribute access functions
      //@{

      //! Return time integrator name
      [[nodiscard]] enum Inpar::STR::DynamicType method_name() const override
      {
        return Inpar::STR::dyna_centrdiff;
      }

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a m-multistep method returns m
      [[nodiscard]] int method_steps() const override { return 1; }

      //! Give local order of accuracy of displacement part
      [[nodiscard]] int method_order_of_accuracy_dis() const override { return 2; }

      //! Give local order of accuracy of velocity part
      [[nodiscard]] int method_order_of_accuracy_vel() const override { return 2; }

      /*! Return linear error coefficient of displacements
       *
       *  The local discretization error reads
       *  \f[
       *  e \approx \frac{1}{6}\Delta t_n^3 \dddot{d_n} + HOT(\Delta t_n^4)
       *  \f]
       */
      [[nodiscard]] double method_lin_err_coeff_dis() const override { return 1. / 6.; }

      /*! Return linear error coefficient of velocities
       *
       *  The local discretization error reads
       *  \f[
       *  e \approx -\frac{1}{12}\Delta t_n^3 d_n^{(4)} + HOT(\Delta t_n^4)
       *  \f]
       */
      [[nodiscard]] double method_lin_err_coeff_vel() const override { return -1. / 12.; }

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
