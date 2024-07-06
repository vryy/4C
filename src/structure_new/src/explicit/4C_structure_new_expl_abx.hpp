/*-----------------------------------------------------------*/
/*! \file

\brief High order Adams-Bashforth time integration for solid dynamics

The Adams-Bashforth scheme will update the displacement and velocities as follows:
+ For second order Adams-Bashforth:

  u_{n+1} = u_{n} + 3/2 dt v_{n} - 1/2 dt v_{n-1}
  v_{n+1} = v_{n} + 3/2 dt a_{n} - 1/2 dt a_{n-1}

+ For fourth order Adams-Bashforth:

  u_{n+1} = u_{n} + 55/24 dt v_{n} - 59/24 dt v_{n-1} + 37/24 dt v_{n-2} - 9/24 dt v_{n-3}
  v_{n+1} = v_{n} + 55/24 dt a_{n} - 59/24 dt a_{n-1} + 37/24 dt a_{n-2} - 9/24 dt a_{n-3}

Adams-Bashforth scheme of order p typically needs previous p-1 steps for explicit extrapolation.
This is generally not available in a new analysis, thus modified forward Euler scheme is employed in
the first p-1 steps. For a restart analysis, those values are restored from the history.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_EXPL_ABX_HPP
#define FOUR_C_STRUCTURE_NEW_EXPL_ABX_HPP

#include "4C_config.hpp"

#include "4C_structure_new_expl_generic.hpp"

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace EXPLICIT
  {
    /*! \brief Helper for high order Adams-Bashforth
     *
     */
    template <int t_order>
    struct AdamsBashforthHelper;

    /*! \brief High order Adams-Bashforth time integration for solid dynamics
     *
     */
    template <int t_order>
    class AdamsBashforthX : public Generic
    {
     public:
      //! constructor
      AdamsBashforthX();

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
      [[nodiscard]] enum Inpar::Solid::DynamicType method_name() const override
      {
        return Inpar::Solid::dyna_ab4;
      }

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a m-multistep method returns m
      [[nodiscard]] int method_steps() const override { return t_order; }

      //! Give local order of accuracy of displacement part
      [[nodiscard]] int method_order_of_accuracy_dis() const override { return t_order; }

      //! Give local order of accuracy of velocity part
      [[nodiscard]] int method_order_of_accuracy_vel() const override { return t_order; }

      //! Return linear error coefficient of displacements
      [[nodiscard]] double method_lin_err_coeff_dis() const override
      {
        FOUR_C_THROW("not yet implemented");
      }

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

      /*! \brief Flag to record the computing phase
       * compute_phase_ = 0..order: compute the initial values using modified forward Euler
       * compute_phase_ > order: compute the value using Adams-Bashforth multistep extrapolation
       * (normal mode) By default, the initial value is 0, and is incremented until the first
       * [order] steps are completed.
       */
      int compute_phase_;
    };

    //! @name Specialization for different type of Adams-Bashforth scheme
    //@{

    template <>
    struct AdamsBashforthHelper<2>
    {
      static constexpr std::array<double, 2> exc{{1.5, -0.5}};  // extrapolation coefficients
      static enum Inpar::Solid::DynamicType method_name() { return Inpar::Solid::dyna_ab2; }
    };

    template <>
    struct AdamsBashforthHelper<4>
    {
      static constexpr std::array<double, 4> exc{
          {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -9.0 / 24.0}};  // extrapolation coefficients
      static enum Inpar::Solid::DynamicType method_name() { return Inpar::Solid::dyna_ab4; }
    };

    //@}
  }  // namespace EXPLICIT
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
