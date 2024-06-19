/*-----------------------------------------------------------*/
/*! \file

\brief BDF-2 time-integration scheme


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_TIMINT_BDF2_HPP
#define FOUR_C_FLUID_TIMINT_BDF2_HPP


#include "4C_config.hpp"

#include "4C_fluid_implicit_integration.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntBDF2 : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntBDF2(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;

    /*!
    \brief Print information about current time step to screen

    */
    void print_time_step_info() override;

    /*!
    \brief Set theta_ to its value, dependent on integration method for GenAlpha and BDF2

    */
    void SetTheta() override;

    /*!
    \brief Set the part of the righthandside belonging to the last
           timestep for incompressible or low-Mach-number flow

       for low-Mach-number flow: distinguish momentum and continuity part
       (continuity part only meaningful for low-Mach-number flow)

       Stationary/af-generalized-alpha:

                     mom: hist_ = 0.0
                    (con: hist_ = 0.0)

       One-step-Theta:

                     mom: hist_ = veln_  + dt*(1-Theta)*accn_
                    (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

       BDF2: for constant time step:

                     mom: hist_ = 4/3 veln_  - 1/3 velnm_
                    (con: hist_ = 4/3 densn_ - 1/3 densnm_)


    */
    void set_old_part_of_righthandside() override;

    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    void SetStateTimInt() override;

    /*!
    \brief Calculate time derivatives for
           stationary/one-step-theta/BDF2/af-generalized-alpha time integration
           for incompressible and low-Mach-number flow
    */
    void calculate_acceleration(const Teuchos::RCP<const Epetra_Vector> velnp,  ///< velocity at n+1
        const Teuchos::RCP<const Epetra_Vector> veln,   ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector> velnm,  ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector> accn,   ///< acceleration at n-1
        const Teuchos::RCP<Epetra_Vector> accnp         ///< acceleration at n+1
        ) override;

    /*!
    \brief Set gamma to a value

    */
    void SetGamma(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Scale separation

    */
    void Sep_Multiply() override;

    /*!
    \brief Output of filtered velocity

    */
    void OutputofFilteredVel(
        Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec) override;

    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void set_element_time_parameter() override;

    /*!
    \brief return scheme-specific time integration parameter

    */
    double TimIntParam() const override { return 0.0; }

    /*!
    \brief return scaling factor for the residual

    */
    double residual_scaling() const override { return 1.0 / (theta_ * dta_); }

    //! @name Time Step Size Adaptivity
    //@{

    //! Give local order of accuracy of velocity part
    int method_order_of_accuracy_vel() const override { return 2; };

    //! Give local order of accuracy of pressure part
    int method_order_of_accuracy_pres() const override { return 2; };

    /*! \brief Return linear error coefficient of velocity
     *
     *  The linear discretization error reads
     *  \f[
     *  e \approx -\frac{\left(\Delta t_n + \Delta t_{n-1}\right)^2}
     *                  {6\Delta t_n\left(2\Delta t_n+\Delta t_{n-1}\right)}
     *             \Delta t_n^3 \dddot{u}(t_n)
     *    + HOT\left(\Delta t_n^4\right)
     *  \f]
     *
     *  \author mayr.mt \date 04/2015
     */
    double method_lin_err_coeff_vel() const override;

    //@}

   protected:
    /// copy constructor
    TimIntBDF2(const TimIntBDF2& old);

    /*!
    \brief velocity required for evaluation of related quantities required on element level

    */
    Teuchos::RCP<const Epetra_Vector> evaluation_vel() override { return velnp_; };


   private:
    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

  };  // class TimIntBDF2

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
