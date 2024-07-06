/*-----------------------------------------------------------*/
/*! \file

\brief Setting of time-integration parameters in fluid element evaluation

This class provides the fluid element parameter for the time integration,
which are unique and equal for every fluid in the problem. Time integration
with different parameters in more than one fluid field is not yet supported.


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_PARAMETER_TIMINT_HPP
#define FOUR_C_FLUID_ELE_PARAMETER_TIMINT_HPP

#include "4C_config.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class FluidEleParameterTimInt
    {
     public:
      //! Singleton access method
      static FluidEleParameterTimInt* instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      //! time parameter are set
      void set_element_time_parameter(Teuchos::ParameterList& params  //> parameter list
      );

      //! print parameter to screen
      void print_fluid_time_parameter();

      //! flag to (de)activate generalized-alpha time-integration scheme
      bool is_genalpha() { return is_genalpha_; };
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate generalized-alpha-np1 time-integration scheme
      bool is_genalpha_np() { return is_genalpha_np_; };
      //! flag to (de)activate stationary formulation
      bool is_stationary() { return is_stationary_; };
      //! flag to (de)activate stationary formulation
      bool is_one_step_theta() { return is_one_step_theta_; };
      //! flag to (de)activate full implicit pressure in momentum equation (and adjust continuity eq
      //! accordingly)
      bool is_full_impl_pressure_and_cont() const { return is_cont_impl_press_impl_; };
      //! flag to (de)activate full implicit pressure in momentum equation (and adjust continuity eq
      //! accordingly)
      bool is_impl_pressure() const { return is_cont_impl_press_normal_; };

      //! flag for new or old implementation. WILL BE REMOVED!
      bool is_new_ost_implementation() const { return ostnew_; };

      /// parameters for the time integration
      //! time algorithm
      Inpar::FLUID::TimeIntegrationScheme time_algo() { return timealgo_; };
      //! actual time to evaluate the body BC
      double time() { return time_; };
      //! time-step length
      double dt() { return dt_; };
      //! timefac = dt_ * ("pseudo"-)theta_
      double time_fac() { return timefac_; };
      //! factor for left-hand side due to one-step-theta time-integration scheme
      double theta() { return theta_; };
      //! factor for right-hand side due to one-step-theta time-integration scheme
      double om_theta() { return omtheta_; };
      //! generalised-alpha parameter (connecting velocity and acceleration)
      double gamma() { return gamma_; };
      //! generalised-alpha parameter (velocity)
      double alpha_f() { return alpha_f_; };
      //! generalised-alpha parameter (acceleration)
      double alpha_m() { return alpha_m_; };
      //! generalised-alpha parameter, alphaF_*gamma_*dt_
      double afgdt() { return afgdt_; };
      //! generalised-alpha parameter, gamma_/alphaM_*dt_
      //! time integration factor for the right hand side (boundary elements)
      double time_fac_rhs() { return timefacrhs_; };
      //! time integration factor for the left hand side (pressure)
      double time_fac_pre() { return timefacpre_; };

     private:
      //! Flag SetGeneralParameter was called
      bool set_general_fluid_timeparameter_;
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate generalized-alpha time-integration scheme
      bool is_genalpha_;
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate generalized-alpha-np1 time-integration scheme
      bool is_genalpha_np_;
      //! flag to (de)activate stationary formulation
      bool is_stationary_;
      //! flag to (de)activate one step theta
      bool is_one_step_theta_;
      //! flag to (de)activate full implicit pressure and continuity eq.
      bool is_cont_impl_press_impl_;
      //! flag to (de)activate implicit handling of continuity eq, but pressure at theta.
      bool is_cont_impl_press_normal_;

      /// parameters for the time integration
      //! time algorithm
      Inpar::FLUID::TimeIntegrationScheme timealgo_;
      //! One Step Theta discretization for the continuity and pressure.
      Inpar::FLUID::OstContAndPress ostalgo_;
      //! New One Step Theta implementation flag
      bool ostnew_;
      //! actual time to evaluate the body BC
      double time_;
      //! time-step length
      double dt_;
      //! timefac = dt_ * ("pseudo"-)theta_
      double timefac_;
      //! factor for left-hand side due to one-step-theta time-integration scheme
      double theta_;
      //! factor for right-hand side due to one-step-theta time-integration scheme
      double omtheta_;
      //! generalised-alpha parameter (connecting velocity and acceleration)
      double gamma_;
      //! generalised-alpha parameter (velocity)
      double alpha_f_;
      //! generalised-alpha parameter (acceleration)
      double alpha_m_;
      //! generalised-alpha parameter, alphaF_*gamma_*dt_
      double afgdt_;
      //! generalised-alpha parameter, gamma_/alphaM_*dt_
      //! time integration factor for the right hand side (boundary elements)
      double timefacrhs_;
      //! time integration factor for the left hand side (pressure)
      double timefacpre_;

      // private constructor
      FluidEleParameterTimInt();


    };  // class FluidEleParameterTimInt

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
