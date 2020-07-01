/*-----------------------------------------------------------*/
/*! \file

\brief Setting of time-integration parameters in fluid element evaluation

This class provides the fluid element parameter for the time integration,
which are unique and equal for every fluid in the problem. Time integration
with different parameters in more than one fluid field is not yet supported.


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_TIMINT_H
#define FLUID_ELE_PARAMETER_TIMINT_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class FluidEleParameterTimInt : public DRT::SingletonDestruction
    {
     public:
      //! Singleton access method
      static FluidEleParameterTimInt* Instance(bool create = true);

      //! Singleton destruction method called from global problem
      void Done();

      //! dtor
      ~FluidEleParameterTimInt() {}

      //! time parameter are set
      void SetElementTimeParameter(Teuchos::ParameterList& params  //> parameter list
      );

      //! print parameter to screen
      void PrintFluidTimeParameter();

      //! flag to (de)activate generalized-alpha time-integration scheme
      bool IsGenalpha() { return is_genalpha_; };
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate generalized-alpha-np1 time-integration scheme
      bool IsGenalphaNP() { return is_genalpha_np_; };
      //! flag to (de)activate stationary formulation
      bool IsStationary() { return is_stationary_; };
      //! flag to (de)activate stationary formulation
      bool IsOneStepTheta() { return is_one_step_theta_; };
      //! flag to (de)activate full implicit pressure in momentum equation (and adjust continuity eq
      //! accordingly)
      bool IsFullImplPressureAndCont() const { return is_cont_impl_press_impl_; };
      //! flag to (de)activate full implicit pressure in momentum equation (and adjust continuity eq
      //! accordingly)
      bool IsImplPressure() const { return is_cont_impl_press_normal_; };

      //! flag for new or old implementation. WILL BE REMOVED!
      bool IsNewOSTImplementation() const { return ostnew_; };

      /// parameters for the time integration
      //! time algorithm
      INPAR::FLUID::TimeIntegrationScheme TimeAlgo() { return timealgo_; };
      //! actual time to evaluate the body BC
      double Time() { return time_; };
      //! time-step length
      double Dt() { return dt_; };
      //! timefac = dt_ * ("pseudo"-)theta_
      double TimeFac() { return timefac_; };
      //! factor for left-hand side due to one-step-theta time-integration scheme
      double Theta() { return theta_; };
      //! factor for right-hand side due to one-step-theta time-integration scheme
      double OmTheta() { return omtheta_; };
      //! generalised-alpha parameter (connecting velocity and acceleration)
      double Gamma() { return gamma_; };
      //! generalised-alpha parameter (velocity)
      double AlphaF() { return alphaF_; };
      //! generalised-alpha parameter (acceleration)
      double AlphaM() { return alphaM_; };
      //! generalised-alpha parameter, alphaF_*gamma_*dt_
      double Afgdt() { return afgdt_; };
      //! generalised-alpha parameter, gamma_/alphaM_*dt_
      //! time integration factor for the right hand side (boundary elements)
      double TimeFacRhs() { return timefacrhs_; };
      //! time integration factor for the left hand side (pressure)
      double TimeFacPre() { return timefacpre_; };

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
      INPAR::FLUID::TimeIntegrationScheme timealgo_;
      //! One Step Theta discretization for the continuity and pressure.
      INPAR::FLUID::OST_Cont_and_Press ostalgo_;
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
      double alphaF_;
      //! generalised-alpha parameter (acceleration)
      double alphaM_;
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
}  // namespace DRT

#endif
