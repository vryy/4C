/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_timint.cpp

\brief Evaluation of general fluid parameter for time integration

This class provides the fluid element parameter for the time integration,
which are unique and equal for every fluid in the problem. Time integration
with different parameters in more than one fluid field is not supported so
far.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_TIMINT_H
#define FLUID_ELE_PARAMETER_TIMINT_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "fluid_ele_parameter_interface.H"
#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters (constant over time)
class FluidEleParameterTimInt : public FluidEleParameterInterface
{
public:

  //! Singleton access method
  static FluidEleParameterTimInt* Instance(bool create=true);

  //! Singleton destruction method called from global problem
  virtual void Done();

  //! dtor
  ~FluidEleParameterTimInt() {}

  //! time parameter are set
  void SetElementTimeParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! print parameter to screen
  void PrintFluidTimeParameter();

  //! flag to (de)activate generalized-alpha time-integration scheme
  bool IsGenalpha(){return is_genalpha_;};
  /// Flags to switch on/off the different fluid formulations
  //! flag to (de)activate generalized-alpha-np1 time-integration scheme
  bool IsGenalphaNP(){return is_genalpha_np_;};
  //! flag to (de)activate stationary formulation
  bool IsStationary(){return is_stationary_;};

  /// parameters for the time integration
  //! time algorithm
  INPAR::FLUID::TimeIntegrationScheme TimeAlgo(){return timealgo_;};
  //! actual time to evaluate the body BC
  double Time(){return time_;};
  //! time-step length
  double Dt(){return dt_;};
  //! timefac = dt_ * ("pseudo"-)theta_
  double TimeFac(){return timefac_;};
  //! factor for left-hand side due to one-step-theta time-integration scheme
  double Theta(){return theta_;};
  //! factor for right-hand side due to one-step-theta time-integration scheme
  double OmTheta(){return omtheta_;};
  //! generalised-alpha parameter (connecting velocity and acceleration)
  double Gamma(){return gamma_;};
  //! generalised-alpha parameter (velocity)
  double AlphaF(){return alphaF_;};
  //! generalised-alpha parameter (acceleration)
  double AlphaM(){return alphaM_;};
  //! generalised-alpha parameter, alphaF_*gamma_*dt_
  double Afgdt(){return afgdt_;};
  //! generalised-alpha parameter, gamma_/alphaM_*dt_
  //! time integration factor for the right hand side (boundary elements)
  double TimeFacRhs(){return timefacrhs_;};
  //! time integration factor for the left hand side (pressure)
  double TimeFacPre(){return timefacpre_;};

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

  /// parameters for the time integration
  //! time algorithm
  INPAR::FLUID::TimeIntegrationScheme timealgo_;
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

private:

  // private constructor
  FluidEleParameterTimInt();
  /// no copying
  FluidEleParameterTimInt(const FluidEleParameterTimInt&);
  /// no copying
  FluidEleParameterTimInt& operator=(const FluidEleParameterTimInt&);

}; // class FluidEleParameterTimInt

} // namespace ELEMENTS
} // namespace DRT

#endif
