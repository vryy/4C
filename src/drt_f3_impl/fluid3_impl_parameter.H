#ifndef FLUID3_IMPL_PARAMETER_H
#define FLUID3_IMPL_PARAMETER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
namespace ELEMENTS
{

/// General fluid evaluation parameters.
class Fluid3ImplParameter
{
public:

  static Fluid3ImplParameter* Instance();

  //! general parameter are set
  void SetParameter(Teuchos::ParameterList& params);

  /// Flags to switch on/off the different fluid formulations
  //! flag to (de)activate generalized-alpha time-integration scheme
  bool is_genalpha_;
  //! flag to (de)activate conservative formulation
  bool is_conservative_;
  //! flag to (de)activate stationary formulation
  bool is_stationary_;
  //! flag to (de)activate Newton linearization
  bool is_newton_;
  //! flag to (de)activate second derivatives
  bool is_inconsistent_;
  //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq)
  INPAR::FLUID::PhysicalType physicaltype_;
  //! Flag to (de)activate time-dependent subgrid stabilization
  INPAR::FLUID::SubscalesTD tds_;
  //! Flag to (de)activate time-dependent term in large-scale momentum equation
  INPAR::FLUID::Transient transient_;
  //! Flag to (de)activate PSPG stabilization
  INPAR::FLUID::PSPG pspg_;
  //! Flag to (de)activate SUPG stabilization
  INPAR::FLUID::SUPG supg_;
  //! Flag to (de)activate viscous term in residual-based stabilization
  INPAR::FLUID::VStab vstab_;
  //! Flag to (de)activate least-squares stabilization of continuity equation
  INPAR::FLUID::CStab cstab_;
  //! Flag to (de)activate cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress cross_;
  //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress reynolds_;
  //! Flag to define tau
  INPAR::FLUID::TauType whichtau_;
  //! flag to (de)activate fine-scale subgrid viscosity
  INPAR::FLUID::FineSubgridVisc fssgv_;
  //! flag to define turbulence model
  INPAR::FLUID::TurbModelAction turb_mod_action_;
  //! flag for material evaluation at Gaussian integration points
  bool mat_gp_;
  //! flag for stabilization parameter evaluation at Gaussian integration points
  bool tau_gp_;

  /// parameters for the time integration
  //! actual time to evaluate the body BC
  double time_;
  //! time-step length
  double dt_;
  //! timefac = dt_ * theta_
  double timefac_;
  //! factor for right-hand side due to time-integration scheme
  double omtheta_;
  //! generalised-alpha parameter (connecting velocity and acceleration)
  double gamma_;
  //! generalised-alpha parameter (velocity)
  double alphaF_;
  //! generalised-alpha parameter (acceleration)
  double alphaM_;
  //! generalised-alpha parameter, alphaF_*gamma_*dt_
  double afgdt_;

  /// constant parameters for the turbulence formulation
  //! smagorinsky constant
  double Cs_;
  //! channel length to normalize the normal wall distance
  double l_tau_;

private:

  // private constructor
  Fluid3ImplParameter();

  /// no copying
  Fluid3ImplParameter(const Fluid3ImplParameter&);

  /// no copying
  Fluid3ImplParameter& operator=(const Fluid3ImplParameter&);

  /// the single instance
  static Fluid3ImplParameter* instance_;
};

}
}

#endif
