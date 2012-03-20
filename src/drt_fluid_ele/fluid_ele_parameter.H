/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter.H

\brief Evaluation of general fluid parameter

FluidEleParameter::SetParameter(Teuchos::ParameterList& params)
set all general fluid parameter once for all elements.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef FLUID_ELE_PARAMETER_H
#define FLUID_ELE_PARAMETER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"

namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters (constant over time)
class FluidEleParameter
{
public:

  //! Singleton access method
  static FluidEleParameter* Instance();

  //! general fluid parameter are set
  void SetElementGeneralFluidParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! time parameter are set
  void SetElementTimeParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! turbulence parameters are set
  void SetElementTurbulenceParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! loma parameters are set
  void SetElementLomaParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! print parameter to screen
  void PrintFluidParameter();

  //! Flag SetGeneralParameter was called
  bool set_general_fluid_parameter_;
  /// Flags to switch on/off the different fluid formulations
  //! flag to (de)activate generalized-alpha time-integration scheme
  bool is_genalpha_;
  /// Flags to switch on/off the different fluid formulations
  //! flag to (de)activate generalized-alpha-np1 time-integration scheme
  bool is_genalpha_np_;
  //! flag to (de)activate conservative formulation
  bool is_conservative_;
  //! flag to (de)activate stationary formulation
  bool is_stationary_;
  //! flag to (de)activate Newton linearization
  bool is_newton_;
  //! flag to (de)activate second derivatives
  bool is_inconsistent_;
  //! flag to (de)activate potential reactive terms
  bool reaction_;
  //! flag to (de)aktivate porous darcy flow
  bool darcy_;
  //! flag to (de)aktivate reaction due to topology optimization
  bool reaction_topopt_;
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
  //! Flag to (de)activate reactive term in residual-based stabilization
  INPAR::FLUID::RStab rstab_;
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

  /// parameter for stabilization
  //! (sign) factor for viscous and reactive stabilization terms
  double viscreastabfac_;

  /// constant parameters for the turbulence formulation
  /// subgrid viscosity models
  //! smagorinsky constant
  double Cs_;
  bool Cs_averaged_;
  //! channel length to normalize the normal wall distance
  double l_tau_;
  //! scale similarity constant
  double Cl_;
  /// multifractal subgrid-scales
  double Csgs_;
  double Csgs_phi_;
  double alpha_;
  bool CalcN_;
  double N_;
  enum INPAR::FLUID::RefVelocity refvel_;
  enum INPAR::FLUID::RefLength reflength_;
  double c_nu_;
  double c_diff_;
  bool near_wall_limit_;
  bool B_gp_;
  double beta_;

  //! some loma specific parameters
  //! flag for material update
  bool update_mat_;
  //! flag to (de)activate continuity SUPG term
  INPAR::FLUID::SUPG conti_supg_;
  //! flag to (de)activate continuity cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress conti_cross_;
  //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress conti_reynolds_;
  //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
  bool multifrac_loma_conti_;

private:

  // private constructor
  FluidEleParameter();

  /// no copying
  FluidEleParameter(const FluidEleParameter&);

  /// no copying
  FluidEleParameter& operator=(const FluidEleParameter&);

  /// the single instance
  static FluidEleParameter* instance_;
};

}
}

#endif
