/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_impl_parameter.H

\brief general element parameter for fluid adjoint equations for topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_FLUIDADJOINT3_IMPL_PARAMETER_H_
#define TOPOPT_FLUIDADJOINT3_IMPL_PARAMETER_H_


#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_topopt.H"


namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters (constant over time)
class FluidAdjoint3ImplParameter
{
public:

  //! Singleton access method
  static Teuchos::RCP<FluidAdjoint3ImplParameter> Instance();

  //! general fluid parameter are set
  void SetElementGeneralAdjointParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! time parameter are set
  void SetElementAdjointTimeParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! print parameter to screen
  void PrintAdjointParameter() const;

  bool GeneralParamsSet() {return set_general_adjoint_parameter_;}
  bool ObjDissipationTerm() {return dissipation_;}
  bool ObjPresDropTerm() {return pressure_drop_;}
  double ObjDissipationFac() {return dissipation_fac_;}
  double ObjPresDropFac() {return pressure_drop_fac_;}

  INPAR::TOPOPT::AdjointType AdjointType() {return adjoint_type_;}
  bool IsStationary() {return is_stationary_;}
  bool IsInconsistent() {return is_inconsistent_;}
  INPAR::FLUID::PSPG PSPG() {return pspg_;}
  INPAR::FLUID::SUPG SUPG() {return supg_;}
  INPAR::FLUID::CStab CStab() {return cstab_;}
  INPAR::FLUID::TauType TauType() {return whichtau_;}
  bool EvalMatAtGP() {return mat_gp_;}
  bool EvalTauAtGP() {return tau_gp_;}

  INPAR::FLUID::TimeIntegrationScheme TimeIntScheme() {return timealgo_;}
  double Time() {return time_;}
  double Dt() {return dt_;}
  double Theta() {return theta_;}
  double OmTheta() {return omtheta_;}
  double ThetaPre() {return theta_pre_;}
  double OmThetaPre() {return omtheta_pre_;}
  double ThetaDiv() {return theta_div_;}
  double OmThetaDiv() {return omtheta_div_;}
  double Timefac() {return timefac_;}
  double TimefacRhs() {return timefacrhs_;}
  double TimefacPre() {return timefacpre_;}
  double TimefacPreRhs() {return timefacprerhs_;}
  double TimefacDiv() {return timefacdiv_;}
  double TimefacDivRhs() {return timefacdivrhs_;}

  double Density() {return dens_;}
  double Viscosity() {return visc_;}

  INPAR::TOPOPT::AdjointCase TestCase() {return testcase_;}

private:

  //! Flag SetGeneralParameter was called
  bool set_general_adjoint_parameter_;

  // set if objective contains dissipation
  bool dissipation_;
  // set if objective contains pressure drop
  bool pressure_drop_;
  // set objective's dissipation factor
  double dissipation_fac_;
  // set objective's pressure drop factor
  double pressure_drop_fac_;

  //! type of adjoint equations
  INPAR::TOPOPT::AdjointType adjoint_type_;
  //! flag to (de)activate stationary formulation
  bool is_stationary_;
  //! flag to (de)activate second derivatives
  bool is_inconsistent_;
  //! Flag to (de)activate PSPG stabilization
  INPAR::FLUID::PSPG pspg_;
  //! Flag to (de)activate SUPG stabilization
  INPAR::FLUID::SUPG supg_;
  //! Flag to (de)activate least-squares stabilization of continuity equation
  INPAR::FLUID::CStab cstab_;
  //! Flag to define tau
  INPAR::FLUID::TauType whichtau_;
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
  //! factor for left-hand side due to one-step-theta time-integration scheme
  double theta_;
  //! factor for right-hand side due to one-step-theta time-integration scheme
  double omtheta_;
  //! factor for pressure left-hand side due to one-step-theta time-integration scheme
  double theta_pre_;
  //! factor for pressure right-hand side due to one-step-theta time-integration scheme
  double omtheta_pre_;
  //! factor for continuity left-hand side due to one-step-theta time-integration scheme
  double theta_div_;
  //! factor for continuity right-hand side due to one-step-theta time-integration scheme
  double omtheta_div_;
  //! timefac = dt_ * theta_
  double timefac_;
  //! timefacrhs = dt_ * (1-theta_)
  double timefacrhs_;
  //! timefac = dt_ * theta_pre_
  double timefacpre_;
  //! timefacrhs = dt_ * (1-theta_pre_)
  double timefacprerhs_;
  //! timefac = dt_ * theta_div_
  double timefacdiv_;
  //! timefacrhs = dt_ * (1-theta_div_)
  double timefacdivrhs_;

  /// physical density
  double dens_;
  /// physical viscosity
  double visc_;

  //! special test case
  INPAR::TOPOPT::AdjointCase testcase_;

  // private constructor
  FluidAdjoint3ImplParameter();

  /// no copying
  FluidAdjoint3ImplParameter(const FluidAdjoint3ImplParameter&);

  /// no copying
  FluidAdjoint3ImplParameter& operator=(const FluidAdjoint3ImplParameter&);

  /// the single instance
  static Teuchos::RCP<FluidAdjoint3ImplParameter> instance_;
};

}
}

#endif /* TOPOPT_FLUIDADJOINT3_IMPL_PARAMETER_H_ */
