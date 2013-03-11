/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter.H

\brief Evaluation of general fluid parameter

FluidEleParameter::SetParameter(Teuchos::ParameterList& params)
set all general fluid parameter once for all elements.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef FLUID_ELE_PARAMETER_H
#define FLUID_ELE_PARAMETER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>

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
  static Teuchos::RCP<FluidEleParameter> Instance();

  //! general fluid parameter are set
  void SetElementGeneralFluidParameter(Teuchos::ParameterList& params,   //> parameter list
                                       int myrank                        //> proc id
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

  //! poro parameters are set
  void SetElementPoroParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! print parameter to screen
  void PrintFluidParameter();

  //! flag to (de)activate generalized-alpha time-integration scheme
  bool IsGenalpha(){return is_genalpha_;};
  /// Flags to switch on/off the different fluid formulations
  //! flag to (de)activate generalized-alpha-np1 time-integration scheme
  bool IsGenalphaNP(){return is_genalpha_np_;};
  //! flag to (de)activate conservative formulation
  bool IsConservative(){return is_conservative_;};
  //! flag to (de)activate stationary formulation
  bool IsStationary(){return is_stationary_;};
  //! flag to (de)activate Newton linearization
  bool IsNewton(){return is_newton_;};
  //! flag to (de)activate second derivatives
  bool IsInconsistent(){return is_inconsistent_;};
  //! flag to (de)activate potential reactive terms
  bool Reaction(){return reaction_;};
  //! flag to (de)aktivate porous darcy flow
  bool Darcy(){return darcy_;};
  //! flag to (de)aktivate reaction due to topology optimization
  bool ReactionTopopt(){return reaction_topopt_;};
  //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq)
  INPAR::FLUID::PhysicalType PhysicalType(){return physicaltype_;};
  //! get the stabtype
  INPAR::FLUID::StabType StabType(){return stabtype_;};

  /// parameter for residual stabilization
  //! Flag to (de)activate time-dependent subgrid stabilization
  INPAR::FLUID::SubscalesTD Tds(){return tds_;};
  //! Flag to (de)activate time-dependent term in large-scale momentum equation
  INPAR::FLUID::Transient Transient(){return transient_;};
  //! Flag to (de)activate PSPG stabilization
  INPAR::FLUID::PSPG PSPG(){return pspg_;};
  //! Flag to (de)activate SUPG stabilization
  INPAR::FLUID::SUPG SUPG(){return supg_;};
  //! Flag to (de)activate viscous term in residual-based stabilization
  INPAR::FLUID::VStab VStab(){return vstab_;};
  //! Flag to (de)activate reactive term in residual-based stabilization
  INPAR::FLUID::RStab RStab(){return rstab_;};
  //! Flag to (de)activate least-squares stabilization of continuity equation
  INPAR::FLUID::CStab CStab(){return cstab_;};
  //! Flag to (de)activate cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress Cross(){return cross_;};
  //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress Reynolds(){return reynolds_;};
  //! Flag to define tau
  INPAR::FLUID::TauType WhichTau(){return whichtau_;};
  //! flag to (de)activate fine-scale subgrid viscosity
  INPAR::FLUID::FineSubgridVisc Fssgv(){return fssgv_;};
  //! (sign) factor for viscous and reactive stabilization terms
  double ViscReaStabFac(){return viscreastabfac_;};


  /// parameter for edge-based stabilizations
  //! Flag to (de)activate pressure stabilization
  INPAR::FLUID::EOS_Pres EOS_Pres(){return EOS_pres_;};
  //! Flag to (de)activate convective streamline stabilization
  INPAR::FLUID::EOS_Conv_Stream EOS_Conv_Stream(){return EOS_conv_stream_;};
  //! Flag to (de)activate convective crosswind stabilization
  INPAR::FLUID::EOS_Conv_Cross EOS_Conv_Cross(){return EOS_conv_cross_;};
  //! Flag to (de)activate divergence stabilization
  INPAR::FLUID::EOS_Div EOS_Div(){return EOS_div_;};
  //! Flag to define element length
  INPAR::FLUID::EOS_ElementLength EOS_element_length(){return EOS_element_lenght_;};
  //! Flag to define tau for edge-based stabilization
  INPAR::FLUID::EOS_TauType EOS_WhichTau(){return EOS_whichtau_;};


  //! flag for material evaluation at Gaussian integration points
  bool MatGp(){return mat_gp_;};
  //! flag for stabilization parameter evaluation at Gaussian integration points
  bool TauGp(){return tau_gp_;};


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

  /// constant parameters for the turbulence formulation
  /// subgrid viscosity models
  //! flag to define turbulence model
  INPAR::FLUID::TurbModelAction TurbModAction(){return turb_mod_action_;};
  double Cs(){return Cs_;};
  bool CsAveraged(){return Cs_averaged_;};
  double Ci(){return Ci_;};
  void SetvanDriestdamping(double damping){van_Driest_damping_=damping; return;};
  double VanDriestdamping(){return van_Driest_damping_;};
  bool IncludeCi(){return include_Ci_;};
  double ltau(){return l_tau_;};
  //! scale similarity constant
  double Cl(){return Cl_;};
  /// multifractal subgrid-scales
  double Csgs(){return Csgs_;};
  double CsgsPhi(){double tmp = 0.0; if(not adapt_Csgs_phi_) tmp = Csgs_phi_; else tmp = Csgs_phi_*meanCai_; return tmp;};
  double Alpha(){return alpha_;};
  bool CalcN(){return CalcN_;};
  double N(){return N_;};
  enum INPAR::FLUID::RefVelocity RefVel(){return refvel_;};
  enum INPAR::FLUID::RefLength RefLength(){return reflength_;};
  double CNu(){return c_nu_;};
  double CDiff(){return c_diff_;};
  bool NearWallLimit(){return near_wall_limit_;};
  bool NearWallLimitScatra(){return near_wall_limit_scatra_;};
  bool BGp(){return B_gp_;};
  double Beta(){return beta_;};
  double MfsIsConservative(){return mfs_is_conservative_;};
  double AdaptCsgsPhi(){return adapt_Csgs_phi_;};
  void SetCsgsPhi(double meanCai){meanCai_ = meanCai; return;};

  //! some loma specific parameters
  //! flag for material update
  bool UpdateMat(){return update_mat_;};
  //! flag to (de)activate continuity SUPG term
  INPAR::FLUID::SUPG ContiSUPG(){return conti_supg_;};
  //! flag to (de)activate continuity cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress ContiCross(){return conti_cross_;};
  //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress ContiReynolds(){return conti_reynolds_;};
  //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
  bool MultiFracLomaConti(){return multifrac_loma_conti_;};

  //! some poro specific parameters
  //! partial integration of porosity gradient in continuity equation
  bool PoroContiPartInt(){return poro_conti_partint_;};

private:

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

  /// parameter to switch the stabilization
  INPAR::FLUID::StabType stabtype_;

  /// parameter for residual based stabilizations
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
  //! (sign) factor for viscous and reactive stabilization terms
  double viscreastabfac_;

  /// parameter for edge-based (EOS,CIP) stabilizations
  //! Flag to (de)activate pressure stabilization
  INPAR::FLUID::EOS_Pres EOS_pres_;
  //! Flag to (de)activate convective streamline stabilization
  INPAR::FLUID::EOS_Conv_Stream EOS_conv_stream_;
  //! Flag to (de)activate convective crosswind stabilization
  INPAR::FLUID::EOS_Conv_Cross EOS_conv_cross_;
  //! Flag to (de)activate divergence stabilization
  INPAR::FLUID::EOS_Div EOS_div_;
  //! Flag to define tau
  INPAR::FLUID::EOS_TauType EOS_whichtau_;
  //! element length definition for edge-base stabilization
  INPAR::FLUID::EOS_ElementLength EOS_element_lenght_;


  /// parameter for evaluation of material and stabilization parameter
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

  /// constant parameters for the turbulence formulation
  /// subgrid viscosity models
  //! flag to define turbulence model
  INPAR::FLUID::TurbModelAction turb_mod_action_;
  //! smagorinsky constant
  double Cs_;
  bool Cs_averaged_;
  //! loma: constant for isotropic part of subgrid-stress tensor (according to Yoshizawa 1986)
  double Ci_;
  bool include_Ci_;
  //! van Driest damping for channel flow
  double van_Driest_damping_;
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
  bool near_wall_limit_scatra_;
  bool B_gp_;
  double beta_;
  bool mfs_is_conservative_;
  bool adapt_Csgs_phi_;
  double meanCai_;

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

  //! some poro specific parameters
  //! partial integration of porosity gradient in continuity equation
  bool poro_conti_partint_;

  // private constructor
  FluidEleParameter();

  /// no copying
  FluidEleParameter(const FluidEleParameter&);

  /// no copying
  FluidEleParameter& operator=(const FluidEleParameter&);

  /// the single instance
  static Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> instance_;
};

}
}

#endif
