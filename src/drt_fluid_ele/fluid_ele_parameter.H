/*----------------------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for element evaluation

This file has to contain all parameters called in fluid_ele_calc.cpp.
Additional parameters required in derived classes of FluidEleCalc have to
be set in problem specific parameter lists derived from this class.

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_H
#define FLUID_ELE_PARAMETER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_twophase.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "fluid_ele_parameter_timint.H"

#include "../drt_lib/drt_dserror.H"



namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameter : public DRT::SingletonDestruction
    {
     public:
      virtual ~FluidEleParameter() {}

      /// called upon destruction
      virtual void Done() = 0;

      /*========================================================================*/
      //! @name set-routines
      /*========================================================================*/

      //! general fluid parameter are set
      void SetElementGeneralFluidParameter(Teuchos::ParameterList& params,  //> parameter list
          int myrank);                                                      //> proc id

      //! turbulence parameters are set
      void SetElementTurbulenceParameters(Teuchos::ParameterList& params);  //> parameter list

      /// set loma parameters
      void SetElementLomaParameter(Teuchos::ParameterList& params);  //> parameter list

      //! set two-phase parameters
      void SetElementTwoPhaseParameter(Teuchos::ParameterList& params);  //> parameter list

      //! topopt parameters are set
      void SetElementTopoptParameter(Teuchos::ParameterList& params);  //> parameter list

      /*========================================================================*/
      //! @name access-routines
      /*========================================================================*/

      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, poro)
      INPAR::FLUID::PhysicalType PhysicalType() const { return physicaltype_; };
      //! flag to (de)activate conservative formulation
      bool IsConservative() const { return is_conservative_; };
      //! flag to (de)activate Newton linearization
      bool IsNewton() const { return is_newton_; };
      //! flag to (de)activate second derivatives
      bool IsInconsistent() const { return is_inconsistent_; };
      //! flag to (de)activate potential reactive terms
      bool Reaction() const { return reaction_; };
      //! Return function number of Oseen advective field
      int OseenFieldFuncNo() const { return oseenfieldfuncno_; };
      //! flag to activate consistent reconstruction of second derivatives
      bool IsReconstructDer() const { return is_reconstructder_; };

      /*----------------------------------------------------*/
      //! @name stabilization parameters
      /*----------------------------------------------------*/

      //! get the stabtype
      INPAR::FLUID::StabType StabType() const { return stabtype_; };
      /// parameter for residual stabilization
      //! Flag to (de)activate time-dependent subgrid stabilization
      INPAR::FLUID::SubscalesTD Tds() const { return tds_; };
      //! Flag to (de)activate time-dependent term in large-scale momentum equation
      INPAR::FLUID::Transient Transient() const { return transient_; };
      //! Flag to (de)activate PSPG stabilization
      bool PSPG() const { return pspg_; };
      //! Flag to (de)activate SUPG stabilization
      bool SUPG() const { return supg_; };
      //! Flag to (de)activate viscous term in residual-based stabilization
      INPAR::FLUID::VStab VStab() const { return vstab_; };
      //! Flag to (de)activate reactive term in residual-based stabilization
      INPAR::FLUID::RStab RStab() const { return rstab_; };
      //! Flag to (de)activate least-squares stabilization of continuity equation
      bool CStab() const { return graddiv_; };
      //! Flag to (de)activate cross-stress term -> residual-based VMM
      INPAR::FLUID::CrossStress Cross() const { return cross_; };
      //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
      INPAR::FLUID::ReynoldsStress Reynolds() const { return reynolds_; };
      //! Flag to define tau
      INPAR::FLUID::TauType WhichTau() const { return whichtau_; };
      //! Flag to define characteristic element length for tau_Mu
      INPAR::FLUID::CharEleLengthU CharEleLengthU() const { return charelelengthu_; };
      //! Flag to define characteristic element length for tau_Mp and tau_C
      INPAR::FLUID::CharEleLengthPC CharEleLengthPC() const { return charelelengthpc_; };
      //! (sign) factor for viscous and reactive stabilization terms
      double ViscReaStabFac() const { return viscreastabfac_; };

      //! Flag to (de)activate polynomial pressure projection stabilization
      bool PPP() const { return ppp_; };

      //! flag for material evaluation at Gaussian integration points
      bool MatGp() const { return mat_gp_; };
      //! flag for stabilization parameter evaluation at Gaussian integration points
      bool TauGp() const { return tau_gp_; };

      /*----------------------------------------------------*/
      //! @name two phase parameters
      /*----------------------------------------------------*/
      double GetInterfaceThickness() const { return interface_thickness_; };
      bool GetEnhancedGaussRuleInInterface() const { return enhanced_gaussrule_; }
      bool GetIncludeSurfaceTension() const { return include_surface_tension_; };
      INPAR::TWOPHASE::SurfaceTensionApprox GetSurfaceTensionApprox() const
      {
        return eval_surfacetension_;
      };


      /*----------------------------------------------------*/
      //! @name turbulence model
      /*----------------------------------------------------*/

      /// constant parameters for the turbulence formulation
      /// subgrid viscosity models
      //! flag to define turbulence model
      INPAR::FLUID::TurbModelAction TurbModAction() const { return turb_mod_action_; };
      double Cs() const { return Cs_; };
      bool CsAveraged() const { return Cs_averaged_; };
      double Ci() const { return Ci_; };
      void SetvanDriestdamping(double damping)
      {
        van_Driest_damping_ = damping;
        return;
      };
      double VanDriestdamping() const { return van_Driest_damping_; };
      bool IncludeCi() const { return include_Ci_; };
      double ltau() const { return l_tau_; };
      //! flag to (de)activate fine-scale subgrid viscosity
      INPAR::FLUID::FineSubgridVisc Fssgv() const { return fssgv_; };
      // Flag to Vreman filter method
      INPAR::FLUID::VremanFiMethod Vrfi() const { return vrfi_; };
      /// multifractal subgrid-scales
      double Csgs() const { return Csgs_; };
      double CsgsPhi() const
      {
        double tmp = 0.0;
        if (not adapt_Csgs_phi_)
          tmp = Csgs_phi_;
        else
          tmp = Csgs_phi_ * meanCai_;
        return tmp;
      };
      double Alpha() const { return alpha_; };
      bool CalcN() const { return CalcN_; };
      double N() const { return N_; };
      enum INPAR::FLUID::RefVelocity RefVel() const { return refvel_; };
      enum INPAR::FLUID::RefLength RefLength() const { return reflength_; };
      double CNu() const { return c_nu_; };
      double CDiff() const { return c_diff_; };
      bool NearWallLimit() const { return near_wall_limit_; };
      bool NearWallLimitScatra() const { return near_wall_limit_scatra_; };
      bool BGp() const { return B_gp_; };
      double Beta() const { return beta_; };
      double MfsIsConservative() const { return mfs_is_conservative_; };
      double AdaptCsgsPhi() const { return adapt_Csgs_phi_; };
      void SetCsgsPhi(double meanCai)
      {
        meanCai_ = meanCai;
        return;
      };
      bool ConsistentMFSResidual() const { return consistent_mfs_residual_; };

      /*----------------------------------------------------*/
      //! @name loma parameters
      /*----------------------------------------------------*/

      //! flag for material update
      virtual bool UpdateMat() const { return update_mat_; };
      //! flag to (de)activate continuity SUPG term
      virtual bool ContiSUPG() const { return conti_supg_; };
      //! flag to (de)activate continuity cross-stress term -> residual-based VMM
      virtual INPAR::FLUID::CrossStress ContiCross() const { return conti_cross_; };
      //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
      virtual INPAR::FLUID::ReynoldsStress ContiReynolds() const { return conti_reynolds_; };
      //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
      virtual bool MultiFracLomaConti() const { return multifrac_loma_conti_; };

      /*----------------------------------------------------*/
      //! @name topopt parameters
      /*----------------------------------------------------*/

      //! flag to (de)aktivate reaction due to topology optimization
      bool ReactionTopopt() const { return reaction_topopt_; };
      //! flag to (de)aktivate reaction due to topology optimization
      const double* TopoptParams() const { return (double*)topopt_params_; };

     protected:
      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag SetGeneralParameter was called
      bool set_general_fluid_parameter_;

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, Poro)
      INPAR::FLUID::PhysicalType physicaltype_;
      //! parameter to switch the stabilization
      INPAR::FLUID::StabType stabtype_;
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate conservative formulation
      bool is_conservative_;
      //! flag to (de)activate Newton linearization
      bool is_newton_;
      //! flag to (de)activate second derivatives
      bool is_inconsistent_;
      //! flag to (de)activate potential reactive terms
      bool reaction_;
      //! function number of Oseen advective field
      int oseenfieldfuncno_;
      //! flag to activate consistent reconstruction of second derivatives
      bool is_reconstructder_;

      /*----------------------------------------------------*/
      //! @name stabilization parameters
      /*----------------------------------------------------*/

      /// parameter for residual based stabilizations
      //! Flag to (de)activate time-dependent subgrid stabilization
      INPAR::FLUID::SubscalesTD tds_;
      //! Flag to (de)activate time-dependent term in large-scale momentum equation
      INPAR::FLUID::Transient transient_;
      //! Flag to (de)activate PSPG stabilization
      bool pspg_;
      //! Flag to (de)activate SUPG stabilization
      bool supg_;
      //! Flag to (de)activate viscous term in residual-based stabilization
      INPAR::FLUID::VStab vstab_;
      //! Flag to (de)activate reactive term in residual-based stabilization
      INPAR::FLUID::RStab rstab_;
      //! Flag to (de)activate least-squares stabilization of continuity equation
      bool graddiv_;
      //! Flag to (de)activate cross-stress term -> residual-based VMM
      INPAR::FLUID::CrossStress cross_;
      //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
      INPAR::FLUID::ReynoldsStress reynolds_;
      //! Flag to define tau
      INPAR::FLUID::TauType whichtau_;
      //! Flag to define characteristic element length for tau_Mu
      INPAR::FLUID::CharEleLengthU charelelengthu_;
      //! Flag to define characteristic element length for tau_Mp and tau_C
      INPAR::FLUID::CharEleLengthPC charelelengthpc_;
      //! (sign) factor for viscous and reactive stabilization terms
      double viscreastabfac_;

      //! Flag to (de)activate PPP (polynomial pressure projection) stabilization
      bool ppp_;

      /// parameter for evaluation of material and stabilization parameter
      //! flag for material evaluation at Gaussian integration points
      bool mat_gp_;
      //! flag for stabilization parameter evaluation at Gaussian integration points
      bool tau_gp_;

      /*----------------------------------------------------*/
      //! @name two phase parameters
      /*----------------------------------------------------*/

      double interface_thickness_;
      bool enhanced_gaussrule_;
      bool include_surface_tension_;
      INPAR::TWOPHASE::SurfaceTensionApprox eval_surfacetension_;

      /*----------------------------------------------------*/
      //! @name turbulence model
      /*----------------------------------------------------*/

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
      //! flag to (de)activate fine-scale subgrid viscosity
      INPAR::FLUID::FineSubgridVisc fssgv_;
      // flag to Vreman filter method
      INPAR::FLUID::VremanFiMethod vrfi_;
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
      bool consistent_mfs_residual_;

      /*----------------------------------------------------*/
      //! @name loma parameters
      /*----------------------------------------------------*/

      //! flag for material update
      bool update_mat_;
      //! flag to (de)activate continuity SUPG term
      bool conti_supg_;
      //! flag to (de)activate continuity cross-stress term -> residual-based VMM
      INPAR::FLUID::CrossStress conti_cross_;
      //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
      INPAR::FLUID::ReynoldsStress conti_reynolds_;
      //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
      bool multifrac_loma_conti_;

      /*----------------------------------------------------*/
      //! @name topopt parameters
      /*----------------------------------------------------*/

      //! flag to (de)activate reaction due to topology optimization
      bool reaction_topopt_;
      //! matrix with values for computation of porosity with respect to topopt density
      double topopt_params_[3];

      /*-----------------------------------------------------*/

      /// private Constructor since we are a Singleton.
      FluidEleParameter();

     private:
      //! access time-integration parameters
      DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
