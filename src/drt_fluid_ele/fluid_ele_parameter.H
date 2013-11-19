/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_std.cpp

\brief Evaluation of general fluid parameter for standard fluid

       Poro specific parameters are defined in a derived class.
       Since there are a couple of terms from the stdfluid called
       by the porofluid this is  inevitable.

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

#include "fluid_ele_parameter_base.H"


namespace DRT
{
namespace ELEMENTS
{

class FluidEleParameter : public FluidEleParameterBase
{

public:

  virtual ~FluidEleParameter() {}

  //! general fluid parameter are set
  void SetElementStdFluidParameter(Teuchos::ParameterList& params,   //> parameter list
                                   int myrank                        //> proc id
                                  );

  //! turbulence parameters are set
  void SetElementTurbulenceParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! topopt parameters are set
  void SetElementTopoptParameter(Teuchos::ParameterList& params    //> parameter list
      );

  /// parameter for residual stabilization
  //! Flag to (de)activate time-dependent subgrid stabilization
  INPAR::FLUID::SubscalesTD Tds(){return tds_;};
  //! Flag to (de)activate time-dependent term in large-scale momentum equation
  INPAR::FLUID::Transient Transient(){return transient_;};
  //! Flag to (de)activate PSPG stabilization
  bool PSPG(){return pspg_;};
  //! Flag to (de)activate SUPG stabilization
  bool SUPG(){return supg_;};
  //! Flag to (de)activate viscous term in residual-based stabilization
  INPAR::FLUID::VStab VStab(){return vstab_;};
  //! Flag to (de)activate reactive term in residual-based stabilization
  INPAR::FLUID::RStab RStab(){return rstab_;};
  //! Flag to (de)activate least-squares stabilization of continuity equation
  bool CStab(){return graddiv_;};
  //! Flag to (de)activate cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress Cross(){return cross_;};
  //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress Reynolds(){return reynolds_;};
  //! Flag to define tau
  INPAR::FLUID::TauType WhichTau(){return whichtau_;};
  //! Flag to define characteristic element length for tau_Mu
  INPAR::FLUID::CharEleLengthU CharEleLengthU(){return charelelengthu_;};
  //! Flag to define characteristic element length for tau_Mp and tau_C
  INPAR::FLUID::CharEleLengthPC CharEleLengthPC(){return charelelengthpc_;};
  //! flag to (de)activate fine-scale subgrid viscosity
  INPAR::FLUID::FineSubgridVisc Fssgv(){return fssgv_;};
  // Flag to Vreman filter method
  INPAR::FLUID::VremanFiMethod Vrfi(){return vrfi_;};
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


protected:

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
  //! flag to (de)activate fine-scale subgrid viscosity
  INPAR::FLUID::FineSubgridVisc fssgv_;
  // flag to Vreman filter method
  INPAR::FLUID::VremanFiMethod vrfi_;
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

  /// private Constructor since we are a Singleton.
  FluidEleParameter();

public:

  //! some loma specific parameters
  //! flag for material update
  virtual bool UpdateMat(){dserror("must be overridden by subclass"); return false;};
  //! flag to (de)activate continuity SUPG term
  virtual bool ContiSUPG(){dserror("must be overridden by subclass"); return false;};
  //! flag to (de)activate continuity cross-stress term -> residual-based VMM
  virtual INPAR::FLUID::CrossStress ContiCross(){dserror("must be overridden by subclass"); return (INPAR::FLUID::CrossStress)0;};
  //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
  virtual INPAR::FLUID::ReynoldsStress ContiReynolds(){dserror("must be overridden by subclass"); return (INPAR::FLUID::ReynoldsStress)0;};
  //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
  virtual bool MultiFracLomaConti(){dserror("must be overridden by subclass"); return false;};

};

}
}

#endif
