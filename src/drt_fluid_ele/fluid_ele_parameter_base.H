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


#ifndef FLUID_ELE_PARAMETER_BASE_H
#define FLUID_ELE_PARAMETER_BASE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"

#include "../drt_lib/drt_dserror.H"
#include "fluid_ele_parameter_interface.H"
#include "fluid_ele_parameter_timint.H"

namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters (constant over time)
class FluidEleParameterBase : public FluidEleParameterInterface
{
public:

  /// destructor
  virtual ~FluidEleParameterBase() {}

  //! general fluid parameter are set
  void SetElementGeneralFluidParameter(Teuchos::ParameterList& params,   //> parameter list
                                       int myrank                        //> proc id
                                      );

  //! print parameter to screen
  void PrintFluidGeneralParameter();

  //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq, poro)
  INPAR::FLUID::PhysicalType PhysicalType(){return physicaltype_;};
  //! get the stabtype
  INPAR::FLUID::StabType StabType(){return stabtype_;};
  //! flag to (de)activate conservative formulation
  bool IsConservative(){return is_conservative_;};
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
  //! flag to (de)aktivate reaction due to topology optimization
  const double* TopoptParams(){return (double*)topopt_params_;};

protected:

  /// protected Constructor since we are a Singleton.
  FluidEleParameterBase();

  DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;

  //! Flag SetGeneralParameter was called
  bool set_general_fluid_parameter_;

  //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq, Poro)
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
  //! flag to (de)aktivate porous darcy flow
  bool darcy_;
  //! flag to (de)aktivate reaction due to topology optimization
  bool reaction_topopt_;
  //! matrix with values for computation of porosity with respect to topopt density
  double topopt_params_[3];

  /// dummy call method section to make std parameters available for return type FluidEleParameter
public:
  /// parameter for residual stabilization
  //! Flag to (de)activate time-dependent subgrid stabilization
  virtual INPAR::FLUID::SubscalesTD Tds()=0;
  //! Flag to (de)activate time-dependent term in large-scale momentum equation
  virtual INPAR::FLUID::Transient Transient()=0;
  //! Flag to (de)activate PSPG stabilization
  virtual bool PSPG()=0;
  //! Flag to (de)activate SUPG stabilization
  virtual bool SUPG()=0;
  //! Flag to (de)activate viscous term in residual-based stabilization
  virtual INPAR::FLUID::VStab VStab()=0;
  //! Flag to (de)activate reactive term in residual-based stabilization
  virtual INPAR::FLUID::RStab RStab()=0;
  //! Flag to (de)activate least-squares stabilization of continuity equation
  virtual bool CStab()=0;
  //! Flag to (de)activate cross-stress term -> residual-based VMM
  virtual INPAR::FLUID::CrossStress Cross()=0;
  //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
  virtual INPAR::FLUID::ReynoldsStress Reynolds()=0;
  //! Flag to define tau
  virtual INPAR::FLUID::TauType WhichTau()=0;
  //! Flag to define characteristic element length for tau_Mu
  virtual INPAR::FLUID::CharEleLengthU CharEleLengthU()=0;
  //! Flag to define characteristic element length for tau_Mp and tau_C
  virtual INPAR::FLUID::CharEleLengthPC CharEleLengthPC()=0;
  //! flag to (de)activate fine-scale subgrid viscosity
  virtual INPAR::FLUID::FineSubgridVisc Fssgv()=0;
  // Flag to Vreman filter method
  virtual INPAR::FLUID::VremanFiMethod Vrfi()=0;
  //! (sign) factor for viscous and reactive stabilization terms
  virtual double ViscReaStabFac()=0;

  /// parameter for edge-based stabilizations
  //! Flag to (de)activate pressure stabilization
  virtual INPAR::FLUID::EOS_Pres EOS_Pres()=0;
  //! Flag to (de)activate convective streamline stabilization
  virtual INPAR::FLUID::EOS_Conv_Stream EOS_Conv_Stream()=0;
  //! Flag to (de)activate convective crosswind stabilization
  virtual INPAR::FLUID::EOS_Conv_Cross EOS_Conv_Cross()=0;
  //! Flag to (de)activate divergence stabilization
  virtual INPAR::FLUID::EOS_Div EOS_Div()=0;
  //! Flag to define element length
  virtual INPAR::FLUID::EOS_ElementLength EOS_element_length()=0;
  //! Flag to define tau for edge-based stabilization
  virtual INPAR::FLUID::EOS_TauType EOS_WhichTau()=0;

  //! flag for material evaluation at Gaussian integration points
  virtual bool MatGp()=0;
  //! flag for stabilization parameter evaluation at Gaussian integration points
  virtual bool TauGp()=0;

  /// constant parameters for the turbulence formulation
  /// subgrid viscosity models
  //! flag to define turbulence model
  virtual INPAR::FLUID::TurbModelAction TurbModAction()=0;
  virtual double Cs()=0;
  virtual bool CsAveraged()=0;
  virtual double Ci()=0;
  virtual void SetvanDriestdamping(double damping)=0;
  virtual double VanDriestdamping()=0;
  virtual bool IncludeCi()=0;
  virtual double ltau()=0;
  //! scale similarity constant
  virtual double Cl()=0;
  /// multifractal subgrid-scales
  virtual double Csgs()=0;
  virtual double CsgsPhi()=0;
  virtual double Alpha()=0;
  virtual bool CalcN()=0;
  virtual double N()=0;
  virtual enum INPAR::FLUID::RefVelocity RefVel()=0;
  virtual enum INPAR::FLUID::RefLength RefLength()=0;
  virtual double CNu()=0;
  virtual double CDiff()=0;
  virtual bool NearWallLimit()=0;
  virtual bool NearWallLimitScatra()=0;
  virtual bool BGp()=0;
  virtual double Beta()=0;
  virtual double MfsIsConservative()=0;
  virtual double AdaptCsgsPhi()=0;
  virtual void SetCsgsPhi(double meanCai)=0;

  /// some loma specific parameters
  virtual bool UpdateMat()=0;
  virtual bool ContiSUPG()=0;
  virtual INPAR::FLUID::CrossStress ContiCross()=0;
  virtual INPAR::FLUID::ReynoldsStress ContiReynolds()=0;
  virtual bool MultiFracLomaConti()=0;

};

}
}

#endif
