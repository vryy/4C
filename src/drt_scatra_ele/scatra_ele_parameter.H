/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter.H

\brief Setting of general scatra parameter for element evaluation

This file has to contain all parameters called in scatra_ele_calc.cpp.
Additional parameters required in derived classes of ScaTraEleCalc have to
be set in problem specific parameter lists derived from this class.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_H
#define SCATRA_ELE_PARAMETER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"


namespace DRT
{
namespace ELEMENTS
{

class ScaTraEleParameter : public DRT::SingletonDestruction
{
protected:

  /// (private) protected constructor, since we are a Singleton.
  /// this constructor is called from a derived class
  /// -> therefore, it has to be protected instead of private
  ScaTraEleParameter(const std::string& disname);

public:

  virtual ~ScaTraEleParameter() {}

  /// called upon destruction
  virtual void Done() = 0;

  /// In this class we do not define a static ScaTraEle...* Instance
  /// since only derived child classes are free to be allocated!!

  /*========================================================================*/
  //! @name set-routines
  /*========================================================================*/

  //! general fluid parameter are set
  virtual void SetElementGeneralParameters(Teuchos::ParameterList& params);

  //! turbulence parameters are set
  void SetElementTurbulenceParameters(Teuchos::ParameterList& params);   //> parameter list

  /*========================================================================*/
  //! @name access-routines
  /*========================================================================*/

  /*--------------------------*/
  //! @name general parameters
  /*--------------------------*/
  const bool IsAle(){return is_ale_;};
  const bool IsConservative(){return is_conservative_;};

  /*--------------------------*/
  //! @name query and output
  /*--------------------------*/

  const INPAR::SCATRA::FluxType WriteFlux(){return writeflux_;};
  Teuchos::RCP<std::vector<int> > WriteFluxIds(){return writefluxids_;};
  const INPAR::SCATRA::FDCheck FDCheck(){return fdcheck_;};
  const double FDCheckEps(){return fdcheckeps_;};
  const double FDCheckTol(){return fdchecktol_;};

  /*--------------------------*/
  //! @name stabilization parameters
  /*--------------------------*/

  const INPAR::SCATRA::StabType StabType(){return stabtype_;};
  const INPAR::SCATRA::TauType TauDef(){return whichtau_;};
  const INPAR::SCATRA::CharEleLength CharEleLength(){return charelelength_;};
  const double USFEMGLSFac(){return diffreastafac_;};
  const bool RBSubGrVel(){return sgvel_;};
  const bool ASSGD(){return assgd_;};
  const INPAR::SCATRA::AssgdType ASSGDType(){return whichassgd_;};
  const bool TauGP(){return tau_gp_;};
  const bool MatGP(){return mat_gp_;};

  /*--------------------------*/
  //! @name turbulence model
  /*--------------------------*/

  const INPAR::FLUID::TurbModelAction TurbModel(){return turbmodel_;};
  const INPAR::FLUID::ScalarForcing ScalarForcing(){return scalarforcing_;};
  const INPAR::SCATRA::FSSUGRDIFF WhichFssgd(){return whichfssgd_;};
  const bool FSSGD(){return fssgd_;};
  const double Cs(){return Cs_;};
  const double TPN(){return tpn_;};
  const bool CsAv(){return Cs_av_;};
  const double Csgs_SgVel(){return Csgs_sgvel_;};
  const double Alpha(){return alpha_;}
  const bool Calc_N(){return calc_N_;};
  const double N_Vel(){return N_vel_;};
  const INPAR::FLUID::RefVelocity RefVel(){return refvel_;};
  const INPAR::FLUID::RefLength RefLength(){return reflength_;};
  const double C_Nu(){return c_nu_;};
  const bool Nwl(){return nwl_;};
  const bool Nwl_ScaTra(){return nwl_scatra_;};
  const bool Beta(){return beta_;};
  const bool BD_Gp(){return BD_gp_;};
  const double Csgs_SgPhi(){double tmp = 0.0; if(adapt_Csgs_phi_ and nwl_) tmp = Csgs_sgvel_*meanCai_; else tmp = Csgs_sgphi_; return tmp;};
  const double C_Diff(){return c_diff_;};
  const bool MfsConservative(){return mfs_conservative_;};
  const void SetCsgsPhi(double meanCai){meanCai_ = meanCai; return;};
  const bool AdaptCsgsPhi(){return adapt_Csgs_phi_;};
  const bool TurbInflow(){return turbinflow_;};

protected:

  /*--------------------------*/
  //! @name general parameters
  /*--------------------------*/

  //! flag for ale
  bool is_ale_;
  //! flag for conservative form
  bool is_conservative_;

  /*--------------------------*/
  //! @name query and output
  /*--------------------------*/

  //! flag for writing the flux vector fields
  INPAR::SCATRA::FluxType writeflux_;
  //! ids of scalars for which flux vectors are written (starting with 1)
  Teuchos::RCP<std::vector<int> > writefluxids_;
  //! flag for finite difference check
  INPAR::SCATRA::FDCheck fdcheck_;
  //! perturbation magnitude for finite difference check
  double fdcheckeps_;
  //! relative tolerance for finite difference check
  double fdchecktol_;

  /*--------------------------*/
  //! @name stabilization parameters
  /*--------------------------*/

  //! type of stabilization
  INPAR::SCATRA::StabType stabtype_;
  //! definition of stabilization parameter
  INPAR::SCATRA::TauType whichtau_;
  //! definition of characteristic element length
  INPAR::SCATRA::CharEleLength charelelength_;
  //! parameter to switch between SUPG, GLS and USFEM
  double diffreastafac_;
  //! flag to include residual-based subgrid-scale velocity
  bool sgvel_;
  //! flag to active artificial diffusion
  bool assgd_;
  //! definition of artificial diffusion
  INPAR::SCATRA::AssgdType whichassgd_;
  //! flag for evaluation of tau at Gauss point
  bool tau_gp_;
  //! flag for evaluation of material at Gauss point
  bool mat_gp_;

  /*--------------------------*/
  //! @name turbulence model
  /*--------------------------*/

  //! definition of turbulence model
  INPAR::FLUID::TurbModelAction turbmodel_;
  //! define forcing for scalar field
  INPAR::FLUID::ScalarForcing scalarforcing_;
  //! flag to activate AVM3
  bool fssgd_;
  INPAR::SCATRA::FSSUGRDIFF whichfssgd_;
  // parameters for subgrid-diffusivity models
  double Cs_;
  double tpn_;
  bool Cs_av_;
  // parameters for multifractal subgrid-scale modeling
  double Csgs_sgvel_;
  double alpha_;
  bool calc_N_;
  double N_vel_;
  INPAR::FLUID::RefVelocity refvel_;
  INPAR::FLUID::RefLength reflength_;
  double c_nu_;
  bool nwl_;
  bool nwl_scatra_;
  bool beta_;
  bool BD_gp_;
  double Csgs_sgphi_;
  double c_diff_;
  bool mfs_conservative_;
  double meanCai_;
  bool adapt_Csgs_phi_;
  // further parameter
  bool turbinflow_;

private:

  //! access time-integration parameters
  DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;

};

}
}

#endif
