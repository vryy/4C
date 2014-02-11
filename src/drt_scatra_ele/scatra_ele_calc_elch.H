/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.H

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_H
#define SCATRA_ELE_CALC_ELCH_H

#include "scatra_ele_calc.H"

#include "../drt_inpar/inpar_elch.H"

// TODO: remove when ScaTraEleInternalVariableManagerElch moved to own file
#include "scatra_ele_parameter_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;
  class ScaTraEleDiffManagerElch;
  template<int NSD, int NEN>
  class ScaTraEleInternalVariableManagerElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElch : public ScaTraEleCalc<distype>
{
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElch(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcElch() {}

  /// Singleton access method
  static ScaTraEleCalcElch<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

  //! evaluate the element
  int Evaluate(DRT::ELEMENTS::Transport*     ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               const std::vector<int> &      lm,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                      Teuchos::ParameterList&       params,
                      DRT::Discretization &         discretization,
                      const std::vector<int> &      lm,
                      Epetra_SerialDenseMatrix&     elemat1_epetra,
                      Epetra_SerialDenseMatrix&     elemat2_epetra,
                      Epetra_SerialDenseVector&     elevec1_epetra,
                      Epetra_SerialDenseVector&     elevec2_epetra,
                      Epetra_SerialDenseVector&     elevec3_epetra);

private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! calculate matrix and rhs. Here the whole thing is hidden.
  void Sysmat(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
  );

  void CalMatAndRhsNP(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    const int&                             k
  );

  void CalMatAndRhsDiffCond(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&               erhs,         //!< element rhs to calculate+
    const int&                               k,            //!< index of current scalar
    const double&                            fac,          //!< domain-integration factor
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,           //!< diffusion manager
    const double&                            rhsint        //!< rhs at Gauss point
  );

  void CalMatAndRhsOutsideScalarLoop(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&               erhs,         //!< element rhs to calculate
    const double&                            fac,          //!< domain-integration factor
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme          //!< diffusion manager
  );

  void CorrectionForFluxAccrosDC(
    DRT::Discretization&        discretization,
    const std::vector<int>&     lm,
    Epetra_SerialDenseMatrix&   emat,
    Epetra_SerialDenseVector&   erhs
  );

  /*========================================================================*/
  //! @name material and related and related functions
  /*========================================================================*/

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );


  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatPhase(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               iphase,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatNewman(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  void MatIon(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,     //!< fluid viscosity
    const int                               iquad     //!< id of current gauss point
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/

  //! CalcMat: Linearization of diffusion coefficient in diffusion term
  void CalcMatDiffCoeffLin(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int&                               k,            //!< index of current scalar
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi       //!< gradient of concentration at GP
  );

  //! CalcMat: Conduction term with inserted current - ohmic overpotential
  void CalcMatCondOhm(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int&                               k,            //!< index of current scalar
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot       //!< gradient of potenial at GP
  );

  //! CalcMat: Conduction term with inserted current - conc. overpotential
  void CalcMatCondConc(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int&                               k,              //!< index of current scalar
    const double&                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double&                            rtffcval,       //!< RT/F^2/Newman_const_c/z_k
    const double&                            newman_const_a, //!< Newman constant a
    const double&                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,            //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,        //!< gradient of concentration at GP
    const std::vector<double>&               conintinv       //!< inverted concentration at GP
  );

  //! CalcMat: Conduction term without inserted current
  void CalcMatCond(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int&                               k,            //!< index of current scalar
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        curint        //!< current at GP
  );

  //! CalcMat: Additional diffusion term without inserted current
  void CalcMatCondDiff(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int&                               k,            //!< index of current scalar
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const std::vector<LINALG::Matrix<my::nsd_,1> >&  gradphi       //!< gradient of concentration at GP
  );

  //! CalcMat: Potential equation div i inserted current - ohmic overpotential
  void CalcMatPotEquDiviOhm(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const double&                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double&                            invf,         //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot       //!< gradient of potenial at GP
  );

  //! Potential equation div i inserted current - concentration overpotential
  void CalcMatPotEquDiviConc(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int&                               k,              //!< index of current scalar
    const double&                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double&                            rtffc,          //!< RT/(F^2 Newman_const_c)
    const double&                            rtf,            //!< RT/F
    const double&                            invf,           //!< 1/F
    const double&                            newman_const_a, //!< Newman constant a
    const double&                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,            //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,        //!< gradient of concentration at GP
    const double&                            conintinv       //!< inverted concentration at GP
  );

  //! Potential equation ENC
  void CalcMatPotEquENC(
    Epetra_SerialDenseMatrix&               emat,   //!< element matrix to be filled
    const int&                               k,      //!< index of current scalar
    const double&                            fac,    //!< domain-integration factor
    const double&                            alphaf, //!< time factor for ENC
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme     //!< diffusion manager
  );

  //! Potential equation div i without inserted current
  void CalcMatPotEquDivi(
    Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
    const double&                timefacfac, //!< domain-integration factor times time-integration factor
    const double&                invf        //!< 1/F
  );

  //! CalcMat: Current equation current
  void CalcMatCurEquCur(
    Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
    const double&                timefacfac, //!< domain-integration factor times time-integration factor
    const double&                invf        //!< 1/F
  );

  //! CalcMat: Current equation ohmic overpotential
  void CalcMatCurEquOhm(
      Epetra_SerialDenseMatrix&               emat,       //!< element matrix to be filled
      const double&                            timefacfac, //!< domain-integration factor times time-integration factor
      const double&                            invf,       //!< 1/F
      Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,        //!< diffusion manager
      const LINALG::Matrix<my::nsd_,1>&              gradpot     //!< gradient of potenial at GP
  );

  //! CalcMat: Current equation concentration overpotential
  void CalcMatCurEquConc(
    Epetra_SerialDenseMatrix&               emat,             //!< element matrix to be filled
    const double&                            timefacfac,       //!< domain-integration factor times time-integration factor
    const double&                            rtf,              //!< RT/F
    const double&                            rtffc,            //!< RT/(F^2 Newman_const_c)
    const std::vector <double>&              invfval,          //!< 1/(F z_k)
    const double&                            newman_const_a,   //!< Newman constant a
    const double&                            newman_const_b,   //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,              //!< diffusion manager
    const std::vector< LINALG::Matrix<my::nsd_,1> >& gradphi,  //!< gradient of concentration at GP
    const std::vector<double>&                     conintinv   //!< inverted concentration at GP
  );

  //! CalcRhs: Conduction term with inserted current - ohmic overpotential
  void CalcRhsCondOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            invfval,  //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot   //!< gradient of potenial at GP
  );

  //! CalcRhs: Conduction term with inserted current - conc. overpotential
  void CalcRhsCondConc(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            rtffcval, //!< RT/(F^2 Newman_const_c z_k)
    const double&                            newman_const_a, //!< Newman constant a
    const double&                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,  //!< gradient of concentration at GP
    const std::vector<double>&               conintinv //!< inverted concentration at GP
  );

  //! CalcRhs: Conduction term without inserted current
  void CalcRhsCond(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&       curint        //!< current at GP
  );

  //! CalcRhs: Additional diffusion term without inserted current
  void CalcRhsCondDiff(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    const std::vector< LINALG::Matrix<my::nsd_,1> >&  gradphi   //!< gradient of concentration at GP
  );

  //! CalcRhs: Potential equation div i inserted current - ohmic overpotential
  void CalcRhsPotEquDiviOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&              gradpot   //!< gradient of potenial at GP
  );

  //! CalcRhs: Potential equation div i inserted current - conc. overpotential
  void CalcRhsPotEquDiviConc(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            rtf,      //!< RT/F
    const std::vector<double>&               invfval,  //!< 1/(F z_k)
    const double&                            rtffc, //!< RT/(F^2 Newman_const_c)
    const double&                            newman_const_a, //!< Newman constant a
    const double&                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,  //!< gradient of concentration at GP
    const double&                            conintinv //!< inverted concentration at GP
  );

  //! CalcRhs: Potential equation ENC
  void CalcRhsPotEquENC(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int&                               k,        //!< index of current scalar
    const double&                            fac,      //!< domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const double&                            conint    //!< concentration at GP
  );

  //! CalcRhs: Potential equation divi without inserted current
  void CalcRhsPotEquDivi(
    Epetra_SerialDenseVector&   erhs,     //!< element vector to be filled
    const double&                rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                invf,     //!< 1/F
    const LINALG::Matrix<my::nsd_,1> &  curint    //!< current at GP
  );

  //! CalcRhs: Current equation - current
  void CalcRhsCurEquCur(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1> &             curint    //!< current at GP
  );

  //! CalcRhs: Current equation - ohmic overpotential
  void CalcRhsCurEquOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double&                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&              gradpot   //!< gradient of potenial at GP
  );

  //! Current equation - concentration overpotential
  void CalcRhsCurEquConc(
    Epetra_SerialDenseVector&                 erhs,     //!< element vector to be filled
    const double&                              rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double&                              rtf,      //!< RT/F
    const std::vector<double>&                 invfval,  //!< 1/(F z_k)
    const double&                              rtffc,    //!< RT/(F^2 Newman_const_c)
    const double&                              newman_const_a, //!< Newman constant a
    const double&                              newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElch>&    dme,      //!< diffusion manager
    const std::vector<LINALG::Matrix<my::nsd_,1> >&  gradphi,  //!< vector of gradient of concentration at GP
    const std::vector<double>&                       conintinv //!< inverted concentration at GP
   );

  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  void CheckElchElementParameter(
    DRT::ELEMENTS::Transport*  ele
    );

  void CalculateFlux(
    LINALG::Matrix<3,my::nen_>&     flux,
    const DRT::Element*             ele,
    INPAR::ELCH::ElchType           elchtype,
    const INPAR::SCATRA::FluxType   fluxtype,
    const int                       k
    );

  void CalErrorComparedToAnalytSolution(
    const DRT::Element*                   ele,
    const enum INPAR::ELCH::ElchType      elchtype,
    Teuchos::ParameterList&               params,
    Epetra_SerialDenseVector&             errors
  );

  void CalculateConductivity(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseVector&         sigma
  );

  void CalculateElectricPotentialField(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  );

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! pointer to parameter class elch
  DRT::ELEMENTS::ScaTraEleParameterElch* elchpara_;

  //! flag for used element formulation (material based)
  INPAR::ELCH::DiffCondMat diffcondmat_;

  //! flag if current is used as a solution variable
  bool cursolvar_;

  // flag if elch system is closed by ENC or div i
  bool equpot_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! electric potential
  LINALG::Matrix<my::nen_,1> epotnp_;

  //! nodal current density values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<my::nsd_,my::nen_> ecurnp_;

  //TODO: ELCH: cleaning
  //! flag for migration velocity in the stabilization operator (ELCH)
  bool migrationstab_;
  //! flag for migration velocity in formula for stabilization parameter (ELCH)
  bool migrationintau_;
  //! flag for migration velocity in residual of strong form (ELCH)
  bool migrationinresidual_;
  //! transference numbers
  double transelim_;
  std::vector<double> diffuselimderiv_;
  double diffuselim_;
  std::vector<double> eps_;
  std::vector<double> tort_;
  std::vector<double> epstort_;
};


/// ScaTraEleDiffManagerElch implementation
/*!
  This class keeps all elch-specific transport parameter needed for the evaluation of an element. The
  ScaTraEleDiffManagerElch is derived from the standard ScaTraEleDiffManager.
  for
*/
class ScaTraEleDiffManagerElch : public ScaTraEleDiffManager
{
  public:

    ScaTraEleDiffManagerElch(int numscal):
    ScaTraEleDiffManager(numscal),
    valence_(numscal,0.0),
    derivdiffcoef_(numscal,std::vector<double>(numscal,0.0)),
    transnum_(numscal,0.0),
    derivtransnum_(numscal,std::vector<double>(numscal,0.0)),
    cond_(0.0),
    derivcond_(numscal,0.0),
    thermfac_(1.0),
    derivthermfac_(numscal,0.0)
    {return;}

    ~ScaTraEleDiffManagerElch(){};

    //! Set valence of the single ionic species
    void SetValence(const double valence,const int k)
    {
      valence_[k]=valence;
      return;
    };

    //! Access routine for valence of the single ionic species
    const double GetValence(const int k){return valence_[k];};

    //! Set derivative of diffusion coefficients with respect to concentrations
    void SetDerivIsoDiffCoef(const double derivdiffcoef,const int k,const int iscal)
    {
      (derivdiffcoef_[k])[iscal] = derivdiffcoef;
      return;
    };

    //! Access routine for derivative of diffusion coefficients with respect to concentrations
    const double GetDerivIsoDiffCoef(const int k,const int iscal){return (derivdiffcoef_[k])[iscal];};

    //! Set transference numbers with respect to single ionic species
    void SetTransNum(const double transnum,const int k)
    {
      transnum_[k] = transnum;
      return;
    };

    //! Access routine for transference numbers with respect to single ionic species
    const double GetTransNum(const int k){return transnum_[k];};

    //! Set derivative of transference numbers with respect to concentrations
    void SetDerivTransNum(const double derivtransnum,const int k,const int iscal)
    {
      (derivtransnum_[k])[iscal] = derivtransnum;
      return;
    };

    //! Access routine for derivative of transference numbers with respect to concentrations
    const double GetDerivTransNum(const int k,const int iscal){return (derivtransnum_[k])[iscal];};

    //! Set conductivity of the electrolyte solution
    void SetCond(const double cond)
    {
      cond_ = cond;
      return;
    };

    //! Access routine for conductivity of the electrolyte solution
    const double GetCond(){return cond_;};

    //! Set derivative of the conductivity with respect to concentrations
    void SetDerivCond(const double derivcond,const int k)
    {
      derivcond_[k] = derivcond;
      return;
    };

    //! Access routine for derivative of the conductivity with respect to concentrations
    const double GetDerivCond(const int k){return derivcond_[k];};

    //! Set thermodynamic factor for a specific electrolyte solution
    void SetThermFac(const double thermfac)
    {
      thermfac_ = thermfac;
      return;
    };

    //! Access routine for the thermodynamic factor for a specific electrolyte solution
    const double GetThermFac(){return thermfac_;};

    //! Set derivative of thermodynamic factor with respect to concentrations
    void SetDerivThermFac(const double derivthermfac,const int k)
    {
      derivthermfac_[k] = derivthermfac;
      return;
    };

    //! Access routine for derivative of thermodynamic factor with respect to concentrations
    const double GetDerivThermFac(const int k){return derivthermfac_[k];};

    void CalcConductivity(const int numscal, const double ffrt, const std::vector<double>& conint)
    {
      double cond = 0.0;
      for(int ispec = 0;ispec<numscal;++ispec)
      {
        // conductivity
        cond += ffrt*valence_[ispec]*valence_[ispec]*diff_[ispec]*conint[ispec];
        // derivation of conductivity wrt concentrations
        derivcond_[ispec]=ffrt*valence_[ispec]*valence_[ispec]*diff_[ispec];
      }

      cond_=cond;

      return;
    };

    void CalcTransNum(const int numscal, const std::vector<double>& conint)
    {
      // conductivity without ffrt
      double sum=0.0;
      for(int k=0; k<numscal;++k)
      {
        sum += valence_[k]*valence_[k]*diff_[k]*conint[k];
      }
      double denomin = sum*sum;

      for(int k = 0;k<numscal;++k)
      {
        transnum_[k]=valence_[k]*valence_[k]*diff_[k]*conint[k]/sum;

        for(int iscal=0; iscal<numscal;++iscal)
        {
          if(k==iscal)
            (derivtransnum_[k])[iscal]=(valence_[k]*valence_[k]*diff_[k]*(sum-valence_[k]*valence_[k]*diff_[k]*conint[k]))/denomin;
          else
            (derivtransnum_[k])[iscal]=(-valence_[k]*valence_[k]*diff_[k]*conint[k]*valence_[iscal]*valence_[iscal]*diff_[iscal])/denomin;
        }
      }
      return;
    };

    //! Output of transport parameter (to screen)
    void OutputTransportParams(const int numscal)
    {
      for(int k=0;k<numscal;++k)
        std::cout << "valence " << k << ":   " << valence_[k] << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "diffusion coefficient " << k << ":   " << diff_[k] << std::endl;

      for(int k=0;k<numscal;++k)
      {
        for(int iscal=0;iscal<numscal;++iscal)
        {
          std::cout << "derivation diffusion coefficient (" << k << "," << iscal << "):  " << (derivdiffcoef_[k])[iscal] << std::endl;
        }
      }
      std::cout << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "transference number " << k << ":   " << transnum_[k] << std::endl;

      for(int k=0;k<numscal;++k)
      {
        for(int iscal=0;iscal<numscal;++iscal)
        {
          std::cout << "derivation transference number (" << k << "," << iscal << "):  " << (derivtransnum_[k])[iscal] << std::endl;
        }
      }
      std::cout << std::endl;

      std::cout << "conductivity:   " << cond_ << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "derivation of conductivity  " << k << ":   " << derivcond_[k] << std::endl;
      std::cout << std::endl;

      std::cout << "thermodynamic factor:   " << thermfac_ << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "derivation of thermodynamic factor " << k << ":   " << derivthermfac_[k] << std::endl;
      std::cout << std::endl;

      return;
    };

    // return
    unsigned int LengthVector(){return diff_.size();};

    void IncreaseLengthVector(const int k, const int numscal)
    {
      if(diff_.size()== (unsigned) numscal)
      {
        valence_.push_back(0.0);
        diff_.push_back(0.0);
      }

      return;
    };


  private:

    //! valence of the single ionic species
    std::vector<double> valence_;

    //! derivative of diffusion coefficients with respect to concentrations
    std::vector< std::vector<double> > derivdiffcoef_;

    //! transference numbers for single ionic species
    std::vector<double> transnum_;

    //! derivative of transference numbers with respect to concentrations
    std::vector< std::vector<double> > derivtransnum_;

    //! conductivity of the electrolyte solution
    double cond_;

    //! derivative of the conductivity with respect to concentrations
    std::vector<double> derivcond_;

    //! thermodynamic factor for a specific electrolyte solution
    //! transport parameter used only in the Newman model
    double thermfac_;

    //! derivative of thermodynamic factor with respect to concentrations
    std::vector<double> derivthermfac_;

};

/// ScaTraEleInternalVariableManagerElch implementation
/*!
  This class keeps all internal variables needed for the evaluation of an element.
*/
template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerElch
{
  public:

    ScaTraEleInternalVariableManagerElch(int numscal, int nsd,DRT::ELEMENTS::ScaTraEleParameterElch* elchpara):
      numscal_(numscal),
      nsd_(nsd),
      // constant internal variables
      faraday_(INPAR::ELCH::faraday_const),
      epsilon_(INPAR::ELCH::epsilon),
      frt_(elchpara->FRT()),
      // TODO: diffcond variables
      invf_(1/faraday_),
      rtf_(1/frt_),
      rtffc_(rtf_*invf_/elchpara->NewmanConstC()),
      invval_(numscal_),
      invfval_(numscal_),
      rtffcval_(numscal_),
      //internal variables evaluated at the Gauss point
      gradphi_(numscal_),
      conint_(numscal_),
      conintinv_(numscal_),
      gradpot_(true),
      //velint_(true),
      convelint_(true),
      conv_(true),
      convphi_(numscal),
      sgconv_(true),
      // TODO: diffcond variables
      curint_(true),
      curdiv_(0.0)
    {return;}

    ~ScaTraEleInternalVariableManagerElch(){};

    void SetInternalVariables(
      Teuchos::RCP<ScaTraEleDiffManagerElch>  dme,      //!< diffusion manager
      LINALG::Matrix<NEN,1>                 funct,      //!< array for shape functions
      LINALG::Matrix<NSD,NEN>               derxy,      //!< global derivatives of shape functions w.r.t x,y,z
      std::vector<LINALG::Matrix<NEN,1> >   ephinp,     //!< nodal concentration values at time n+1
      LINALG::Matrix<NEN,1>                 epotnp,     //!< nodal electric potential values at time n+1
      LINALG::Matrix<NSD,NEN>               ecurnp,     //!< nodal current density values at t_(n+1) or t_(n+alpha_F)
      LINALG::Matrix<NSD,NEN>               econvelnp   //!< nodal velocity values at t_(n+1) or t_(n+alpha_F)
    )
    {
      // constant material-depending variables
      for(int k=0; k<numscal_; ++k)
      {
        invval_[k]=1.0/dme->GetValence(k);
        invfval_[k]=invf_/dme->GetValence(k);
        rtffcval_[k]=rtffc_/dme->GetValence(k);
      }

      // loop over all transported scalars
      // get concentration of transported scalar k at integration point
      // evaluation of all concentrations is necessary at this point since
      // -> homogeneous reactions of scalar k may depend on all concentrations
      // -> concentration depending material parameters for the diffusion-confection formulation
      // -> avoiding of possible errors (concentration was always defined as a vector where only on
      //    entry was filled)
      for (int k = 0; k < numscal_;++k)
      {
        // calculate concentration at GP
        conint_[k] = funct.Dot(ephinp[k]);
        // calculate 1/concentration at GP
        conintinv_[k] = 1/conint_[k];
        // calculatd gradient of concentration at GP
        gradphi_[k].Multiply(derxy,ephinp[k]);
      }
      // calculate gradient of electric potential at GP
      gradpot_.Multiply(derxy,epotnp);

      // calculate velocity at GP
      //velint_.Multiply(evelnp,funct);
      // calculate convective velocity at GP
      convelint_.Multiply(econvelnp,funct);
      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      conv_.MultiplyTN(derxy,convelint_);

      // TODO: Elch variables
      // convective term using current scalar value
      for (int k = 0; k < numscal_;++k)
      {
        //convphi_[k] = convelint_.Dot(gradphi_[k]);
        convphi_[k] = conv_.Dot(ephinp[k]);
      }

      // TODO: diffcond variables
      curint_.Multiply(ecurnp,funct);
      // Calculate divergence of current desity
      // since it is a summation it is important to set curdiv to zero
      curdiv_=0.0;
      for (int idim = 0; idim <nsd_; ++idim)
      {
        LINALG::Matrix<NSD,NSD> curderxy(true);
        curderxy.MultiplyNT(ecurnp,derxy);
        curdiv_ += curderxy(idim, idim);
      }

      return;
    };

  /*========================================================================*/
  //! @name return constant internal variables
  /*========================================================================*/

    const double Faraday(){return faraday_;};
    const double Epsilon(){return epsilon_;};
    const double FRT(){return frt_;};
    // TODO: diffcond variables
    const double InvF(){return invf_;};
    const double RTF(){return rtf_;};
    const double RTFFC(){return rtffc_;};
    const std::vector <double> InvVal(){return invval_;};
    const double InvFVal(const int k){return invfval_[k];};
    const std::vector <double> InvFVal(){return invfval_;};
    const double RTFFCVal(const int k){return rtffcval_[k];};
    const std::vector <double> RTFFCVal(){return rtffcval_;};

  /*========================================================================*/
  //! @name return methods for GP values
  /*========================================================================*/
    // TODO: general variables
    //! return gradient of phi for species k at GP
    const LINALG::Matrix<NSD,1> GradPhi(const int k){return gradphi_[k];};
    //! return gradient of phi for all species in a vector at GP
    const std::vector<LINALG::Matrix<NSD,1> > GradPhi(){return gradphi_;};
    //! return concentration of species k at GP
    const double ConInt(const int k){return conint_[k];};
    //! return 1/concentration of species k at GP
    const double ConIntInv(const int k){return conintinv_[k];};
    //! return 1/concentration of all species in a vector at GP
    const std::vector<double> ConIntInv(){return conintinv_;};
    //! return gradient of electric potential at GP
    const LINALG::Matrix<NSD,1> GradPot(){return gradpot_;};
    //! return convective part in convective form
    const LINALG::Matrix<NEN,1> Conv(){return conv_;};
    //! return convective contribution
    const double ConvPhi(const int k ){return convphi_[k];};
    //! return subgrid velocity
    const LINALG::Matrix<NEN,1> SGConv(){return sgconv_;};

    // TODO: diffcond variables
    //! return current density at GP
    const LINALG::Matrix<NSD,1> CurInt(){return curint_;};
    //! return divergence of current density
    const double CurDiv(){return curdiv_;};

  private:

    //! number of scalars
    int numscal_;
    //! number of space dimensions
    int nsd_;

  /*========================================================================*/
  //! @name constant internal variables
  /*========================================================================*/

    // TODO: diffcond variables
    //! Faraday constant
    double faraday_;
    //! di
    double epsilon_;
    //! pre-calculation of divisions which are regularly used in following
    //! (a division is much more expensive than a multiplication)
    //! constant parameter 1/F
    double frt_;
    //! constant parameter F/RT
    double invf_;
    //! constant parameter RT/F
    double rtf_;
    //! constant parameter RT/(F^2*Newman_const_C)
    //! attention: this is a newman specific parameter
    double rtffc_;
    //! constant parameter 1/(z_k)
    std::vector <double> invval_;
    //! constant parameter 1/(F z_k)
    std::vector <double> invfval_;
    //! constant parameter RT/(F^2*Newman_const_C*z_k)
    //! attention: this is a newman specific parameter
    std::vector <double> rtffcval_;

  /*========================================================================*/
  //! @name internal variables evaluated at the Gauss point
  /*========================================================================*/

    // TODO: general variables
    //! gradient of concentration at GP
    std::vector<LINALG::Matrix<NSD,1> > gradphi_;
    //! concentration at GP
    std::vector<double> conint_;
    //! 1/concentration at GP
    std::vector<double> conintinv_;
    //! gradient of electric potential at integration point
    LINALG::Matrix<NSD,1> gradpot_;
    //! velocity at GP
    //LINALG::Matrix<NSD,1> velint_;
    //! convective velocity at GP
    LINALG::Matrix<NSD,1> convelint_;
    //! convective part in convective form: u_x*N,x+ u_y*N,y
    LINALG::Matrix<NEN,1> conv_;
    //! convective contribution
    std::vector<double> convphi_;
    // subgrind velocity
    LINALG::Matrix<NEN,1> sgconv_;


    // TODO: diffcond variables
    //! current density at Gauss point
    LINALG::Matrix<NSD,1> curint_;
    //! divergence of current density at GP
    double curdiv_;

};

}
}

#endif

