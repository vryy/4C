/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_diffcond.H

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_DIFFCOND_H
#define SCATRA_ELE_CALC_ELCH_DIFFCOND_H

#include "scatra_ele_calc_elch.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;
  class ScaTraEleDiffManagerElchDiffCond;
  template<int NSD, int NEN>
  class ScaTraEleInternalVariableManagerElchDiffCond;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElchDiffCond : public ScaTraEleCalcElch<distype>
{
private:

  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElchDiffCond(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcElch<distype> myelch;

public:

  virtual ~ScaTraEleCalcElchDiffCond() {}

  /// Singleton access method
  static ScaTraEleCalcElchDiffCond<distype> * Instance(
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


private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! calculate contributions to matrix and rhs (inside of loop over all scalars)
  virtual void CalcMatAndRhs(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&               emat,        //!< element matrix to calculate
    Epetra_SerialDenseVector&               erhs,        //!< element rhs to calculate+
    const int                               k,           //!< index of current scalar
    const double                            fac,         //!< domain-integration factor
    const double                            timefacfac,  //!< domain-integration factor times time-integration factor
    const double                            rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    const double                            taufac,       //!< tau times domain-integration factor
    const double                            timetaufac,   //!< domain-integration factor times tau times time-integration factor
    const double                            rhstaufac,    //!< time-integration factor for rhs times tau times domain-integration factor
    LINALG::Matrix<my::nen_,1>&             tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
    Teuchos::RCP<ScaTraEleDiffManagerElch>& dme,         //!< diffusion manager
    double&                                 rhsint,      //!< rhs at Gauss point
    const double                            hist         //!< history
  );

  //! calculate contributions to matrix and rhs (outside of loop over all scalars)
  virtual void CalcMatAndRhsOutsideScalarLoop(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&                emat,        //!< element matrix to calculate
    Epetra_SerialDenseVector&                erhs,        //!< element rhs to calculate
    const double                             fac,         //!< domain-integration factor
    const double                             timefacfac,  //!< domain-integration factor times time-integration factor
    const double                             rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme          //!< diffusion manager
  );

  //! Correction for additional flux terms / currents across Dirichlet boundaries
  void CorrectionForFluxAcrossDC(
      DRT::Discretization&        discretization,   //!< discretization
      const std::vector<int>&     lm,               //!< location vector
      Epetra_SerialDenseMatrix&   emat,             //!< element matrix to calculate
      Epetra_SerialDenseVector&   erhs              //!< element rhs to calculate
  );

  //! Formulation-specific internal variables are set in the variable manager
  virtual void SetFormulationSpecificInternalVariables(
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                   dme,  //!< rcp to diffusion manager Elch
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >&  vm    //!< rcp to variable manager Elch
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
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager, //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,  //!< reaction manager
    double&             visc,                       //!< fluid viscosity
    const int           iquad=-1                    //!< id of current gauss point (default = -1)
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
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad=-1      //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatPhase(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               iphase,   //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>  dmedc,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,     //!< reaction manager
    double&                                 visc,           //!< fluid viscosity
    const int                               iquad=-1        //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatNewman(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>  dmedc,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,     //!< reaction manager
    double&                                 visc,           //!< fluid viscosity
    const int                               iquad=-1        //!< id of current gauss point (default = -1)
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/

  //! CalcMat: Linearization of diffusion coefficient in diffusion term
  void CalcMatDiffCoeffLin(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int                               k,              //!< index of current scalar
    const double                            timefacfac,     //!< domain-integration factor times time-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi        //!< gradient of concentration at GP
  );

  //! CalcMat: Conduction term with inserted current - ohmic overpotential
  void CalcMatCondOhm(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int                               k,              //!< index of current scalar
    const double                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double                            invfval,        //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot        //!< gradient of potential at GP
  );

  //! CalcMat: Conduction term with inserted current - conc. overpotential
  void CalcMatCondConc(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int                               k,              //!< index of current scalar
    const double                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double                            rtffcval,       //!< RT/F^2/Newman_const_c/z_k
    const double                            newman_const_a, //!< Newman constant a
    const double                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,       //!< gradient of concentration at GP
    const std::vector<double>&               conintinv      //!< inverted concentration at GP
  );

  //! CalcMat: Conduction term without inserted current
  void CalcMatCond(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int                               k,            //!< index of current scalar
    const double                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        curint         //!< current at GP
  );

  //! CalcMat: Additional diffusion term without inserted current
  void CalcMatCondDiff(
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
    const int                               k,            //!< index of current scalar
    const double                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double                            invfval,      //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,  //!< diffusion manager
    const std::vector<LINALG::Matrix<my::nsd_,1> >&  gradphi //!< gradient of concentration at GP
  );

  //! CalcMat: Potential equation div i inserted current - ohmic overpotential
  void CalcMatPotEquDiviOhm(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const double                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double                            invf,           //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot        //!< gradient of potenial at GP
  );

  //! Potential equation div i inserted current - concentration overpotential
  void CalcMatPotEquDiviConc(
    Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
    const int                               k,              //!< index of current scalar
    const double                            timefacfac,     //!< domain-integration factor times time-integration factor
    const double                            rtffc,          //!< RT/(F^2 Newman_const_c)
    const double                            rtf,            //!< RT/F
    const double                            invf,           //!< 1/F
    const double                            newman_const_a, //!< Newman constant a
    const double                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,       //!< gradient of concentration at GP
    const double                            conintinv       //!< inverted concentration at GP
  );

  //! Potential equation div i without inserted current
  void CalcMatPotEquDivi(
    Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
    const double                timefacfac, //!< domain-integration factor times time-integration factor
    const double                invf        //!< 1/F
  );

  //! CalcMat: Current equation current
  void CalcMatCurEquCur(
    Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
    const double                timefacfac, //!< domain-integration factor times time-integration factor
    const double                invf        //!< 1/F
  );

  //! CalcMat: Current equation ohmic overpotential
  void CalcMatCurEquOhm(
      Epetra_SerialDenseMatrix&               emat,       //!< element matrix to be filled
      const double                            timefacfac, //!< domain-integration factor times time-integration factor
      const double                            invf,       //!< 1/F
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,   //!< diffusion manager
      const LINALG::Matrix<my::nsd_,1>&               gradpot  //!< gradient of potenial at GP
  );

  //! CalcMat: Current equation concentration overpotential
  void CalcMatCurEquConc(
    Epetra_SerialDenseMatrix&               emat,             //!< element matrix to be filled
    const double                            timefacfac,       //!< domain-integration factor times time-integration factor
    const double                            rtf,              //!< RT/F
    const double                            rtffc,            //!< RT/(F^2 Newman_const_c)
    const std::vector <double>&             invfval,          //!< 1/(F z_k)
    const double                            newman_const_a,   //!< Newman constant a
    const double                            newman_const_b,   //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,   //!< diffusion manager
    const std::vector< LINALG::Matrix<my::nsd_,1> >& gradphi, //!< gradient of concentration at GP
    const std::vector<double>&              conintinv         //!< inverted concentration at GP
  );

  //! CalcRhs: Conduction term with inserted current - ohmic overpotential
  void CalcRhsCondOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int                               k,        //!< index of current scalar
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invfval,  //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradpot        //!< gradient of potenial at GP
  );

  //! CalcRhs: Conduction term with inserted current - conc. overpotential
  void CalcRhsCondConc(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int                               k,        //!< index of current scalar
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            rtffcval, //!< RT/(F^2 Newman_const_c z_k)
    const double                            newman_const_a, //!< Newman constant a
    const double                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc, //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&        gradphi,       //!< gradient of concentration at GP
    const std::vector<double>&               conintinv      //!< inverted concentration at GP
  );

  //! CalcRhs: Conduction term without inserted current
  void CalcRhsCond(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int                               k,        //!< index of current scalar
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invfval,  //!< 1/(F z_k)
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,  //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&                curint  //!< current at GP
  );

  //! CalcRhs: Additional diffusion term without inserted current
  void CalcRhsCondDiff(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int                               k,        //!< index of current scalar
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,  //!< diffusion manager
    const std::vector< LINALG::Matrix<my::nsd_,1> >& gradphi //!< gradient of concentration at GP
  );

  //! CalcRhs: Potential equation div i inserted current - ohmic overpotential
  void CalcRhsPotEquDiviOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,    //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&               gradpot   //!< gradient of potenial at GP
  );

  //! CalcRhs: Potential equation div i inserted current - conc. overpotential
  void CalcRhsPotEquDiviConc(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const int                               k,        //!< index of current scalar
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            rtf,      //!< RT/F
    const std::vector<double>&              invfval,  //!< 1/(F z_k)
    const double                            rtffc,    //!< RT/(F^2 Newman_const_c)
    const double                            newman_const_a, //!< Newman constant a
    const double                            newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,  //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&      gradphi,   //!< gradient of concentration at GP
    const double                           conintinv  //!< inverted concentration at GP
  );

  //! CalcRhs: Potential equation divi without inserted current
  void CalcRhsPotEquDivi(
    Epetra_SerialDenseVector&           erhs,     //!< element vector to be filled
    const double                        rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                        invf,     //!< 1/F
    const LINALG::Matrix<my::nsd_,1> &  curint    //!< current at GP
  );

  //! CalcRhs: Current equation - current
  void CalcRhsCurEquCur(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,    //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1> &              curint    //!< current at GP
  );

  //! CalcRhs: Current equation - ohmic overpotential
  void CalcRhsCurEquOhm(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,  //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&              gradpot   //!< gradient of potenial at GP
  );

  //! Current equation - concentration overpotential
  void CalcRhsCurEquConc(
    Epetra_SerialDenseVector&                   erhs,     //!< element vector to be filled
    const double                                rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                                rtf,      //!< RT/F
    const std::vector<double>&                  invfval,  //!< 1/(F z_k)
    const double                                rtffc,    //!< RT/(F^2 Newman_const_c)
    const double                                newman_const_a, //!< Newman constant a
    const double                                newman_const_b, //!< Newman constant b
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,     //!< diffusion manager
    const std::vector<LINALG::Matrix<my::nsd_,1> >&  gradphi,   //!< vector of gradient of concentration at GP
    const std::vector<double>&                       conintinv  //!< inverted concentration at GP
   );


  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  //! validity check for all elements with respect to
  //! formulation-specific parameter, degree's of freedom, number of scalars, ...
  virtual void CheckElchElementParameter(
    DRT::ELEMENTS::Transport*  ele  //!< the element we are dealing with
  );

  //! Solving for intial time derivative: finish preparation of matrix and rhs
  virtual void PrepMatAndRhsInitialTimeDerivative(
    Epetra_SerialDenseMatrix&  elemat1_epetra,  //!< element matrix
    Epetra_SerialDenseVector&  elevec1_epetra   //!< element rhs
  );

  //! Correct RHS calculated from CalcRHSLinMass() for the linearized mass term
  virtual void CorrectRHSFromCalcRHSLinMass(
    Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
    const int                     k,           //!< index of current scalar
    const double                  fac,         //!< domain-integration factor
    const double                  densnp,      //!< density at time_(n+1)
    const double                  phinp        //!< scalar at time_(n+1)
  );

  // Get conductivity from material
  virtual void GetConductivity(
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    double&                           sigma_all,  //!< conductivity of electrolyte solution
    Epetra_SerialDenseVector&         sigma       //!< conductivity or a single ion + overall electrolyte solution
  );

  // Solving for electric potential field : Calculate Matrix and Rhs
  virtual void CalcMatAndRhsElectricPotentialField(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    Epetra_SerialDenseMatrix&         emat,       //!< element matrix
    Epetra_SerialDenseVector&         erhs,       //!< element rhs
    const double                      fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>& dme   //!< diffusion manager
  );

  //!  calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
  virtual void CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,      //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme                                 //!< diffusion manager
  );

  //!  calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
  virtual void CalculateCurrent(
    LINALG::Matrix<my::nsd_,1>&     q,      //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const double                    fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                  dme  //!< diffusion manager
  );

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! flag for used element formulation (material based)
  INPAR::ELCH::DiffCondMat diffcondmat_;

  //! flag if current is used as a solution variable
  bool cursolvar_;

  // flag if elch system is closed by ENC or div i
  enum INPAR::ELCH::EquPot equpot_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! nodal current density values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<my::nsd_,my::nen_> ecurnp_;
};


/// ScaTraEleDiffManagerElch implementation
/*!
  This class keeps all Diffusion-Conduction-specific transport parameter needed for the evaluation of an element.
  The ScaTraEleDiffManagerElch is derived from the standard ScaTraEleDiffManager.
  for
*/
class ScaTraEleDiffManagerElchDiffCond : public ScaTraEleDiffManagerElch
{
  public:

    ScaTraEleDiffManagerElchDiffCond(int numscal):
    ScaTraEleDiffManagerElch(numscal),
    derivdiffcoef_(numscal,std::vector<double>(numscal,0.0)),
    transnum_(numscal,0.0),
    derivtransnum_(numscal,std::vector<double>(numscal,0.0)),
    cond_(0.0),
    derivcond_(numscal,0.0),
    thermfac_(1.0),
    derivthermfac_(numscal,0.0),
    // In the moment, we have only one phase but the framework is flexible
    eps_(1,1.0),
    tort_(1,1.0)
    {return;}

    ~ScaTraEleDiffManagerElchDiffCond(){};

    /*========================================================================*/
    //! @name access methods for transport parameter
    /*========================================================================*/

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

    //! Calculate conductivity based on diffusion coefficients
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

    //! Calculate transference numbers based on diffusion coefficients
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

    /*========================================================================*/
    //! @name access methods for geometrical parameter
    /*========================================================================*/

    //! Set transference numbers with respect to single ionic species
    void SetPhasePoro(const double eps,const int phase)
    {
      eps_[phase]= eps;
      return;
    }

    const double GetPhasePoro(const int phase){return eps_[phase];};

    //! Set transference numbers with respect to single ionic species
    void SetPhaseTort(const double tort,const int phase)
    {
      tort_[phase]= tort;
      return;
    }

    // get geometrical parameter: tortuosity
    const double GetPhaseTort(const int phase){return tort_[phase];};

    // get geometrical parameter: porosity*tortuosity
    const double GetPhasePoroTort(const int phase){
      double epstort =eps_[phase]*tort_[phase];
      return epstort;};

    /*========================================================================*/
    //! @name output
    /*========================================================================*/

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

  private:

    /*========================================================================*/
    //! @name transport parameter
    /*========================================================================*/

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


    /*========================================================================*/
    //! @name geometrical parameter of the porous medium
    /*========================================================================*/

    //! porosity
    std::vector<double> eps_;
    //! tortuosity
    std::vector<double> tort_;

};


/// ScaTraEleInternalVariableManagerElchDiffCond implementation
/*!
  This class keeps all internal (diffusion-conduction-specific) variables needed for the evaluation of an element.
*/
template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerElchDiffCond : public ScaTraEleInternalVariableManagerElch<NSD,NEN>
{
  public:

    typedef ScaTraEleInternalVariableManagerElch<NSD,NEN> vm;

    ScaTraEleInternalVariableManagerElchDiffCond(int numscal, int nsd,DRT::ELEMENTS::ScaTraEleParameterElch* elchpara):
      ScaTraEleInternalVariableManagerElch<NSD,NEN>(numscal,nsd,elchpara),
      // constant internal variables
      invf_(1/vm::faraday_),
      rtf_(1/vm::frt_),
      rtffc_(rtf_*invf_/elchpara->NewmanConstC()),
      invval_(numscal),
      invfval_(numscal),
      rtffcval_(numscal),
      //internal variables evaluated at the Gauss point
      curint_(true),
      curdiv_(0.0)
    {return;}

    virtual ~ScaTraEleInternalVariableManagerElchDiffCond(){};

    // compute and set internal variables only used by the Diffusion-Conduction formulation
    void SetInternalVariablesElchDiffCond(
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,  //!< diffusion manager
      LINALG::Matrix<NEN,1>&                          funct,  //!< array for shape functions
      LINALG::Matrix<NSD,NEN>&                        derxy,  //!< global derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<NSD,NEN>&                        ecurnp  //!< nodal current density values at t_(n+1) or t_(n+alpha_F)
    )
    {
      // constant material-depending variables
      for(int k=0; k<vm::numscal_; ++k)
      {
        invval_[k]=1.0/dmedc->GetValence(k);
        invfval_[k]=invf_/dmedc->GetValence(k);
        rtffcval_[k]=rtffc_/dmedc->GetValence(k);
      }

      curint_.Multiply(ecurnp,funct);
      // Calculate divergence of current desity
      // since it is a summation it is important to set curdiv to zero
      curdiv_=0.0;
      for (int idim = 0; idim <vm::nsd_; ++idim)
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

    //! return current density at GP
    const LINALG::Matrix<NSD,1> CurInt(){return curint_;};
    //! return divergence of current density
    const double CurDiv(){return curdiv_;};

  protected:

  /*========================================================================*/
  //! @name constant internal variables
  /*========================================================================*/

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

    //! current density at Gauss point
    LINALG::Matrix<NSD,1> curint_;
    //! divergence of current density at GP
    double curdiv_;
};


}
}

#endif

