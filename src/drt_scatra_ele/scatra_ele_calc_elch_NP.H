/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_NP.H

\brief evaluation of ScaTra elements for Nernst-Planck ion-transport equations

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_NP_H
#define SCATRA_ELE_CALC_ELCH_NP_H

#include "scatra_ele_calc_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElchNP : public ScaTraEleCalcElch<distype>
{
private:

  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElchNP(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcElch<distype> myelch;

public:

  virtual ~ScaTraEleCalcElchNP() {}

  /// Singleton access method
  static ScaTraEleCalcElchNP<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();


private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! calculate contributions to matrix and rhs (inside of loop over all scalars)
  virtual void CalcMatAndRhs(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >&   vm,
    Epetra_SerialDenseMatrix&                                                  emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&                                                  erhs,         //!< element rhs to calculate+
    const int                                                                  k,            //!< index of current scalar
    const double                                                               fac,          //!< domain-integration factor
    const double                                                               timefacfac,   //!< domain-integration factor times time-integration factor
    const double                                                               rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    const double                                                               taufac,       //!< tau times domain-integration factor
    const double                                                               timetaufac,   //!< domain-integration factor times tau times time-integration factor
    const double                                                               rhstaufac,    //!< time-integration factor for rhs times tau times domain-integration factor
    LINALG::Matrix<my::nen_,1>&                                                tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                    dme,          //!< diffusion manager
    double&                                                                    rhsint,       //!< rhs at Gauss point
    const double                                                               hist          //!< history
  );

  //! Correction for additional flux terms / currents across Dirichlet boundaries
  virtual void CorrectionForFluxAcrossDC(
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
  //! @name material and related functions
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
    const int                               k,        //!< index of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  /*========================================================================*/
  //! @name stabilization and related functions
  /*========================================================================*/

  //! Calculate derivative of tau w.r.t. electric potential according to Taylor, Hughes and Zarins
  virtual void CalcTauDerPotTaylorHughesZarins(
    LINALG::Matrix<my::nen_,1>&         tauderpot,       //!< derivatives of stabilization parameter w.r.t. electric potential
    double&                             tau,             //!< stabilization parameter
    const double                        densnp,          //!< density at t_(n+1)
    const double                        frt,             //!< F/(RT)
    const double                        diffusvalence,   //!< diffusion coefficient times valence
    const LINALG::Matrix<my::nsd_,1>&   veleff           //!< effective convective velocity (fluid velocity plus migration velocity if applicable)
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/

  //! CalcRes: Residual of Nernst-Planck equation in strong form
  double CalcRes(
    const int                                 k,            //!< index of current scalar
    const double                              conint,       //!< concentration at GP
    const double                              hist,         //!< history value at GP
    const double                              convphi,      //!< convective term (without convective part of migration term)
    const double                              frt,          //!< F/(RT)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         migconv,      //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
    const double                              rhsint        //!< rhs of Nernst-Planck equation (not of Newton-Raphson scheme) at Gauss point
  );

  //! CalcMat: SUPG Stabilization of convective term due to fluid flow and migration
  void CalcMatConvStab(
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to calculate
    const int                                 k,            //!< index of current scalar
    const double                              timefacfac,   //!< domain-integration factor times time-integration factor
    const double                              taufac,       //!< stabilization parameter tau times domain-integration factor
    const double                              timetaufac,   //!< domain-integration factor times tau times time-integration factor
    LINALG::Matrix<my::nen_,1>&               tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
    const double                              frt,          //!< F/(RT)
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         conv,         //!< convection operator: u_x*N,x + u_y*N,y + u_z*N,z
    const LINALG::Matrix<my::nen_,1>&         migconv,      //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
    const double                              conint,       //!< concentration at GP
    const LINALG::Matrix<my::nsd_,1>&         gradphi,      //!< gradient of concentration at GP
    const double                              residual      //!< residual of Nernst-Planck equation in strong form
  );

  //! CalcMat: Migration term
  void CalcMatMigr(
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to be filled
    const int                                 k,            //!< index of current scalar
    const double                              timefacfac,   //!< domain-integration factor times time-integration factor
    const double                              frt,          //!< F/RT
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         migconv,      //!< migration operator
    const double                              conint        //!< concentration at GP
  );

  //! CalcMat: Electroneutrality condition in PDE form
  void CalcMatPotEquENCPDE(
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to be filled
    const int                                 k,            //!< index of current scalar
    const double                              timefacfac,   //!< domain-integration factor times time-integration factor
    const double                              frt,          //!< F/RT
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         migconv,      //!< migration operator
    const double                              conint        //!< concentration at GP
  );

  //! CalcMat: Electroneutrality condition in PDE form with Nernst-Planck equation for species m eliminated
  void CalcMatPotEquENCPDEElim(
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to be filled
    const int                                 k,            //!< index of current scalar
    const double                              timefacfac,   //!< domain-integration factor times time-integration factor
    const double                              frt,          //!< F/RT
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         migconv,      //!< migration operator
    const double                              conint        //!< concentration at GP
  );

  //! CalcMat: Poisson equation for electric potential
  void CalcMatPotEquPoisson(
    Epetra_SerialDenseMatrix&                 emat,      //!< element matrix to be filled
    const int                                 k,         //!< index of current scalar
    const double                              fac,       //!< domain-integration factor
    const double                              epsilon,   //!< dielectric constant
    const double                              faraday,   //!< Faraday constant
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme        //!< diffusion manager
  );

  //! CalcMat: Laplace equation for electric potential
  void CalcMatPotEquLaplace(
    Epetra_SerialDenseMatrix&                 emat,      //!< element matrix to be filled
    const int                                 k,         //!< index of current scalar
    const double                              fac        //!< domain-integration factor
  );

  //! CalcRhs: Additional contributions from conservative formulation of Nernst-Planck equations
  void CalcRhsConvAddCons(
    Epetra_SerialDenseVector&                 erhs,     //!< element vector to be filled
    const int                                 k,        //!< index of current scalar
    const double                              rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                              conint,   //!< concentration at GP
    const double                              vdiv      //!< velocity divergence
  );

  //! CalcRhs: SUPG Stabilization of convective term due to fluid flow and migration
  void CalcRhsConvStab(
    Epetra_SerialDenseVector&                 erhs,        //!< element vector to be filled
    const int                                 k,           //!< index of current scalar
    const double                              rhstaufac,   //!< time-integration factor for rhs times tau times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,         //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         conv,        //!< convection operator: u_x*N,x + u_y*N,y + u_z*N,z
    const LINALG::Matrix<my::nen_,1>&         migconv,     //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
    const double                              residual     //!< residual of Nernst-Planck equation in strong form
  );

  //! CalcRhs: Migration term
  void CalcRhsMigr(
    Epetra_SerialDenseVector&                 erhs,      //!< element vector to be filled
    const int                                 k,         //!< index of current scalar
    const double                              rhsfac,    //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,       //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&         migconv,   //!< migration operator
    const double                              conint     //!< concentration at GP
  );

  //! CalcRhs: Electroneutrality condition in PDE form
  void CalcRhsPotEquENCPDE(
    Epetra_SerialDenseVector&                         erhs,      //!< element vector to be filled
    const int                                         k,         //!< index of current scalar
    const double                                      fac,       //!< domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&           dme,       //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&                 migconv,   //!< migration operator
    const double                                      conint,    //!< concentration at GP
    const LINALG::Matrix<my::nsd_,1>&                 gradphi    //!< gradient of concentration at GP
  );

  //! CalcRhs: Electroneutrality condition in PDE form with Nernst-Planck equation for species m eliminated
  void CalcRhsPotEquENCPDEElim(
    Epetra_SerialDenseVector&                         erhs,      //!< element vector to be filled
    const int                                         k,         //!< index of current scalar
    const double                                      fac,       //!< domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&           dme,       //!< diffusion manager
    const LINALG::Matrix<my::nen_,1>&                 migconv,   //!< migration operator
    const double                                      conint,    //!< concentration at GP
    const LINALG::Matrix<my::nsd_,1>&                 gradphi    //!< gradient of concentration at GP
  );

  //! CalcRhs: Poisson equation for electric potential
  void CalcRhsPotEquPoisson(
    Epetra_SerialDenseVector&                 erhs,      //!< element vector to be filled
    const int                                 k,         //!< index of current scalar
    const double                              fac,       //!< domain-integration factor
    const double                              epsilon,   //!< dielectric constant
    const double                              faraday,   //!< Faraday constant
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,       //!< diffusion manager
    const double                              conint,    //!< concentration at GP
    const LINALG::Matrix<my::nsd_,1>&         gradpot    //!< gradient of potential at GP
  );

  //! CalcRhs: Laplace equation for electric potential
  void CalcRhsPotEquLaplace(
    Epetra_SerialDenseVector&                 erhs,      //!< element vector to be filled
    const int                                 k,         //!< index of current scalar
    const double                              fac,       //!< domain-integration factor
    const LINALG::Matrix<my::nsd_,1>&         gradpot    //!< gradient of potential at GP
  );

  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  //! Solving for intial time derivative: finish preparation of matrix and rhs
  virtual void PrepMatAndRhsInitialTimeDerivative(
    Epetra_SerialDenseMatrix&  elemat1_epetra,  //!< element matrix
    Epetra_SerialDenseVector&  elevec1_epetra   //!< element rhs
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
    const double                    fac,         //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme                                 //!< diffusion manager
  );

  //!  calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
  //!  This would be a nice feature, e.g., for equpot::equpot_enc_pde
  virtual void CalculateCurrent(
    LINALG::Matrix<my::nsd_,1>&     q,      //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const double                    fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                  dme  //!< diffusion manager
  ){return;};

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/



  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

};


/// ScaTraEleInternalVariableManagerElchNP implementation
/*!
  This class keeps all internal (Nernst-Planck-specific) variables needed for the evaluation of an element.
*/
template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerElchNP : public ScaTraEleInternalVariableManagerElch<NSD,NEN>
{
  public:

    typedef ScaTraEleInternalVariableManagerElch<NSD,NEN> vm;

    ScaTraEleInternalVariableManagerElchNP(int numscal, int nsd, DRT::ELEMENTS::ScaTraEleParameterElch* elchpara):
      ScaTraEleInternalVariableManagerElch<NSD,NEN>(numscal,nsd,elchpara),
      // constant internal variables
      // empty

      // internal variables evaluated at the Gauss point
      migconv_(true)
    {return;}

    virtual ~ScaTraEleInternalVariableManagerElchNP(){};

    // compute and set internal variables only used by the Nernst-Planck formulation
    void SetInternalVariablesElchNP(
      LINALG::Matrix<NSD,NEN>&                 derxy   //!< global derivatives of shape functions w.r.t x,y,z
    )
    {
      //! convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi}*\grad{N}
      migconv_.MultiplyTN(-vm::frt_,derxy,vm::gradpot_);

      return;
    };

  /*========================================================================*/
  //! @name return constant internal variables
  /*========================================================================*/

    // empty

  /*========================================================================*/
  //! @name return methods for internal variables
  /*========================================================================*/

    //! return convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi} * \grad{N}
    const LINALG::Matrix<NEN,1> MigConv(){return migconv_;};

  protected:

  /*========================================================================*/
  //! @name constant internal variables
  /*========================================================================*/


  /*========================================================================*/
  //! @name internal variables evaluated at element center or Gauss point
  /*========================================================================*/

    //! convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi} * \grad{N}
    LINALG::Matrix<NEN,1> migconv_;

}; // class ScaTraEleInternalVariableManagerElchNP

} // namespace ELEMENTS

} // namespace DRT

#endif
