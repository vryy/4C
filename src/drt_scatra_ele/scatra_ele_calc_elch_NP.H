/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_NP.H

\brief evalution of ScaTra elements for ion-transport equation

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

//#include "scatra_ele_calc.H"
#include "scatra_ele_calc_elch.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;
//  class ScaTraEleDiffManagerElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElchNP : public ScaTraEleCalcElch<distype>
{
public:

  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElchNP(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcElch<distype> myelch;

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

  virtual void CalMatAndRhs(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&               emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&               erhs,         //!< element rhs to calculate+
    const int                               k,            //!< index of current scalar
    const double                            fac,          //!< domain-integration factor
    const double                            timefacfac,   //!< domain-integration factor times time-integration factor
    const double                            rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>& dme,          //!< diffusion manager
    double&                                 rhsint,       //!< rhs at Gauss point
    const double                            hist          //!< history
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

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/


  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  virtual void PrepMatAndRhsInitialTimeDerivative(
    Epetra_SerialDenseMatrix&  elemat1_epetra,
    Epetra_SerialDenseVector&  elevec1_epetra
  );

  virtual void GetConductivity(
    const enum INPAR::ELCH::ElchType  elchtype,
    double&                           sigma_all,
    Epetra_SerialDenseVector&         sigma
  );

  virtual void CalMatAndRhsElectricPotentialField(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs,
    const double                      fac,
    Teuchos::RCP<ScaTraEleDiffManagerElch>& dme          //!< diffusion manager
  );

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/



  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

};

}
}

#endif

