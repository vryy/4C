/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.H

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_ADVANCED_REACTION_H_
#define SCATRA_ELE_CALC_ADVANCED_REACTION_H_

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcAdvReac : public ScaTraEleCalc<distype>
{
    /// private constructor, since we are a Singleton.
  ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcAdvReac() {}

  /// Singleton access method
  static ScaTraEleCalcAdvReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

protected:

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

  //! material BioFilm
  void MatBioFilm(
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

  //! calculation of reactive element matrix
  virtual void CalcMatReact(
    Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
    const int                          k,           //!< index of current scalar
    const double                       timefacfac,  //!< domain-integration factor times time-integration factor
    const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
    const double                       taufac,      //!< domain-integration factor times tau
    const double                       densnp,      //!< density at time_(n+1)
    const double                       phinp,       //!< scalar at time_(n+1)
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,  //!< reaction manager
    const LINALG::Matrix<my::nen_,1>&      conv,        //!< convective operator
    const LINALG::Matrix<my::nen_,1>&      sgconv,      //!< subgrid-scale convective operator
    const LINALG::Matrix<my::nen_,1>&      diff         //!< laplace term
  );

};
}
}


#endif /* SCATRA_ELE_CALC_ADVANCED_REACTION_H_ */
