/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_H_
#define SCATRA_ELE_CALC_PORO_H_

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcPoro : public ScaTraEleCalc<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcPoro(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalcPoro<distype> my;

public:
  virtual ~ScaTraEleCalcPoro() {}

  /// Singleton access method
  static ScaTraEleCalcPoro<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

protected:

  //! material ScaTra
  virtual void MatScaTra(
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
//! @name routines for additional element evaluations (called from EvaluateService)
/*========================================================================*/

};
}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */
