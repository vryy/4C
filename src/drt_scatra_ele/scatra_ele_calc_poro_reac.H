/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro_reac.H

 \brief

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_REAC_H_
#define SCATRA_ELE_CALC_PORO_REAC_H_

#include "scatra_ele_calc.H"
#include "scatra_ele_calc_poro.H"
#include "scatra_ele_calc_advanced_reaction.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcPoroReac : public ScaTraEleCalcPoro<distype>, public ScaTraEleCalcAdvReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcPoroReac() {}

  /// Singleton access method
  static ScaTraEleCalcPoroReac<distype> * Instance(
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
