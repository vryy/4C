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

namespace DRT
{
namespace ELEMENTS
{

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

private:

/*========================================================================*/
//! @name routines for additional element evaluations (called from EvaluateService)
/*========================================================================*/

};
}
}

#endif

