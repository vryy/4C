/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_std.H

\brief Factory of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_STD_H
#define SCATRA_ELE_CALC_STD_H

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcStd : public ScaTraEleCalc<distype>
{
    /// private constructor, since we are a Singleton.
    ScaTraEleCalcStd(const int numdofpernode,const int numscal);

    typedef ScaTraEleCalc<distype> my;

public:
    virtual ~ScaTraEleCalcStd() {}

    /// Singleton access method
    static ScaTraEleCalcStd<distype> * Instance(
      const int numdofpernode,
      const int numscal,
      bool create=true );

    /// called upon destruction
    virtual void Done();

private:

};
}
}

#endif

