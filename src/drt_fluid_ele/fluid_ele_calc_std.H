/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_std.H

\brief standard routines for calculation of fluid element

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_STD_H
#define FLUID_ELE_CALC_STD_H

#include "fluid_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcStd : public FluidEleCalc<distype>
{
    /// private constructor, since we are a Singleton.
    FluidEleCalcStd(int num=0);

    typedef FluidEleCalc<distype> my;
    typedef FluidEleInterface MY;

public:
    virtual ~FluidEleCalcStd() {}

    /// Singleton access method
    static FluidEleCalcStd<distype>* Instance( bool create=true, int num=0 );

    /// called upon destruction
    virtual void Done();

    int num_;

private:

};
}
}

#endif

