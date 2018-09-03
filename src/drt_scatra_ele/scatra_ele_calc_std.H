/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_std.H

\brief evaluation of scalar transport elements for standard scalar transport problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
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
    // class implementation
    template <DRT::Element::DiscretizationType distype, int probdim>
    class ScaTraEleCalcStd : public ScaTraEleCalc<distype, probdim>
    {
     public:
      //! abbreviation
      typedef ScaTraEleCalc<distype, probdim> my;

      //! destructor
      virtual ~ScaTraEleCalcStd() { return; };

      //! singleton access method
      static ScaTraEleCalcStd<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname, const ScaTraEleCalcStd* delete_me = NULL);

      //! called upon destruction
      virtual void Done();

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcStd(const int numdofpernode, const int numscal, const std::string& disname);
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
