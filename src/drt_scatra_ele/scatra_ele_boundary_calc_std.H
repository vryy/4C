/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.H

\brief evaluation of scatra boundary terms at integration points

\level 1

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_STD_H
#define SCATRA_ELE_BOUNDARY_CALC_STD_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcStd<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleBoundaryCalcStd* delete_me = NULL);

      /// Empty destructor
      virtual ~ScaTraEleBoundaryCalcStd() { return; };

      /// called upon destruction
      virtual void Done();

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcStd(
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraEleBoundaryCalcStd
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
