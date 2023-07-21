/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_BOUNDARY_CALC_STD_H
#define BACI_SCATRA_ELE_BOUNDARY_CALC_STD_H

#include "baci_scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype,
        int probdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim + 1>
    class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcStd<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// Empty destructor
      virtual ~ScaTraEleBoundaryCalcStd() { return; };

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcStd(
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraEleBoundaryCalcStd
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
