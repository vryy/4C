/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems


\level 1
*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_SCATRA_ELE_CALC_STD_H
#define BACI_SCATRA_ELE_CALC_STD_H

#include "baci_scatra_ele_calc.H"

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


      //! singleton access method
      static ScaTraEleCalcStd<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcStd(const int numdofpernode, const int numscal, const std::string& disname);
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
