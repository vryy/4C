/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems


\level 1
*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_SCATRA_ELE_CALC_STD_HPP
#define BACI_SCATRA_ELE_CALC_STD_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <CORE::FE::CellType distype, int probdim>
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
BACI_NAMESPACE_CLOSE

#endif
