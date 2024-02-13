/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_BOUNDARY_CALC_STD_HPP
#define BACI_SCATRA_ELE_BOUNDARY_CALC_STD_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_boundary_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcStd<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcStd(
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraEleBoundaryCalcStd
  }     // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
