/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STD_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STD_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;

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
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
