// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STD_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STD_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::Elements::ScaTraEleBoundaryCalc<distype, probdim> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcStd<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcStd(
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraEleBoundaryCalcStd
  }  // namespace Elements
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
