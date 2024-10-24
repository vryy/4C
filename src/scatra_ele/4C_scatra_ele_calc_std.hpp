// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_STD_HPP
#define FOUR_C_SCATRA_ELE_CALC_STD_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // class implementation
    template <Core::FE::CellType distype, int probdim>
    class ScaTraEleCalcStd : public ScaTraEleCalc<distype, probdim>
    {
     public:
      //! abbreviation
      typedef ScaTraEleCalc<distype, probdim> my;


      //! singleton access method
      static ScaTraEleCalcStd<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcStd(const int numdofpernode, const int numscal, const std::string& disname);
    };
  }  // namespace Elements
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
