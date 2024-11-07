// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_RESULTTEST_ELCH_HPP
#define FOUR_C_SCATRA_RESULTTEST_ELCH_HPP

#include "4C_config.hpp"

#include "4C_scatra_resulttest.hpp"
#include "4C_scatra_timint_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{

  // class implementation
  class ElchResultTest : public ScaTraResultTest
  {
   public:
    //! constructor
    ElchResultTest(std::shared_ptr<ScaTraTimIntElch> elchtimint);

   private:
    //! return pointer to elch time integrator after cast
    std::shared_ptr<const ScaTra::ScaTraTimIntElch> elch_tim_int() const
    {
      return std::dynamic_pointer_cast<const ScaTra::ScaTraTimIntElch>(scatratimint_);
    };

    //! get special result to be tested
    double result_special(const std::string quantity  //! name of quantity to be tested
    ) const override;
  };  // class ElchResultTest : public ScaTraResultTest
}  // namespace ScaTra
FOUR_C_NAMESPACE_CLOSE

#endif
