// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_RESULTTEST_HDG_HPP
#define FOUR_C_SCATRA_RESULTTEST_HDG_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_scatra_resulttest.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  // forward declaration
  class TimIntHDG;

  // class implementation
  class HDGResultTest : public ScaTraResultTest
  {
   public:
    //! constructor
    HDGResultTest(std::shared_ptr<ScaTraTimIntImpl> timint);

   private:
    //! get nodal result to be tested
    double result_node(const std::string quantity,  //! name of quantity to be tested
        Core::Nodes::Node* node                     //! node carrying the result to be tested
    ) const override;

    //! time integrator
    std::shared_ptr<const TimIntHDG> scatratiminthdg_;

    std::shared_ptr<Core::LinAlg::SerialDenseVector> errors_;

  };  // class HDGResultTest : public ScaTraResultTest
}  // namespace ScaTra
FOUR_C_NAMESPACE_CLOSE

#endif
