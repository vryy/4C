// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ELEMAG_RESULTTEST_HPP
#define FOUR_C_ELEMAG_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace EleMag
{
  class ElemagTimeInt;

  class ElemagResultTest : public Core::Utils::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ElemagResultTest(ElemagTimeInt& elemagalgo);


    /// nodal value tests
    /*!
      Possible position flags is only "pressure"
     */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to elemagstical discretization
    Teuchos::RCP<Core::FE::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Core::LinAlg::Vector<double>> mysol_;
    /// Error vector
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> error_;

  };  // class ElemagResultTest

}  // namespace EleMag

FOUR_C_NAMESPACE_CLOSE

#endif