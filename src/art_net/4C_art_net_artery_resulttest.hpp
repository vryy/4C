// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP
#define FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Arteries
{
  // forward declaration
  class ArtNetExplicitTimeInt;
  class ArtNetImplStationary;

  /*!
    \brief artnet specific result test class

    \author Mahmoud Ismail
    \date 11/11
  */
  class ArteryResultTest : public Core::Utils::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetExplicitTimeInt& art_net);

    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetImplStationary& art_net);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "phi"
     */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

    /// our version of element value tests
    void test_element(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// std::shared_ptr to scalar transport discretization
    std::shared_ptr<Core::FE::Discretization> dis_;
    /// std::shared_ptr to solution vector
    std::shared_ptr<const Core::LinAlg::Vector<double>> mysol_;
    /// std::shared_ptr to element volumetric flow
    std::shared_ptr<const Core::LinAlg::Vector<double>> myelevolflow_;
    /// std::shared_ptr to element radius
    std::shared_ptr<const Core::LinAlg::Vector<double>> myeleradius_;
  };

}  // namespace Arteries

FOUR_C_NAMESPACE_CLOSE

#endif
