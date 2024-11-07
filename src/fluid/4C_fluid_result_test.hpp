// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_RESULT_TEST_HPP
#define FOUR_C_FLUID_RESULT_TEST_HPP


#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FLD
{
  // forward declarations
  class FluidImplicitTimeInt;
  class FluidGenAlphaIntegration;


  /*!
  \brief Fluid specific result test class

  \author u.kue
  */
  class FluidResultTest : public Core::Utils::ResultTest
  {
   public:
    /*!
    \brief not documented yet
    */
    FluidResultTest(FluidImplicitTimeInt& fluid);

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz", "pressure",
    "tractionx", "tractiony", "tractionz", "errvel", "errpre" and "divu".
    With the obvious meaning.
    */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// pointer to fluid discretization
    std::shared_ptr<Core::FE::Discretization> fluiddis_;
    /// pointer to unknown vector with nodal values
    std::shared_ptr<Core::LinAlg::Vector<double>> mysol_;
    /// pointer to traction vector with values
    std::shared_ptr<Core::LinAlg::Vector<double>> mytraction_;
    /// pointer to traction vector with values
    std::shared_ptr<Core::LinAlg::Vector<double>> mywss_;
    /// pointer to error evaluation
    std::shared_ptr<std::vector<double>> myerror_;
    /// pointer to div u evaluation
    std::shared_ptr<double> mydivu_;
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
