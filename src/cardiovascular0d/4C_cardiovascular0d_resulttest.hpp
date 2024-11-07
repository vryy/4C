// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CARDIOVASCULAR0D_RESULTTEST_HPP
#define FOUR_C_CARDIOVASCULAR0D_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Utils
{
  class Cardiovascular0D;
  class Cardiovascular0DManager;
}  // namespace Utils

namespace Core::IO
{
  class DiscretizationWriter;
}

/*!
  \brief Structure specific result test class
*/
class Cardiovascular0DResultTest : public Core::Utils::ResultTest
{
 public:
  Cardiovascular0DResultTest(Utils::Cardiovascular0DManager& cardvasc0dman,
      std::shared_ptr<Core::FE::Discretization> discr);

  void test_special(
      const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;



 private:
  std::shared_ptr<Core::FE::Discretization> actdisc_;  ///< standard discretization

  const std::shared_ptr<Core::LinAlg::Vector<double>> cardvasc0d_dof_;

  const bool havecardio_4elementwindkessel_;
  const bool havecardio_arterialproxdist_;
  const bool havecardio_syspulcirculation_;
  const bool havecardiorespir_syspulperiphcirculation_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
