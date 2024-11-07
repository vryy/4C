// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include <gmock/gmock.h>

#include "4C_io_input_file_utils.hpp"

#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

namespace
{
  using namespace FourC;

  TEST(YamlMetadata, Dump)
  {
    Teuchos::ParameterList pl;
    auto& pl_A = pl.sublist("A");
    Core::Utils::int_parameter("a", 1, "The parameter 'a'.", &pl_A);
    pl_A.set("b", true);
    auto& pl_AB = pl_A.sublist("B");
    Core::Utils::string_parameter(
        "c", "string1", "The parameter 'c'.", &pl_AB, {"string1", "string2"});
    pl_AB.set("d", 2.0);
    pl_AB.set("e", false);
    pl_AB.set("f", "string2");

    std::stringstream out;
    Core::IO::InputFileUtils::print_metadata_yaml(out, pl);

    const std::string expected_output =
        R"(parameters:
  A:
    a:
      type: int
      default: 1
      description: The parameter 'a'.
    b:
      type: bool
      default: 1
  A/B:
    c:
      type: string
      default: string1
      description: The parameter 'c'.
      valid options:
        - string1
        - string2
    d:
      type: double
      default: 2
    e:
      type: bool
      default: 0
    f:
      type: string
      default: string2)";

    EXPECT_THAT(out.str(), ::testing::HasSubstr(expected_output));
  }
}  // namespace