// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_demangle.hpp"

#include <Teuchos_RCPDecl.hpp>

namespace DemangleTest
{
  struct TestStruct
  {
  };
}  // namespace DemangleTest

namespace
{

  using namespace FourC;

  TEST(DemangleTest, Struct)
  {
    EXPECT_EQ(Core::Utils::try_demangle(typeid(DemangleTest::TestStruct).name()),
        "DemangleTest::TestStruct");
  }

  TEST(DemangleTest, RCP)
  {
    EXPECT_EQ(Core::Utils::try_demangle(typeid(Teuchos::RCP<DemangleTest::TestStruct>).name()),
        "Teuchos::RCP<DemangleTest::TestStruct>");
  }
}  // namespace
