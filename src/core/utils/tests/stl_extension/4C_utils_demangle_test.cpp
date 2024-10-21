// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_demangle.hpp"

#include <Teuchos_RCPDecl.hpp>

namespace DEMANGLE_TEST
{
  struct TestStruct
  {
  };
}  // namespace DEMANGLE_TEST

namespace
{

  using namespace FourC;

  TEST(DemangleTest, Struct)
  {
    EXPECT_EQ(Core::Utils::try_demangle(typeid(DEMANGLE_TEST::TestStruct).name()),
        "DEMANGLE_TEST::TestStruct");
  }

  TEST(DemangleTest, RCP)
  {
    EXPECT_EQ(Core::Utils::try_demangle(typeid(Teuchos::RCP<DEMANGLE_TEST::TestStruct>).name()),
        "Teuchos::RCP<DEMANGLE_TEST::TestStruct>");
  }
}  // namespace
