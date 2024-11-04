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

  struct Base
  {
    virtual ~Base() = default;
  };

  struct Derived : Base
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

  TEST(DemangleTest, PointerOfRCP)
  {
    Teuchos::RCP<DemangleTest::TestStruct>* ptr = nullptr;
    EXPECT_EQ(Core::Utils::get_dynamic_type_name(ptr), "Teuchos::RCP<DemangleTest::TestStruct>*");
  }

  TEST(DemangleTest, BaseDerivedRef)
  {
    DemangleTest::Derived d;
    DemangleTest::Base& base_ref = d;
    EXPECT_EQ(Core::Utils::get_dynamic_type_name(base_ref), "DemangleTest::Derived");
  }

  TEST(DemangleTest, BaseDerivedPtr)
  {
    DemangleTest::Derived d;
    DemangleTest::Base* base_ptr = &d;
    EXPECT_EQ(Core::Utils::get_dynamic_type_name(base_ptr), "DemangleTest::Base*");
    EXPECT_EQ(Core::Utils::get_dynamic_type_name(*base_ptr), "DemangleTest::Derived");
  }
}  // namespace
