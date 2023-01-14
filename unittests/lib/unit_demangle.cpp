/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for demangle utility
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include "lib_demangle.H"
#include "Teuchos_RCPDecl.hpp"

namespace DEMANGLE_TEST
{
  struct TestStruct
  {
  };
}  // namespace DEMANGLE_TEST

namespace
{

  TEST(DemangleTest, Struct)
  {
    EXPECT_EQ(DRT::UTILS::TryDemangle(typeid(DEMANGLE_TEST::TestStruct).name()),
        "DEMANGLE_TEST::TestStruct");
  }

  TEST(DemangleTest, RCP)
  {
    EXPECT_EQ(DRT::UTILS::TryDemangle(typeid(Teuchos::RCP<DEMANGLE_TEST::TestStruct>).name()),
        "Teuchos::RCP<DEMANGLE_TEST::TestStruct>");
  }
}  // namespace
