// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_discretization.hpp"
#include "4C_io_domainreader.hpp"

#include <Epetra_SerialComm.h>

namespace
{
  using namespace FourC;

  class DatFileReaderStub : public Core::IO::DatFileReader
  {
   public:
    DatFileReaderStub() = default;
    Teuchos::RCP<Epetra_Comm> comm() const { return Teuchos::make_rcp<Epetra_SerialComm>(); }
  };

  class DomainReaderTest : public ::testing::Test
  {
   public:
    DomainReaderTest()
    {
      testdis_ = Teuchos::make_rcp<Core::FE::Discretization>("dummy", testreader_.comm(), 3);
    }

   protected:
    Teuchos::RCP<Core::FE::Discretization> testdis_;
    DatFileReaderStub testreader_;
  };

  TEST_F(DomainReaderTest, TestMyDis0)
  {
    Core::IO::DomainReader dr(testdis_, testreader_, "unittestsection");
    EXPECT_EQ(testdis_, dr.my_dis());
  }

}  // namespace
