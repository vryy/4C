/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the DomainReader class

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_lib_discret.H"
#include "baci_lib_domainreader.H"

#include <Epetra_SerialComm.h>

namespace
{
  using namespace BACI;

  class DatFileReaderStub : public INPUT::DatFileReader
  {
   public:
    DatFileReaderStub() = default;
    Teuchos::RCP<Epetra_Comm> Comm() const { return Teuchos::rcp(new Epetra_SerialComm); }
  };

  class DomainReaderTest : public ::testing::Test
  {
   public:
    DomainReaderTest()
    {
      testdis_ = Teuchos::rcp(new DRT::Discretization("dummy", testreader_.Comm()));
    }

   protected:
    Teuchos::RCP<DRT::Discretization> testdis_;
    DatFileReaderStub testreader_;
  };

  TEST_F(DomainReaderTest, TestMyDis0)
  {
    INPUT::DomainReader dr(testdis_, testreader_, "unittestsection");
    EXPECT_EQ(testdis_, dr.MyDis());
  }

}  // namespace
