/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the DomainReader class

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <Epetra_SerialComm.h>
#include "drt_discret.H"
#include "drt_domainreader.H"

namespace
{
  class DatFileReaderStub : public DRT::INPUT::DatFileReader
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
    DRT::INPUT::DomainReader dr(testdis_, testreader_, "unittestsection");
    EXPECT_EQ(testdis_, dr.MyDis());
  }

  TEST_F(DomainReaderTest, TestMyDis1)
  {
    DRT::INPUT::DomainReader dr(testdis_, testreader_, "unittestsection", "unittesttype");
    EXPECT_EQ(testdis_, dr.MyDis());
  }

  TEST_F(DomainReaderTest, TestMyDis2)
  {
    std::set<std::string> dummyset;
    dummyset.insert("type0");
    dummyset.insert("type1");
    DRT::INPUT::DomainReader dr(testdis_, testreader_, "unittestsection", "unittesttype");
    EXPECT_EQ(testdis_, dr.MyDis());
  }
}  // namespace
