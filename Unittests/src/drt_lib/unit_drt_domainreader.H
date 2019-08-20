/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the DomainReader class

\level 1

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/

#ifndef UNIT_DRT_DOMAINREADER_H_
#define UNIT_DRT_DOMAINREADER_H_

#include <cxxtest/TestSuite.h>

#include "src/drt_lib/drt_domainreader.H"
#include "stubs/DatFileReader_stub.H"

#include "src/drt_lib/drt_discret.H"

namespace DRT
{
  namespace INPUT
  {
    class DomainReader_TestSuite;
  }
}  // namespace DRT

class DRT::INPUT::DomainReader_TestSuite : public CxxTest::TestSuite
{
 private:
  Teuchos::RCP<Discretization> testdis_;
  DatFileReaderStub testreader_;

 public:
  void setUp() { testdis_ = Teuchos::rcp(new Discretization("dummy", testreader_.Comm())); }

  void tearDown() { testdis_ = Teuchos::null; }

  void test_MyDis0()
  {
    DomainReader dr(testdis_, testreader_, "unittestsection");
    TS_ASSERT_EQUALS(testdis_, dr.MyDis());
  }

  void test_MyDis1()
  {
    DomainReader dr(testdis_, testreader_, "unittestsection", "unittesttype");
    TS_ASSERT_EQUALS(testdis_, dr.MyDis());
  }

  void test_MyDis2()
  {
    std::set<std::string> dummyset;
    dummyset.insert("type0");
    dummyset.insert("type1");
    DomainReader dr(testdis_, testreader_, "unittestsection", dummyset);
    TS_ASSERT_EQUALS(testdis_, dr.MyDis());
  }

};  // class DomainReader_TestSuite

#endif  // UNIT_DRT_DOMAINREADER_H_
