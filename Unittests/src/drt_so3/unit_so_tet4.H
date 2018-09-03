/*!---------------------------------------------------------------------*
\file unit_so_tet4.H

\brief Unittests for the So_tet4 class

\level 3

\maintainer Michael Gee
*-----------------------------------------------------------------------*/

#ifndef UNIT_SO_TET4_H_
#define UNIT_SO_TET4_H_

#include <cxxtest/TestSuite.h>
#include <Epetra_SerialComm.h>
#include "src/drt_so3/so_tet4.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class So_tet4_TestSuite;
  }
}  // namespace DRT

class DRT::ELEMENTS::So_tet4_TestSuite : public CxxTest::TestSuite
{
 private:
  Teuchos::RCP<Discretization> testdis_;
  Teuchos::RCP<So_tet4> testele_;
  Teuchos::RCP<So_tet4> copytestele_;

 public:
  void setUp()
  {
    // create a discretization, that creates node to element pointers and keeps the nodes alive
    testdis_ = Teuchos::rcp(new DRT::Discretization("dummy", Teuchos::rcp(new Epetra_SerialComm)));

    // create 4 nodes
    const int nodeids[] = {0, 1, 2, 3};
    double coords[] = {-0.1, -0.2, -0.5, 1.25, 0.23, 0.66, 1.20, 0.99, 0.5, -0.10, -0.2, 1.96};
    for (int lid = 0; lid < 4; ++lid)
      testdis_->AddNode(Teuchos::rcp(new DRT::Node(lid, &coords[3 * lid], 0)));

    // create 1 element
    testele_ = Teuchos::rcp(new DRT::ELEMENTS::So_tet4(0, 0));
    testele_->SetNodeIds(4, nodeids);
    testdis_->AddElement(testele_);
    testdis_->FillComplete(false, false, false);

    copytestele_ = Teuchos::rcp(new DRT::ELEMENTS::So_tet4(*testele_));
  }

  void tearDown()
  {
    copytestele_ = Teuchos::null;
    testele_ = Teuchos::null;
    testdis_ = Teuchos::null;
  }

  void test_NumDofPerNode()
  {
    double pd[3] = {1, 2, 3};
    DRT::Node node_dummy(0, pd, false);
    TS_ASSERT_EQUALS(testele_->NumDofPerNode(node_dummy), 3);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerNode(node_dummy), 3);
  }

  void test_NumDofPerElement()
  {
    TS_ASSERT_EQUALS(testele_->NumDofPerElement(), 0);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerElement(), 0);
  }

  void test_Degree()
  {
    TS_ASSERT_EQUALS(testele_->Degree(), 1);
    TS_ASSERT_EQUALS(copytestele_->Degree(), 1);
  }

  void test_NumVolume()
  {
    TS_ASSERT_EQUALS(testele_->NumVolume(), 1);
    TS_ASSERT_EQUALS(copytestele_->NumVolume(), 1);
  }

  void test_NumSurface()
  {
    TS_ASSERT_EQUALS(testele_->NumSurface(), 4);
    TS_ASSERT_EQUALS(copytestele_->NumSurface(), 4);
  }

  void test_NumLine()
  {
    TS_ASSERT_EQUALS(testele_->NumLine(), 6);
    TS_ASSERT_EQUALS(copytestele_->NumLine(), 6);
  }

  void test_ElementCenterRefeCoords()
  {
    double midpoint[3] = {0.5625, 0.2050, 0.6550};
    for (int i = 0; i < 3; ++i)
    {
      TS_ASSERT_DELTA(testele_->ElementCenterRefeCoords()[i], midpoint[i], 1e-14);
      TS_ASSERT_DELTA(copytestele_->ElementCenterRefeCoords()[i], midpoint[i], 1e-14);
    }
  }
};

#endif /* UNIT_SO_TET4_H_ */
