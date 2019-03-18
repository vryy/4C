/*!---------------------------------------------------------------------*
\file unit_so_hex8.H

\brief Unittests for the So_hex8 class

\level 3

\maintainer Anh-Tu Vuong
*-----------------------------------------------------------------------*/

#ifndef UNIT_SO_HEX8_H_
#define UNIT_SO_HEX8_H_

#include <cxxtest/TestSuite.h>
#include <Epetra_SerialComm.h>
#include "src/drt_so3/so_hex8.H"
#include "src/drt_lib/drt_globalproblem.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class So_hex8_TestSuite;
  }
}  // namespace DRT

/**
 * Class to test main functionality of hex8 element
 */
class DRT::ELEMENTS::So_hex8_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * Set up the testing environment.
   */
  void setUp()
  {
    // create a discretization, that creates node to element pointers and keeps the nodes alive
    testdis_ = Teuchos::rcp(new Discretization("dummy", Teuchos::rcp(new Epetra_SerialComm)));

    // create 8 nodes
    const int nodeids[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    double coords[24] = {-0.1, -0.2, -0.5, 1.25, 0.23, 0.66, 1.20, 0.99, 0.5, -0.11, 1.20, 0.66,
        -0.10, -0.2, 1.9, 1.00, 0.00, 1.90, 1.20, 0.99, 1.50, -0.11, -0.20, 1.66};
    for (int lid = 0; lid < 8; ++lid)
      testdis_->AddNode(Teuchos::rcp(new Node(lid, &coords[3 * lid], 0)));

    // create 1 element
    testele_ = Teuchos::rcp(new So_hex8(0, 0));
    testele_->SetNodeIds(8, nodeids);
    testdis_->AddElement(testele_);
    testdis_->FillComplete(false, false, false);

    copytestele_ = Teuchos::rcp(new So_hex8(*testele_));
  }

  /**
   * Delete pointers.
   */
  void tearDown()
  {
    copytestele_ = Teuchos::null;
    testele_ = Teuchos::null;
    testdis_ = Teuchos::null;

    // We need to make sure the DRT::Problem instance created in setUp is deleted again. If this is
    // not done, some troubles arise where unit tests influence each other on some configurations.
    // We suspect that missing singleton destruction might be the reason for that.
    DRT::Problem::Done();
  }

  /**
   * Test Number of DOFs per element
   */
  void test_NumDofPerElement()
  {
    TS_ASSERT_EQUALS(testele_->NumDofPerElement(), 0);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerElement(), 0);
  }

  /**
   * Test Number of DOFs per node function
   */
  void test_NumDofPerNode()
  {
    double pd[3] = {1, 2, 3};
    DRT::Node node_dummy(0, pd, false);
    TS_ASSERT_EQUALS(testele_->NumDofPerNode(node_dummy), 3);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerNode(node_dummy), 3);
  }

  /**
   * Test the polynomial degree
   */
  void test_Degree()
  {
    TS_ASSERT_EQUALS(testele_->Degree(), 1);
    TS_ASSERT_EQUALS(copytestele_->Degree(), 1);
  }

  /**
   * Test the number of volumes the element is composed of
   */
  void test_NumVolume()
  {
    TS_ASSERT_EQUALS(testele_->NumVolume(), 1);
    TS_ASSERT_EQUALS(copytestele_->NumVolume(), 1);
  }

  /**
   * Test the number of surfaces the element is composed of
   */
  void test_NumSurface()
  {
    TS_ASSERT_EQUALS(testele_->NumSurface(), 6);
    TS_ASSERT_EQUALS(copytestele_->NumSurface(), 6);
  }

  /**
   * Test the number of lines the element is composed of
   */
  void test_NumLine()
  {
    TS_ASSERT_EQUALS(testele_->NumLine(), 12);
    TS_ASSERT_EQUALS(copytestele_->NumLine(), 12);
  }

  /**
   * Test the calculation of the element center coordinates
   */
  void test_ElementCenterRefeCoords()
  {
    double midpoint[3] = {0.528750000000000, 0.351250000000000, 1.035000000000000};
    for (int i = 0; i < 3; ++i)
    {
      TS_ASSERT_DELTA(testele_->ElementCenterRefeCoords()[i], midpoint[i], 1e-14);
      TS_ASSERT_DELTA(copytestele_->ElementCenterRefeCoords()[i], midpoint[i], 1e-14);
    }
  }

 private:
  //! dummy discretization for holding element and node pointers
  Teuchos::RCP<Discretization> testdis_;
  //! the hex8 element to be tested
  Teuchos::RCP<So_hex8> testele_;
  //! a copy of the hex8 element to test the copy constructor
  Teuchos::RCP<So_hex8> copytestele_;

};  // class So_hex8_TestSuite


#endif /* UNIT_SO_Hex8_H_ */
