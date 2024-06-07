/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the So_hex8 class

\level 3

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_so3_hex8.hpp"

#include <Epetra_SerialComm.h>

#include <array>

namespace
{
  using namespace FourC;

  class SoHex8Test : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // create a discretization, that creates node to element pointers and keeps the nodes alive
      testdis_ = Teuchos::rcp(
          new Discret::Discretization("dummy", Teuchos::rcp(new Epetra_SerialComm), 3));

      // create 8 nodes
      const std::array<int, 8> nodeids = {0, 1, 2, 3, 4, 5, 6, 7};
      std::vector<std::vector<double>> coords = {{-0.1, -0.2, -0.5}, {1.25, 0.23, 0.66},
          {1.20, 0.99, 0.5}, {-0.11, 1.20, 0.66}, {-0.10, -0.2, 1.9}, {1.00, 0.00, 1.90},
          {1.20, 0.99, 1.50}, {-0.11, -0.20, 1.66}};
      for (int lid = 0; lid < 8; ++lid)
        testdis_->AddNode(Teuchos::rcp(new Core::Nodes::Node(lid, coords[lid], 0)));

      // create 1 element
      testele_ = Teuchos::rcp(new Discret::ELEMENTS::SoHex8(0, 0));
      testele_->SetNodeIds(8, nodeids.data());
      testdis_->add_element(testele_);
      testdis_->fill_complete(false, false, false);

      copytestele_ = Teuchos::rcp(new Discret::ELEMENTS::SoHex8(*testele_));
    }

    // Delete pointers.
    void TearDown() override
    {
      copytestele_ = Teuchos::null;
      testele_ = Teuchos::null;
      testdis_ = Teuchos::null;

      // We need to make sure the Global::Problem instance created in setUp is deleted again. If
      // this is not done, some troubles arise where unit tests influence each other on some
      // configurations. We suspect that missing singleton destruction might be the reason for that.
      Global::Problem::Done();
    }
    //! dummy discretization for holding element and node pointers
    Teuchos::RCP<Discret::Discretization> testdis_;
    //! the hex8 element to be tested
    Teuchos::RCP<Discret::ELEMENTS::SoHex8> testele_;
    //! a copy of the hex8 element to test the copy constructor
    Teuchos::RCP<Discret::ELEMENTS::SoHex8> copytestele_;
  };

  /**
   * Test Number of DOFs per element
   */
  TEST_F(SoHex8Test, num_dof_per_element)
  {
    EXPECT_EQ(testele_->num_dof_per_element(), 0);
    EXPECT_EQ(copytestele_->num_dof_per_element(), 0);
  }

  /**
   * Test Number of DOFs per node function
   */
  TEST_F(SoHex8Test, TestNumDofPerNode)
  {
    std::vector<double> pd = {1, 2, 3};
    Core::Nodes::Node node_dummy(0, pd, false);
    EXPECT_EQ(testele_->NumDofPerNode(node_dummy), 3);
    EXPECT_EQ(copytestele_->NumDofPerNode(node_dummy), 3);
  }

  /**
   * Test the polynomial degree
   */
  TEST_F(SoHex8Test, TestDegree)
  {
    EXPECT_EQ(testele_->Degree(), 1);
    EXPECT_EQ(copytestele_->Degree(), 1);
  }

  /**
   * Test the number of volumes the element is composed of
   */
  TEST_F(SoHex8Test, TestNumVolume)
  {
    EXPECT_EQ(testele_->NumVolume(), 1);
    EXPECT_EQ(copytestele_->NumVolume(), 1);
  }

  /**
   * Test the number of surfaces the element is composed of
   */
  TEST_F(SoHex8Test, TestNumSurface)
  {
    EXPECT_EQ(testele_->NumSurface(), 6);
    EXPECT_EQ(copytestele_->NumSurface(), 6);
  }

  /**
   * Test the number of lines the element is composed of
   */
  TEST_F(SoHex8Test, TestNumLine)
  {
    EXPECT_EQ(testele_->NumLine(), 12);
    EXPECT_EQ(copytestele_->NumLine(), 12);
  }

  /**
   * Test the calculation of the element center coordinates
   */
  TEST_F(SoHex8Test, TestElementCenterRefeCoords)
  {
    std::array<double, 3> midpoint = {0.528750000000000, 0.351250000000000, 1.035000000000000};
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(testele_->element_center_refe_coords()[i], midpoint[i], 1e-14);
      EXPECT_NEAR(copytestele_->element_center_refe_coords()[i], midpoint[i], 1e-14);
    }
  }

}  // namespace
