// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_global_data.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_utils_singleton_owner.hpp"

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
      testdis_ = std::make_shared<Core::FE::Discretization>("dummy", MPI_COMM_WORLD, 3);

      // create 8 nodes
      const std::array<int, 8> nodeids = {0, 1, 2, 3, 4, 5, 6, 7};
      std::vector<std::vector<double>> coords = {{-0.1, -0.2, -0.5}, {1.25, 0.23, 0.66},
          {1.20, 0.99, 0.5}, {-0.11, 1.20, 0.66}, {-0.10, -0.2, 1.9}, {1.00, 0.00, 1.90},
          {1.20, 0.99, 1.50}, {-0.11, -0.20, 1.66}};
      for (int lid = 0; lid < 8; ++lid)
        testdis_->add_node(std::make_shared<Core::Nodes::Node>(lid, coords[lid], 0));

      // create 1 element
      testele_ = std::make_shared<Discret::Elements::SoHex8>(0, 0);
      testele_->set_node_ids(8, nodeids.data());
      testdis_->add_element(testele_);
      testdis_->fill_complete(false, false, false);

      copytestele_ = std::make_shared<Discret::Elements::SoHex8>(*testele_);
    }

    //! dummy discretization for holding element and node pointers
    std::shared_ptr<Core::FE::Discretization> testdis_;
    //! the hex8 element to be tested
    std::shared_ptr<Discret::Elements::SoHex8> testele_;
    //! a copy of the hex8 element to test the copy constructor
    std::shared_ptr<Discret::Elements::SoHex8> copytestele_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
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
    EXPECT_EQ(testele_->num_dof_per_node(node_dummy), 3);
    EXPECT_EQ(copytestele_->num_dof_per_node(node_dummy), 3);
  }

  /**
   * Test the polynomial degree
   */
  TEST_F(SoHex8Test, TestDegree)
  {
    EXPECT_EQ(testele_->degree(), 1);
    EXPECT_EQ(copytestele_->degree(), 1);
  }

  /**
   * Test the number of volumes the element is composed of
   */
  TEST_F(SoHex8Test, TestNumVolume)
  {
    EXPECT_EQ(testele_->num_volume(), 1);
    EXPECT_EQ(copytestele_->num_volume(), 1);
  }

  /**
   * Test the number of surfaces the element is composed of
   */
  TEST_F(SoHex8Test, TestNumSurface)
  {
    EXPECT_EQ(testele_->num_surface(), 6);
    EXPECT_EQ(copytestele_->num_surface(), 6);
  }

  /**
   * Test the number of lines the element is composed of
   */
  TEST_F(SoHex8Test, TestNumLine)
  {
    EXPECT_EQ(testele_->num_line(), 12);
    EXPECT_EQ(copytestele_->num_line(), 12);
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
