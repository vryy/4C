/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the beam3_kirchhoff class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_beam3_kirchhoff.hpp"

#include "4C_fem_general_element.hpp"

#include <Epetra_SerialComm.h>

#include <array>

const double testTolerance = 1e-14;

namespace
{
  using namespace FourC;

  class Beam3k : public ::testing::Test
  {
   public:
    Beam3k()
    {
      testdis_ = Teuchos::rcp(
          new Core::FE::Discretization("Beam3k", Teuchos::rcp(new Epetra_SerialComm), 3));

      std::vector<std::vector<double>> xrefe{{-0.05, 0.05, 0.3}, {0.45, -0.05, 0.1}};

      for (int lid = 0; lid < 2; ++lid)
        testdis_->AddNode(Teuchos::rcp(new Core::Nodes::Node(lid, xrefe[lid], 0)));

      testele_ = Teuchos::rcp(new Discret::ELEMENTS::Beam3k(0, 0));
      std::array<int, 2> node_ids{0, 1};
      testele_->SetNodeIds(2, node_ids.data());

      // create 1 element discretization
      testdis_->add_element(testele_);
      testdis_->fill_complete(false, false, false);

      // setup internal beam element parameters
      // different data layout is necessary to call this method
      Core::LinAlg::Matrix<3, 1> coord1(true);
      coord1(0) = xrefe[0][0];
      coord1(1) = xrefe[0][1];
      coord1(2) = xrefe[0][2];
      Core::LinAlg::Matrix<3, 1> coord2(true);
      coord2(0) = xrefe[1][0];
      coord2(1) = xrefe[1][1];
      coord2(2) = xrefe[1][2];
      std::vector<Core::LinAlg::Matrix<3, 1>> xrefe_setup{coord1, coord2};

      // setup internal beam element parameters
      std::vector<double> rotrefe(9);
      rotrefe[0] = -2.135698785951414;
      rotrefe[1] = -1.1055190408131161;
      rotrefe[2] = -0.45792098016648797;
      rotrefe[3] = 0.09071600605476587;
      rotrefe[4] = -0.31314870676006484;
      rotrefe[5] = -0.5590172175309829;
      rotrefe[6] = -0.44757433200569813;
      rotrefe[7] = -0.14845112617443665;
      rotrefe[8] = -0.628849061811312;

      testele_->set_up_initial_rotations(rotrefe);
      testele_->set_up_reference_geometry(xrefe_setup);
    }

   protected:
    //! dummy discretization for holding element and node pointers
    Teuchos::RCP<Core::FE::Discretization> testdis_;
    //! the beam3k element to be tested
    Teuchos::RCP<Discret::ELEMENTS::Beam3k> testele_;
  };

  /**
   * Test reference length calculation of Kirchhoff-Love beam
   */
  TEST_F(Beam3k, RefLength)
  {
    EXPECT_NEAR(testele_->RefLength(), 0.61920435714496047, testTolerance);
  }

}  // namespace
