/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for building the nodal coordinate vector of a discretization
       based on a nodal rowmap

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_io_gridgenerator.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_par_material.hpp"

#include <Epetra_MpiComm.h>

namespace
{
  using namespace FourC;

  void CreateMaterialInGlobalProblem()
  {
    const auto mat_stvenant = Teuchos::rcp(new MAT::PAR::Material(
        1, CORE::Materials::MaterialType::m_stvenant, "MAT_Struct_StVenantKirchhoff"));

    mat_stvenant->Add("YOUNG", 1.0);
    mat_stvenant->Add("NUE", 0.1);
    mat_stvenant->Add("DENS", 2.0);

    GLOBAL::Problem::Instance()->Materials()->Insert(1, mat_stvenant);
  }

  // Serial discretization nodal method tests
  class BuildNodeCoordinatesTest : public testing::Test
  {
   public:
    BuildNodeCoordinatesTest()
    {
      CreateMaterialInGlobalProblem();

      comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      test_discretization_ = Teuchos::rcp(new DRT::Discretization("dummy", comm_));

      IO::cout.setup(false, false, false, IO::standard, comm_, 0, 0, "dummyFilePrefix");

      // results in 27 nodes
      inputData_.bottom_corner_point_ = std::array<double, 3>{0.0, 0.0, 0.0};
      inputData_.top_corner_point_ = std::array<double, 3>{1.0, 1.0, 1.0};
      inputData_.interval_ = std::array<int, 3>{2, 2, 4};
      inputData_.node_gid_of_first_new_node_ = 0;

      inputData_.elementtype_ = "SOLIDH8";
      inputData_.distype_ = "HEX8";
      inputData_.elearguments_ = "MAT 1 KINEM nonlinear EAS none";

      IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(
          *test_discretization_, inputData_, true);

      test_discretization_->FillComplete(false, false, false);
    }

    void TearDown() override { IO::cout.close(); }

   protected:
    IO::GRIDGENERATOR::RectangularCuboidInputs inputData_{};
    Teuchos::RCP<DRT::Discretization> test_discretization_;
    Teuchos::RCP<Epetra_Comm> comm_;
  };

  TEST_F(BuildNodeCoordinatesTest, NodalCoordinatesDefault)
  {
    // build node coordinates based on the node row map of the whole discretization
    Teuchos::RCP<Epetra_MultiVector> nodal_test_coordinates =
        test_discretization_->BuildNodeCoordinates();

    if (comm_->MyPID() == 0)
    {
      EXPECT_EQ(nodal_test_coordinates->MyLength(), test_discretization_->NumMyRowNodes());
      EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

      std::array<double, 54> coords;
      nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

      // first coordinate
      EXPECT_NEAR(coords[0], 0.0, 1e-14);
      EXPECT_NEAR(coords[18], 0.0, 1e-14);
      EXPECT_NEAR(coords[36], 0.0, 1e-14);

      // last coordinate
      EXPECT_NEAR(coords[17], 1.0, 1e-14);
      EXPECT_NEAR(coords[35], 1.0, 1e-14);
      EXPECT_NEAR(coords[53], 0.25, 1e-14);
    }
    else if (comm_->MyPID() == 1)
    {
      EXPECT_EQ(nodal_test_coordinates->MyLength(), test_discretization_->NumMyRowNodes());
      EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

      std::array<double, 54> coords;
      nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

      // first coordinate
      EXPECT_NEAR(coords[0], 0.0, 1e-14);
      EXPECT_NEAR(coords[18], 0.0, 1e-14);
      EXPECT_NEAR(coords[36], 0.5, 1e-14);

      // last coordinate
      EXPECT_NEAR(coords[17], 1.0, 1e-14);
      EXPECT_NEAR(coords[35], 1.0, 1e-14);
      EXPECT_NEAR(coords[53], 0.75, 1e-14);
    }
    else if (comm_->MyPID() == 2)
    {
      EXPECT_EQ(nodal_test_coordinates->MyLength(), test_discretization_->NumMyRowNodes());
      EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

      std::array<double, 27> coords;
      nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

      // first coordinate
      EXPECT_NEAR(coords[0], 0.0, 1e-14);
      EXPECT_NEAR(coords[9], 0.0, 1e-14);
      EXPECT_NEAR(coords[18], 1.0, 1e-14);

      // last coordinate
      EXPECT_NEAR(coords[8], 1.0, 1e-14);
      EXPECT_NEAR(coords[17], 1.0, 1e-14);
      EXPECT_NEAR(coords[26], 1.0, 1e-14);
    }
  }

  TEST_F(BuildNodeCoordinatesTest, NodalCoordinatesPartialMap)
  {
    // build node coordinates based on the node row map of first partial discretization
    {
      std::array<int, 4> nodeList{0, 2, 4, 10};  // GID list of first 4 elements
      Teuchos::RCP<Epetra_Map> node_row_map =
          Teuchos::rcp(new Epetra_Map(-1, nodeList.size(), nodeList.data(), 0, *comm_));
      Teuchos::RCP<Epetra_MultiVector> nodal_test_coordinates =
          test_discretization_->BuildNodeCoordinates(node_row_map);

      if (comm_->MyPID() == 0)
      {
        EXPECT_EQ(nodal_test_coordinates->MyLength(), 4);
        EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

        std::array<double, 12> coords;
        nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

        // first coordinate
        EXPECT_NEAR(coords[0], 0.0, 1e-14);
        EXPECT_NEAR(coords[4], 0.0, 1e-14);
        EXPECT_NEAR(coords[8], 0.0, 1e-14);

        // last coordinate
        EXPECT_NEAR(coords[3], 0.0, 1e-14);
        EXPECT_NEAR(coords[7], 0.5, 1e-14);
        EXPECT_NEAR(coords[11], 0.0, 1e-14);
      }
    }

    // build node coordinates based on the node row map of second partial discretization
    {
      Teuchos::RCP<Epetra_Map> node_row_map = Teuchos::null;
      if (comm_->MyPID() == 0)
      {
        std::array<int, 2> nodeList{50, 62};
        node_row_map =
            Teuchos::rcp(new Epetra_Map(-1, nodeList.size(), nodeList.data(), 0, *comm_));
      }
      else if (comm_->MyPID() == 1)
      {
        std::array<int, 1> nodeList{114};
        node_row_map =
            Teuchos::rcp(new Epetra_Map(-1, nodeList.size(), nodeList.data(), 0, *comm_));
      }
      else if (comm_->MyPID() == 2)
      {
        std::array<int, 1> nodeList{212};
        node_row_map =
            Teuchos::rcp(new Epetra_Map(-1, nodeList.size(), nodeList.data(), 0, *comm_));
      }

      Teuchos::RCP<Epetra_MultiVector> nodal_test_coordinates =
          test_discretization_->BuildNodeCoordinates(node_row_map);

      if (comm_->MyPID() == 0)
      {
        EXPECT_EQ(nodal_test_coordinates->MyLength(), 2);
        EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

        std::array<double, 6> coords;
        nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

        EXPECT_NEAR(coords[0], 0.0, 1e-14);
        EXPECT_NEAR(coords[2], 0.0, 1e-14);
        EXPECT_NEAR(coords[4], 0.25, 1e-14);

        EXPECT_NEAR(coords[1], 0.5, 1e-14);
        EXPECT_NEAR(coords[3], 0.5, 1e-14);
        EXPECT_NEAR(coords[5], 0.25, 1e-14);
      }
      else if (comm_->MyPID() == 1)
      {
        EXPECT_EQ(nodal_test_coordinates->MyLength(), 1);
        EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

        std::array<double, 3> coords;
        nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

        EXPECT_NEAR(coords[0], 1.0, 1e-14);
        EXPECT_NEAR(coords[1], 0.5, 1e-14);
        EXPECT_NEAR(coords[2], 0.5, 1e-14);
      }
      else if (comm_->MyPID() == 2)
      {
        EXPECT_EQ(nodal_test_coordinates->MyLength(), 1);
        EXPECT_EQ(nodal_test_coordinates->NumVectors(), 3);

        std::array<double, 3> coords;
        nodal_test_coordinates->ExtractCopy(coords.data(), nodal_test_coordinates->MyLength());

        EXPECT_NEAR(coords[0], 0.5, 1e-14);
        EXPECT_NEAR(coords[1], 0.5, 1e-14);
        EXPECT_NEAR(coords[2], 1.0, 1e-14);
      }
    }
  }
}  // namespace
