/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for grid generator functionality

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_io_gridgenerator.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_par_material.hpp"

#include <Epetra_SerialComm.h>

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

  class GridGeneratorTest : public ::testing::Test
  {
   public:
    GridGeneratorTest()
    {
      inputData_.bottom_corner_point_ = std::array<double, 3>{-1.0, -2.0, -3.0};
      inputData_.top_corner_point_ = std::array<double, 3>{2.5, 3.5, 4.5};
      inputData_.interval_ = std::array<int, 3>{5, 10, 15};
      inputData_.node_gid_of_first_new_node_ = 17;
    };

   protected:
    void SetUp() override
    {
      CreateMaterialInGlobalProblem();
      comm_ = Teuchos::rcp(new Epetra_SerialComm);
      IO::cout.setup(false, false, false, IO::standard, comm_, 0, 0, "dummyFilePrefix");
      testdis_ = Teuchos::rcp(new DRT::Discretization("dummy", comm_));
    }

    void TearDown() override { IO::cout.close(); }

   public:
    IO::GRIDGENERATOR::RectangularCuboidInputs inputData_{};
    Teuchos::RCP<DRT::Discretization> testdis_;
    Teuchos::RCP<Epetra_Comm> comm_;
  };

  TEST_F(GridGeneratorTest, TestGridGeneratorWithHex8Elements)
  {
    inputData_.elementtype_ = "SOLIDH8";
    inputData_.distype_ = "HEX8";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear EAS none";

    IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const auto nodePosition = lastNode->X();

    EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
    EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
    EXPECT_NEAR(nodePosition[2], 4.5, 1e-14);
    EXPECT_EQ(testdis_->NumMyRowNodes(), 1056);
    EXPECT_EQ(testdis_->NumMyRowElements(), 750);
    EXPECT_EQ(lastNode->Id(), 7177);
  }

  TEST_F(GridGeneratorTest, TestGridGeneratorWithRotatedHex8Elements)
  {
    inputData_.elementtype_ = "SOLIDH8";
    inputData_.distype_ = "HEX8";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear EAS none";
    inputData_.rotation_angle_ = std::array<double, 3>{30.0, 10.0, 7.0};

    IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const auto nodePosition = lastNode->X();

    EXPECT_NEAR(nodePosition[0], 2.6565639116964181, 1e-14);
    EXPECT_NEAR(nodePosition[1], 4.8044393443812901, 1e-14);
    EXPECT_NEAR(nodePosition[2], 2.8980306453470042, 1e-14);
    EXPECT_EQ(testdis_->NumMyRowNodes(), 1056);
    EXPECT_EQ(testdis_->NumMyRowElements(), 750);
    EXPECT_EQ(lastNode->Id(), 7177);
  }

  TEST_F(GridGeneratorTest, TestGridGeneratorWithHex27Elements)
  {
    inputData_.elementtype_ = "SOLIDH27";
    inputData_.distype_ = "HEX27";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear";

    IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const auto nodePosition = lastNode->X();

    EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
    EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
    EXPECT_NEAR(nodePosition[2], 4.5, 1e-14);
    EXPECT_EQ(testdis_->NumMyRowNodes(), 7161);
    EXPECT_EQ(testdis_->NumMyRowElements(), 750);
    EXPECT_EQ(lastNode->Id(), 7177);
  }

  TEST_F(GridGeneratorTest, TestGridGeneratorWithWedge6Elements)
  {
    inputData_.elementtype_ = "SOLIDW6";
    inputData_.distype_ = "WEDGE6";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear";
    inputData_.autopartition_ = true;

    IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const auto nodePosition = lastNode->X();

    EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
    EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
    EXPECT_NEAR(nodePosition[2], 4.5, 1e-14);
    EXPECT_EQ(testdis_->NumMyRowNodes(), 1056);
    EXPECT_EQ(testdis_->NumMyRowElements(), 1500);
    EXPECT_EQ(lastNode->Id(), 7177);
  }

}  // namespace
