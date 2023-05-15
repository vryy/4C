/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the beam3_euler_bernoulli class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <Epetra_SerialComm.h>
#include <array>

#include "lib_element.H"
#include "beam3_euler_bernoulli.H"

const double testTolerance = 1e-14;

namespace
{
  class Beam3eb : public ::testing::Test
  {
   public:
    Beam3eb()
    {
      testdis_ =
          Teuchos::rcp(new DRT::Discretization("Beam3eb", Teuchos::rcp(new Epetra_SerialComm)));

      std::vector<double> xrefe{-0.05, 0.05, 0.3, 0.45, -0.05, 0.1};

      for (int lid = 0; lid < 2; ++lid)
        testdis_->AddNode(Teuchos::rcp(new DRT::Node(lid, &xrefe[3 * lid], 0)));

      testele_ = Teuchos::rcp(new DRT::ELEMENTS::Beam3eb(0, 0));
      std::array<int, 2> node_ids{0, 1};
      testele_->SetNodeIds(2, node_ids.data());

      // create 1 element discretization
      testdis_->AddElement(testele_);
      testdis_->FillComplete(false, false, false);

      testele_->SetUpReferenceGeometry(xrefe);
    }

   protected:
    //! dummy discretization for holding element and node pointers
    Teuchos::RCP<DRT::Discretization> testdis_;
    //! the beam3eb element to be tested
    Teuchos::RCP<DRT::ELEMENTS::Beam3eb> testele_;
  };

  /**
   * Test reference length calculation of Euler-Bernoulli beam
   */
  TEST_F(Beam3eb, RefLength)
  {
    EXPECT_NEAR(testele_->RefLength(), 0.5477225575051661, testTolerance);
  }

}  // namespace
