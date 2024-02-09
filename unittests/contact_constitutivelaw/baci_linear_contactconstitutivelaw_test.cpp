/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the linear constact constitutivelaw

\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"

namespace
{

  using namespace BACI;

  class LinearConstitutiveLawTest : public ::testing::Test
  {
   public:
    LinearConstitutiveLawTest()
    {
      // initialize container for material parameters
      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
          Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
              1, INPAR::CONTACT::ConstitutiveLawType::colaw_linear, "Linear Constitutivelaw"));

      // add parameters to container
      container->Add("A", 1.5);
      container->Add("B", 0.0);
      container->Add("Offset", 0.5);

      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> linearcoconstlaw =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
      coconstlaw_ = linearcoconstlaw;
    }

    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
  };

  //! test member function Evaluate
  TEST_F(LinearConstitutiveLawTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(1.0));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(-0.25));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->Evaluate(-0.75), -0.375, 1.e-15);
  }

  //! test member function EvaluateDeriv
  TEST_F(LinearConstitutiveLawTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->EvaluateDeriv(-0.75), 1.5, 1.e-15);
    EXPECT_ANY_THROW(coconstlaw_->EvaluateDeriv(-0.25));
  }
}  // namespace