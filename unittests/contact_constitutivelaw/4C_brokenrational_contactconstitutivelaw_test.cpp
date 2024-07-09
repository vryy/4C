/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the brokenrational contact constitutivelaw

\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_contact_constitutivelaw_brokenrational_contactconstitutivelaw.hpp"
#include "4C_contact_node.hpp"

namespace
{
  using namespace FourC;

  class BrokenrationalConstitutiveLawTest : public ::testing::Test
  {
   public:
    BrokenrationalConstitutiveLawTest()
    {
      // initialize container for material parameters
      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
          Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(1,
              Inpar::CONTACT::ConstitutiveLawType::colaw_brokenrational,
              "Brokenrational Constitutivelaw"));

      // add parameters to container
      container->add("A", -2.);
      container->add("B", 4.);
      container->add("C", -0.5);
      container->add("Offset", 0.5);

      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> brokenrationalcoconstlaw =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::factory(container);
      coconstlaw_ = brokenrationalcoconstlaw;
    }

    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;

    Teuchos::RCP<CONTACT::Node> cnode;
  };

  //! test member function Evaluate
  TEST_F(BrokenrationalConstitutiveLawTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->evaluate(1.0, cnode.get()));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->evaluate(-0.25, cnode.get()));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->evaluate(-2.5, cnode.get()), -0.5, 1.e-15);
  }

  //! test member function EvaluateDeriv
  TEST_F(BrokenrationalConstitutiveLawTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->evaluate_deriv(-2.5, cnode.get()), 0.5, 1.e-15);
    EXPECT_ANY_THROW(coconstlaw_->evaluate_deriv(-0.25, cnode.get()));
  }
}  // namespace