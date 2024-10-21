// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_contact_constitutivelaw_cubic_contactconstitutivelaw.hpp"
#include "4C_contact_node.hpp"

namespace
{
  using namespace FourC;
  // class implementation
  class CubicConstitutiveLawTest : public ::testing::Test
  {
   public:
    CubicConstitutiveLawTest()
    {
      // initialize container for material parameters
      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
          Teuchos::make_rcp<CONTACT::CONSTITUTIVELAW::Container>(
              1, Inpar::CONTACT::ConstitutiveLawType::colaw_cubic, "Cubic Constitutivelaw");

      // add parameters to container
      container->add("A", 1.5);
      container->add("B", 2.0);
      container->add("C", 3.0);
      container->add("D", 0.0);
      container->add("Offset", 0.5);

      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> cubiccoconstlaw =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::factory(container);
      coconstlaw_ = cubiccoconstlaw;
    }

    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;

    Teuchos::RCP<CONTACT::Node> cnode;
  };


  //! test member function Evaluate
  TEST_F(CubicConstitutiveLawTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->evaluate(1.0, cnode.get()));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->evaluate(-0.25, cnode.get()));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->evaluate(-0.75, cnode.get()), -0.8984375, 1.e-15);
  }

  //! test member function EvaluateDeriv
  TEST_F(CubicConstitutiveLawTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->evaluate_deriv(-0.75, cnode.get()), 4.28125, 1.e-15);
    EXPECT_ANY_THROW(coconstlaw_->evaluate_deriv(-0.25, cnode.get()));
  }
}  // namespace
