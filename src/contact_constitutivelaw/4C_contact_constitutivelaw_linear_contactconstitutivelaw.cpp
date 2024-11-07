// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::LinearConstitutiveLawParams(
    const std::shared_ptr<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->get<double>("A")),
      b_(container->get<double>("B"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::create_constitutive_law()
{
  return std::make_shared<CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw>(this);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::LinearConstitutiveLaw(
    CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params)
    : params_(params)
{
}
/*----------------------------------------------------------------------*
 |  Evaluate the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::evaluate(double gap, CONTACT::Node* cnode)
{
  if (gap + params_->get_offset() > 0)
  {
    FOUR_C_THROW("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  return params_->getdata() * (gap + params_->get_offset()) + params_->get_b();
}  // end of linear_coconstlaw evaluate
/*----------------------------------------------------------------------*
 |  Calculate the derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::evaluate_deriv(
    double gap, CONTACT::Node* cnode)
{
  if (gap + params_->get_offset() > 0)
  {
    FOUR_C_THROW("You should not be here. The Evaluate function is only tested for active nodes.");
  }
  return params_->getdata();
}

FOUR_C_NAMESPACE_CLOSE
