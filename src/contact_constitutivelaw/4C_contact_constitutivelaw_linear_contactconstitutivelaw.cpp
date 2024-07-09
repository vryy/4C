/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a linear contact constitutive law
\level 3

*----------------------------------------------------------------------*/


#include "4C_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::LinearConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->get<double>("A")),
      b_(container->get<double>("B"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::create_constitutive_law()
{
  return Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw(this));
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
