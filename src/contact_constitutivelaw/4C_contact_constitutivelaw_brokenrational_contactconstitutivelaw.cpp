/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a contact constitutive law following a broken rational function

\level 1

*----------------------------------------------------------------------*/


#include "4C_contact_constitutivelaw_brokenrational_contactconstitutivelaw.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams::BrokenRationalConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->Get<double>("A")),
      b_(container->Get<double>("B")),
      c_(container->Get<double>("C"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams::CreateConstitutiveLaw()
{
  return Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw::BrokenRationalConstitutiveLaw(
    CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams* params)
    : params_(params)
{
}
/*----------------------------------------------------------------------*
 |  Evaluate the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw::Evaluate(
    double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0)
  {
    FOUR_C_THROW("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  double result = -1.0;
  gap = -gap;
  result *=
      (params_->GetA() * 1. / (gap - params_->GetOffset() - params_->GetB()) + params_->GetC());
  if (result > 0)
    FOUR_C_THROW(
        "The constitutive function you are using seems to be positive, even though the gap is "
        "negative. Please check your coefficients!");
  return result;
}
/*----------------------------------------------------------------------*
 |  Calculate the derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw::EvaluateDeriv(
    double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0)
  {
    FOUR_C_THROW("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  gap = -gap;
  return (-params_->GetA() * 1. /
          ((gap - params_->GetOffset() - params_->GetB()) *
              (gap - params_->GetOffset() - params_->GetB())));
}

FOUR_C_NAMESPACE_CLOSE
