/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a linear contact constitutive law
\level 3

*----------------------------------------------------------------------*/


#include "baci_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"

#include "baci_global_data.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <vector>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::LinearConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(*container->Get<double>("A")),
      b_(*container->Get<double>("B"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::CreateConstitutiveLaw()
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
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::Evaluate(double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  return params_->GetA() * (gap + params_->GetOffset()) + params_->GetB();
}  // end of linear_coconstlaw evaluate
/*----------------------------------------------------------------------*
 |  Calculate the derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::EvaluateDeriv(
    double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes.");
  }
  return params_->GetA();
}

BACI_NAMESPACE_CLOSE
