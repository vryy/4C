/*----------------------------------------------------------------------------*/
/*! \file
 *
\brief implements a power law as contact constitutive law

\level 3


*/
/*----------------------------------------------------------------------*/


#include "baci_contact_constitutivelaw_power_contactconstitutivelaw.hpp"

#include "baci_global_data.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <math.h>

#include <vector>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams::PowerConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(*container->Get<double>("A")),
      b_(*container->Get<double>("B"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams::CreateConstitutiveLaw()
{
  return Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw(this));
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::PowerConstitutiveLaw(
    CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params)
    : params_(params)
{
}
/*----------------------------------------------------------------------*
 |  Evaluate the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::Evaluate(double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  double result = 1;
  gap *= -1;
  result = -1;
  result *= (params_->GetA() * pow(gap - params_->GetOffset(), params_->GetB()));
  if (result > 0)
    dserror(
        "The constitutive function you are using seems to be positive, even though the gap is "
        "negative. Please check your coefficients!");
  return result;
}  // end of Power_coconstlaw evaluate
/*----------------------------------------------------------------------*
 |  Calculate the derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::EvaluateDeriv(
    double gap, CONTACT::Node* cnode)
{
  if (gap + params_->GetOffset() > 0.0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  gap = -gap;
  return params_->GetA() * params_->GetB() * pow(gap - params_->GetOffset(), params_->GetB() - 1);
}

BACI_NAMESPACE_CLOSE
