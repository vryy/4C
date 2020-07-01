/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a contact constitutive law following a broken rational function

\level 1

*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "brokenrational_contactconstitutivelaw.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams::BrokenRationalConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->GetDouble("A")),
      b_(container->GetDouble("B")),
      c_(container->GetDouble("C"))
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
double CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw::Evaluate(double gap)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  double result = -1.0;
  gap = -gap;
  result *=
      (params_->GetA() * 1. / (gap - params_->GetOffset() - params_->GetB()) + params_->GetC());
  if (result > 0)
    dserror(
        "The constitutive function you are using seems to be positive, even though the gap is "
        "negative. Please check your coefficients!");
  return result;
}
/*----------------------------------------------------------------------*
 |  Calculate the derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLaw::EvaluateDeriv(double gap)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  gap = -gap;
  return (-params_->GetA() * 1. /
          ((gap - params_->GetOffset() - params_->GetB()) *
              (gap - params_->GetOffset() - params_->GetB())));
}
