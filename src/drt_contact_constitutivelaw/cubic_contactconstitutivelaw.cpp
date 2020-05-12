/*----------------------------------------------------------------------------*/
/*! \file
 *
\brief implements a Cubic contact constitutive law

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/


#include "cubic_contactconstitutivelaw.H"

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams::CubicConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->GetDouble("A")),
      b_(container->GetDouble("B")),
      c_(container->GetDouble("C")),
      d_(container->GetDouble("D"))
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams::CreateConstitutiveLaw()
{
  return Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::CubicConstitutiveLaw(this));
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::CubicConstitutiveLaw::CubicConstitutiveLaw(
    CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params)
    : params_(params)
{
}
/*----------------------------------------------------------------------*
 |  Evaluate Contact Constitutive Law
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::CubicConstitutiveLaw::Evaluate(double gap)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  double result = 1.0;
  gap = -gap;
  result = -1.0;
  result *= params_->GetA() * (gap - params_->GetOffset()) * (gap - params_->GetOffset()) *
                (gap - params_->GetOffset()) +
            params_->GetB() * (gap - params_->GetOffset()) * (gap - params_->GetOffset()) +
            params_->GetC() * (gap - params_->GetOffset()) + params_->GetD();
  if (result > 0)
    dserror(
        "The constitutive function you are using seems to be positive, even though the gap is "
        "negative. Please check your coefficients!");
  return result;
}  // end of LinearConstitutiveLaw evaluate
/*----------------------------------------------------------------------*
 |  Evaluate derivative of the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::CubicConstitutiveLaw::EvaluateDeriv(double gap)
{
  if (gap + params_->GetOffset() > 0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  gap = -gap;
  return 3 * params_->GetA() * (gap - params_->GetOffset()) * (gap - params_->GetOffset()) +
         2 * params_->GetB() * (gap - params_->GetOffset()) + params_->GetC();
}
