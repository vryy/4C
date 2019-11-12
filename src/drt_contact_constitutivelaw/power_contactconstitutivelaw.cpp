/*----------------------------------------------------------------------------*/
/*! \file
 *
\brief implements a power law as contact constitutive law

\level 3


\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/


#include "power_contactconstitutivelaw.H"

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../drt_lib/drt_globalproblem.H"
#include <math.h>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams::PowerConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->GetDouble("A")),
      b_(container->GetDouble("B"))
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
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawType
    CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawType::instance_;
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // matid
  int CONSTITUTIVELAWid = -1;
  if (params_ != NULL) CONSTITUTIVELAWid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, CONSTITUTIVELAWid);
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  int CONSTITUTIVELAWid;
  ExtractfromPack(position, data, CONSTITUTIVELAWid);
}
/*----------------------------------------------------------------------*
 |  Evaluate the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::Evaluate(double gap)
{
  if (gap + params_->GetOffset() >= 0)
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
double CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw::EvaluateDeriv(double gap)
{
  if (gap + params_->GetOffset() >= 0.0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  gap = -gap;
  return params_->GetA() * params_->GetB() * pow(gap - params_->GetOffset(), params_->GetB() - 1);
}
