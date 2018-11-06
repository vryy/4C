/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a linear contact constitutive law
\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "linear_coconstlaw.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams::LinearConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      a_(container->GetDouble("A")),
      b_(container->GetDouble("B"))
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
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawType
    CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawType::instance_;
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::Pack(DRT::PackBuffer& data) const
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
void CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // matid
  int CONSTITUTIVELAWid = -1;
  ExtractfromPack(position, data, CONSTITUTIVELAWid);
}
/*----------------------------------------------------------------------*
 |  Evaluate the contact constitutive law|
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::Evaluate(double gap)
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
double CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw::EvaluateDeriv(double gap)
{
  if (gap + params_->GetOffset() > 0)
  {
    printf("You should not be here. The Evaluate function is only tested for active nodes.");
  }
  return params_->GetA();
}
