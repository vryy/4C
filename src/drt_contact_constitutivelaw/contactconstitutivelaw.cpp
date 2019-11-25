/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for contact constitutive laws

\level 1

\maintainer Nora Hagmeyer

*/
/*----------------------------------------------------------------------*/


#include "contactconstitutivelaw.H"

#include "../drt_lib/drt_globalproblem.H"
#include "brokenrational_contactconstitutivelaw.H"
#include "contact_constitutivelaw_bundle.H"
#include "contactconstitutivelaw_parameter.H"
#include "cubic_contactconstitutivelaw.H"
#include "linear_contactconstitutivelaw.H"
#include "power_contactconstitutivelaw.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(const int id)
{
  const int probinst = DRT::Problem::Instance()->ContactConstitutiveLaws()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->ContactConstitutiveLaws() == Teuchos::null)
    dserror("Cannot work out problem instance!");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->ContactConstitutiveLaws()->Num() == 0)
    dserror("Cannot find any contact constitutive law!");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> coconstlawdata =
      DRT::Problem::Instance(probinst)->ContactConstitutiveLaws()->ById(id);
  return CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(coconstlawdata);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> coconstlawdata)
{
  switch (coconstlawdata->Type())
  {
    case INPAR::CONTACT::ConstitutiveLawType::colaw_cubic:
    {
      CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams(coconstlawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_brokenrational:
    {
      CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams(coconstlawdata);
      return params->CreateConstitutiveLaw();
    }

    case INPAR::CONTACT::ConstitutiveLawType::colaw_linear:
    {
      CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams(coconstlawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_power:
    {
      CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams(coconstlawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_none:
    {
      dserror("No contact constitutive law found\n");
      break;
    }
    default:
      dserror("unknown type of contact constitutive law %d\n", coconstlawdata->Type());
      break;
  }

  return Teuchos::null;
}
