/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for contact constitutive laws

\level 1

\maintainer Nora Hagmeyer

*/
/*----------------------------------------------------------------------*/


#include "coconstitutivelaw.H"

#include "../drt_lib/drt_globalproblem.H"
#include "brokenrational_coconstlaw.H"
#include "cubic_coconstlaw.H"
#include "linear_coconstlaw.H"
#include "power_coconstlaw.H"
#include "coconstlaw_parameter.H"



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
