/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for contact constitutive laws

\level 1


*/
/*----------------------------------------------------------------------*/


#include "baci_contact_constitutivelaw_contactconstitutivelaw.hpp"

#include "baci_contact_constitutivelaw_brokenrational_contactconstitutivelaw.hpp"
#include "baci_contact_constitutivelaw_bundle.hpp"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "baci_contact_constitutivelaw_cubic_contactconstitutivelaw.hpp"
#include "baci_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"
#include "baci_contact_constitutivelaw_power_contactconstitutivelaw.hpp"
#include "baci_global_data.hpp"

#ifdef BACI_WITH_MIRCO
#include "baci_contact_constitutivelaw_mirco_contactconstitutivelaw.hpp"
#endif

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(const int id)
{
  const int probinst = GLOBAL::Problem::Instance()->ContactConstitutiveLaws()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->ContactConstitutiveLaws() == Teuchos::null)
    dserror("Cannot work out problem instance!");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->ContactConstitutiveLaws()->Num() == 0)
    dserror("Cannot find any contact constitutive law!");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> coconstlawdata =
      GLOBAL::Problem::Instance(probinst)->ContactConstitutiveLaws()->ById(id);
  return CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(coconstlawdata);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> contactconstitutivelawdata)
{
  switch (contactconstitutivelawdata->Type())
  {
    case INPAR::CONTACT::ConstitutiveLawType::colaw_cubic:
    {
      CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams(contactconstitutivelawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_brokenrational:
    {
      CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams(
              contactconstitutivelawdata);
      return params->CreateConstitutiveLaw();
    }

    case INPAR::CONTACT::ConstitutiveLawType::colaw_linear:
    {
      CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams(contactconstitutivelawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_power:
    {
      CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams(contactconstitutivelawdata);
      return params->CreateConstitutiveLaw();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_mirco:
    {
#ifdef BACI_WITH_MIRCO
      CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams(contactconstitutivelawdata);
      return params->CreateConstitutiveLaw();
#else
      dserror(
          "You are trying to use MIRCO contact consitutive law with BACI_WITH_MIRCO flag turned "
          "off. Please enable this flag and build BACI again");
#endif
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_none:
    {
      dserror("No contact constitutive law found\n");
      break;
    }
    default:
      dserror("unknown type of contact constitutive law %d\n", contactconstitutivelawdata->Type());
      break;
  }

  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
