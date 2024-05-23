/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for contact constitutive laws

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_contact_constitutivelaw_contactconstitutivelaw.hpp"

#include "4C_contact_constitutivelaw_brokenrational_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_bundle.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_contact_constitutivelaw_cubic_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_linear_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_power_contactconstitutivelaw.hpp"
#include "4C_global_data.hpp"

#ifdef FOUR_C_WITH_MIRCO
#include "4C_contact_constitutivelaw_mirco_contactconstitutivelaw.hpp"
#endif

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(const int id)
{
  const int probinst =
      GLOBAL::Problem::Instance()->contact_constitutive_laws()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->contact_constitutive_laws() == Teuchos::null)
    FOUR_C_THROW("Cannot work out problem instance!");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->contact_constitutive_laws()->Num() == 0)
    FOUR_C_THROW("Cannot find any contact constitutive law!");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> coconstlawdata =
      GLOBAL::Problem::Instance(probinst)->contact_constitutive_laws()->ById(id);
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
      return params->create_constitutive_law();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_brokenrational:
    {
      CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::BrokenRationalConstitutiveLawParams(
              contactconstitutivelawdata);
      return params->create_constitutive_law();
    }

    case INPAR::CONTACT::ConstitutiveLawType::colaw_linear:
    {
      CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams(contactconstitutivelawdata);
      return params->create_constitutive_law();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_power:
    {
      CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams(contactconstitutivelawdata);
      return params->create_constitutive_law();
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_mirco:
    {
#ifdef FOUR_C_WITH_MIRCO
      CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams* params =
          new CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams(contactconstitutivelawdata);
      return params->create_constitutive_law();
#else
      FOUR_C_THROW(
          "You are trying to use MIRCO contact consitutive law with FOUR_C_WITH_MIRCO flag turned "
          "off. Please enable this flag and build 4C again");
#endif
    }
    case INPAR::CONTACT::ConstitutiveLawType::colaw_none:
    {
      FOUR_C_THROW("No contact constitutive law found\n");
      break;
    }
    default:
      FOUR_C_THROW(
          "unknown type of contact constitutive law %d\n", contactconstitutivelawdata->Type());
      break;
  }

  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
