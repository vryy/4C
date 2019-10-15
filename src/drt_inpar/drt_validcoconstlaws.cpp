/*! \file
\maintainer Nora Hagmeyer

\brief Setup of the list of valid contact constitutive laws for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "drt_validcoconstlaws.H"
#include "../drt_contact_constitutivelaw/constitutivelaw_definition.H"
#include "inpar_material.H"
#include "../drt_lib/drt_colors.H"
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyContactConstitutiveLawDefinitions(std::ostream& stream,
    std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist, bool color)
{
  std::string blue2light = "";
  std::string bluelight = "";
  std::string redlight = "";
  std::string yellowlight = "";
  std::string greenlight = "";
  std::string magentalight = "";
  std::string endcolor = "";

  if (color)
  {
    blue2light = BLUE2_LIGHT;
    bluelight = BLUE_LIGHT;
    redlight = RED_LIGHT;
    yellowlight = YELLOW_LIGHT;
    greenlight = GREEN_LIGHT;
    magentalight = MAGENTA_LIGHT;
    endcolor = END_COLOR;
  }

  const std::string sectionname = "Contact Constitutive Law";
  const unsigned l = sectionname.length();
  stream << redlight << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << greenlight << sectionname << endcolor << '\n';

  for (unsigned i = 0; i < coconstlawlist.size(); ++i)
  {
    coconstlawlist[i]->Print(stream, NULL, color);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintContactConstitutiveLawDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> coconstlawlist =
      DRT::INPUT::ValidContactConstitutiveLaws();
  DRT::INPUT::PrintEmptyContactConstitutiveLawDefinitions(std::cout, *coconstlawlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
DRT::INPUT::ValidContactConstitutiveLaws()
{
  // a list containing all valid contact constitutivelaw definitions
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> vm =
      Teuchos::rcp(new std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>());

  // convenience
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist = *vm;


  /*----------------------------------------------------------------------*/
  // broken rational function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_brokenrational",
            "Broken rational law", INPAR::CONTACT::colaw_brokenrational));

    AddNamedReal(m, "A", "scaling factor");
    AddNamedReal(m, "B", "asymptote");
    AddNamedReal(m, "C", "y intercept");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawDefinition(coconstlawlist, m);
  }
  // power law function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_power", "Power law", INPAR::CONTACT::colaw_power));

    AddNamedReal(m, "A", "scaling factor");
    AddNamedReal(m, "B", "power coefficient");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawDefinition(coconstlawlist, m);
  }

  // cubic function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_cubic", "Cubic function", INPAR::CONTACT::colaw_cubic));

    AddNamedReal(m, "A", "A");
    AddNamedReal(m, "B", "B");
    AddNamedReal(m, "C", "C");
    AddNamedReal(m, "D", "D");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawDefinition(coconstlawlist, m);
  }

  // linear function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_linear", "Linear function", INPAR::CONTACT::colaw_linear));

    AddNamedReal(m, "A", "slope");
    AddNamedReal(m, "B", "y intercept");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawDefinition(coconstlawlist, m);
  }

  // deliver
  return vm;
}
