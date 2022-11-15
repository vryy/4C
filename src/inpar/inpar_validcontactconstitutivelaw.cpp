/*! \file

\brief Setup of the list of valid contact constitutive laws for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "inpar_validcontactconstitutivelaw.H"

#include "contact_constitutivelaw_constitutivelaw_definition.H"
#include "contact_constitutivelaw_contactconstitutivelaw_parameter.H"
#include "lib_colors.H"
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
            "Broken rational law", INPAR::CONTACT::ConstitutiveLawType::colaw_brokenrational));

    AddNamedReal(m, "A", "scaling factor");
    AddNamedReal(m, "B", "asymptote");
    AddNamedReal(m, "C", "y intercept");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }
  // power law function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_power", "Power law", INPAR::CONTACT::ConstitutiveLawType::colaw_power));

    AddNamedReal(m, "A", "scaling factor");
    AddNamedReal(m, "B", "power coefficient");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // cubic function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_cubic",
            "Cubic function", INPAR::CONTACT::ConstitutiveLawType::colaw_cubic));

    AddNamedReal(m, "A", "A");
    AddNamedReal(m, "B", "B");
    AddNamedReal(m, "C", "C");
    AddNamedReal(m, "D", "D");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // linear function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_linear",
            "Linear function", INPAR::CONTACT::ConstitutiveLawType::colaw_linear));

    AddNamedReal(m, "A", "slope");
    AddNamedReal(m, "B", "y intercept");
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // mirco function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_mirco",
            "Mirco function", INPAR::CONTACT::ConstitutiveLawType::colaw_mirco));

    AddNamedReal(m, "FirstMatID", "First material ID");
    AddNamedReal(m, "SecondMatID", "Second material ID");
    AddNamedReal(m, "LateralLength", "length of lateral side of the BEM patch");
    AddNamedReal(m, "Resolution", "resolution of the surface");
    AddNamedReal(m, "InitialTopologyStdDeviation",
        "Initial Standard deviation for the random-midpoint generator", 20, true);
    AddNamedReal(m, "HurstExponent", "Hurst exponent of the surface", 0.7, true);
    AddNamedReal(m, "RandomTopologyFlag", "Use random midpoint generator flag", 1, true);
    AddNamedReal(m, "RandomSeedFlag", "Random seed flag", 0, true);
    AddNamedReal(m, "RandomGeneratorSeed", "Use random seed to reproduce results", 95, true);
    AddNamedReal(m, "Tolerance", "Tolerance for the convergence of force", 0.01, true);
    AddNamedReal(m, "MaxIteration", "Maximum iteration of NNLS", 1000, true);
    AddNamedReal(m, "WarmStartingFlag", "Warm-starting flag, solution accelerator", 1, true);
    AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);
    AddNamedReal(m, "FiniteDifferenceFraction",
        "Fraction of pertubation difference compared to the actual gap", 0.001, true);
    AddNamedReal(m, "ActiveGapTolerance", "Minimum gap to consider a node as active", 1e-6, true);
    AddNamedString(
        m, "TopologyFilePath", "Absolute path to file with micro-topology data", "", true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // deliver
  return vm;
}
