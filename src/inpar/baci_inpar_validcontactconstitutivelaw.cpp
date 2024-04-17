/*! \file

\brief Setup of the list of valid contact constitutive laws for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_validcontactconstitutivelaw.hpp"

#include "baci_contact_constitutivelaw_constitutivelaw_definition.hpp"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::PrintEmptyContactConstitutiveLawDefinitions(std::ostream& stream,
    std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& contactconstitutivlawlist)
{
  const std::string sectionname = "Contact Constitutive Law";
  const unsigned l = sectionname.length();
  stream << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << sectionname << '\n';

  for (unsigned i = 0; i < contactconstitutivlawlist.size(); ++i)
  {
    contactconstitutivlawlist[i]->Print(stream, nullptr);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintContactConstitutiveLawDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> coconstlawlist =
      INPUT::ValidContactConstitutiveLaws();
  INPUT::PrintEmptyContactConstitutiveLawDefinitions(std::cout, *coconstlawlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
INPUT::ValidContactConstitutiveLaws()
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

    AddNamedInt(m, "FirstMatID", "First material ID");
    AddNamedInt(m, "SecondMatID", "Second material ID");
    AddNamedReal(m, "LateralLength", "length of lateral side of the BEM patch");
    AddNamedInt(m, "Resolution", "resolution of the surface");
    AddNamedBool(m, "PressureGreenFunFlag",
        "Use pressure-based Green function instead of a point-force-based", true, true);
    AddNamedReal(m, "InitialTopologyStdDeviation",
        "Initial Standard deviation for the random-midpoint generator", 20, true);
    AddNamedReal(m, "HurstExponent", "Hurst exponent of the surface", 0.7, true);
    AddNamedBool(m, "RandomTopologyFlag", "Use random midpoint generator flag", true, true);
    AddNamedBool(m, "RandomSeedFlag", "Random seed flag", false, true);
    AddNamedInt(m, "RandomGeneratorSeed", "Use random seed to reproduce results", 95, true);
    AddNamedReal(m, "Tolerance", "Tolerance for the convergence of force", 0.01, true);
    AddNamedInt(m, "MaxIteration", "Maximum iteration of NNLS", 1000, true);
    AddNamedBool(m, "WarmStartingFlag", "Warm-starting flag, solution accelerator", true, true);
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

FOUR_C_NAMESPACE_CLOSE
