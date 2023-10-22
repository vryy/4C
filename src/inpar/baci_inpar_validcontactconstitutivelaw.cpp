/*! \file

\brief Setup of the list of valid contact constitutive laws for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_validcontactconstitutivelaw.H"

#include "baci_contact_constitutivelaw_constitutivelaw_definition.H"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.H"
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyContactConstitutiveLawDefinitions(std::ostream& stream,
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

    ::INPUT::AddNamedReal(m, "A", "scaling factor");
    ::INPUT::AddNamedReal(m, "B", "asymptote");
    ::INPUT::AddNamedReal(m, "C", "y intercept");
    ::INPUT::AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }
  // power law function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_power", "Power law", INPAR::CONTACT::ConstitutiveLawType::colaw_power));

    ::INPUT::AddNamedReal(m, "A", "scaling factor");
    ::INPUT::AddNamedReal(m, "B", "power coefficient");
    ::INPUT::AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // cubic function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_cubic",
            "Cubic function", INPAR::CONTACT::ConstitutiveLawType::colaw_cubic));

    ::INPUT::AddNamedReal(m, "A", "A");
    ::INPUT::AddNamedReal(m, "B", "B");
    ::INPUT::AddNamedReal(m, "C", "C");
    ::INPUT::AddNamedReal(m, "D", "D");
    ::INPUT::AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // linear function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_linear",
            "Linear function", INPAR::CONTACT::ConstitutiveLawType::colaw_linear));

    ::INPUT::AddNamedReal(m, "A", "slope");
    ::INPUT::AddNamedReal(m, "B", "y intercept");
    ::INPUT::AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // mirco function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_mirco",
            "Mirco function", INPAR::CONTACT::ConstitutiveLawType::colaw_mirco));

    ::INPUT::AddNamedReal(m, "FirstMatID", "First material ID");
    ::INPUT::AddNamedReal(m, "SecondMatID", "Second material ID");
    ::INPUT::AddNamedReal(m, "LateralLength", "length of lateral side of the BEM patch");
    ::INPUT::AddNamedReal(m, "Resolution", "resolution of the surface");
    ::INPUT::AddNamedReal(m, "InitialTopologyStdDeviation",
        "Initial Standard deviation for the random-midpoint generator", 20, true);
    ::INPUT::AddNamedReal(m, "HurstExponent", "Hurst exponent of the surface", 0.7, true);
    ::INPUT::AddNamedReal(m, "RandomTopologyFlag", "Use random midpoint generator flag", 1, true);
    ::INPUT::AddNamedReal(m, "RandomSeedFlag", "Random seed flag", 0, true);
    ::INPUT::AddNamedReal(
        m, "RandomGeneratorSeed", "Use random seed to reproduce results", 95, true);
    ::INPUT::AddNamedReal(m, "Tolerance", "Tolerance for the convergence of force", 0.01, true);
    ::INPUT::AddNamedReal(m, "MaxIteration", "Maximum iteration of NNLS", 1000, true);
    ::INPUT::AddNamedReal(
        m, "WarmStartingFlag", "Warm-starting flag, solution accelerator", 1, true);
    ::INPUT::AddNamedReal(m, "Offset", "offset for contact to start", 0.0, true);
    ::INPUT::AddNamedReal(m, "FiniteDifferenceFraction",
        "Fraction of pertubation difference compared to the actual gap", 0.001, true);
    ::INPUT::AddNamedReal(
        m, "ActiveGapTolerance", "Minimum gap to consider a node as active", 1e-6, true);
    ::INPUT::AddNamedString(
        m, "TopologyFilePath", "Absolute path to file with micro-topology data", "", true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // deliver
  return vm;
}
