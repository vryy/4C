/*! \file

\brief Setup of the list of valid contact constitutive laws for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_validcontactconstitutivelaw.hpp"

#include "4C_contact_constitutivelaw_constitutivelaw_definition.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_io_dat_file_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Input::PrintEmptyContactConstitutiveLawDefinitions(std::ostream& stream,
    std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& contactconstitutivlawlist)
{
  const std::string sectionname = "Contact Constitutive Law";
  Core::IO::DatFileUtils::print_section_header(stream, sectionname);

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
      Input::ValidContactConstitutiveLaws();
  Input::PrintEmptyContactConstitutiveLawDefinitions(std::cout, *coconstlawlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
Input::ValidContactConstitutiveLaws()
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
            "Broken rational law", Inpar::CONTACT::ConstitutiveLawType::colaw_brokenrational));

    add_named_real(m, "A", "scaling factor");
    add_named_real(m, "B", "asymptote");
    add_named_real(m, "C", "y intercept");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }
  // power law function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition(
            "CoConstLaw_power", "Power law", Inpar::CONTACT::ConstitutiveLawType::colaw_power));

    add_named_real(m, "A", "scaling factor");
    add_named_real(m, "B", "power coefficient");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // cubic function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_cubic",
            "Cubic function", Inpar::CONTACT::ConstitutiveLawType::colaw_cubic));

    add_named_real(m, "A", "A");
    add_named_real(m, "B", "B");
    add_named_real(m, "C", "C");
    add_named_real(m, "D", "D");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // linear function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_linear",
            "Linear function", Inpar::CONTACT::ConstitutiveLawType::colaw_linear));

    add_named_real(m, "A", "slope");
    add_named_real(m, "B", "y intercept");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // mirco function
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::LawDefinition("CoConstLaw_mirco",
            "Mirco function", Inpar::CONTACT::ConstitutiveLawType::colaw_mirco));

    add_named_int(m, "FirstMatID", "First material ID");
    add_named_int(m, "SecondMatID", "Second material ID");
    add_named_real(m, "LateralLength", "length of lateral side of the BEM patch");
    add_named_int(m, "Resolution", "resolution of the surface");
    add_named_bool(m, "PressureGreenFunFlag",
        "Use pressure-based Green function instead of a point-force-based", true, true);
    add_named_int(m, "InitialTopologyStdDeviationFunct",
        "Function id for Initial Standard deviation for the random-midpoint generator", 100, true);
    add_named_int(m, "HurstExponentFunct", "Function for Hurst exponent of the surface", 100, true);
    add_named_bool(m, "RandomTopologyFlag", "Use random midpoint generator flag", true, true);
    add_named_bool(m, "RandomSeedFlag", "Random seed flag", false, true);
    add_named_int(m, "RandomGeneratorSeed", "Use random seed to reproduce results", 95, true);
    add_named_real(m, "Tolerance", "Tolerance for the convergence of force", 0.01, true);
    add_named_int(m, "MaxIteration", "Maximum iteration of NNLS", 1000, true);
    add_named_bool(m, "WarmStartingFlag", "Warm-starting flag, solution accelerator", true, true);
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);
    add_named_real(m, "FiniteDifferenceFraction",
        "Fraction of pertubation difference compared to the actual gap", 0.001, true);
    add_named_real(m, "ActiveGapTolerance", "Minimum gap to consider a node as active", 1e-6, true);
    add_named_string(
        m, "TopologyFilePath", "Absolute path to file with micro-topology data", "", true);

    AppendCoConstLawComponentDefinition(coconstlawlist, m);
  }

  // deliver
  return vm;
}

FOUR_C_NAMESPACE_CLOSE
