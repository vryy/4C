// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_validcontactconstitutivelaw.hpp"

#include "4C_contact_constitutivelaw_constitutivelaw_definition.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_io_input_file_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Input::print_empty_contact_constitutive_law_definitions(std::ostream& stream,
    std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>&
        contactconstitutivlawlist)
{
  const std::string sectionname = "Contact Constitutive Law";
  Core::IO::InputFileUtils::print_section_header(stream, sectionname);

  for (unsigned i = 0; i < contactconstitutivlawlist.size(); ++i)
  {
    contactconstitutivlawlist[i]->print(stream, nullptr);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void print_contact_constitutive_law_dat_header()
{
  std::shared_ptr<std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
      coconstlawlist = Input::valid_contact_constitutive_laws();
  Input::print_empty_contact_constitutive_law_definitions(std::cout, *coconstlawlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
Input::valid_contact_constitutive_laws()
{
  // a list containing all valid contact constitutivelaw definitions
  std::shared_ptr<std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>> vm =
      std::make_shared<std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>>();

  // convenience
  std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist = *vm;

  /*----------------------------------------------------------------------*/
  // broken rational function
  {
    std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        std::make_shared<CONTACT::CONSTITUTIVELAW::LawDefinition>("CoConstLaw_brokenrational",
            "Broken rational law", Inpar::CONTACT::ConstitutiveLawType::colaw_brokenrational);

    add_named_real(m, "A", "scaling factor");
    add_named_real(m, "B", "asymptote");
    add_named_real(m, "C", "y intercept");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    append_co_const_law_component_definition(coconstlawlist, m);
  }
  // power law function
  {
    std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        std::make_shared<CONTACT::CONSTITUTIVELAW::LawDefinition>(
            "CoConstLaw_power", "Power law", Inpar::CONTACT::ConstitutiveLawType::colaw_power);

    add_named_real(m, "A", "scaling factor");
    add_named_real(m, "B", "power coefficient");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    append_co_const_law_component_definition(coconstlawlist, m);
  }

  // cubic function
  {
    std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        std::make_shared<CONTACT::CONSTITUTIVELAW::LawDefinition>(
            "CoConstLaw_cubic", "Cubic function", Inpar::CONTACT::ConstitutiveLawType::colaw_cubic);

    add_named_real(m, "A", "A");
    add_named_real(m, "B", "B");
    add_named_real(m, "C", "C");
    add_named_real(m, "D", "D");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    append_co_const_law_component_definition(coconstlawlist, m);
  }

  // linear function
  {
    std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        std::make_shared<CONTACT::CONSTITUTIVELAW::LawDefinition>("CoConstLaw_linear",
            "Linear function", Inpar::CONTACT::ConstitutiveLawType::colaw_linear);

    add_named_real(m, "A", "slope");
    add_named_real(m, "B", "y intercept");
    add_named_real(m, "Offset", "offset for contact to start", 0.0, true);

    append_co_const_law_component_definition(coconstlawlist, m);
  }

  // mirco function
  {
    std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition> m =
        std::make_shared<CONTACT::CONSTITUTIVELAW::LawDefinition>(
            "CoConstLaw_mirco", "Mirco function", Inpar::CONTACT::ConstitutiveLawType::colaw_mirco);

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

    append_co_const_law_component_definition(coconstlawlist, m);
  }

  // deliver
  return vm;
}

FOUR_C_NAMESPACE_CLOSE
