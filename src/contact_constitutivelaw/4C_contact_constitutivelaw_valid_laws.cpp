// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_constitutivelaw_valid_laws.hpp"

#include "4C_contact_constitutivelaw_bundle.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::InputSpec CONTACT::CONSTITUTIVELAW::valid_contact_constitutive_laws()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  std::vector<Inpar::CONTACT::ConstitutiveLawType> group_index_to_type;

  /*----------------------------------------------------------------------*/
  // broken rational function
  {
    specs.emplace_back(group("CoConstLaw_brokenrational",
        {
            entry<double>("A", {.description = "scaling factor"}),
            entry<double>("B", {.description = "asymptote"}),
            entry<double>("C", {.description = "y intercept"}),
            entry<double>(
                "Offset", {.description = "offset for contact to start", .default_value = 0.0}),
        },
        {
            .description = "Broken rational law",
        }));
    group_index_to_type.push_back(Inpar::CONTACT::ConstitutiveLawType::colaw_brokenrational);
  }
  // power law function
  {
    specs.emplace_back(group("CoConstLaw_power",
        {
            entry<double>("A", {.description = "scaling factor"}),
            entry<double>("B", {.description = "power coefficient"}),
            entry<double>(
                "Offset", {.description = "offset for contact to start", .default_value = 0.0}),
        },
        {
            .description = "Power law",
        }));
    group_index_to_type.push_back(Inpar::CONTACT::ConstitutiveLawType::colaw_power);
  }

  // cubic function
  {
    specs.emplace_back(group("CoConstLaw_cubic",
        {
            entry<double>("A", {.description = "A"}),
            entry<double>("B", {.description = "B"}),
            entry<double>("C", {.description = "C"}),
            entry<double>("D", {.description = "D"}),
            entry<double>(
                "Offset", {.description = "offset for contact to start", .default_value = 0.0}),
        },
        {
            .description = "Cubic function",
        }));
    group_index_to_type.push_back(Inpar::CONTACT::ConstitutiveLawType::colaw_cubic);
  }

  // linear function
  {
    specs.emplace_back(group("CoConstLaw_linear",
        {
            entry<double>("A", {.description = "slope"}),
            entry<double>("B", {.description = "y intercept"}),
            entry<double>(
                "Offset", {.description = "offset for contact to start", .default_value = 0.0}),
        },
        {
            .description = "Linear function",
        }));
    group_index_to_type.push_back(Inpar::CONTACT::ConstitutiveLawType::colaw_linear);
  }

  // mirco function
  {
    specs.emplace_back(group("CoConstLaw_mirco",
        {
            entry<int>("FirstMatID", {.description = "First material ID"}),
            entry<int>("SecondMatID", {.description = "Second material ID"}),
            entry<double>(
                "LateralLength", {.description = "length of lateral side of the BEM patch"}),
            entry<int>("Resolution", {.description = "resolution of the surface"}),
            entry<bool>("PressureGreenFunFlag",
                {.description = "Use pressure-based Green function instead of a point-force-based",
                    .default_value = true}),
            entry<int>("InitialTopologyStdDeviationFunct",
                {.description = "Function id for Initial Standard deviation for the "
                                "random-midpoint generator",
                    .default_value = 100}),
            entry<int>(
                "HurstExponentFunct", {.description = "Function for Hurst exponent of the surface",
                                          .default_value = 100}),
            entry<bool>("RandomTopologyFlag",
                {.description = "Use random midpoint generator flag", .default_value = true}),
            entry<bool>(
                "RandomSeedFlag", {.description = "Random seed flag", .default_value = false}),
            entry<int>("RandomGeneratorSeed",
                {.description = "Use random seed to reproduce results", .default_value = 95}),
            entry<double>("Tolerance",
                {.description = "Tolerance for the convergence of force", .default_value = 0.01}),
            entry<int>("MaxIteration",
                {.description = "Maximum iteration of NNLS", .default_value = 1000}),
            entry<bool>("WarmStartingFlag",
                {.description = "Warm-starting flag, solution accelerator", .default_value = true}),
            entry<double>(
                "Offset", {.description = "offset for contact to start", .default_value = 0.0}),
            entry<double>("FiniteDifferenceFraction",
                {.description = "Fraction of perturbation difference compared to the actual gap",
                    .default_value = 0.001}),
            entry<double>("ActiveGapTolerance",
                {.description = "Minimum gap to consider a node as active", .default_value = 1e-6}),
            entry<std::string>("TopologyFilePath",
                {.description = "Path to file with micro-topology data", .default_value = ""}),
        },
        {
            .description = "Mirco function",
        }));
    group_index_to_type.push_back(Inpar::CONTACT::ConstitutiveLawType::colaw_mirco);
  }

  auto valid_law = group({
      entry<int>("LAW"),
      one_of(specs, store_index_as("LAW_TYPE", group_index_to_type)),
  });

  return valid_law;
}

void CONTACT::CONSTITUTIVELAW::create_contact_constitutive_law_from_input(
    const Core::IO::InputParameterContainer& container)
{
  const auto id = container.get<int>("LAW");

  Global::Problem::instance()->contact_constitutive_laws()->insert(id, container);
}

FOUR_C_NAMESPACE_CLOSE
