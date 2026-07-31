// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_loma_input.hpp"

#include "4C_io_input_spec_builders.hpp"
FOUR_C_NAMESPACE_OPEN

Core::IO::InputSpec LowMach::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("LOMA CONTROL",
      {

          parameter<bool>(
              "MONOLITHIC", {.description = "monolithic solver", .default_value = false}),
          parameter<int>(
              "NUMSTEP", {.description = "Total number of time steps", .default_value = 24}),

          parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}),
          parameter<double>(
              "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}),
          parameter<int>(
              "ITEMAX", {.description = "Maximum number of outer iterations", .default_value = 10}),
          parameter<int>(
              "ITEMAX_BEFORE_SAMPLING", {.description = "Maximum number of outer iterations before "
                                                        "sampling (for turbulent flows only)",
                                            .default_value = 1}),
          parameter<double>(
              "CONVTOL", {.description = "Tolerance for convergence check", .default_value = 1e-6}),
          parameter<int>("RESULTSEVERY",
              {.description = "Increment for writing solution", .default_value = 1}),
          parameter<int>(
              "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}),


          deprecated_selection<std::string>("CONSTHERMPRESS", {"No_energy", "No_mass", "Yes"},
              {.description = "treatment of thermodynamic pressure in time",
                  .default_value = "Yes"}),

          parameter<bool>("SGS_MATERIAL_UPDATE",
              {.description = "update material by adding subgrid-scale scalar field",
                  .default_value = false}),

          // number of linear solver used for LOMA solver
          parameter<int>(
              "LINEAR_SOLVER", {.description = "number of linear solver used for LOMA problem",
                                   .default_value = -1})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE