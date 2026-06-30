// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_io.hpp"

#include "4C_io_input_spec_builders.hpp"
#include "4C_io_pstream.hpp"
#include "4C_structure_new_input.hpp"
#include "4C_thermo_input.hpp"
FOUR_C_NAMESPACE_OPEN

std::vector<Core::IO::InputSpec> Inpar::IO::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("IO",
      {

          parameter<bool>("OUTPUT_GMSH", {.description = "", .default_value = false}),
          parameter<bool>("OUTPUT_ROT", {.description = "", .default_value = false}),

          parameter<bool>("OUTPUT_SPRING",
              {.description = "Flag indicating whether quantities related to the spring dashpot "
                              "conditions should be written to output.",
                  .default_value = false}),
          parameter<bool>("WRITE_TIMINGS",
              {.description =
                      "Flag indicating whether timings of the time monitor should be written to a "
                      "yaml file. The timing information is written to "
                      "<outputprefix>-timings.yaml.",
                  .default_value = true}),
          parameter<bool>("OUTPUT_BIN",
              {.description = "Do you want to have binary output?", .default_value = true}),

          parameter<bool>("ELEMENT_MAT_ID",
              {.description = "Output of the material id of each element", .default_value = false}),

          // Structural output
          parameter<bool>(
              "STRUCT_ELE", {.description = "Output of element properties", .default_value = true}),
          parameter<bool>(
              "STRUCT_DISP", {.description = "Output of displacements", .default_value = true}),
          deprecated_selection<Solid::StressType>("STRUCT_STRESS",
              {
                  {"No", Solid::stress_none},
                  {"no", Solid::stress_none},
                  {"NO", Solid::stress_none},
                  {"Yes", Solid::stress_2pk},
                  {"yes", Solid::stress_2pk},
                  {"YES", Solid::stress_2pk},
                  {"Cauchy", Solid::stress_cauchy},
                  {"cauchy", Solid::stress_cauchy},
                  {"2PK", Solid::stress_2pk},
                  {"2pk", Solid::stress_2pk},
              },
              {.description = "Output of stress. This is only possible in the new vtk-based "
                              "output. Make sure to set INTERVAL_STEPS accordingly.",
                  .default_value = Solid::stress_none}),
          deprecated_selection<Solid::StrainType>("STRUCT_STRAIN",
              {
                  {"No", Solid::strain_none},
                  {"no", Solid::strain_none},
                  {"NO", Solid::strain_none},
                  {"Yes", Solid::strain_gl},
                  {"yes", Solid::strain_gl},
                  {"YES", Solid::strain_gl},
                  {"EA", Solid::strain_ea},
                  {"ea", Solid::strain_ea},
                  {"GL", Solid::strain_gl},
                  {"gl", Solid::strain_gl},
                  {"LOG", Solid::strain_log},
                  {"log", Solid::strain_log},
              },
              {.description = "Output of strains. This is only possible in the new vtk-based "
                              "output. Make sure to set INTERVAL_STEPS accordingly.",
                  .default_value = Solid::strain_none}),
          deprecated_selection<Solid::StrainType>("STRUCT_PLASTIC_STRAIN",
              {
                  {"No", Solid::strain_none},
                  {"no", Solid::strain_none},
                  {"NO", Solid::strain_none},
                  {"Yes", Solid::strain_gl},
                  {"yes", Solid::strain_gl},
                  {"YES", Solid::strain_gl},
                  {"EA", Solid::strain_ea},
                  {"ea", Solid::strain_ea},
                  {"GL", Solid::strain_gl},
                  {"gl", Solid::strain_gl},
              },
              {.description = "", .default_value = Solid::strain_none}),

          parameter<bool>("STRUCT_SURFACTANT", {.description = "", .default_value = false}),

          parameter<bool>("STRUCT_JACOBIAN_MATLAB", {.description = "", .default_value = false}),

          deprecated_selection<Solid::ConditionNumber>("STRUCT_CONDITION_NUMBER",
              {
                  {"gmres_estimate", Solid::ConditionNumber::gmres_estimate},
                  {"max_min_ev_ratio", Solid::ConditionNumber::max_min_ev_ratio},
                  {"one-norm", Solid::ConditionNumber::one_norm},
                  {"inf-norm", Solid::ConditionNumber::inf_norm},
                  {"none", Solid::ConditionNumber::none},
              },
              {.description =
                      "Compute the condition number of the structural system matrix and write "
                      "it to a text file.",
                  .default_value = Solid::ConditionNumber::none}),

          parameter<bool>("FLUID_STRESS", {.description = "", .default_value = false}),

          parameter<bool>("FLUID_WALL_SHEAR_STRESS", {.description = "", .default_value = false}),

          parameter<bool>("FLUID_ELEDATA_EVERY_STEP", {.description = "", .default_value = false}),

          parameter<bool>("FLUID_NODEDATA_FIRST_STEP", {.description = "", .default_value = false}),

          parameter<int>(
              "FILESTEPS", {.description = "Amount of timesteps written to a single result file",
                               .default_value = 1000}),
          parameter<int>(
              "STDOUTEVERY", {.description = "Print to screen every n step", .default_value = 1}),

          parameter<bool>(
              "WRITE_TO_SCREEN", {.description = "Write screen output", .default_value = true}),
          parameter<bool>("WRITE_TO_FILE",
              {.description = "Write the output into a file", .default_value = false}),

          parameter<bool>("WRITE_INITIAL_STATE",
              {.description = "Do you want to write output for initial state ?",
                  .default_value = true}),
          parameter<bool>("WRITE_FINAL_STATE",
              {.description = "Enforce to write restart data at the final "
                              "state regardless of the other output or restart intervals",
                  .default_value = false}),

          parameter<bool>("PREFIX_GROUP_ID",
              {.description = "Put a <GroupID>: in front of every line", .default_value = false}),
          parameter<int>("LIMIT_OUTP_TO_PROC",
              {.description = "Only the specified procs will write output", .default_value = -1}),
          deprecated_selection<Core::IO::Verbositylevel>("VERBOSITY",
              {
                  {"minimal", Core::IO::minimal},
                  {"Minimal", Core::IO::minimal},
                  {"standard", Core::IO::standard},
                  {"Standard", Core::IO::standard},
                  {"verbose", Core::IO::verbose},
                  {"Verbose", Core::IO::verbose},
                  {"debug", Core::IO::debug},
                  {"Debug", Core::IO::debug},
              },
              {.description = "", .default_value = Core::IO::verbose}),

          parameter<double>("RESTARTWALLTIMEINTERVAL",
              {.description = "Enforce restart after this walltime interval "
                              "(in seconds), smaller zero to disable",
                  .default_value = -1.0}),
          parameter<int>("RESTARTEVERY",
              {.description = "write restart every RESTARTEVERY steps", .default_value = -1}),
          parameter<bool>("PER_RANK_EVAL_TIME",
              {.description = "write per rank wall time for element evaluation to csv file",
                  .default_value = false})},
      {.required = false}));
  return specs;
}

FOUR_C_NAMESPACE_CLOSE