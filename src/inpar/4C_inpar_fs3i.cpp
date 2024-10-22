// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_fs3i.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::FS3I::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& fs3idyn =
      list.sublist("FS3I DYNAMIC", false, "control parameters for FS3I problems\n");

  Core::Utils::double_parameter("TIMESTEP", 0.1, "Time increment dt", &fs3idyn);
  Core::Utils::int_parameter("NUMSTEP", 20, "Total number of time steps", &fs3idyn);
  Core::Utils::double_parameter("MAXTIME", 1000.0, "Total simulation time", &fs3idyn);
  Core::Utils::int_parameter("RESULTSEVRY", 1, "Increment for writing solution", &fs3idyn);
  Core::Utils::int_parameter("RESTARTEVRY", 1, "Increment for writing restart", &fs3idyn);
  setStringToIntegralParameter<Inpar::ScaTra::SolverType>("SCATRA_SOLVERTYPE", "nonlinear",
      "type of scalar transport solver", tuple<std::string>("linear", "nonlinear"),
      tuple<Inpar::ScaTra::SolverType>(
          Inpar::ScaTra::solvertype_linear_incremental, Inpar::ScaTra::solvertype_nonlinear),
      &fs3idyn);
  Core::Utils::bool_parameter("INF_PERM", "yes", "Flag for infinite permeability", &fs3idyn);
  std::vector<std::string> consthermpress_valid_input = {"No_energy", "No_mass", "Yes"};
  Core::Utils::string_parameter("CONSTHERMPRESS", "Yes",
      "treatment of thermodynamic pressure in time", &fs3idyn, consthermpress_valid_input);

  // number of linear solver used for fs3i problems
  Core::Utils::int_parameter(
      "COUPLED_LINEAR_SOLVER", -1, "number of linear solver used for fs3i problem", &fs3idyn);
  Core::Utils::int_parameter(
      "LINEAR_SOLVER1", -1, "number of linear solver used for fluid problem", &fs3idyn);
  Core::Utils::int_parameter(
      "LINEAR_SOLVER2", -1, "number of linear solver used for structural problem", &fs3idyn);

  setStringToIntegralParameter<Inpar::ScaTra::ConvForm>("STRUCTSCAL_CONVFORM", "conservative",
      "form of convective term of structure scalar",
      tuple<std::string>("convective", "conservative"),
      tuple<Inpar::ScaTra::ConvForm>(
          Inpar::ScaTra::convform_convective, Inpar::ScaTra::convform_conservative),
      &fs3idyn);

  setStringToIntegralParameter<Inpar::ScaTra::InitialField>("STRUCTSCAL_INITIALFIELD", "zero_field",
      "Initial Field for structure scalar transport problem",
      tuple<std::string>("zero_field", "field_by_function"),
      tuple<Inpar::ScaTra::InitialField>(
          Inpar::ScaTra::initfield_zero_field, Inpar::ScaTra::initfield_field_by_function),
      &fs3idyn);

  Core::Utils::int_parameter("STRUCTSCAL_INITFUNCNO", -1,
      "function number for structure scalar transport initial field", &fs3idyn);

  // Type of coupling strategy between structure and structure-scalar field
  setStringToIntegralParameter<VolumeCoupling>("STRUCTSCAL_FIELDCOUPLING", "volume_matching",
      "Type of coupling strategy between structure and structure-scalar field",
      tuple<std::string>("volume_matching", "volume_nonmatching"),
      tuple<VolumeCoupling>(coupling_match, coupling_nonmatch), &fs3idyn);

  // Type of coupling strategy between fluid and fluid-scalar field
  setStringToIntegralParameter<VolumeCoupling>("FLUIDSCAL_FIELDCOUPLING", "volume_matching",
      "Type of coupling strategy between fluid and fluid-scalar field",
      tuple<std::string>("volume_matching", "volume_nonmatching"),
      tuple<VolumeCoupling>(coupling_match, coupling_nonmatch), &fs3idyn);

  // type of scalar transport
  setStringToIntegralParameter<Inpar::ScaTra::ImplType>("FLUIDSCAL_SCATRATYPE",
      "ConvectionDiffusion", "Type of scalar transport problem",
      tuple<std::string>("Undefined", "ConvectionDiffusion", "Loma", "Advanced_Reaction",
          "Chemotaxis", "Chemo_Reac"),
      tuple<Inpar::ScaTra::ImplType>(Inpar::ScaTra::impltype_undefined, Inpar::ScaTra::impltype_std,
          Inpar::ScaTra::impltype_loma, Inpar::ScaTra::impltype_advreac,
          Inpar::ScaTra::impltype_chemo, Inpar::ScaTra::impltype_chemoreac),
      &fs3idyn);

  // Restart from FSI instead of FS3I
  Core::Utils::bool_parameter("RESTART_FROM_PART_FSI", "No",
      "restart from partitioned fsi problem (e.g. from prestress calculations) instead of fs3i",
      &fs3idyn);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned FS3I */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fs3idynpart = fs3idyn.sublist(
      "PARTITIONED", false, "partioned fluid-structure-scalar-scalar interaction control section");

  // Coupling strategy for partitioned FS3I
  setStringToIntegralParameter<SolutionSchemeOverFields>("COUPALGO", "fs3i_IterStagg",
      "Coupling strategies for FS3I solvers",
      tuple<std::string>("fs3i_SequStagg", "fs3i_IterStagg"),
      tuple<SolutionSchemeOverFields>(fs3i_SequStagg, fs3i_IterStagg), &fs3idynpart);

  // convergence tolerance of outer iteration loop
  Core::Utils::double_parameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration within partitioned FS3I", &fs3idynpart);

  Core::Utils::int_parameter("ITEMAX", 10, "Maximum number of outer iterations", &fs3idynpart);

  /*----------------------------------------------------------------------  */
  /* parameters for stabilization of the structure-scalar field             */
  /*----------------------------------------------------------------------  */
  Teuchos::ParameterList& fs3idynstructscalstab = fs3idyn.sublist("STRUCTURE SCALAR STABILIZATION",
      false, "parameters for stabilization of the structure-scalar field");

  Teuchos::ParameterList& scatradyn = list.sublist(
      "SCALAR TRANSPORT DYNAMIC", true, "control parameters for scalar transport problems\n");
  fs3idynstructscalstab = scatradyn.sublist("STABILIZATION", true,
      "control parameters for the stabilization of scalar transport problems");
}

FOUR_C_NAMESPACE_CLOSE
