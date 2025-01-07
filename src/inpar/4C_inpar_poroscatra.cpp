// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_poroscatra.hpp"

#include "4C_inpar_poroelast.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::PoroScaTra::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& poroscatradyn = list.sublist(
      "POROSCATRA CONTROL", false, "Control parameters for scatra porous media coupling");

  // Output type
  Core::Utils::int_parameter(
      "RESTARTEVERY", 1, "write restart possibility every RESTARTEVERY steps", &poroscatradyn);
  // Time loop control
  Core::Utils::int_parameter("NUMSTEP", 200, "maximum number of Timesteps", &poroscatradyn);
  Core::Utils::double_parameter("MAXTIME", 1000.0, "total simulation time", &poroscatradyn);
  Core::Utils::double_parameter("TIMESTEP", 0.05, "time step size dt", &poroscatradyn);
  Core::Utils::int_parameter("RESULTSEVERY", 1, "increment for writing solution", &poroscatradyn);
  Core::Utils::int_parameter(
      "ITEMAX", 10, "maximum number of iterations over fields", &poroscatradyn);
  Core::Utils::int_parameter(
      "ITEMIN", 1, "minimal number of iterations over fields", &poroscatradyn);

  // Iterationparameters
  Core::Utils::double_parameter("TOLRES_GLOBAL", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLINC_GLOBAL", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLRES_DISP", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLINC_DISP", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLRES_VEL", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLINC_VEL", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLRES_PRES", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLINC_PRES", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLRES_SCALAR", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poroscatradyn);
  Core::Utils::double_parameter("TOLINC_SCALAR", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poroscatradyn);

  setStringToIntegralParameter<Inpar::PoroElast::ConvNorm>("NORM_INC", "AbsSingleFields",
      "type of norm for primary variables convergence check",
      tuple<std::string>("AbsGlobal", "AbsSingleFields"),
      tuple<Inpar::PoroElast::ConvNorm>(
          Inpar::PoroElast::convnorm_abs_global, Inpar::PoroElast::convnorm_abs_singlefields),
      &poroscatradyn);

  setStringToIntegralParameter<Inpar::PoroElast::ConvNorm>("NORM_RESF", "AbsSingleFields",
      "type of norm for residual convergence check",
      tuple<std::string>("AbsGlobal", "AbsSingleFields"),
      tuple<Inpar::PoroElast::ConvNorm>(
          Inpar::PoroElast::convnorm_abs_global, Inpar::PoroElast::convnorm_abs_singlefields),
      &poroscatradyn);

  setStringToIntegralParameter<Inpar::PoroElast::BinaryOp>("NORMCOMBI_RESFINC", "And",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or"),
      tuple<Inpar::PoroElast::BinaryOp>(Inpar::PoroElast::bop_and, Inpar::PoroElast::bop_or),
      &poroscatradyn);

  setStringToIntegralParameter<Inpar::PoroElast::VectorNorm>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<Inpar::PoroElast::VectorNorm>(Inpar::PoroElast::norm_l1,
          Inpar::PoroElast::norm_l1_scaled, Inpar::PoroElast::norm_l2, Inpar::PoroElast::norm_rms,
          Inpar::PoroElast::norm_inf),
      &poroscatradyn);

  setStringToIntegralParameter<Inpar::PoroElast::VectorNorm>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<Inpar::PoroElast::VectorNorm>(Inpar::PoroElast::norm_l1,
          Inpar::PoroElast::norm_l1_scaled, Inpar::PoroElast::norm_l2, Inpar::PoroElast::norm_rms,
          Inpar::PoroElast::norm_inf),
      &poroscatradyn);

  // number of linear solver used for poroelasticity
  Core::Utils::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for monolithic poroscatra problems", &poroscatradyn);

  // Coupling strategy for poroscatra solvers
  setStringToIntegralParameter<SolutionSchemeOverFields>("COUPALGO", "solid_to_scatra",
      "Coupling strategies for poroscatra solvers",
      tuple<std::string>("monolithic", "scatra_to_solid", "solid_to_scatra", "two_way"),
      tuple<SolutionSchemeOverFields>(
          Monolithic, Part_ScatraToPoro, Part_PoroToScatra, Part_TwoWay),
      &poroscatradyn);

  Core::Utils::bool_parameter("MATCHINGGRID", "Yes", "is matching grid", &poroscatradyn);
}

FOUR_C_NAMESPACE_CLOSE
