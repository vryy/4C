// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_lubrication.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::Lubrication::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& lubricationdyn =
      list.sublist("LUBRICATION DYNAMIC", false, "control parameters for Lubrication problems\n");

  Core::Utils::double_parameter("MAXTIME", 1000.0, "Total simulation time", &lubricationdyn);
  Core::Utils::int_parameter("NUMSTEP", 20, "Total number of time steps", &lubricationdyn);
  Core::Utils::double_parameter("TIMESTEP", 0.1, "Time increment dt", &lubricationdyn);
  Core::Utils::int_parameter("RESULTSEVRY", 1, "Increment for writing solution", &lubricationdyn);
  Core::Utils::int_parameter("RESTARTEVRY", 1, "Increment for writing restart", &lubricationdyn);

  setStringToIntegralParameter<Inpar::Lubrication::CalcError>("CALCERROR", "No",
      "compute error compared to analytical solution",
      tuple<std::string>("No", "error_by_function"),
      tuple<Inpar::Lubrication::CalcError>(calcerror_no, calcerror_byfunction), &lubricationdyn);

  Core::Utils::int_parameter(
      "CALCERRORNO", -1, "function number for lubrication error computation", &lubricationdyn);

  setStringToIntegralParameter<Inpar::Lubrication::VelocityField>("VELOCITYFIELD", "zero",
      "type of velocity field used for lubrication problems",
      tuple<std::string>("zero", "function", "EHL"),
      tuple<Inpar::Lubrication::VelocityField>(velocity_zero, velocity_function, velocity_EHL),
      &lubricationdyn);

  Core::Utils::int_parameter(
      "VELFUNCNO", -1, "function number for lubrication velocity field", &lubricationdyn);

  setStringToIntegralParameter<Inpar::Lubrication::HeightField>("HEIGHTFEILD", "zero",
      "type of height field used for lubrication problems",
      tuple<std::string>("zero", "function", "EHL"),
      tuple<Inpar::Lubrication::HeightField>(height_zero, height_function, height_EHL),
      &lubricationdyn);

  Core::Utils::int_parameter(
      "HFUNCNO", -1, "function number for lubrication height field", &lubricationdyn);

  Core::Utils::bool_parameter(
      "OUTMEAN", "No", "Output of mean values for scalars and density", &lubricationdyn);

  Core::Utils::bool_parameter(
      "OUTPUT_GMSH", "No", "Do you want to write Gmsh postprocessing files?", &lubricationdyn);

  Core::Utils::bool_parameter("MATLAB_STATE_OUTPUT", "No",
      "Do you want to write the state solution to Matlab file?", &lubricationdyn);

  /// linear solver id used for lubrication problems
  Core::Utils::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for the Lubrication problem", &lubricationdyn);

  Core::Utils::int_parameter("ITEMAX", 10, "max. number of nonlin. iterations", &lubricationdyn);
  Core::Utils::double_parameter("ABSTOLRES", 1e-14,
      "Absolute tolerance for deciding if residual of nonlinear problem is already zero",
      &lubricationdyn);
  Core::Utils::double_parameter(
      "CONVTOL", 1e-13, "Tolerance for convergence check", &lubricationdyn);

  // convergence criteria adaptivity
  Core::Utils::bool_parameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &lubricationdyn);
  Core::Utils::double_parameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &lubricationdyn);

  setStringToIntegralParameter<ConvNorm>("NORM_PRE", "Abs",
      "type of norm for temperature convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<ConvNorm>(convnorm_abs, convnorm_rel, convnorm_mix), &lubricationdyn);

  setStringToIntegralParameter<ConvNorm>("NORM_RESF", "Abs",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<ConvNorm>(convnorm_abs, convnorm_rel, convnorm_mix), &lubricationdyn);

  setStringToIntegralParameter<VectorNorm>("ITERNORM", "L2",
      "type of norm to be applied to residuals", tuple<std::string>("L1", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(norm_l1, norm_l2, norm_rms, norm_inf), &lubricationdyn);

  /// Iterationparameters
  Core::Utils::double_parameter("TOLPRE", 1.0E-06,
      "tolerance in the temperature norm of the Newton iteration", &lubricationdyn);

  Core::Utils::double_parameter("TOLRES", 1.0E-06,
      "tolerance in the residual norm for the Newton iteration", &lubricationdyn);

  Core::Utils::double_parameter(
      "PENALTY_CAVITATION", 0., "penalty parameter for regularized cavitation", &lubricationdyn);

  Core::Utils::double_parameter(
      "GAP_OFFSET", 0., "Additional offset to the fluid gap", &lubricationdyn);

  Core::Utils::double_parameter(
      "ROUGHNESS_STD_DEVIATION", 0., "standard deviation of surface roughness", &lubricationdyn);

  /// use modified reynolds equ.
  Core::Utils::bool_parameter("MODIFIED_REYNOLDS_EQU", "No",
      "the lubrication problem will use the modified reynolds equ. in order to consider surface"
      " roughness",
      &lubricationdyn);

  /// Flag for considering the Squeeze term in Reynolds Equation
  Core::Utils::bool_parameter("ADD_SQUEEZE_TERM", "No",
      "the squeeze term will also be considered in the Reynolds Equation", &lubricationdyn);

  /// Flag for considering the pure Reynolds Equation
  Core::Utils::bool_parameter("PURE_LUB", "No", "the problem is pure lubrication", &lubricationdyn);
}

FOUR_C_NAMESPACE_CLOSE
