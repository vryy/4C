// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_levelset_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

std::vector<Core::IO::InputSpec> LevelSet::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("LEVEL-SET CONTROL",
      {

          parameter<int>(
              "NUMSTEP", {.description = "Total number of time steps", .default_value = 24}),

          parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}),
          parameter<double>(
              "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}),
          parameter<int>("RESULTSEVERY",
              {.description = "Increment for writing solution", .default_value = 1}),
          parameter<int>(
              "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}),


          deprecated_selection<LevelSet::CalcErrorLevelSet>("CALCERROR",
              {
                  {"No", LevelSet::calcerror_no_ls},
                  {"InitialField", LevelSet::calcerror_initial_field},
              },
              {.description = "compute error compared to analytical solution",
                  .default_value = LevelSet::calcerror_no_ls}),

          parameter<bool>("EXTRACT_INTERFACE_VEL",
              {.description =
                      "replace computed velocity at nodes of given distance of interface by "
                      "approximated interface velocity",
                  .default_value = false}),
          parameter<int>("NUM_CONVEL_LAYERS",
              {.description = "number of layers around the interface which keep their computed "
                              "convective velocity",
                  .default_value = -1})},
      {.required = false}));
  /*----------------------------------------------------------------------*/
  specs.push_back(group("LEVEL-SET CONTROL/REINITIALIZATION",
      {

          deprecated_selection<LevelSet::ReInitialAction>("REINITIALIZATION",
              {
                  {"None", LevelSet::reinitaction_none},
                  {"Signed_Distance_Function", LevelSet::reinitaction_signeddistancefunction},
                  {"Sussman", LevelSet::reinitaction_sussman},
                  {"EllipticEq", LevelSet::reinitaction_ellipticeq},
              },
              {.description = "Type of reinitialization strategy for level set function",
                  .default_value = LevelSet::reinitaction_none}),

          parameter<bool>("REINIT_INITIAL",
              {.description = "Has level set field to be reinitialized before first time step?",
                  .default_value = false}),
          parameter<int>(
              "REINITINTERVAL", {.description = "reinitialization interval", .default_value = 1}),

          // parameters for signed distance reinitialization
          parameter<bool>("REINITBAND",
              {.description =
                      "reinitialization only within a band around the interface, or entire domain?",
                  .default_value = false}),
          parameter<double>("REINITBANDWIDTH",
              {.description = "level-set value defining band width for reinitialization",
                  .default_value = 1.0}),

          // parameters for reinitialization equation
          parameter<int>("NUMSTEPSREINIT",
              {.description = "(maximal) number of pseudo-time steps", .default_value = 1}),
          parameter<double>("TIMESTEPREINIT",
              {.description = "pseudo-time step length (usually a * characteristic "
                              "element length of discretization with a>0)",
                  .default_value = 1.0}),
          parameter<double>("THETAREINIT",
              {.description = "theta for time discretization of reinitialization equation",
                  .default_value = 1.0}),
          deprecated_selection<Inpar::ScaTra::StabType>("STABTYPEREINIT",
              {
                  {"no_stabilization", Inpar::ScaTra::stabtype_no_stabilization},
                  {"SUPG", Inpar::ScaTra::stabtype_SUPG},
                  {"GLS", Inpar::ScaTra::stabtype_GLS},
                  {"USFEM", Inpar::ScaTra::stabtype_USFEM},
              },
              {.description =
                      "Type of stabilization (if any). No stabilization is only reasonable for "
                      "low-Peclet-number.",
                  .default_value = Inpar::ScaTra::stabtype_SUPG}),
          // this parameter selects the tau definition applied
          deprecated_selection<Inpar::ScaTra::TauType>("DEFINITION_TAU_REINIT",
              {
                  {"Taylor_Hughes_Zarins", Inpar::ScaTra::tau_taylor_hughes_zarins},
                  {"Taylor_Hughes_Zarins_wo_dt", Inpar::ScaTra::tau_taylor_hughes_zarins_wo_dt},
                  {"Franca_Valentin", Inpar::ScaTra::tau_franca_valentin},
                  {"Franca_Valentin_wo_dt", Inpar::ScaTra::tau_franca_valentin_wo_dt},
                  {"Shakib_Hughes_Codina", Inpar::ScaTra::tau_shakib_hughes_codina},
                  {"Shakib_Hughes_Codina_wo_dt", Inpar::ScaTra::tau_shakib_hughes_codina_wo_dt},
                  {"Codina", Inpar::ScaTra::tau_codina},
                  {"Codina_wo_dt", Inpar::ScaTra::tau_codina_wo_dt},
                  {"Exact_1D", Inpar::ScaTra::tau_exact_1d},
                  {"Zero", Inpar::ScaTra::tau_zero},
              },
              {.description = "Definition of tau",
                  .default_value = Inpar::ScaTra::tau_taylor_hughes_zarins}),
          // this parameter governs whether all-scale subgrid diffusivity is included
          deprecated_selection<LevelSet::ArtDiff>("ARTDIFFREINIT",
              {
                  {"no", LevelSet::artdiff_none},
                  {"isotropic", LevelSet::artdiff_isotropic},
                  {"crosswind", LevelSet::artdiff_crosswind},
              },
              {.description = "potential incorporation of all-scale subgrid diffusivity (a.k.a. "
                              "discontinuity-capturing) term",
                  .default_value = LevelSet::artdiff_none}),
          // this parameter selects the all-scale subgrid-diffusivity definition applied
          deprecated_selection<Inpar::ScaTra::AssgdType>("DEFINITION_ARTDIFFREINIT",
              {
                  {"artificial_linear", Inpar::ScaTra::assgd_artificial},
                  {"artificial_linear_reinit", Inpar::ScaTra::assgd_lin_reinit},
                  {"Hughes_etal_86_nonlinear", Inpar::ScaTra::assgd_hughes},
                  {"Tezduyar_Park_86_nonlinear", Inpar::ScaTra::assgd_tezduyar},
                  {"Tezduyar_Park_86_nonlinear_wo_phizero",
                      Inpar::ScaTra::assgd_tezduyar_wo_phizero},
                  {"doCarmo_Galeao_91_nonlinear", Inpar::ScaTra::assgd_docarmo},
                  {"Almeida_Silva_97_nonlinear", Inpar::ScaTra::assgd_almeida},
                  {"YZbeta_nonlinear", Inpar::ScaTra::assgd_yzbeta},
                  {"Codina_nonlinear", Inpar::ScaTra::assgd_codina},
              },
              {.description = "Definition of (all-scale) subgrid diffusivity",
                  .default_value = Inpar::ScaTra::assgd_artificial}),


          deprecated_selection<LevelSet::SmoothedSignType>("SMOOTHED_SIGN_TYPE",
              {
                  {"NonSmoothed", LevelSet::signtype_nonsmoothed},
                  {"SussmanFatemi1999", LevelSet::signtype_SussmanFatemi1999},
                  {"SussmanSmerekaOsher1994", LevelSet::signtype_SussmanSmerekaOsher1994},
                  {"PengEtAl1999", LevelSet::signtype_PengEtAl1999},
              },
              {.description = "sign function for reinitialization equation",
                  .default_value = LevelSet::signtype_SussmanSmerekaOsher1994}),

          deprecated_selection<LevelSet::CharEleLengthReinit>("CHARELELENGTHREINIT",
              {
                  {"root_of_volume", LevelSet::root_of_volume_reinit},
                  {"streamlength", LevelSet::streamlength_reinit},
              },
              {.description = "characteristic element length for sign function",
                  .default_value = LevelSet::root_of_volume_reinit}),
          parameter<double>("INTERFACE_THICKNESS",
              {.description = "factor for interface thickness (multiplied by element length)",
                  .default_value = 1.0}),
          deprecated_selection<LevelSet::VelReinit>("VELREINIT",
              {
                  {"integration_point_based", LevelSet::vel_reinit_integration_point_based},
                  {"node_based", LevelSet::vel_reinit_node_based},
              },
              {.description =
                      "evaluate velocity at integration point or compute node-based velocity",
                  .default_value = LevelSet::vel_reinit_integration_point_based}),
          parameter<LevelSet::LinReinit>("LINEARIZATIONREINIT",
              {.description = "linearization scheme for nonlinear convective term of "
                              "reinitialization equation",
                  .default_value = LevelSet::newton}),
          parameter<bool>("CORRECTOR_STEP",
              {.description = "correction of interface position via volume constraint "
                              "according to Sussman & Fatemi",
                  .default_value = true}),
          parameter<double>("CONVTOL_REINIT",
              {.description = "tolerance for convergence check according to Sussman et "
                              "al. 1994 (turned off negative)",
                  .default_value = -1.0}),

          parameter<bool>("REINITVOLCORRECTION",
              {.description = "volume correction after reinitialization", .default_value = false}),

          parameter<double>(
              "PENALTY_PARA", {.description = "penalty parameter for elliptic reinitialization",
                                  .default_value = -1.0}),

          deprecated_selection<LevelSet::LSDim>("DIMENSION",
              {
                  {"3D", LevelSet::ls_3D},
                  {"2Dx", LevelSet::ls_2Dx},
                  {"2Dy", LevelSet::ls_2Dy},
                  {"2Dz", LevelSet::ls_2Dz},
              },
              {.description = "number of space dimensions for handling of quasi-2D problems with "
                              "3D elements",
                  .default_value = LevelSet::ls_3D}),

          parameter<bool>(
              "PROJECTION", {.description = "use L2-projection for grad phi and related quantities",
                                .default_value = true}),
          parameter<double>("PROJECTION_DIFF",
              {.description = "use diffusive term for L2-projection", .default_value = 0.0}),
          parameter<bool>("LUMPING",
              {.description = "use lumped mass matrix for L2-projection", .default_value = false}),

          parameter<LevelSet::DiffFunc>("DIFF_FUNC",
              {.description = "function for diffusivity", .default_value = LevelSet::hyperbolic})},
      {.required = false}));
  return specs;
}



void LevelSet::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  /*--------------------------------------------------------------------*/
  // Taylor Galerkin outflow Boundaries for level set transport equation

  Core::Conditions::ConditionDefinition surfOutflowTaylorGalerkin(
      "TAYLOR GALERKIN OUTFLOW SURF CONDITIONS", "TaylorGalerkinOutflow",
      "Surface Taylor Galerkin Outflow", Core::Conditions::TaylorGalerkinOutflow, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(surfOutflowTaylorGalerkin);

  /*--------------------------------------------------------------------*/

  Core::Conditions::ConditionDefinition surfneumanninflowTaylorGalerkin(
      "TAYLOR GALERKIN NEUMANN INFLOW SURF CONDITIONS", "TaylorGalerkinNeumannInflow",
      "Surface Taylor Galerkin Neumann Inflow", Core::Conditions::TaylorGalerkinNeumannInflow, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(surfneumanninflowTaylorGalerkin);


  /*--------------------------------------------------------------------*/
  // Characteristic Galerkin Boundaries for LevelSet-reinitialization

  Core::Conditions::ConditionDefinition surfreinitializationtaylorgalerkin(
      "REINITIALIZATION TAYLOR GALERKIN SURF CONDITIONS", "ReinitializationTaylorGalerkin",
      "Surface reinitialization Taylor Galerkin", Core::Conditions::ReinitializationTaylorGalerkin,
      true, Core::Conditions::geometry_type_surface);

  condlist.push_back(surfreinitializationtaylorgalerkin);

  /*--------------------------------------------------------------------*/
  // level-set condition for contact points

  Core::Conditions::ConditionDefinition linelscontact("DESIGN LINE LEVEL SET CONTACT CONDITION",
      "LsContact", "level-set condition for contact points", Core::Conditions::LsContact, false,
      Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition pointlscontact("DESIGN POINT LEVEL SET CONTACT CONDITION",
      "LsContact", "level-set condition for contact points", Core::Conditions::LsContact, false,
      Core::Conditions::geometry_type_point);

  condlist.push_back(linelscontact);
  condlist.push_back(pointlscontact);
}

FOUR_C_NAMESPACE_CLOSE