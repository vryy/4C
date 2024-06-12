/*-----------------------------------------------------------*/
/*! \file

\brief Input parameter for beam-to-solid interaction.


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_inpar_beam_to_solid.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_inpar_geometry_pair.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
void Inpar::BeamToSolid::BeamToSolidInteractionGetString(
    const Inpar::BEAMINTERACTION::BeamInteractionConditions& interaction,
    std::array<std::string, 2>& condition_names)
{
  if (interaction ==
      Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying)
  {
    condition_names[0] = "BeamToSolidVolumeMeshtyingLine";
    condition_names[1] = "BeamToSolidVolumeMeshtyingVolume";
  }
  else if (interaction ==
           Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying)
  {
    condition_names[0] = "BeamToSolidSurfaceMeshtyingLine";
    condition_names[1] = "BeamToSolidSurfaceMeshtyingSurface";
  }
  else if (interaction ==
           Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_contact)
  {
    condition_names[0] = "BeamToSolidSurfaceContactLine";
    condition_names[1] = "BeamToSolidSurfaceContactSurface";
  }
  else
    FOUR_C_THROW("Got unexpected beam-to-solid interaction type.");
}

/**
 *
 */
void Inpar::BeamToSolid::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAM INTERACTION", false, "");

  // Beam to solid volume mesh tying parameters.
  Teuchos::ParameterList& beam_to_solid_volume_mestying =
      beaminteraction.sublist("BEAM TO SOLID VOLUME MESHTYING", false, "");
  {
    setStringToIntegralParameter<BeamToSolidContactDiscretization>("CONTACT_DISCRETIZATION", "none",
        "Type of employed contact discretization",
        tuple<std::string>("none", "gauss_point_to_segment", "mortar", "gauss_point_cross_section"),
        tuple<BeamToSolidContactDiscretization>(BeamToSolidContactDiscretization::none,
            BeamToSolidContactDiscretization::gauss_point_to_segment,
            BeamToSolidContactDiscretization::mortar,
            BeamToSolidContactDiscretization::gauss_point_cross_section),
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidConstraintEnforcement>("CONSTRAINT_STRATEGY", "none",
        "Type of employed constraint enforcement strategy", tuple<std::string>("none", "penalty"),
        tuple<BeamToSolidConstraintEnforcement>(
            BeamToSolidConstraintEnforcement::none, BeamToSolidConstraintEnforcement::penalty),
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidMortarShapefunctions>("MORTAR_SHAPE_FUNCTION", "none",
        "Shape function for the mortar Lagrange-multipliers",
        tuple<std::string>("none", "line2", "line3", "line4"),
        tuple<BeamToSolidMortarShapefunctions>(BeamToSolidMortarShapefunctions::none,
            BeamToSolidMortarShapefunctions::line2, BeamToSolidMortarShapefunctions::line3,
            BeamToSolidMortarShapefunctions::line4),
        &beam_to_solid_volume_mestying);

    Core::UTILS::DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid volume meshtying", &beam_to_solid_volume_mestying);

    // Add the geometry pair input parameters.
    Inpar::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_volume_mestying);

    // This option only has an effect during a restart simulation.
    // - No:  (default) The coupling is treated the same way as during a non restart simulation,
    //        i.e. the initial configurations (zero displacement) of the beams and solids are
    //        coupled.
    // - Yes: The beam and solid states at the restart configuration are coupled. This allows to
    //        pre-deform the structures and then couple them.
    Core::UTILS::BoolParameter("COUPLE_RESTART_STATE", "No",
        "Enable / disable the coupling of the restart configuration.",
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidRotationCoupling>("ROTATION_COUPLING", "none",
        "Type of rotational coupling",
        tuple<std::string>("none", "deformation_gradient_3d_general_in_cross_section_plane",
            "polar_decomposition_2d", "deformation_gradient_y_2d", "deformation_gradient_z_2d",
            "deformation_gradient_average_2d", "fix_triad_2d", "deformation_gradient_3d_local_1",
            "deformation_gradient_3d_local_2", "deformation_gradient_3d_local_3",
            "deformation_gradient_3d_general",

            "deformation_gradient_3d_base_1"),
        tuple<BeamToSolidRotationCoupling>(BeamToSolidRotationCoupling::none,
            BeamToSolidRotationCoupling::deformation_gradient_3d_general_in_cross_section_plane,
            BeamToSolidRotationCoupling::polar_decomposition_2d,
            BeamToSolidRotationCoupling::deformation_gradient_y_2d,
            BeamToSolidRotationCoupling::deformation_gradient_z_2d,
            BeamToSolidRotationCoupling::deformation_gradient_average_2d,
            BeamToSolidRotationCoupling::fix_triad_2d,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_1,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_2,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_3,
            BeamToSolidRotationCoupling::deformation_gradient_3d_general,
            BeamToSolidRotationCoupling::deformation_gradient_3d_base_1),
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidMortarShapefunctions>(
        "ROTATION_COUPLING_MORTAR_SHAPE_FUNCTION", "none",
        "Shape function for the mortar Lagrange-multipliers",
        tuple<std::string>("none", "line2", "line3", "line4"),
        tuple<BeamToSolidMortarShapefunctions>(BeamToSolidMortarShapefunctions::none,
            BeamToSolidMortarShapefunctions::line2, BeamToSolidMortarShapefunctions::line3,
            BeamToSolidMortarShapefunctions::line4),
        &beam_to_solid_volume_mestying);

    Core::UTILS::DoubleParameter("ROTATION_COUPLING_PENALTY_PARAMETER", 0.0,
        "Penalty parameter for rotational coupling in beam-to-solid volume mesh tying",
        &beam_to_solid_volume_mestying);
  }

  // Beam to solid volume mesh tying output parameters.
  Teuchos::ParameterList& beam_to_solid_volume_mestying_output =
      beam_to_solid_volume_mestying.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write visualization output at all for btsvmt.
    Core::UTILS::BoolParameter("WRITE_OUTPUT", "No",
        "Enable / disable beam-to-solid volume mesh tying output.",
        &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        &beam_to_solid_volume_mestying_output);

    Core::UTILS::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("SEGMENTATION", "No",
        "Enable / disable output of segmentation points.", &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        &beam_to_solid_volume_mestying_output);

    Core::UTILS::BoolParameter("UNIQUE_IDS", "No",
        "Enable / disable output of unique IDs (mainly for testing of created VTK files).",
        &beam_to_solid_volume_mestying_output);
  }

  // Beam to solid surface mesh tying parameters.
  Teuchos::ParameterList& beam_to_solid_surface_mestying =
      beaminteraction.sublist("BEAM TO SOLID SURFACE MESHTYING", false, "");
  {
    setStringToIntegralParameter<BeamToSolidContactDiscretization>("CONTACT_DISCRETIZATION", "none",
        "Type of employed contact discretization",
        tuple<std::string>("none", "gauss_point_to_segment", "mortar"),
        tuple<BeamToSolidContactDiscretization>(BeamToSolidContactDiscretization::none,
            BeamToSolidContactDiscretization::gauss_point_to_segment,
            BeamToSolidContactDiscretization::mortar),
        &beam_to_solid_surface_mestying);

    setStringToIntegralParameter<BeamToSolidConstraintEnforcement>("CONSTRAINT_STRATEGY", "none",
        "Type of employed constraint enforcement strategy", tuple<std::string>("none", "penalty"),
        tuple<BeamToSolidConstraintEnforcement>(
            BeamToSolidConstraintEnforcement::none, BeamToSolidConstraintEnforcement::penalty),
        &beam_to_solid_surface_mestying);

    setStringToIntegralParameter<BeamToSolidSurfaceCoupling>("COUPLING_TYPE", "none",
        "How the coupling constraints are formulated/",
        tuple<std::string>("none", "reference_configuration_forced_to_zero",
            "reference_configuration_forced_to_zero_fad", "displacement", "displacement_fad",
            "consistent_fad"),
        tuple<BeamToSolidSurfaceCoupling>(BeamToSolidSurfaceCoupling::none,
            BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero,
            BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero_fad,
            BeamToSolidSurfaceCoupling::displacement, BeamToSolidSurfaceCoupling::displacement_fad,
            BeamToSolidSurfaceCoupling::consistent_fad),
        &beam_to_solid_surface_mestying);

    setStringToIntegralParameter<BeamToSolidMortarShapefunctions>("MORTAR_SHAPE_FUNCTION", "none",
        "Shape function for the mortar Lagrange-multipliers",
        tuple<std::string>("none", "line2", "line3", "line4"),
        tuple<BeamToSolidMortarShapefunctions>(BeamToSolidMortarShapefunctions::none,
            BeamToSolidMortarShapefunctions::line2, BeamToSolidMortarShapefunctions::line3,
            BeamToSolidMortarShapefunctions::line4),
        &beam_to_solid_surface_mestying);

    Core::UTILS::DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid surface meshtying", &beam_to_solid_surface_mestying);

    // Parameters for rotational coupling.
    Core::UTILS::BoolParameter("ROTATIONAL_COUPLING", "No", "Enable / disable rotational coupling",
        &beam_to_solid_surface_mestying);
    Core::UTILS::DoubleParameter("ROTATIONAL_COUPLING_PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid surface rotational meshtying",
        &beam_to_solid_surface_mestying);
    setStringToIntegralParameter<BeamToSolidSurfaceRotationCoupling>(
        "ROTATIONAL_COUPLING_SURFACE_TRIAD", "none", "Construction method for surface triad",
        tuple<std::string>("none", "surface_cross_section_director", "averaged"),
        tuple<BeamToSolidSurfaceRotationCoupling>(BeamToSolidSurfaceRotationCoupling::none,
            BeamToSolidSurfaceRotationCoupling::surface_cross_section_director,
            BeamToSolidSurfaceRotationCoupling::averaged),
        &beam_to_solid_surface_mestying);

    // Add the geometry pair input parameters.
    Inpar::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_surface_mestying);

    // Add the surface options.
    Inpar::GEOMETRYPAIR::SetValidParametersLineToSurface(beam_to_solid_surface_mestying);
  }

  // Beam to solid surface contact parameters.
  Teuchos::ParameterList& beam_to_solid_surface_contact =
      beaminteraction.sublist("BEAM TO SOLID SURFACE CONTACT", false, "");
  {
    setStringToIntegralParameter<BeamToSolidContactDiscretization>("CONTACT_DISCRETIZATION", "none",
        "Type of employed contact discretization",
        tuple<std::string>("none", "gauss_point_to_segment", "mortar"),
        tuple<BeamToSolidContactDiscretization>(BeamToSolidContactDiscretization::none,
            BeamToSolidContactDiscretization::gauss_point_to_segment,
            BeamToSolidContactDiscretization::mortar),
        &beam_to_solid_surface_contact);

    setStringToIntegralParameter<BeamToSolidConstraintEnforcement>("CONSTRAINT_STRATEGY", "none",
        "Type of employed constraint enforcement strategy", tuple<std::string>("none", "penalty"),
        tuple<BeamToSolidConstraintEnforcement>(
            BeamToSolidConstraintEnforcement::none, BeamToSolidConstraintEnforcement::penalty),
        &beam_to_solid_surface_contact);

    Core::UTILS::DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid surface contact", &beam_to_solid_surface_contact);

    setStringToIntegralParameter<BeamToSolidSurfaceContact>("CONTACT_TYPE", "none",
        "How the contact constraints are formulated",
        tuple<std::string>("none", "gap_variation", "potential"),
        tuple<BeamToSolidSurfaceContact>(BeamToSolidSurfaceContact::none,
            BeamToSolidSurfaceContact::gap_variation, BeamToSolidSurfaceContact::potential),
        &beam_to_solid_surface_contact);

    setStringToIntegralParameter<BeamToSolidSurfaceContactPenaltyLaw>("PENALTY_LAW", "none",
        "Type of penalty law", tuple<std::string>("none", "linear", "linear_quadratic"),
        tuple<BeamToSolidSurfaceContactPenaltyLaw>(BeamToSolidSurfaceContactPenaltyLaw::none,
            BeamToSolidSurfaceContactPenaltyLaw::linear,
            BeamToSolidSurfaceContactPenaltyLaw::linear_quadratic),
        &beam_to_solid_surface_contact);

    Core::UTILS::DoubleParameter("PENALTY_PARAMETER_G0", 0.0,
        "First penalty regularization parameter G0 >=0: For gap<G0 contact is active",
        &beam_to_solid_surface_contact);

    setStringToIntegralParameter<BeamToSolidSurfaceContactMortarDefinedIn>(
        "MORTAR_CONTACT_DEFINED_IN", "none", "Configuration where the mortar contact is defined",
        tuple<std::string>("none", "reference_configuration", "current_configuration"),
        tuple<BeamToSolidSurfaceContactMortarDefinedIn>(
            BeamToSolidSurfaceContactMortarDefinedIn::none,
            BeamToSolidSurfaceContactMortarDefinedIn::reference_configuration,
            BeamToSolidSurfaceContactMortarDefinedIn::current_configuration),
        &beam_to_solid_surface_contact);

    // Add the geometry pair input parameters.
    Inpar::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_surface_contact);

    // Add the surface options.
    Inpar::GEOMETRYPAIR::SetValidParametersLineToSurface(beam_to_solid_surface_contact);

    // Define the mortar shape functions for contact
    setStringToIntegralParameter<BeamToSolidMortarShapefunctions>("MORTAR_SHAPE_FUNCTION", "none",
        "Shape function for the mortar Lagrange-multipliers", tuple<std::string>("none", "line2"),
        tuple<BeamToSolidMortarShapefunctions>(
            BeamToSolidMortarShapefunctions::none, BeamToSolidMortarShapefunctions::line2),
        &beam_to_solid_surface_contact);
  }

  // Beam to solid surface parameters.
  Teuchos::ParameterList& beam_to_solid_surface =
      beaminteraction.sublist("BEAM TO SOLID SURFACE", false, "");

  // Beam to solid surface output parameters.
  Teuchos::ParameterList& beam_to_solid_surface_output =
      beam_to_solid_surface.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write visualization output at all.
    Core::UTILS::BoolParameter("WRITE_OUTPUT", "No",
        "Enable / disable beam-to-solid volume mesh tying output.", &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("AVERAGED_NORMALS", "No",
        "Enable / disable output of averaged nodal normals on the surface.",
        &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        &beam_to_solid_surface_output);

    Core::UTILS::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("SEGMENTATION", "No",
        "Enable / disable output of segmentation points.", &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        &beam_to_solid_surface_output);

    Core::UTILS::BoolParameter("UNIQUE_IDS", "No",
        "Enable / disable output of unique IDs (mainly for testing of created VTK files).",
        &beam_to_solid_surface_output);
  }
}

/**
 *
 */
void Inpar::BeamToSolid::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  // Beam-to-volume mesh tying conditions.
  {
    std::array<std::string, 2> condition_names;
    BeamToSolidInteractionGetString(
        Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying,
        condition_names);

    Teuchos::RCP<Core::Conditions::ConditionDefinition> beam_to_solid_volume_meshtying_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING VOLUME", condition_names[1],
            "Beam-to-volume mesh tying conditions - volume definition",
            Core::Conditions::BeamToSolidVolumeMeshtyingVolume, true,
            Core::Conditions::geometry_type_volume));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_volume_meshtying_condition);

    beam_to_solid_volume_meshtying_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING LINE", condition_names[0],
            "Beam-to-volume mesh tying conditions - line definition",
            Core::Conditions::BeamToSolidVolumeMeshtyingLine, true,
            Core::Conditions::geometry_type_line));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_volume_meshtying_condition);
  }

  // Beam-to-surface mesh tying conditions.
  {
    std::array<std::string, 2> condition_names;
    BeamToSolidInteractionGetString(
        Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying,
        condition_names);

    Teuchos::RCP<Core::Conditions::ConditionDefinition> beam_to_solid_surface_meshtying_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING SURFACE", condition_names[1],
            "Beam-to-surface mesh tying conditions - surface definition",
            Core::Conditions::BeamToSolidSurfaceMeshtyingSurface, true,
            Core::Conditions::geometry_type_surface));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_surface_meshtying_condition);

    beam_to_solid_surface_meshtying_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING LINE", condition_names[0],
            "Beam-to-surface mesh tying conditions - line definition",
            Core::Conditions::BeamToSolidSurfaceMeshtyingLine, true,
            Core::Conditions::geometry_type_line));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_surface_meshtying_condition);
  }

  // Beam-to-surface contact conditions.
  {
    std::array<std::string, 2> condition_names;
    BeamToSolidInteractionGetString(
        Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_contact,
        condition_names);

    Teuchos::RCP<Core::Conditions::ConditionDefinition> beam_to_solid_surface_contact_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID SURFACE CONTACT SURFACE", condition_names[1],
            "Beam-to-surface contact conditions - surface definition",
            Core::Conditions::BeamToSolidSurfaceContactSurface, true,
            Core::Conditions::geometry_type_surface));
    beam_to_solid_surface_contact_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_surface_contact_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_surface_contact_condition);

    beam_to_solid_surface_contact_condition =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "BEAM INTERACTION/BEAM TO SOLID SURFACE CONTACT LINE", condition_names[0],
            "Beam-to-surface contact conditions - line definition",
            Core::Conditions::BeamToSolidSurfaceContactLine, true,
            Core::Conditions::geometry_type_line));
    beam_to_solid_surface_contact_condition->AddComponent(
        Teuchos::rcp(new Input::SeparatorComponent("COUPLING_ID")));
    beam_to_solid_surface_contact_condition->AddComponent(
        Teuchos::rcp(new Input::IntComponent("COUPLING_ID")));
    condlist.push_back(beam_to_solid_surface_contact_condition);
  }
}

FOUR_C_NAMESPACE_CLOSE
