/*-----------------------------------------------------------*/
/*! \file

\brief Input parameter for beam-to-solid interaction.


\level 2

*/
/*-----------------------------------------------------------*/


#include "inpar_beam_to_solid.H"

#include "drt_validparameters.H"
#include "inpar_beaminteraction.H"
#include "inpar_geometry_pair.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_conditiondefinition.H"


/**
 *
 */
void INPAR::BEAMTOSOLID::BeamToSolidInteractionGetString(
    const INPAR::BEAMINTERACTION::BeamInteractionConditions& interaction,
    std::array<std::string, 2>& condition_names)
{
  if (interaction ==
      INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying)
  {
    condition_names[0] = "BeamToSolidVolumeMeshtyingLine";
    condition_names[1] = "BeamToSolidVolumeMeshtyingVolume";
  }
  else if (interaction ==
           INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying)
  {
    condition_names[0] = "BeamToSolidSurfaceMeshtyingLine";
    condition_names[1] = "BeamToSolidSurfaceMeshtyingSurface";
  }
  else
    dserror("Got unexpected beam-to-solid interaction type.");
}

/**
 *
 */
void INPAR::BEAMTOSOLID::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAM INTERACTION", false, "");

  // Beam to solid volume mesh tying parameters.
  Teuchos::ParameterList& beam_to_solid_volume_mestying =
      beaminteraction.sublist("BEAM TO SOLID VOLUME MESHTYING", false, "");
  {
    setStringToIntegralParameter<BeamToSolidContactDiscretization>("CONTACT_DISCRETIZATION", "none",
        "Type of employed contact discretization",
        tuple<std::string>("none", "gauss_point_to_segment", "mortar", "gauss_point_cross_section",
            "gauss_point_cross_section_rotation"),
        tuple<BeamToSolidContactDiscretization>(BeamToSolidContactDiscretization::none,
            BeamToSolidContactDiscretization::gauss_point_to_segment,
            BeamToSolidContactDiscretization::mortar,
            BeamToSolidContactDiscretization::gauss_point_cross_section,
            BeamToSolidContactDiscretization::gauss_point_cross_section_rotation),
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

    DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid volume meshtying", &beam_to_solid_volume_mestying);

    // Add the geometry pair input parameters.
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_volume_mestying);

    // This option only has an effect during a restart simulation.
    // - No:  (default) The coupling is treated the same way as during a non restart simulation,
    //        i.e. the initial configurations (zero displacement) of the beams and solids are
    //        coupled.
    // - Yes: The beam and solid states at the restart configuration are coupled. This allows to
    //        pre-deform the structures and then couple them.
    BoolParameter("COUPLE_RESTART_STATE", "No",
        "Enable / disable the coupling of the restart configuration.",
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidRotationCoupling>("ROTATION_COUPLING", "none",
        "Type of rotational coupling",
        tuple<std::string>("none", "polar_decomposition_2d", "deformation_gradient_y_2d",
            "deformation_gradient_z_2d", "deformation_gradient_average_2d", "fix_triad_2d",
            "deformation_gradient_3d_local_1", "deformation_gradient_3d_local_2",
            "deformation_gradient_3d_local_3", "deformation_gradient_3d_general",
            "deformation_gradient_3d_general_in_cross_section_plane",
            "deformation_gradient_3d_base_1"),
        tuple<BeamToSolidRotationCoupling>(BeamToSolidRotationCoupling::none,
            BeamToSolidRotationCoupling::polar_decomposition_2d,
            BeamToSolidRotationCoupling::deformation_gradient_y_2d,
            BeamToSolidRotationCoupling::deformation_gradient_z_2d,
            BeamToSolidRotationCoupling::deformation_gradient_average_2d,
            BeamToSolidRotationCoupling::fix_triad_2d,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_1,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_2,
            BeamToSolidRotationCoupling::deformation_gradient_3d_local_3,
            BeamToSolidRotationCoupling::deformation_gradient_3d_general,
            BeamToSolidRotationCoupling::deformation_gradient_3d_general_in_cross_section_plane,
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

    DoubleParameter("ROTATION_COUPLING_PENALTY_PARAMETER", 0.0,
        "Penalty parameter for rotational coupling in beam-to-solid volume mesh tying",
        &beam_to_solid_volume_mestying);
  }

  // Beam to solid volume mesh tying output parameters.
  Teuchos::ParameterList& beam_to_solid_volume_mestying_vtk =
      beam_to_solid_volume_mestying.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write vtp output at all for btsvmt.
    BoolParameter("WRITE_OUTPUT", "No", "Enable / disable beam-to-solid volume mesh tying output.",
        &beam_to_solid_volume_mestying_vtk);

    BoolParameter("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        &beam_to_solid_volume_mestying_vtk);

    BoolParameter("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        &beam_to_solid_volume_mestying_vtk);

    BoolParameter("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        &beam_to_solid_volume_mestying_vtk);

    DRT::INPUT::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_volume_mestying_vtk);

    BoolParameter("SEGMENTATION", "No", "Enable / disable output of segmentation points.",
        &beam_to_solid_volume_mestying_vtk);

    BoolParameter("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        &beam_to_solid_volume_mestying_vtk);

    BoolParameter("UNIQUE_IDS", "No",
        "Enable / disable output of unique IDs (mainly for testing of created VTK files).",
        &beam_to_solid_volume_mestying_vtk);
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

    DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid surface meshtying", &beam_to_solid_surface_mestying);

    // Parameters for rotational coupling.
    BoolParameter("ROTATIONAL_COUPLING", "No", "Enable / disable rotational coupling",
        &beam_to_solid_surface_mestying);
    DoubleParameter("ROTATIONAL_COUPLING_PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid surface rotational meshtying",
        &beam_to_solid_surface_mestying);

    // Add the geometry pair input parameters.
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_surface_mestying);

    // Add the surface options.
    INPAR::GEOMETRYPAIR::SetValidParametersLineToSurface(beam_to_solid_surface_mestying);
  }

  // Beam to solid surface parameters.
  Teuchos::ParameterList& beam_to_solid_surface =
      beaminteraction.sublist("BEAM TO SOLID SURFACE", false, "");

  // Beam to solid surface output parameters.
  Teuchos::ParameterList& beam_to_solid_surface_vtk =
      beam_to_solid_surface.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write vtp output at all.
    BoolParameter("WRITE_OUTPUT", "No", "Enable / disable beam-to-solid volume mesh tying output.",
        &beam_to_solid_surface_vtk);

    BoolParameter("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        &beam_to_solid_surface_vtk);

    BoolParameter("AVERAGED_NORMALS", "No",
        "Enable / disable output of averaged nodal normals on the surface.",
        &beam_to_solid_surface_vtk);

    BoolParameter("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        &beam_to_solid_surface_vtk);

    BoolParameter("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        &beam_to_solid_surface_vtk);

    DRT::INPUT::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_surface_vtk);

    BoolParameter("SEGMENTATION", "No", "Enable / disable output of segmentation points.",
        &beam_to_solid_surface_vtk);

    BoolParameter("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        &beam_to_solid_surface_vtk);

    BoolParameter("UNIQUE_IDS", "No",
        "Enable / disable output of unique IDs (mainly for testing of created VTK files).",
        &beam_to_solid_surface_vtk);
  }
}

/**
 *
 */
void INPAR::BEAMTOSOLID::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  // Beam-to-volume mesh tying conditions.
  {
    std::array<std::string, 2> condition_names;
    BeamToSolidInteractionGetString(
        INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying,
        condition_names);

    Teuchos::RCP<ConditionDefinition> beam_to_solid_volume_meshtying_condition = Teuchos::rcp(
        new ConditionDefinition("BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING VOLUME",
            condition_names[1], "Beam-to-volume mesh tying conditions - volume definition",
            DRT::Condition::BeamToSolidVolumeMeshtyingVolume, true, DRT::Condition::Volume));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("COUPLING_ID")));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new IntConditionComponent("COUPLING_ID", false, false)));
    condlist.push_back(beam_to_solid_volume_meshtying_condition);

    beam_to_solid_volume_meshtying_condition =
        Teuchos::rcp(new ConditionDefinition("BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING LINE",
            condition_names[0], "Beam-to-volume mesh tying conditions - line definition",
            DRT::Condition::BeamToSolidVolumeMeshtyingLine, true, DRT::Condition::Line));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("COUPLING_ID")));
    beam_to_solid_volume_meshtying_condition->AddComponent(
        Teuchos::rcp(new IntConditionComponent("COUPLING_ID", false, false)));
    condlist.push_back(beam_to_solid_volume_meshtying_condition);
  }

  // Beam-to-surface mesh tying conditions.
  {
    std::array<std::string, 2> condition_names;
    BeamToSolidInteractionGetString(
        INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying,
        condition_names);

    Teuchos::RCP<ConditionDefinition> beam_to_solid_surface_meshtying_condition = Teuchos::rcp(
        new ConditionDefinition("BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING SURFACE",
            condition_names[1], "Beam-to-surface mesh tying conditions - surface definition",
            DRT::Condition::BeamToSolidSurfaceMeshtyingSurface, true, DRT::Condition::Surface));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("COUPLING_ID")));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new IntConditionComponent("COUPLING_ID", false, false)));
    condlist.push_back(beam_to_solid_surface_meshtying_condition);

    beam_to_solid_surface_meshtying_condition = Teuchos::rcp(
        new ConditionDefinition("BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING LINE",
            condition_names[0], "Beam-to-surface mesh tying conditions - line definition",
            DRT::Condition::BeamToSolidSurfaceMeshtyingLine, true, DRT::Condition::Line));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("COUPLING_ID")));
    beam_to_solid_surface_meshtying_condition->AddComponent(
        Teuchos::rcp(new IntConditionComponent("COUPLING_ID", false, false)));
    condlist.push_back(beam_to_solid_surface_meshtying_condition);
  }
}
