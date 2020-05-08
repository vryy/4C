/*-----------------------------------------------------------*/
/*! \file

\brief Input parameter for beam-to-solid interaction.

\maintainer Ivo Steinbrecher

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

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

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
    setStringToIntegralParameter<int>("COUPLE_RESTART_STATE", "No",
        "Enable / disable the coupling of the restart configuration.", yesnotuple, yesnovalue,
        &beam_to_solid_volume_mestying);
  }

  // Beam to solid volume mesh tying output parameters.
  Teuchos::ParameterList& beam_to_solid_volume_mestying_vtk =
      beam_to_solid_volume_mestying.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write vtp output at all for btsvmt.
    setStringToIntegralParameter<int>("WRITE_OUTPUT", "No",
        "Enable / disable beam-to-solid volume mesh tying output.", yesnotuple, yesnovalue,
        &beam_to_solid_volume_mestying_vtk);

    setStringToIntegralParameter<int>("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        yesnotuple, yesnovalue, &beam_to_solid_volume_mestying_vtk);

    setStringToIntegralParameter<int>("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        yesnotuple, yesnovalue, &beam_to_solid_volume_mestying_vtk);

    setStringToIntegralParameter<int>("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        yesnotuple, yesnovalue, &beam_to_solid_volume_mestying_vtk);

    DRT::INPUT::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_volume_mestying_vtk);

    setStringToIntegralParameter<int>("SEGMENTATION", "No",
        "Enable / disable output of segmentation points.", yesnotuple, yesnovalue,
        &beam_to_solid_volume_mestying_vtk);

    setStringToIntegralParameter<int>("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        yesnotuple, yesnovalue, &beam_to_solid_volume_mestying_vtk);
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
        tuple<std::string>("none", "configurations_forced_to_zero",
            "configurations_forced_to_zero_fad", "displacements", "displacements_fad"),
        tuple<BeamToSolidSurfaceCoupling>(BeamToSolidSurfaceCoupling::none,
            BeamToSolidSurfaceCoupling::configurations_forced_to_zero,
            BeamToSolidSurfaceCoupling::configurations_forced_to_zero_fad,
            BeamToSolidSurfaceCoupling::displacements,
            BeamToSolidSurfaceCoupling::displacements_fad),
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

    // Add the geometry pair input parameters.
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_surface_mestying);
  }

  // Beam to solid surface parameters.
  Teuchos::ParameterList& beam_to_solid_surface =
      beaminteraction.sublist("BEAM TO SOLID SURFACE", false, "");

  // Beam to solid surface output parameters.
  Teuchos::ParameterList& beam_to_solid_surface_vtk =
      beam_to_solid_surface.sublist("RUNTIME VTK OUTPUT", false, "");
  {
    // Whether to write vtp output at all.
    setStringToIntegralParameter<int>("WRITE_OUTPUT", "No",
        "Enable / disable beam-to-solid volume mesh tying output.", yesnotuple, yesnovalue,
        &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("NODAL_FORCES", "No",
        "Enable / disable output of the resulting nodal forces due to beam to solid interaction.",
        yesnotuple, yesnovalue, &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("AVERAGED_NORMALS", "No",
        "Enable / disable output of averaged nodal normals on the surface.", yesnotuple, yesnovalue,
        &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("MORTAR_LAMBDA_DISCRET", "No",
        "Enable / disable output of the discrete Lagrange multipliers at the node of the Lagrange "
        "multiplier shape functions.",
        yesnotuple, yesnovalue, &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("MORTAR_LAMBDA_CONTINUOUS", "No",
        "Enable / disable output of the continuous Lagrange multipliers function along the beam.",
        yesnotuple, yesnovalue, &beam_to_solid_surface_vtk);

    DRT::INPUT::IntParameter("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS", 5,
        "Number of segments for continuous mortar output", &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("SEGMENTATION", "No",
        "Enable / disable output of segmentation points.", yesnotuple, yesnovalue,
        &beam_to_solid_surface_vtk);

    setStringToIntegralParameter<int>("INTEGRATION_POINTS", "No",
        "Enable / disable output of used integration points. If the contact method has 'forces' at "
        "the integration point, they will also be output.",
        yesnotuple, yesnovalue, &beam_to_solid_surface_vtk);
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
