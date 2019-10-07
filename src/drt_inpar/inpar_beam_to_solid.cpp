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
DRT::UTILS::GaussRule1D INPAR::BEAMTOSOLID::IntToGaussRule1D(const int n_gauss_points)
{
  switch (n_gauss_points)
  {
    case 1:
      return DRT::UTILS::GaussRule1D::intrule_line_1point;
    case 2:
      return DRT::UTILS::GaussRule1D::intrule_line_2point;
    case 3:
      return DRT::UTILS::GaussRule1D::intrule_line_3point;
    case 4:
      return DRT::UTILS::GaussRule1D::intrule_line_4point;
    case 5:
      return DRT::UTILS::GaussRule1D::intrule_line_5point;
    case 6:
      return DRT::UTILS::GaussRule1D::intrule_line_6point;
    case 7:
      return DRT::UTILS::GaussRule1D::intrule_line_7point;
    case 8:
      return DRT::UTILS::GaussRule1D::intrule_line_8point;
    case 9:
      return DRT::UTILS::GaussRule1D::intrule_line_9point;
    case 10:
      return DRT::UTILS::GaussRule1D::intrule_line_10point;
    case 16:
      return DRT::UTILS::GaussRule1D::intrule_line_16point;
    case 20:
      return DRT::UTILS::GaussRule1D::intrule_line_20point;
    case 32:
      return DRT::UTILS::GaussRule1D::intrule_line_32point;
    case 50:
      return DRT::UTILS::GaussRule1D::intrule_line_50point;
    default:
    {
      dserror("No Gauss rule defined for %d points", n_gauss_points);
      return DRT::UTILS::GaussRule1D::intrule1D_undefined;
    }
  }
};

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
    setStringToIntegralParameter<BeamToSolidVolumeContactDiscretization>("CONTACT_DISCRETIZATION",
        "none", "Type of employed contact discretization",
        tuple<std::string>("none", "gauss_point_to_segment", "mortar", "gauss_point_cross_section"),
        tuple<BeamToSolidVolumeContactDiscretization>(BeamToSolidVolumeContactDiscretization::none,
            BeamToSolidVolumeContactDiscretization::gauss_point_to_segment,
            BeamToSolidVolumeContactDiscretization::mortar,
            BeamToSolidVolumeContactDiscretization::gauss_point_cross_section),
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidVolumeConstraintEnforcement>("CONSTRAINT_STRATEGY",
        "none", "Type of employed constraint enforcement strategy",
        tuple<std::string>("none", "penalty"),
        tuple<BeamToSolidVolumeConstraintEnforcement>(BeamToSolidVolumeConstraintEnforcement::none,
            BeamToSolidVolumeConstraintEnforcement::penalty),
        &beam_to_solid_volume_mestying);

    setStringToIntegralParameter<BeamToSolidVolumeMortarShapefunctions>("MORTAR_SHAPE_FUNCTION",
        "none", "Shape function for the mortar Lagrange-multiplicators",
        tuple<std::string>("none", "line2", "line3", "line4"),
        tuple<BeamToSolidVolumeMortarShapefunctions>(BeamToSolidVolumeMortarShapefunctions::none,
            BeamToSolidVolumeMortarShapefunctions::line2,
            BeamToSolidVolumeMortarShapefunctions::line3,
            BeamToSolidVolumeMortarShapefunctions::line4),
        &beam_to_solid_volume_mestying);

    DoubleParameter("PENALTY_PARAMETER", 0.0,
        "Penalty parameter for beam-to-solid volume meshtying", &beam_to_solid_volume_mestying);

    IntParameter("GAUSS_POINTS", 6, "Number of Gauss Points for the integral evaluations",
        &beam_to_solid_volume_mestying);

    IntParameter("INTEGRATION_POINTS_CIRCUMFENCE", 6,
        "Number of Integration points along the circumfencial direction of the beam. This is "
        "parameter is only used in beam to cylinder meshtying. No gauss integration is "
        "used along the circumfencial direction, equally spaced integration points are used.",
        &beam_to_solid_volume_mestying);

    // Add the geometry pair input parameters.
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_solid_volume_mestying);
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
