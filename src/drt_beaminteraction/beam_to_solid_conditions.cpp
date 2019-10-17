/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_conditions.H"

#include "beam_contact_params.H"
#include "beam_contact_pair.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point.H"
#include "beam_to_solid_volume_meshtying_pair_mortar.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section.H"
#include "beam_to_solid_volume_meshtying_params.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_so3/so_base.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidCondition::BeamToSolidCondition(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamInteractionConditionBase(condition_line), condition_other_(condition_other)
{
}

/**
 *
 */
bool BEAMINTERACTION::BeamToSolidCondition::IdsInCondition(
    const int id_line, const int id_other) const
{
  if (line_ids_.find(id_line) != line_ids_.end())
    if (IdInOther(id_other)) return true;
  return false;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BeamToSolidConditionVolumeMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
  // Create the geometry evaluation data for this condition.
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData());
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  // Check if the given elements are in this condition.
  if (!IdsInCondition(ele_ptrs[0]->Id(), ele_ptrs[1]->Id())) return Teuchos::null;

  // Cast the solid element.
  DRT::ELEMENTS::So_base const* solidele = dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]);
  DRT::Element::DiscretizationType shape = solidele->Shape();

  // Get the contact discretization method.
  INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization contact_discretization =
      params_ptr->BeamToSolidVolumeMeshtyingParams()->GetContactDiscretization();

  // Check which contact discretization is wanted.
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> new_pair = Teuchos::null;
  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
        break;
      case DRT::Element::hex20:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
        break;
      case DRT::Element::hex27:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
        break;
      case DRT::Element::tet4:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
        break;
      case DRT::Element::tet10:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
        break;
      case DRT::Element::nurbs27:
        new_pair = Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_nurbs27>());
        break;

      default:
        dserror("Wrong element type for solid element.");
    }
  }
  else if (contact_discretization ==
           INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::mortar)
  {
    INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions mortar_shape_function =
        params_ptr->BeamToSolidVolumeMeshtyingParams()->GetMortarShapeFunctionType();

    switch (mortar_shape_function)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line2:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>());
            break;
          case DRT::Element::hex20:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>());
            break;
          case DRT::Element::hex27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>());
            break;
          case DRT::Element::tet4:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>());
            break;
          case DRT::Element::tet10:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>());
            break;
          case DRT::Element::nurbs27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line2>());
            break;
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line3:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>());
            break;
          case DRT::Element::hex20:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>());
            break;
          case DRT::Element::hex27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>());
            break;
          case DRT::Element::tet4:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>());
            break;
          case DRT::Element::tet10:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>());
            break;
          case DRT::Element::nurbs27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line3>());
            break;
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line4:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line4>());
            break;
          case DRT::Element::hex20:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line4>());
            break;
          case DRT::Element::hex27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line4>());
            break;
          case DRT::Element::tet4:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line4>());
            break;
          case DRT::Element::tet10:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line4>());
            break;
          case DRT::Element::nurbs27:
            new_pair = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line4>());
            break;
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      default:
        dserror("Wrong mortar shape function.");
    }
  }
  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::gauss_point_cross_section)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        new_pair =
            Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>());
        break;
      case DRT::Element::hex20:
        new_pair =
            Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>());
        break;
      case DRT::Element::hex27:
        new_pair =
            Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>());
        break;
      case DRT::Element::tet4:
        new_pair =
            Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>());
        break;
      case DRT::Element::tet10:
        new_pair =
            Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>());
        break;
      default:
        dserror("Wrong element type for solid element.");
    }
  }

  // Create the geometry pair on the beam contact pair.
  if (new_pair != Teuchos::null)
    new_pair->CreateGeometryPair(geometry_evaluation_data_);
  else
    dserror("No pair was created!");

  // Return the newly created pair.
  return new_pair;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BuildIdSets()
{
  // Call the parent method to build the line maps.
  BeamToSolidCondition::BuildIdSets();

  // Build the volume map.
  std::vector<int> volume_ids;
  ConditionToElementIds(condition_other_, volume_ids);
  volume_ids_ = std::set<int>(volume_ids.begin(), volume_ids.end());
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::BeamToSolidConditionSurfaceMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
  // Create the geometry evaluation data for this condition.
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData());
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  dserror("Not yet implemented");
  return Teuchos::null;
}
