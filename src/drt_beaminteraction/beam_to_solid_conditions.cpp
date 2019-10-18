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
#include "../drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
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
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidCondition::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  // Check if the given elements are in this condition.
  if (!IdsInCondition(ele_ptrs[0]->Id(), ele_ptrs[1]->Id())) return Teuchos::null;

  // Create the beam contact pair.
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> contact_pair =
      CreateContactPairInternal(ele_ptrs, params_ptr);

  // Create the geometry pair on the beam contact pair.
  if (contact_pair != Teuchos::null)
    contact_pair->CreateGeometryPair(geometry_evaluation_data_);
  else
    dserror(
        "No contact pair was created. This is fatal, since we want to create a geometry pair on "
        "the contact pair here.");

  // Return the newly created pair.
  return contact_pair;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BeamToSolidConditionVolumeMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
  // Get the input parameter list that will be passed to the geometry pair.
  const Teuchos::ParameterList& input_parameter_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID VOLUME MESHTYING");

  // Create the geometry evaluation data for this condition.
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData(input_parameter_list));
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
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::CreateContactPairInternal(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  // Cast the solid element.
  const DRT::ELEMENTS::So_base* solidele = dynamic_cast<const DRT::ELEMENTS::So_base*>(ele_ptrs[1]);
  const DRT::Element::DiscretizationType shape = solidele->Shape();

  // Get the contact discretization method.
  INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization contact_discretization =
      params_ptr->BeamToSolidVolumeMeshtyingParams()->GetContactDiscretization();

  // Check which contact discretization is wanted.
  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
      case DRT::Element::nurbs27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_nurbs27>());

      default:
        dserror("Wrong element type for solid element.");
    }
  }
  else if (contact_discretization == INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::mortar)
  {
    INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function =
        params_ptr->BeamToSolidVolumeMeshtyingParams()->GetMortarShapeFunctionType();

    switch (mortar_shape_function)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line2>());
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line3>());
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line4>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line4>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line4>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line4>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line4>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line4>());
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
      INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_cross_section)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>());
      default:
        dserror("Wrong element type for solid element.");
    }
  }

  // Default return value.
  return Teuchos::null;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::BeamToSolidConditionSurfaceMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
  // Get the input parameter list that will be passed to the geometry pair.
  const Teuchos::ParameterList& input_parameter_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID SURFACE MESHTYING");

  // Create the geometry evaluation data for this condition.
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineToSurfaceEvaluationData<double>>(
      new GEOMETRYPAIR::LineToSurfaceEvaluationData<double>(input_parameter_list));
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::BuildIdSets()
{
  // Call the parent method to build the line maps.
  BeamToSolidCondition::BuildIdSets();

  // Build the surface map.
  surface_ids_.clear();
  for (const auto& map_item : condition_other_->Geometry())
  {
    if (!map_item.second->IsFaceElement()) dserror("Expected FaceElement");
    Teuchos::RCP<const DRT::FaceElement> face_element =
        Teuchos::rcp_dynamic_cast<const DRT::FaceElement>(map_item.second);
    const int solid_id = face_element->ParentElementId();
    surface_ids_[solid_id] = face_element;
  }

  // The size of the surface id set and the geometry in the conditions have to match, otherwise
  // there are two faces connected to the same element, which is not implemented at this point.
  if (surface_ids_.size() != condition_other_->Geometry().size())
    dserror(
        "There are multiple faces connected to one solid element in this condition. This case is "
        "currently not implemented.");
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::CreateContactPairInternal(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  dserror("Not yet implemented");
  return Teuchos::null;
}
