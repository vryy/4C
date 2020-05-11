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
#include "beam_to_solid_surface_meshtying_pair_gauss_point.H"
#include "beam_to_solid_surface_meshtying_pair_gauss_point_FAD.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "beam_to_solid_surface_meshtying_params.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
#include "../drt_so3/so_base.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidCondition::BeamToSolidCondition(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamInteractionConditionBase(condition_line),
      condition_other_(condition_other),
      condition_contact_pairs_()
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
void BEAMINTERACTION::BeamToSolidCondition::Clear()
{
  BeamInteractionConditionBase::Clear();
  condition_contact_pairs_.clear();
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

  if (contact_pair != Teuchos::null)
  {
    // Create the geometry pair on the beam contact pair.
    contact_pair->CreateGeometryPair(geometry_evaluation_data_);

    // Add to the internal vector which keeps track of the created contact pairs.
    condition_contact_pairs_.push_back(contact_pair);
  }
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
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
      new GEOMETRYPAIR::LineToSurfaceEvaluationData(input_parameter_list));
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
void BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  // Call the parent method.
  BeamToSolidCondition::Setup(discret);

  // Pointer to the beam contact parameters.
  Teuchos::RCP<const BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>
      beal_to_solid_surface_params = Teuchos::null;
  if (condition_contact_pairs_.size() > 0)
    beal_to_solid_surface_params =
        condition_contact_pairs_[0]->Params()->BeamToSolidSurfaceMeshtyingParams();

  // Loop over all pairs and add the needed face elements.
  std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> pair_face_elemets;
  pair_face_elemets.clear();
  for (const auto& pair : condition_contact_pairs_)
  {
    const int solid_id = pair->Element2()->Id();
    auto find_in_condition = surface_ids_.find(solid_id);
    if (find_in_condition != surface_ids_.end())
    {
      // Check if the face is already in the pair_face_elemets map.
      auto find_in_pair = pair_face_elemets.find(solid_id);
      if (find_in_pair == pair_face_elemets.end())
      {
        // The face element has to be created and added to the contact pair.
        Teuchos::RCP<GEOMETRYPAIR::FaceElement> new_face_element = GEOMETRYPAIR::FaceElementFactory(
            find_in_condition->second, beal_to_solid_surface_params->GetIsFAD());
        new_face_element->SetPartOfPair(true);
        pair_face_elemets[solid_id] = new_face_element;
        pair->SetFaceElement(new_face_element);
      }
      else
      {
        // Add the existing face element to the contact pair.
        pair->SetFaceElement(find_in_pair->second);
      }
    }
    else
    {
      dserror("The face of the solid element %d is not in the current condition!",
          pair->Element2()->Id());
    }
  }

  // Now all faces of contact pairs are in pair_face_elemets, we still need to add faces that are
  // needed for averaged normal calculation, but are not contained in any pair.
  std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> face_elements_needed;
  face_elements_needed = pair_face_elemets;
  for (const auto& face_element_iterator : pair_face_elemets)
  {
    // Loop over the nodes of the face element.
    const DRT::Node* const* nodes = face_element_iterator.second->GetDrtFaceElement()->Nodes();
    for (int i_node = 0; i_node < face_element_iterator.second->GetDrtFaceElement()->NumNode();
         i_node++)
    {
      // Loop over the elements connected to that node and check if they are in this condition.
      const DRT::Element* const* elements = nodes[i_node]->Elements();
      for (int i_element = 0; i_element < nodes[i_node]->NumElement(); i_element++)
      {
        const int element_id = elements[i_element]->Id();
        auto find_in_condition = surface_ids_.find(element_id);
        if (find_in_condition != surface_ids_.end())
        {
          // The element exists in this condition, check if it is already in the needed faces map.
          auto find_in_needed = face_elements_needed.find(element_id);
          if (find_in_needed == face_elements_needed.end())
          {
            // It is not already in the needed faces -> add it.
            face_elements_needed[element_id] = GEOMETRYPAIR::FaceElementFactory(
                find_in_condition->second, beal_to_solid_surface_params->GetIsFAD());
          }
        }
        else
        {
          // The element is not part of this condition, i.e. it will not be used for the
          // calculation of averaged normals. This allows for 'sharp' corners.
        }
      }
    }
  }

  // Cast the geometry evaluation data to the correct type.
  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
          geometry_evaluation_data_, false);
  if (line_to_surface_evaluation_data == Teuchos::null)
    dserror("Could not cast to GEOMETRYPAIR::LineToSurfaceEvaluationData.");

  // Setup the geometry data for the surface patch.
  line_to_surface_evaluation_data->Setup(discret, face_elements_needed);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::SetState(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  // Cast the geometry evaluation data to the correct type.
  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
          geometry_evaluation_data_, false);
  if (line_to_surface_evaluation_data == Teuchos::null)
    dserror("Could not cast to GEOMETRYPAIR::LineToSurfaceEvaluationData.");

  // Setup the geometry data for the surface patch.
  line_to_surface_evaluation_data->SetState(beaminteraction_data_state);
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::CreateContactPairInternal(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  const Teuchos::RCP<const DRT::FaceElement>& face_element = surface_ids_[ele_ptrs[1]->Id()];
  const DRT::Element::DiscretizationType shape = face_element->Shape();

  INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type =
      params_ptr->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType();
  switch (coupling_type)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::configurations_forced_to_zero:
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacements:
    {
      switch (shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>());
        case DRT::Element::tri6:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>());
        case DRT::Element::quad4:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>());
        case DRT::Element::quad8:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>());
        case DRT::Element::quad9:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }

    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::configurations_forced_to_zero_fad:
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacements_fad:
    {
      switch (shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
              GEOMETRYPAIR::t_tri3>());
        case DRT::Element::tri6:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
              GEOMETRYPAIR::t_tri6>());
        case DRT::Element::quad4:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
              GEOMETRYPAIR::t_quad4>());
        case DRT::Element::quad8:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
              GEOMETRYPAIR::t_quad8>());
        case DRT::Element::quad9:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
              GEOMETRYPAIR::t_quad9>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
              Sacado::ELRFad::SLFad<
                  Sacado::ELRFad::SLFad<double,
                      GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
                  GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
              GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }

    default:
      dserror("Wrong coupling type.");
  }

  return Teuchos::null;
}
