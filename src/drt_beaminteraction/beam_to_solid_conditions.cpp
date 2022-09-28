/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
*/


#include "beam_to_solid_conditions.H"

#include "beam_contact_params.H"
#include "beam_contact_pair.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point.H"
#include "beam_to_solid_volume_meshtying_pair_mortar.H"
#include "beam_to_solid_volume_meshtying_pair_mortar_rotation.H"
#include "beam_to_solid_volume_meshtying_pair_2d-3d_plane.H"
#include "beam_to_solid_volume_meshtying_pair_2d-3d_full.H"
#include "beam_to_solid_surface_contact_pair.H"
#include "beam_to_solid_surface_meshtying_pair_gauss_point.H"
#include "beam_to_solid_surface_meshtying_pair_gauss_point_FAD.H"
#include "beam_to_solid_surface_meshtying_pair_mortar.H"
#include "beam_to_solid_surface_meshtying_pair_mortar_FAD.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "beam_to_solid_surface_contact_params.H"
#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H"
#include "str_model_evaluator_beaminteraction_datastate.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
#include "../drt_so3/so_base.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidCondition::BeamToSolidCondition(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other,
    const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params)
    : BeamInteractionConditionBase(condition_line),
      geometry_evaluation_data_(Teuchos::null),
      condition_other_(condition_other),
      condition_contact_pairs_(),
      beam_to_solid_params_(beam_to_solid_params)
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
  geometry_evaluation_data_->Clear();
  condition_contact_pairs_.clear();
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidCondition::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs)
{
  // Check if the given elements are in this condition.
  if (!IdsInCondition(ele_ptrs[0]->Id(), ele_ptrs[1]->Id())) return Teuchos::null;

  // Create the beam contact pair.
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> contact_pair = CreateContactPairInternal(ele_ptrs);

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
Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManager>
BEAMINTERACTION::BeamToSolidCondition::CreateIndirectAssemblyManager(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  if (beam_to_solid_params_->GetContactDiscretization() ==
      INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::mortar)
    return Teuchos::rcp(new SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect(
        condition_contact_pairs_, discret, beam_to_solid_params_));
  return Teuchos::null;
}


/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BeamToSolidConditionVolumeMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other,
    const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params)
    : BeamToSolidCondition(condition_line, condition_other, beam_to_solid_params)
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
void BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BuildIdSets(
    const Teuchos::RCP<const DRT::Discretization>& discretization)
{
  // Call the parent method to build the line maps.
  BeamToSolidCondition::BuildIdSets(discretization);

  // Build the volume map.
  std::vector<int> volume_ids;
  ConditionToElementIds(condition_other_, volume_ids);
  volume_ids_ = std::set<int>(volume_ids.begin(), volume_ids.end());
}

/**
 *
 */
template <template <typename...> class bts_class, typename... bts_template_arguments>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BEAMINTERACTION::CreateBeamToSolidVolumePairShape(
    const DRT::Element::DiscretizationType shape)
{
  switch (shape)
  {
    case DRT::Element::hex8:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8,
          bts_template_arguments...>());
    case DRT::Element::hex20:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20,
          bts_template_arguments...>());
    case DRT::Element::hex27:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27,
          bts_template_arguments...>());
    case DRT::Element::tet4:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4,
          bts_template_arguments...>());
    case DRT::Element::tet10:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10,
          bts_template_arguments...>());
    case DRT::Element::nurbs27:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27,
          bts_template_arguments...>());
    default:
      dserror("Wrong element type for solid element.");
      return Teuchos::null;
  }
}

/**
 *
 */
template <template <typename...> class bts_class, typename... bts_template_arguments>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::CreateBeamToSolidVolumePairShapeNoNurbs(
    const DRT::Element::DiscretizationType shape)
{
  switch (shape)
  {
    case DRT::Element::hex8:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8,
          bts_template_arguments...>());
    case DRT::Element::hex20:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20,
          bts_template_arguments...>());
    case DRT::Element::hex27:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27,
          bts_template_arguments...>());
    case DRT::Element::tet4:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4,
          bts_template_arguments...>());
    case DRT::Element::tet10:
      return Teuchos::rcp(new bts_class<GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10,
          bts_template_arguments...>());
    default:
      dserror("Wrong element type for solid element.");
      return Teuchos::null;
  }
}

/**
 *
 */
template <template <typename...> class bts_class, typename... bts_mortar_template_arguments,
    typename... bts_mortar_shape>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BEAMINTERACTION::CreateBeamToSolidVolumePairMortar(
    const DRT::Element::DiscretizationType shape,
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function,
    bts_mortar_shape... other_mortar_shape_function)
{
  switch (mortar_shape_function)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
      return CreateBeamToSolidVolumePairMortar<bts_class, bts_mortar_template_arguments...,
          GEOMETRYPAIR::t_line2>(shape, other_mortar_shape_function...);
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
      return CreateBeamToSolidVolumePairMortar<bts_class, bts_mortar_template_arguments...,
          GEOMETRYPAIR::t_line3>(shape, other_mortar_shape_function...);
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
      return CreateBeamToSolidVolumePairMortar<bts_class, bts_mortar_template_arguments...,
          GEOMETRYPAIR::t_line4>(shape, other_mortar_shape_function...);
    default:
      dserror("Wrong mortar shape function.");
      return Teuchos::null;
  }
}

/**
 *
 */
template <template <typename...> class bts_class, typename... bts_mortar_template_arguments>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BEAMINTERACTION::CreateBeamToSolidVolumePairMortar(
    const DRT::Element::DiscretizationType shape)
{
  return CreateBeamToSolidVolumePairShape<bts_class, bts_mortar_template_arguments...>(shape);
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::CreateContactPairInternal(
    const std::vector<DRT::Element const*>& ele_ptrs)
{
  const DRT::ELEMENTS::So_base* solidele = dynamic_cast<const DRT::ELEMENTS::So_base*>(ele_ptrs[1]);
  const DRT::Element::DiscretizationType shape = solidele->Shape();
  const auto beam_to_volume_params =
      Teuchos::rcp_dynamic_cast<const BeamToSolidVolumeMeshtyingParams>(
          beam_to_solid_params_, true);
  const INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization contact_discretization =
      beam_to_volume_params->GetContactDiscretization();

  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_to_segment)
  {
    // Create the Gauss point to segment pairs.
    return CreateBeamToSolidVolumePairShape<BeamToSolidVolumeMeshtyingPairGaussPoint>(shape);
  }
  else if (contact_discretization == INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::mortar)
  {
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function =
        beam_to_volume_params->GetMortarShapeFunctionType();
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shape_function_rotation =
        beam_to_volume_params->GetMortarShapeFunctionRotationType();

    if (mortar_shape_function_rotation == INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::none)
    {
      // Create the positional mortar pairs.
      return CreateBeamToSolidVolumePairMortar<BeamToSolidVolumeMeshtyingPairMortar>(
          shape, mortar_shape_function);
    }
    else
    {
      // Create the rotational mortart pairs.
      return CreateBeamToSolidVolumePairMortar<BeamToSolidVolumeMeshtyingPairMortarRotation>(
          shape, mortar_shape_function, mortar_shape_function_rotation);
    }
  }
  else if (contact_discretization ==
           INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_cross_section)
  {
    // Depending on the type of beam element we create the correct beam-to-solid pair here.
    const auto sr_beam = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(ele_ptrs[0]);
    const auto eb_beam = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(ele_ptrs[0]);
    if (sr_beam != nullptr)
      return CreateBeamToSolidVolumePairShapeNoNurbs<BeamToSolidVolumeMeshtyingPair2D3DFull>(shape);
    else if (eb_beam != nullptr)
      return CreateBeamToSolidVolumePairShapeNoNurbs<BeamToSolidVolumeMeshtyingPair2D3DPlane>(
          shape);
    else
      dserror(
          "2D-3D coupling is only implemented for Simo-Reissner and torsion free beam elements.");
  }

  // Default return value.
  return Teuchos::null;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionSurface::BeamToSolidConditionSurface(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other,
    const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params, const bool is_mesh_tying)
    : BeamToSolidCondition(condition_line, condition_other, beam_to_solid_params),
      is_mesh_tying_(is_mesh_tying)
{
  // Get the input parameter list that will be passed to the geometry pair.
  std::string condition_name;
  if (IsMeshTying())
    condition_name = "BEAM TO SOLID SURFACE MESHTYING";
  else
    condition_name = "BEAM TO SOLID SURFACE CONTACT";

  const Teuchos::ParameterList& input_parameter_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist(condition_name);

  // Create the geometry evaluation data for this condition.
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
      new GEOMETRYPAIR::LineToSurfaceEvaluationData(input_parameter_list));
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionSurface::BuildIdSets(
    const Teuchos::RCP<const DRT::Discretization>& discretization)
{
  // Call the parent method to build the line maps.
  BeamToSolidCondition::BuildIdSets(discretization);

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
void BEAMINTERACTION::BeamToSolidConditionSurface::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  // Call the parent method.
  BeamToSolidCondition::Setup(discret);

  // Cast the geometry evaluation data to the correct type.
  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
          geometry_evaluation_data_, false);
  if (line_to_surface_evaluation_data == Teuchos::null)
    dserror("Could not cast to GEOMETRYPAIR::LineToSurfaceEvaluationData.");

  // If the pairs are FAD, i.e., if the averaged normals have to be evaluated using FAD.
  int fad_order = 0;
  if (condition_contact_pairs_.size() > 0)
  {
    if (IsMeshTying())
    {
      fad_order =
          condition_contact_pairs_[0]->Params()->BeamToSolidSurfaceMeshtyingParams()->GetFADOrder();
    }
    else
    {
      fad_order =
          condition_contact_pairs_[0]->Params()->BeamToSolidSurfaceContactParams()->GetFADOrder();
    }
  }

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
        Teuchos::RCP<GEOMETRYPAIR::FaceElement> new_face_element =
            GEOMETRYPAIR::FaceElementFactory(find_in_condition->second, fad_order,
                line_to_surface_evaluation_data->GetSurfaceNormalStrategy());
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
            face_elements_needed[element_id] =
                GEOMETRYPAIR::FaceElementFactory(find_in_condition->second, fad_order,
                    line_to_surface_evaluation_data->GetSurfaceNormalStrategy());
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

  // Setup the geometry data for the surface patch.
  line_to_surface_evaluation_data->Setup(discret, face_elements_needed);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionSurface::SetState(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  // For contact we reset the evaluation data in each iteration (we don't call Clear() here, since
  // we want to keep the contact pairs).
  if (IsContact())
  {
    auto line_to_other_evaluation_data =
        Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineTo3DEvaluationData>(
            geometry_evaluation_data_, true);
    line_to_other_evaluation_data->ResetTracker();
  }

  // Cast the geometry evaluation data to the correct type.
  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
          geometry_evaluation_data_, false);
  if (line_to_surface_evaluation_data == Teuchos::null)
    dserror("Could not cast to GEOMETRYPAIR::LineToSurfaceEvaluationData.");

  // Setup the geometry data for the surface patch.
  line_to_surface_evaluation_data->SetState(beaminteraction_data_state->GetDisColNp());
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionSurface::CreateContactPairInternal(
    const std::vector<DRT::Element const*>& ele_ptrs)
{
  using namespace GEOMETRYPAIR;

  const Teuchos::RCP<const DRT::FaceElement>& face_element = surface_ids_[ele_ptrs[1]->Id()];
  const DRT::Element::DiscretizationType shape = face_element->Shape();

  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineToSurfaceEvaluationData>(
          geometry_evaluation_data_, false);
  if (line_to_surface_evaluation_data == Teuchos::null)
    dserror("Could not cast to GEOMETRYPAIR::LineToSurfaceEvaluationData.");
  auto surface_normal_strategy = line_to_surface_evaluation_data->GetSurfaceNormalStrategy();

  if (IsMeshTying())
  {
    // Create beam-to-surface pairs for mesh tying.
    auto beam_to_surface_params =
        Teuchos::rcp_dynamic_cast<const BeamToSolidSurfaceMeshtyingParams>(
            beam_to_solid_params_, true);

    INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type =
        beam_to_surface_params->GetCouplingType();

    INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization coupling_discretization =
        beam_to_surface_params->GetContactDiscretization();

    switch (coupling_discretization)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_to_segment:
      {
        switch (coupling_type)
        {
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
              reference_configuration_forced_to_zero:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement:
          {
            switch (shape)
            {
              case DRT::Element::tri3:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_tri3>());
              case DRT::Element::tri6:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_tri6>());
              case DRT::Element::quad4:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad4>());
              case DRT::Element::quad8:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad8>());
              case DRT::Element::quad9:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad9>());
              case DRT::Element::nurbs9:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_nurbs9>());
              default:
                dserror("Wrong element type for surface element.");
            }
            break;
          }

          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
              reference_configuration_forced_to_zero_fad:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad:
          {
            if (surface_normal_strategy == INPAR::GEOMETRYPAIR::SurfaceNormals::standard)
            {
              switch (shape)
              {
                case DRT::Element::tri3:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type, t_hermite, t_tri3>());
                case DRT::Element::tri6:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type, t_hermite, t_tri6>());
                case DRT::Element::quad4:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type, t_hermite, t_quad4>());
                case DRT::Element::quad8:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type, t_hermite, t_quad8>());
                case DRT::Element::quad9:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type, t_hermite, t_quad9>());
                case DRT::Element::nurbs9:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite,
                      t_nurbs9>());
                default:
                  dserror("Wrong element type for surface element.");
              }
              break;
            }
            else if (surface_normal_strategy ==
                     INPAR::GEOMETRYPAIR::SurfaceNormals::extended_volume)
            {
              switch (shape)
              {
                case DRT::Element::quad4:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite,
                      t_quad4>());
                case DRT::Element::quad8:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite,
                      t_quad8>());
                case DRT::Element::quad9:
                  return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
                      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite,
                      t_quad9>());
                default:
                  dserror("Wrong element type for surface element.");
              }
            }
          }

          default:
            dserror("Wrong coupling type.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::mortar:
      {
        INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction =
            beam_to_surface_params->GetMortarShapeFunctionType();
        bool rotational_coupling = beam_to_surface_params->GetIsRotationalCoupling();

        switch (coupling_type)
        {
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
              reference_configuration_forced_to_zero:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement:
            return BeamToSolidSurfaceMeshtyingPairMortarFactory(shape, mortar_shapefunction);
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
              reference_configuration_forced_to_zero_fad:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad:
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad:
            return BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
                shape, mortar_shapefunction, rotational_coupling, surface_normal_strategy);
          default:
            dserror("Wrong coupling type.");
        }
        break;
      }
      default:
        dserror("Wrong coupling discretization.");
    }
  }
  else
  {
    // Create beam-to-surface pairs for contact.
    auto beam_to_surface_contact_params =
        Teuchos::rcp_dynamic_cast<const BeamToSolidSurfaceContactParams>(
            beam_to_solid_params_, true);

    INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization contact_discretization =
        beam_to_surface_contact_params->GetContactDiscretization();

    switch (contact_discretization)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::gauss_point_to_segment:
      {
        INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact contact_type =
            beam_to_surface_contact_params->GetContactType();

        switch (contact_type)
        {
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::gap_variation:
            switch (shape)
            {
              case DRT::Element::tri3:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3>());
              case DRT::Element::tri6:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6>());
              case DRT::Element::quad4:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4>());
              case DRT::Element::quad8:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8>());
              case DRT::Element::quad9:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9>());
              case DRT::Element::nurbs9:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairGapVariation<
                    line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>,
                    t_hermite, t_nurbs9>());
              default:
                dserror("Wrong element type for surface element.");
            }
            break;
          case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::potential:
            switch (shape)
            {
              case DRT::Element::tri3:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
                        t_hermite, t_tri3>());
              case DRT::Element::tri6:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
                        t_hermite, t_tri6>());
              case DRT::Element::quad4:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
                        t_hermite, t_quad4>());
              case DRT::Element::quad8:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
                        t_hermite, t_quad8>());
              case DRT::Element::quad9:
                return Teuchos::rcp(
                    new BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
                        t_hermite, t_quad9>());
              case DRT::Element::nurbs9:
                return Teuchos::rcp(new BeamToSolidSurfaceContactPairPotential<
                    line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite,
                    t_nurbs9>());
              default:
                dserror("Wrong element type for surface element.");
            }
            break;
          default:
            dserror("Wrong contact type.");
        }
        break;
      }
      default:
        dserror("Wrong contact discretization.");
    }
  }

  // Default return value. This will cause an error.
  return Teuchos::null;
}
