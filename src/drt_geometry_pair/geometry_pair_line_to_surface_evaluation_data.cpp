/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_surface_evaluation_data.H"

#include "geometry_pair_utility_classes.H"
#include "geometry_pair_element_faces.H"

#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"


/**
 *
 */
void GEOMETRYPAIR::LineToSurfaceEvaluationData::Clear()
{
  // Call reset on the base method.
  LineTo3DEvaluationData::Clear();
  face_elements_.clear();
}

/**
 *
 */
void GEOMETRYPAIR::LineToSurfaceEvaluationData::Setup(
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  face_elements_ = face_elements;

  // Setup each face element.
  for (auto& face_element_map_item : face_elements) face_element_map_item.second->Setup();
}

/**
 *
 */
void GEOMETRYPAIR::LineToSurfaceEvaluationData::SetState(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  for (const auto& face_element_iterator : face_elements_)
  {
    // Set the configurations for all faces.
    face_element_iterator.second->SetState(discret, beaminteraction_data_state->GetDisColNp());
  }

  for (const auto& face_element_iterator : face_elements_)
  {
    if (face_element_iterator.second->IsPartOfPair())
      // Calculate the averaged normals on the faces that are contained in pairs.
      face_element_iterator.second->CalculateAveragedNormals(face_elements_);
  }
}
