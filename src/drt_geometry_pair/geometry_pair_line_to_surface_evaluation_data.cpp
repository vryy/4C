/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
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
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  face_elements_ = face_elements;

  for (const auto& face_element_iterator : face_elements_)
    face_element_iterator.second->Setup(discret, face_elements_);

  // The averaged reference normals have to be calculated after each face element is set up.
  for (const auto& face_element_iterator : face_elements_)
    if (face_element_iterator.second->IsPartOfPair())
      face_element_iterator.second->CalculateAveragedReferenceNormals(face_elements_);
}

/**
 *
 */
void GEOMETRYPAIR::LineToSurfaceEvaluationData::SetState(
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  for (const auto& face_element_iterator : face_elements_)
    if (face_element_iterator.second->IsPartOfPair())
      face_element_iterator.second->SetState(
          beaminteraction_data_state->GetDisColNp(), face_elements_);
}
