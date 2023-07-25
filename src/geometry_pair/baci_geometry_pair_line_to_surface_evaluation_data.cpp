/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
*/


#include "baci_geometry_pair_line_to_surface_evaluation_data.H"

#include "baci_geometry_pair_utility_classes.H"
#include "baci_geometry_pair_element_faces.H"

#include "baci_beaminteraction_str_model_evaluator_datastate.H"


/**
 *
 */
GEOMETRYPAIR::LineToSurfaceEvaluationData::LineToSurfaceEvaluationData(
    const Teuchos::ParameterList& input_parameter_list)
    : LineTo3DEvaluationData(input_parameter_list),
      face_elements_(),
      surface_normal_strategy_(INPAR::GEOMETRYPAIR::SurfaceNormals::standard)
{
  surface_normal_strategy_ = Teuchos::getIntegralValue<INPAR::GEOMETRYPAIR::SurfaceNormals>(
      input_parameter_list, "GEOMETRY_PAIR_SURFACE_NORMALS");
}

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
    const Teuchos::RCP<const Epetra_Vector>& displacement_col_np)
{
  for (const auto& [id, face_element] : face_elements_)
    if (face_element->IsPartOfPair()) face_element->SetState(displacement_col_np, face_elements_);
}
