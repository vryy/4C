/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_GLOBAL_H
#define GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_GLOBAL_H


#include "geometry_pair_line_to_3D_evaluation_data.H"

#include <Teuchos_RCP.hpp>
#include <unordered_map>


// Forward declarations.
namespace STR
{
  namespace MODELEVALUATOR
  {
    class BeamInteractionDataState;
  }
}  // namespace STR
namespace GEOMETRYPAIR
{
  class FaceElement;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Class to manage input parameters and evaluation data for line to surface interactions.
   */
  class LineToSurfaceEvaluationData : public LineTo3DEvaluationData
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineToSurfaceEvaluationData(const Teuchos::ParameterList& input_parameter_list)
        : LineTo3DEvaluationData(input_parameter_list), face_elements_(){};

    /**
     * \brief Destructor.
     */
    ~LineToSurfaceEvaluationData() {}

    /**
     * \brief Reset the evaluation data (derived).
     */
    void Clear() override;

    /**
     * \brief Setup the surface data.
     *
     * \param discret (in) Pointer to the discretization.
     * \param face_elements (in) Map to all face elements in this condition on this rank.
     */
    void Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements);

    /**
     * \brief Calculate the averaged nodal normals.
     */
    void SetState(const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
            beaminteraction_data_state);

    /**
     * \brief Get a reference to the face element map.
     */
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& GetFaceElements() const
    {
      return face_elements_;
    }

   private:
    //! A map of all face elements needed for this surface.
    std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> face_elements_;
  };
}  // namespace GEOMETRYPAIR

#endif
