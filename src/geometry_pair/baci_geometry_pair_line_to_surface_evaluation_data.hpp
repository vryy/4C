/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_HPP


#include "baci_config.hpp"

#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"

#include <Teuchos_RCP.hpp>

#include <unordered_map>

// Forward declarations.
class Epetra_Vector;

BACI_NAMESPACE_OPEN

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
    LineToSurfaceEvaluationData(const Teuchos::ParameterList& input_parameter_list);

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
    void SetState(const Teuchos::RCP<const Epetra_Vector>& displacement_col_np);

    /**
     * \brief Get a reference to the face element map.
     */
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& GetFaceElements() const
    {
      return face_elements_;
    }

    /**
     * \brief Return the strategy to be used for the surface normals.
     */
    INPAR::GEOMETRYPAIR::SurfaceNormals GetSurfaceNormalStrategy() const
    {
      return surface_normal_strategy_;
    }

   private:
    //! A map of all face elements needed for this surface.
    std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> face_elements_;

    //! Strategy to be used for surface normals.
    INPAR::GEOMETRYPAIR::SurfaceNormals surface_normal_strategy_;
  };
}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE

#endif
