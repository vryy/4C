/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to 3D pairs, as well as global evaluation data.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_3D_EVALUATION_DATA_GLOBAL_H
#define GEOMETRY_PAIR_LINE_TO_3D_EVALUATION_DATA_GLOBAL_H


#include "geometry_pair_evaluation_data_base.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_inpar/inpar_geometry_pair.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class to manage inout parameters and evaluation data for line to 3D interactions.
   */
  class LineTo3DEvaluationData : public GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineTo3DEvaluationData(const Teuchos::ParameterList& input_parameter_list);

    /**
     * \brief Destructor.
     */
    ~LineTo3DEvaluationData(){};

    /**
     * \brief Clear the evaluation data.
     */
    void Clear() override;

    /**
     * \brief Get the segmentation strategy.
     * @return flag for segmentation strategy.
     */
    inline INPAR::GEOMETRYPAIR::LineTo3DStrategy GetStrategy() const { return strategy_; }

    /**
     * \brief Get the number of search points for segmentation search.
     * @return number of search points.
     */
    inline unsigned int GetNumberOfSearchPoints() const { return n_search_points_; }

    /**
     * \brief Get the Gauss rule to be used for Gauss point projection method.
     * @return Gauss rule
     */
    inline DRT::UTILS::GaussRule1D GetGaussRule() const { return gauss_rule_; };

    /**
     * \brief Get the Gauss rule to be used for Gauss point projection method.
     * @return Gauss rule
     */
    inline DRT::UTILS::IntegrationPoints1D GetGaussPoints() const
    {
      return DRT::UTILS::IntegrationPoints1D(gauss_rule_);
    };

    /**
     * \brief Get the number of Gauss points.
     * @return Gauss rule
     */
    inline int GetNumberOfGaussPoints() const { return GetGaussPoints().nquad; };

    /**
     * \brief Returns the number of integration points along the circumference of the line cross
     * section.
     * @return number of integration points in circumferencial direction.
     */
    inline unsigned int GetNumberOfIntegrationPointsCircumference() const
    {
      return integration_points_circumference_;
    }

    /**
     * \brief Return a reference to the gauss point projection tracker.
     * @return Projection tracker.
     */
    inline std::map<int, std::vector<bool>>& GetGaussPointProjectionTrackerMutable()
    {
      return gauss_point_projection_tracker_;
    };

    /**
     * \brief Return a reference to the segment tracker.
     * @return Segment tracker.
     */
    std::map<int, std::set<LineSegment<double>>>& GetSegmentTrackerMutable()
    {
      return segment_tracker_;
    };

   private:
    //! Strategy to be used for contact search.
    INPAR::GEOMETRYPAIR::LineTo3DStrategy strategy_;

    //! Gauss rule for Gauss point projection method.
    DRT::UTILS::GaussRule1D gauss_rule_;

    //! Number of integration points in the circumferencial direction of the line cross section.
    unsigned int integration_points_circumference_;

    //! Gauss point projection tracking vector.
    std::map<int, std::vector<bool>> gauss_point_projection_tracker_;

    //! Number of points for segmentation search.
    unsigned int n_search_points_;

    //! Segment tracking vector for segmentation. We use double in this case, because otherwise the
    //! class would have to be templated on the type of this tracker.
    std::map<int, std::set<LineSegment<double>>> segment_tracker_;
  };
}  // namespace GEOMETRYPAIR

#endif
