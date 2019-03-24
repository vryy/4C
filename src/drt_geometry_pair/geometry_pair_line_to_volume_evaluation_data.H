/*!
\file geometry_pair_line_to_volume_evaluation_data.H

\brief container for parameters for line to volume pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_EVALUATION_DATA_GLOBAL_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_EVALUATION_DATA_GLOBAL_H


#include "geometry_pair_evaluation_data_base.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_inpar/inpar_geometry_pair.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class to manage inout parameters and evaluation data for line to volume interactions.
   */
  class LineToVolumeEvaluationData : public GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor.
     */
    LineToVolumeEvaluationData();

    /**
     * \brief Destructor.
     */
    ~LineToVolumeEvaluationData(){};

    /**
     * \brief Initialize parameters. Load them from the input file.
     */
    void Init();

    /**
     * \brief Setup function.
     */
    void Setup();

    /**
     * \brief Get the segmentation strategy.
     * @return flag for segmentation strategy.
     */
    inline INPAR::GEOMETRYPAIR::LineToVolumeStrategy GetStrategy() const { return strategy_; }

    /**
     * \brief Get the number of search points for segmentation search.
     * @return number of search points.
     */
    inline unsigned int GetNumberOfSearchPoints() const { return n_search_points_; }

    /**
     * \brief Set the Gauss rule to be used for Gauss point projection method.
     * @param gauss_rule
     */
    void SetGaussRule(DRT::UTILS::GaussRule1D gauss_rule) { gauss_rule_ = gauss_rule; };

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
     * \brief Sets the number of integration points along the circumfence of the line cross
     * section.
     */
    inline void SetNumberOfIntegrationPointsCircumfence(
        const unsigned int integration_points_circumfence)
    {
      integration_points_circumfence_ = integration_points_circumfence;
    }

    /**
     * \brief Returns the number of integration points along the circumfence of the line cross
     * section.
     * @return number of integration points in cirumfencial direction.
     */
    inline unsigned int GetNumberOfIntegrationPointsCircumfence() const
    {
      return integration_points_circumfence_;
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
    INPAR::GEOMETRYPAIR::LineToVolumeStrategy strategy_;

    //! Gauss rule for Gauss point projection method.
    DRT::UTILS::GaussRule1D gauss_rule_;

    //! Number of integration points in the circumfencial direction of the line cross section.
    unsigned int integration_points_circumfence_;

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
