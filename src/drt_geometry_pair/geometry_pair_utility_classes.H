/*----------------------------------------------------------------------*/
/*! \file

\brief Utility classes for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_UTILITY_CLASSES_H
#define GEOMETRY_PAIR_UTILITY_CLASSES_H


#include "geometry_pair_constants.H"

#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../drt_lib/drt_dserror.H"
#include "../headers/FAD_utils.H"

#include <vector>


namespace GEOMETRYPAIR
{
  /**
   * \brief Result of a projection with the geometry pairs.
   */
  enum class ProjectionResult
  {
    //! Default value
    none,
    //! System of equations could not be solved.
    projection_not_found,
    //! Projection found, but the parameter coordinates are not all valid.
    projection_found_not_valid,
    //! Projection found and the parameter coordinates are valid.
    projection_found_valid
  };

  /**
   * \brief Class that represents a projection from a 1D structure (usually a line) to a 3D
   * structure (can be volume, as well as surface including normal direction).
   */
  template <typename scalar_type>
  class ProjectionPoint1DTo3D
  {
   public:
    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on line.
     * @param xi Parameter coordinates in volume.
     * @param gauss_weight Gauss weight for this point.
     */
    ProjectionPoint1DTo3D(
        scalar_type eta, LINALG::Matrix<3, 1, scalar_type> xi, double gauss_weight)
        : eta_(eta),
          xi_(xi),
          projection_result_(ProjectionResult::none),
          gauss_weight_(gauss_weight),
          eta_cross_section_(true),
          is_cross_section_point_(false){};

    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on line.
     * @param xi Parameter coordinates in volume.
     */
    ProjectionPoint1DTo3D(scalar_type eta, LINALG::Matrix<3, 1, scalar_type> xi)
        : ProjectionPoint1DTo3D(eta, xi, -1.){};

    /**
     * \brief Destructor.
     */
    virtual ~ProjectionPoint1DTo3D(){};

    /**
     * \brief Set the parameter coordinate on the line.
     */
    inline void SetEta(const scalar_type& eta) { eta_ = eta; };

    /**
     * \brief Get the parameter coordinate on the line.
     */
    inline scalar_type GetEta() const { return eta_; };

    /**
     * \brief Set the parameter coordinates in the volume.
     */
    inline void SetXi(const LINALG::Matrix<3, 1, scalar_type>& xi) { xi_ = xi; };

    /**
     * \brief Get the parameter coordinates in the volume.
     */
    inline LINALG::Matrix<3, 1, scalar_type> GetXi() const { return xi_; };

    /**
     * \brief Get the parameter coordinates in the volume as a reference.
     */
    inline LINALG::Matrix<3, 1, scalar_type>& GetXiMutable() { return xi_; };

    /**
     * \brief Set the parameter coordinates in the cross section.
     */
    inline void SetEtaCrossSection(LINALG::Matrix<2, 1, scalar_type> eta_cross_section)
    {
      eta_cross_section_ = eta_cross_section;
      is_cross_section_point_ = true;
    };

    /**
     * \brief Get the parameter coordinates in the cross section.
     */
    inline LINALG::Matrix<2, 1, scalar_type> GetEtaCrossSection() const
    {
      if (!is_cross_section_point_) dserror("The cross section coordinate has not been set!");
      return eta_cross_section_;
    };

    /**
     * Set the projection result for this projection point.
     * @param projection_result
     */
    inline void SetProjectionResult(ProjectionResult projection_result)
    {
      projection_result_ = projection_result;
    }

    /**
     * Get the projection result for this projection point.
     */
    inline ProjectionResult GetProjectionResult() const { return projection_result_; }

    /**
     * Get the projection result for this projection point as a reference.
     */
    inline ProjectionResult& GetProjectionResultMutable() { return projection_result_; }

    /**
     * \brief Set the Gauss weight for this point.
     * @param gauss_weight
     */
    inline void SetGaussWeight(double gauss_weight) { gauss_weight_ = gauss_weight; };

    /**
     * \brief Get the Gauss weight for this point, if none is defined, an error is thrown.
     */
    inline double GetGaussWeight() const
    {
      if (gauss_weight_ < 0.)
        dserror(
            "Negative Gauss weight not possible. Probably the default value was not overwritten!");
      return gauss_weight_;
    };

    /**
     * \brief Overloaded $<$ operator.
     * @param lhs
     * @param rhs
     * @return True if smaller, false if larger.
     */
    friend bool operator<(const ProjectionPoint1DTo3D<scalar_type>& lhs,
        const ProjectionPoint1DTo3D<scalar_type>& rhs)
    {
      if (lhs.GetEta() < rhs.GetEta() - CONSTANTS::projection_xi_eta_tol)
        return true;
      else
        return false;
    };

    /**
     * \brief Overloaded $>$ operator.
     * @param lhs
     * @param rhs
     * @return False if smaller, true if larger.
     */
    friend bool operator>(const ProjectionPoint1DTo3D<scalar_type>& lhs,
        const ProjectionPoint1DTo3D<scalar_type>& rhs)
    {
      if (lhs.GetEta() > rhs.GetEta() + CONSTANTS::projection_xi_eta_tol)
        return true;
      else
        return false;
    };

   private:
    //! Parameter coordinate on line.
    scalar_type eta_;

    //! Parameter coordinates in volume.
    LINALG::Matrix<3, 1, scalar_type> xi_;

    //! Projection result.
    ProjectionResult projection_result_;

    //! Gauss weight for this point.
    double gauss_weight_;

    //! Parameter coordinates in the cross section.
    LINALG::Matrix<2, 1, scalar_type> eta_cross_section_;

    //! Flag if this is a point on a cross section.
    bool is_cross_section_point_;
  };

  /**
   * \brief Class to manage a segment on a line.
   */
  template <typename scalar_type>
  class LineSegment
  {
   public:
    /**
     * \brief Default constructor.
     */
    LineSegment() : LineSegment(-1., 1.){};

    /**
     * \brief Constructor. Set the range of the segment.
     * @param eta_a
     * @param eta_b
     */
    LineSegment(scalar_type eta_a, scalar_type eta_b)
        : eta_a_(eta_a), eta_b_(eta_b), segment_projection_points_()
    {
      // Sanity check that eta_a is larger than eta_b.
      if (!(eta_a_ < eta_b_))
        dserror(
            "The segment is created with eta_a=%f and eta_b=%f, this is not possible, as eta_a has "
            "to be smaller than eta_b!",
            eta_a_, eta_b);
    }

    /**
     * \brief Get the length of the segment in parameter coordinates.
     * @return Segment length.
     */
    inline scalar_type GetSegmentLength() const { return eta_b_ - eta_a_; }

    /**
     * \brief Return a const reference to eta_a_.
     * @return eta_a_
     */
    inline const scalar_type& GetEtaA() const { return eta_a_; };

    /**
     * \brief Return a const reference to eta_b_.
     * @return eta_a_
     */
    inline const scalar_type& GetEtaB() const { return eta_b_; };

    /**
     * \brief Add a projection point to the projection point vector.
     * @param projection_point projection point to add to the end of the vector.
     */
    inline void AddProjectionPoint(ProjectionPoint1DTo3D<scalar_type> projection_point)
    {
      segment_projection_points_.push_back(projection_point);
    }

    /**
     * \brief Return the number of projection points in this segment.
     * @return Number of projection points.
     */
    inline unsigned int GetNumberOfProjectionPoints() const
    {
      return segment_projection_points_.size();
    }

    /**
     * \brief Return a const reference to the projection points in this segment.
     * @return Reference to projection point vector.
     */
    inline const std::vector<ProjectionPoint1DTo3D<scalar_type>>& GetProjectionPoints() const
    {
      return segment_projection_points_;
    }

    /**
     * \brief Return a mutable reference to the projection points in this segment.
     * @return Reference to projection point vector.
     */
    inline std::vector<ProjectionPoint1DTo3D<scalar_type>>& GetProjectionPointsMutable()
    {
      return segment_projection_points_;
    }

    /**
     * \brief Overloaded $<$ operator.
     * @param lhs
     * @param rhs
     * @return True if smaller, false if larger. Throw error if the segments are overlapping.
     */
    friend bool operator<(const LineSegment<scalar_type>& lhs, const LineSegment<scalar_type>& rhs)
    {
      if (lhs.GetEtaB() < rhs.GetEtaA() + CONSTANTS::projection_xi_eta_tol)
        return true;
      else if (lhs.GetEtaA() > rhs.GetEtaB() - CONSTANTS::projection_xi_eta_tol)
      {
        // The segments do not overlap.
      }
      else if (abs(lhs.GetEtaA() - rhs.GetEtaA()) < CONSTANTS::projection_xi_eta_tol &&
               abs(lhs.GetEtaB() - rhs.GetEtaB()) < CONSTANTS::projection_xi_eta_tol)
      {
        // The segments are equal.
      }
      else
        dserror("The two segments are overlapping. This is fatal!");

      return false;
    };

    /**
     * \brief Overloaded $>$ operator.
     * @param lhs
     * @param rhs
     * @return False if smaller, true if larger. Throw error if the segments are overlapping.
     */
    friend bool operator>(const LineSegment<scalar_type>& lhs, const LineSegment<scalar_type>& rhs)
    {
      if (lhs.GetEtaA() > rhs.GetEtaB() - CONSTANTS::projection_xi_eta_tol)
        return true;
      else if (lhs.GetEtaB() < rhs.GetEtaA() + CONSTANTS::projection_xi_eta_tol)
      {
        // The segments do not overlap.
      }
      else if (abs(lhs.GetEtaA() - rhs.GetEtaA()) < CONSTANTS::projection_xi_eta_tol &&
               abs(lhs.GetEtaB() - rhs.GetEtaB()) < CONSTANTS::projection_xi_eta_tol)
      {
        // The segments are equal.
      }
      else
        dserror("The two segments are overlapping. This is fatal!");

      return false;
    };

   private:
    //! Lower limit of segment.
    scalar_type eta_a_;

    //! Upper limit of segment.
    scalar_type eta_b_;

    //! Vector to store projection points for this segment.
    std::vector<ProjectionPoint1DTo3D<scalar_type>> segment_projection_points_;
  };
}  // namespace GEOMETRYPAIR


#endif
