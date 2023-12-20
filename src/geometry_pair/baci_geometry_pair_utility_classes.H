/*----------------------------------------------------------------------*/
/*! \file

\brief Utility classes for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_UTILITY_CLASSES_H
#define BACI_GEOMETRY_PAIR_UTILITY_CLASSES_H


#include "baci_config.H"

#include "baci_geometry_pair_constants.H"
#include "baci_linalg_utils_densematrix_inverse.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_fad.H"

#include <vector>

BACI_NAMESPACE_OPEN

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
   * @tparam scalar_type Scalar type of the parameter coordinate values.
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
        scalar_type eta, CORE::LINALG::Matrix<3, 1, scalar_type> xi, double gauss_weight)
        : eta_(eta),
          xi_(xi),
          projection_result_(ProjectionResult::none),
          gauss_weight_(gauss_weight),
          intersection_face_(-1),
          eta_cross_section_(true),
          is_cross_section_point_(false){};

    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on line.
     * @param xi Parameter coordinates in volume.
     */
    ProjectionPoint1DTo3D(scalar_type eta, CORE::LINALG::Matrix<3, 1, scalar_type> xi)
        : ProjectionPoint1DTo3D(eta, xi, -1.){};

    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on the line.
     */
    ProjectionPoint1DTo3D(scalar_type eta)
        : ProjectionPoint1DTo3D(eta, CORE::LINALG::Matrix<3, 1, scalar_type>(true), -1.){};

    /**
     * \brief Empty constructor.
     */
    ProjectionPoint1DTo3D() : ProjectionPoint1DTo3D(0.0){};

    /**
     * \brief Destructor.
     */
    virtual ~ProjectionPoint1DTo3D() = default;

    /**
     * \brief Construct the point from another point where all scalar values are cast to double.
     * @param point_double Projection point with scalar type double.
     */
    template <typename scalar_type_other>
    inline void SetFromOtherPointDouble(const ProjectionPoint1DTo3D<scalar_type_other>& point_other)
    {
      eta_ = CORE::FADUTILS::CastToDouble(point_other.GetEta());
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        xi_(i_dim) = CORE::FADUTILS::CastToDouble(point_other.GetXi()(i_dim));
      projection_result_ = point_other.GetProjectionResult();
      gauss_weight_ = point_other.GetGaussWeightNoCheck();
      intersection_face_ = point_other.GetIntersectionFace();
    }

    /**
     * \brief Set the parameter coordinate on the line.
     */
    inline void SetEta(const scalar_type& eta) { eta_ = eta; };

    /**
     * \brief Get the parameter coordinate on the line.
     */
    inline const scalar_type& GetEta() const { return eta_; };

    /**
     * \brief Get a mutable reference to the parameter coordinate on the line.
     */
    inline scalar_type& GetEta() { return eta_; };

    /**
     * \brief Set the parameter coordinates in the volume.
     */
    inline void SetXi(const CORE::LINALG::Matrix<3, 1, scalar_type>& xi) { xi_ = xi; };

    /**
     * \brief Get the parameter coordinates in the volume.
     */
    inline const CORE::LINALG::Matrix<3, 1, scalar_type>& GetXi() const { return xi_; };

    /**
     * \brief Get the parameter coordinates in the volume as a reference.
     */
    inline CORE::LINALG::Matrix<3, 1, scalar_type>& GetXi() { return xi_; };

    /**
     * \brief Set the parameter coordinates in the cross section.
     */
    inline void SetEtaCrossSection(CORE::LINALG::Matrix<2, 1, scalar_type> eta_cross_section)
    {
      eta_cross_section_ = eta_cross_section;
      is_cross_section_point_ = true;
    };

    /**
     * \brief Get the parameter coordinates in the cross section.
     */
    inline CORE::LINALG::Matrix<2, 1, scalar_type> GetEtaCrossSection() const
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
    inline ProjectionResult& GetProjectionResult() { return projection_result_; }

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
    }

    /**
     * \brief Get the Gauss weight for this point.
     */
    inline double GetGaussWeightNoCheck() const { return gauss_weight_; }

    /**
     * \brief Set the index of the intersection face.
     */
    inline void SetIntersectionFace(const int intersection_face)
    {
      intersection_face_ = intersection_face;
    }

    /**
     * \brief Get the index of the intersection face.
     */
    inline int GetIntersectionFace() const { return intersection_face_; }

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
    CORE::LINALG::Matrix<3, 1, scalar_type> xi_;

    //! Projection result.
    ProjectionResult projection_result_;

    //! Gauss weight for this point.
    double gauss_weight_;

    //! If this point is an intersection point, this is the index of the local face that the
    //! intersection occurs on.
    int intersection_face_;

    //! Parameter coordinates in the cross section.
    CORE::LINALG::Matrix<2, 1, scalar_type> eta_cross_section_;

    //! Flag if this is a point on a cross section.
    bool is_cross_section_point_;
  };

  /**
   * \brief Class to manage a segment on a line.
   * @tparam scalar_type Scalar type of the parameter coordinate values.
   */
  template <typename scalar_type>
  class LineSegment
  {
   public:
    /**
     * \brief Default constructor, the segment is from -1 to 1.
     */
    LineSegment() : LineSegment(scalar_type(-1.0), scalar_type(1.0)){};

    /**
     * \brief Constructor. Set the range of the segment.
     * @param start_point
     * @param end_point
     */
    LineSegment(ProjectionPoint1DTo3D<scalar_type> start_point,
        ProjectionPoint1DTo3D<scalar_type> end_point)
        : start_point_(start_point), end_point_(end_point), segment_projection_points_()
    {
      // Sanity check that eta_a is larger than eta_b.
      if (!(GetEtaA() < GetEtaB()))
        dserror(
            "The segment is created with eta_a=%f and eta_b=%f, this is not possible, as eta_a "
            "has "
            "to be smaller than eta_b!",
            CORE::FADUTILS::CastToDouble(GetEtaA()), CORE::FADUTILS::CastToDouble(GetEtaB()));
    }

    /**
     * \brief Get the length of the segment in parameter coordinates.
     * @return Segment length.
     */
    inline scalar_type GetSegmentLength() const { return GetEtaB() - GetEtaA(); }

    /**
     * \brief Return a const reference to eta start.
     */
    inline const scalar_type& GetEtaA() const { return start_point_.GetEta(); };

    /**
     * \brief Return a const reference to eta end.
     */
    inline const scalar_type& GetEtaB() const { return end_point_.GetEta(); };

    /**
     * \brief Return a const reference to the start point.
     */
    inline const ProjectionPoint1DTo3D<scalar_type>& GetStartPoint() const { return start_point_; };

    /**
     * \brief Return a mutable reference to the start point.
     */
    inline ProjectionPoint1DTo3D<scalar_type>& GetStartPoint() { return start_point_; };

    /**
     * \brief Return a const reference to the end point.
     */
    inline const ProjectionPoint1DTo3D<scalar_type>& GetEndPoint() const { return end_point_; };

    /**
     * \brief Return a mutable reference to the end point.
     */
    inline ProjectionPoint1DTo3D<scalar_type>& GetEndPoint() { return end_point_; };

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
    inline std::vector<ProjectionPoint1DTo3D<scalar_type>>& GetProjectionPoints()
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
    //! Start point of the segment.
    ProjectionPoint1DTo3D<scalar_type> start_point_;

    //! Endpoint of the segment.
    ProjectionPoint1DTo3D<scalar_type> end_point_;

    //! Vector to store projection points for this segment.
    std::vector<ProjectionPoint1DTo3D<scalar_type>> segment_projection_points_;
  };
}  // namespace GEOMETRYPAIR


BACI_NAMESPACE_CLOSE

#endif
