/*!
\file geometry_pair_utility_classes.H

\brief Utility classes for the geometry pairs.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#ifndef GEOMETRY_PAIR_UTILITY_CLASSES_H
#define GEOMETRY_PAIR_UTILITY_CLASSES_H


#include "../drt_lib/drt_dserror.H"

#include <vector>


// Forward declarations.
namespace LINALG
{
  template <class value_type, unsigned int rows, unsigned int cols>
  class TMatrix;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Type of solid element.
   */
  enum class DiscretizationTypeVolume
  {
    none,      //!< none
    hexaeder,  //!< hexaeder
    tetraeder  //!< tetraeder
  };


  /**
   * \brief Result of a projection with the geometry pairs.
   */
  enum class ProjectionResult
  {
    none,                        //! Default value
    projection_not_found,        //! System of equations could not be solved.
    projection_found_not_valid,  //! Projection found, but the parameter coordinates are not all
                                 //! valid.
    projection_found_valid       //! Projection found and the parameter coordinates are valid.
  };


  /**
   * \brief Class that represents a projection from a line to a volume.
   */
  template <typename scalar_type>
  class ProjectionPointLineToVolume
  {
   public:
    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on line.
     * @param xi Parameter coordinates in volume.
     * @param gauss_weight_ Gauss weight for this point.
     */
    ProjectionPointLineToVolume(
        scalar_type eta, LINALG::TMatrix<scalar_type, 3, 1> xi, double gauss_weight)
        : eta_(eta), xi_(xi), gauss_weight_(gauss_weight){};

    /**
     * \brief Constructor.
     * @param eta Parameter coordinate on line.
     * @param xi Parameter coordinates in volume.
     */
    ProjectionPointLineToVolume(scalar_type eta, LINALG::TMatrix<scalar_type, 3, 1> xi)
        : eta_(eta), xi_(xi), gauss_weight_(-1.){};

    /**
     * \brief Get the parameter coordinate on the line.
     */
    inline scalar_type GetEta() const { return eta_; };

    /**
     * \brief Get the parameter coordinates in the volume.
     */
    inline LINALG::TMatrix<scalar_type, 3, 1> GetXi() const { return xi_; };

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

   private:
    //! Parameter coordinate on line.
    scalar_type eta_;

    //! Parameter coordinates in volume.
    LINALG::TMatrix<scalar_type, 3, 1> xi_;

    //! Gauss weight for this point.
    double gauss_weight_;
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
        : eta_a_(eta_a), eta_b_(eta_b), segment_gauss_points_(){};

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
     * \brief Add a Gauss point to the Gauss point vector.
     * @param gauss_point Gauss point to add to the end of the vector.
     */
    inline void AddGaussPoint(ProjectionPointLineToVolume<scalar_type> gauss_point)
    {
      segment_gauss_points_.push_back(gauss_point);
    }

    /**
     * \brief Return the number of Gauss points in this segment.
     * @return Number of Gauss points.
     */
    inline unsigned int GetNumerOfGaussPoints() const { return segment_gauss_points_.size(); }

    /**
     * \brief Return a const reference to the Gauss points in this segment.
     * @return Reference to Gauss point vector.
     */
    inline const std::vector<ProjectionPointLineToVolume<scalar_type>>& GetGaussPoints() const
    {
      return segment_gauss_points_;
    }

    /**
     * \brief Return a mutable reference to the Gauss points in this segment.
     * @return Reference to Gauss point vector.
     */
    inline std::vector<ProjectionPointLineToVolume<scalar_type>>& GetGaussPointsMutable()
    {
      return segment_gauss_points_;
    }

   private:
    //! Lower limit of segment.
    scalar_type eta_a_;

    //! Upper limit of segment.
    scalar_type eta_b_;

    //! Vector to store Gauss points for this segment.
    std::vector<ProjectionPointLineToVolume<scalar_type>> segment_gauss_points_;
  };

}  // namespace GEOMETRYPAIR


#endif
