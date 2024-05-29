/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_INTERSECTION_SERVICE_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_INTERSECTION_SERVICE_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_discretization_geometry_geo_utils.hpp"
#include "4C_discretization_geometry_intersection_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CORE::GEO
{
  std::map<int, CORE::LINALG::Matrix<3, 2>> getCurrentXAABBs(const DRT::Discretization& dis,
      const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions);

  /*!
  \brief checks if two 18 Dops are intersecting (note : for efficiency it only checks slabs
         which are not present for XAABBs)

  \param cutterDOP   (in)         : DOP of the cutting element
  \param xfemDOP     (in)         : DOP of the xfem element
  \return true if the DOP's intersect or false otherwise
   */
  bool intersectionOfKDOPs(
      const CORE::LINALG::Matrix<9, 2>& cutterDOP, const CORE::LINALG::Matrix<9, 2>& xfemDOP);

  /*!
  \brief checks the intersection between two bounding volumes (AABB)
  \param currentBV   (in)         : AABB of the current element
  \param queryBV     (in)         : AABB of the query element
  \return true if the AABB's intersect or false otherwise
   */
  bool intersectionOfBVs(
      const CORE::LINALG::Matrix<3, 2>& currentBV, const CORE::LINALG::Matrix<3, 2>& queryBV);

  /*!
  \brief checks the overlap of two intervals in one coordinate
  \param smin     (in)         : minimum value of the current interval
  \param smax     (in)         : maximum value of the current interval
  \param omin     (in)         : minimum value of the query interval
  \param omax     (in)         : maximum value of the query interval
  \return true if the intervals's overlap or false otherwise
   */
  bool overlap(double smin, double smax, double omin, double omax);

  /*!
  \brief checks if an element is Cartesian, linear or higherorder

  \param element        (in)         : element
  \param xyze_element   (in)         : coordinates of the element
  \param eleGeoType     (out)        : element geometric type CARTESIAN LINEAR or HIGHERORDER
   */
  void checkGeoType(const CORE::Elements::Element* element,
      const CORE::LINALG::SerialDenseMatrix& xyze_element, EleGeoType& eleGeoType);

}  // namespace CORE::GEO


FOUR_C_NAMESPACE_CLOSE

#endif
