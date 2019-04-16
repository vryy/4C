/*!
\file geometry_pair_utility_functions.H

\brief Utility functions for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_UTILITY_FUNCTIONS_H
#define GEOMETRY_PAIR_UTILITY_FUNCTIONS_H


#include "../linalg/linalg_utils.H"

#include "geometry_pair_element_types.H"
#include "../drt_beam3/beam3.H"
#include "../headers/FAD_utils.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Get the spatial position in a volume defined by a centerline and a cross section.
   *
   * The parameter coordinates in the cross section are from -1 t 1 in eta2 and eta3 direction. The
   * radius is taken from the element pointer. The current implementation does not take rotations
   * along the centerline into account. Therefore following triad along the centerline is assumed:
   *   - director 1: Tangent along the centerline.
   *   - director 2: Normal to director 1 in x-y plane.
   *   - director 3: Normal to director 1 and director 2.
   * This works fine as long as the line only deforms in the x-y plane. Other deformations will
   * cause a runtime error.
   *
   * @param eta (in) Parameter coordinate along the centerline.
   * @param eta_cross_section (in) Parameter coordinates in the cross section. They will be
   * multiplied with the radius of the cylindrical surface on the centerline.
   * @param q (in) DOF for of the line.
   * @param r (out) Spatial position in the cross section.
   * @param element (in) Pointer to the line element.
   */
  template <typename element_type, typename T1, typename T2, typename V, typename scalar_type>
  void EvaluatePositionLineCrossSection(const T1& eta, const T2& eta_cross_section, const V& q,
      LINALG::TMatrix<scalar_type, 3, 1>& r, const DRT::Element* element = NULL)
  {
    // Get the radius from the beam element.
    const double radius = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element))
                              ->GetCircularCrossSectionRadiusForInteractions();

    // Get the position and tangent along the centerline.
    EvaluatePosition<element_type>(eta, q, r, element);
    LINALG::TMatrix<scalar_type, 3, 1> tangent, cross_section_director_2, cross_section_director_3;
    EvaluatePositionDerivative1<element_type>(eta, q, tangent, element);

    // The current implementation only allows tangential vectors in the x-y plane.
    if (abs(tangent(2)) > CONSTANTS::pos_tol)
      dserror("The tangent vector can not have a component in z direction! The component is %f!",
          FADUTILS::CastToDouble(tangent(2)));

    // Create the director vectors in the cross section.
    // Director 2 is the one in the y-axis (reference configuration).
    // Director 3 is the one in the z-axis (reference configuration).
    tangent.Scale(1. / FADUTILS::VectorNorm(tangent));
    cross_section_director_2.Clear();
    cross_section_director_2(0) = -tangent(1);
    cross_section_director_2(1) = tangent(0);
    cross_section_director_3.Clear();
    cross_section_director_3(2) = 1.;

    // Add the in cross section components of the position.
    for (unsigned int dir = 0; dir < 3; dir++)
      r(dir) += radius * (eta_cross_section(0) * cross_section_director_2(dir) +
                             eta_cross_section(1) * cross_section_director_3(dir));
  }
}  // namespace GEOMETRYPAIR


#endif
