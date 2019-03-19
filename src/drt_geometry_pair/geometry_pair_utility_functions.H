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
   * \brief Get the in a line element but with the line interpreted as a cylinder. The parameter
   * space for this element is [-1, 1] in all directions, but it is a circle in xi2 and xi2
   * direction. The radius is taken from the element pointer.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param r (out) Position on the element.
   * @param element (in) Pointer to the element, needed for beam elements.
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
      dserror("The tangent vector can not have a component in z direction!");

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
