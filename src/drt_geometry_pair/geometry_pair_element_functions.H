/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_FUNCTIONS_H
#define GEOMETRY_PAIR_ELEMENT_FUNCTIONS_H


#include "../drt_beam3/beam3.H"
#include "../headers/FAD_utils.H"
#include "geometry_pair_element.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Get the position in the element.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param r (out) Position on the element.
   * @param element (in) Pointer to the element, needed for beam elements.
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePosition(const T& xi, const V& q, LINALG::Matrix<3, 1, scalar_type>& r,
      const DRT::Element* element = NULL)
  {
    // Matrix for shape function values.
    LINALG::Matrix<1, element_type::n_nodes_ * element_type::n_val_, scalar_type> N(true);

    // Evaluate the shape function values.
    element_type::EvaluateShapeFunction(
        N, xi, std::integral_constant<unsigned int, element_type::dim_>{}, element);

    // Calculate the position.
    r.Clear();
    for (unsigned int node = 0; node < element_type::n_nodes_; node++)
      for (unsigned int dim = 0; dim < 3; dim++)
        for (unsigned int val = 0; val < element_type::n_val_; val++)
          r(dim) += q(3 * element_type::n_val_ * node + 3 * val + dim) *
                    N(element_type::n_val_ * node + val);
  }

  /**
   * \brief Get the derivative of the position w.r.t xi in the element.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param dr (out) Derivative of the position on the element.
   * @param element (in) Pointer to the element, needed for beam elements.
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePositionDerivative1(const T& xi, const V& q,
      LINALG::Matrix<3, element_type::dim_, scalar_type>& dr, const DRT::Element* element = NULL)
  {
    // Matrix for shape function values.
    LINALG::Matrix<element_type::dim_, element_type::n_nodes_ * element_type::n_val_, scalar_type>
        dN(true);

    // Evaluate the shape function values.
    element_type::EvaluateShapeFunctionDeriv1(
        dN, xi, std::integral_constant<unsigned int, element_type::dim_>{}, element);

    // Calculate the derivative of the position.
    dr.Clear();
    for (unsigned int dim = 0; dim < 3; dim++)
      for (unsigned int direction = 0; direction < element_type::dim_; direction++)
        for (unsigned int node = 0; node < element_type::n_nodes_; node++)
          for (unsigned int val = 0; val < element_type::n_val_; val++)
            dr(dim, direction) += q(3 * element_type::n_val_ * node + 3 * val + dim) *
                                  dN(direction, element_type::n_val_ * node + val);
  }

  /**
   * \brief Get a normal on an element. This only makes sense for a 2D element.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param n (out) Normal of the element.
   * @param element (in) Pointer to the element.
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluateNormal(const T& xi, const V& q, LINALG::Matrix<3, 1, scalar_type>& n,
      const DRT::Element* element = NULL)
  {
    // Calculate derivative of the position on the element.
    LINALG::Matrix<3, element_type::dim_, scalar_type>& dr;
    EvaluatePositionDerivative1<element_type>(xi, q, dr, element);

    // Calculate the cross product between the two direction vectors.
    LINALG::Matrix<3, 1, scalar_type>& dr_0;
    LINALG::Matrix<3, 1, scalar_type>& dr_1;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      dr_0[i_dim] = dr[i_dim, 0];
      dr_1[i_dim] = dr[i_dim, 1];
    }
    n.CrossProduct(dr_0, dr_1);
    n.Scale(1. / FADUTILS::VectorNorm(n));
  }

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
      LINALG::Matrix<3, 1, scalar_type>& r, const DRT::Element* element = NULL)
  {
    // Get the radius from the beam element.
    const double radius = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element))
                              ->GetCircularCrossSectionRadiusForInteractions();

    // Get the position and tangent along the centerline.
    EvaluatePosition<element_type>(eta, q, r, element);
    LINALG::Matrix<3, 1, scalar_type> tangent, cross_section_director_2, cross_section_director_3;
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
