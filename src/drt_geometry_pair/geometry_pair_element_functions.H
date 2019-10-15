/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs elements.

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

  /**
   * \brief Evaluate the normal at a position on the surface.
   *
   * @template surface Class of surface element.
   * @template scalar_type Scalar type used for the degrees of freedom of the element.
   * @template scalar_type_evaluate Scalar type used for the parameter coordinates, this will be
   * the result type of this function.
   *
   * @param q_surface (in) Degrees of freedom for the surface.
   * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
   * parameter coordiantes, the third one is in the normal direction).
   * @param normal (out) Normal on the surface.
   * @param element (in) Pointer to the element object.
   * @param nodal_normals (in) Optional - Normals on the nodes.
   */
  template <typename surface, typename scalar_type, typename scalar_type_evaluate>
  void EvaluateSurfaceNormal(const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
      const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
      LINALG::Matrix<3, 1, scalar_type_evaluate>& normal, const DRT::Element* element,
      const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals = NULL)
  {
    // Check at compile time if a surface (2D) element is given.
    static_assert(surface::dim_ == 2, "EvaluateSurfaceNormal can only be called for 2D elements!");

    if (nodal_normals == NULL)
    {
      // Calculate the normal as the geometrical normal on the element.
      LINALG::Matrix<3, 2, scalar_type_evaluate> dr;
      LINALG::Matrix<3, 1, scalar_type_evaluate> dr_0;
      LINALG::Matrix<3, 1, scalar_type_evaluate> dr_1;
      GEOMETRYPAIR::EvaluatePositionDerivative1<surface>(xi, q_surface, dr, element);
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
      {
        dr_0(i_dir) = dr(i_dir, 0);
        dr_1(i_dir) = dr(i_dir, 1);
      }
      normal.CrossProduct(dr_0, dr_1);
      normal.Scale(1.0 / FADUTILS::VectorNorm(normal));
    }
    else
    {
      // Calculate the normal as a interpolation of nodal normals.
      GEOMETRYPAIR::EvaluatePosition<surface>(xi, *nodal_normals, normal, element);
      normal.Scale(1.0 / FADUTILS::VectorNorm(normal));
    }
  }

  /**
   * \brief Evaluate a position on the surface.
   *
   * @template surface Class of surface element.
   * @template scalar_type Scalar type used for the degrees of freedom of the element.
   * @template scalar_type_evaluate Scalar type used for the parameter coordinates, this will be
   * the result type of this function.
   *
   * @param q_surface (in) Degrees of freedom for the surface.
   * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
   * parameter coordiantes, the third one is in the normal direction).
   * @param r (out) Position on the surface.
   * @param element (in) Pointer to the element object.
   * @param nodal_normals (in) Optional - Normals on the nodes.
   */
  template <typename surface, typename scalar_type, typename scalar_type_evaluate>
  void EvaluateSurfacePosition(const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
      const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
      LINALG::Matrix<3, 1, scalar_type_evaluate>& r, const DRT::Element* element,
      const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals = NULL)
  {
    // Check at compile time if a surface (2D) element is given.
    static_assert(surface::dim_ == 2, "EvaluateSurfaceNormal can only be called for 2D elements!");

    // Evaluate the normal.
    LINALG::Matrix<3, 1, scalar_type_evaluate> normal;
    EvaluateSurfaceNormal<surface>(q_surface, xi, normal, element, nodal_normals);

    // Evaluate the position on the surface.
    GEOMETRYPAIR::EvaluatePosition<surface>(xi, q_surface, r, element);

    // Add the normal part to the position.
    normal.Scale(xi(2));
    r += normal;
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 1D element.
   * @param xi (in) Parameter coordinate on the line.
   * @return True if xi is in a valid range for the element.
   */
  template <typename T>
  bool ValidParameter1D(const T& xi)
  {
    const double xi_limit = 1.0 + CONSTANTS::projection_xi_eta_tol;
    if (fabs(xi) < xi_limit)
      return true;
    else
      return false;
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 2D element.
   * @param xi (in) Parameter coordinate on the line.
   * @return True if xi is in a valid range for the element.
   */
  template <typename element_type, typename T>
  bool ValidParameter2D(const T& xi)
  {
    const double xi_limit = 1.0 + CONSTANTS::projection_xi_eta_tol;
    if (element_type::surface_type_ == DiscretizationTypeSurface::quad)
    {
      if (fabs(xi(0)) < xi_limit && fabs(xi(1)) < xi_limit) return true;
    }
    else if (element_type::surface_type_ == DiscretizationTypeSurface::triangle)
    {
      if (xi(0) > -CONSTANTS::projection_xi_eta_tol && xi(1) > -CONSTANTS::projection_xi_eta_tol &&
          xi(0) + xi(1) < 1.0 + CONSTANTS::projection_xi_eta_tol)
        return true;
    }
    else
    {
      dserror("Wrong surface type given!");
    }

    // Default value.
    return false;
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 3D element.
   * @param xi (in) Parameter coordinate on the line.
   * @return True if xi is in a valid range for the element.
   */
  template <typename element_type, typename T>
  bool ValidParameter3D(const T& xi)
  {
    const double xi_limit = 1.0 + CONSTANTS::projection_xi_eta_tol;
    if (element_type::volume_type_ == DiscretizationTypeVolume::hexaeder)
    {
      if (fabs(xi(0)) < xi_limit && fabs(xi(1)) < xi_limit && fabs(xi(2)) < xi_limit) return true;
    }
    else if (element_type::volume_type_ == DiscretizationTypeVolume::tetraeder)
    {
      if (xi(0) > -CONSTANTS::projection_xi_eta_tol && xi(1) > -CONSTANTS::projection_xi_eta_tol &&
          xi(2) > -CONSTANTS::projection_xi_eta_tol &&
          xi(0) + xi(1) + xi(2) < 1.0 + CONSTANTS::projection_xi_eta_tol)
        return true;
    }
    else
    {
      dserror("Wrong volume type given!");
    }

    // Default value.
    return false;
  }

  /**
   * \brief Set the parameter coordinate to the middle of an 1D element.
   * @param xi(in/out) Parameter coordinate.
   */
  template <typename T>
  void SetStartValuesElement1D(T& xi)
  {
    xi = 0.;
  }

  /**
   * \brief Set the parameter coordinate to the middle of a 2D element.
   * @param xi(in/out) Parameter coordinates.
   */
  template <typename element_type, typename T>
  void SetStartValuesElement2D(T& xi)
  {
    if (element_type::surface_type_ == DiscretizationTypeSurface::quad)
      xi.PutScalar(0.0);
    else
      dserror("Wrong volume type given!");
  }

  /**
   * \brief Set the parameter coordinate to the middle of a 2D element.
   * @param xi(in/out) Parameter coordinates.
   */
  template <typename element_type, typename T>
  void SetStartValuesElement3D(T& xi)
  {
    if (element_type::volume_type_ == DiscretizationTypeVolume::hexaeder)
      xi.PutScalar(0.0);
    else if (element_type::volume_type_ == DiscretizationTypeVolume::tetraeder)
      xi.PutScalar(0.25);
    else
      dserror("Wrong element type given!");
  }
}  // namespace GEOMETRYPAIR


#endif
