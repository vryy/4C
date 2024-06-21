/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_fiber_node_utils.hpp"
#include "4C_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_nurbs_discretization_utils.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_mat_so3_material.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype, typename Enable = void>
  struct ElementNodes;
}  // namespace Discret::ELEMENTS

namespace Discret::ELEMENTS::DETAIL
{
  template <Core::FE::CellType celltype>
  inline static constexpr int num_nodes = Core::FE::num_nodes<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dim = Core::FE::dim<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_str = num_dim<celltype>*(num_dim<celltype> + 1) / 2;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dof_per_ele = num_nodes<celltype>* num_dim<celltype>;
}  // namespace Discret::ELEMENTS::DETAIL

namespace Discret::ELEMENTS
{

  /*!
   * @brief Calculate the lumped mass matrix
   */
  inline void LumpMatrix(Core::LinAlg::SerialDenseMatrix& matrix)
  {
    FOUR_C_ASSERT(
        matrix.numRows() == matrix.numCols(), "The provided mass matrix is not a square matrix!");

    // we assume mass is a square matrix
    for (int c = 0; c < matrix.numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < matrix.numRows(); ++r)  // parse rows
      {
        d += matrix(r, c);  // accumulate row entries
        matrix(r, c) = 0.0;
      }
      matrix(c, c) = d;  // apply sum of row entries on diagonal
    }
  }

  /*!
   * @brief A type holding information of the nodes of an element,
   * such as their reference and current position. Additional information
   * is stored for NURBS elements
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype, typename Enable>
  struct ElementNodes
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>
        reference_coordinates_;

    /*!
     * @brief Displacements of the element nodes
     */
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>> displacements_;
  };

  template <Core::FE::CellType celltype>
  struct ElementNodes<celltype, std::enable_if_t<Core::FE::is_nurbs<celltype>>>
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>
        reference_coordinates_;

    /*!
     * @brief Displacements of the element nodes
     */
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>> displacements_;

    /*!
     * @brief Knot span of a NURBS element
     */
    std::vector<Core::LinAlg::SerialDenseVector> knots_;

    /*!
     * @brief Weights of control points
     */
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1, double> weights_;
  };

  /*!
   * @brief Evaluate element node information given the element and a displacement vector
   *
   * @tparam celltype
   * @param ele (in): Element
   * @param disp (in) : Vector of nodal displacements of the element
   * @return ElementNodes<celltype>
   */
  template <Core::FE::CellType celltype>
  ElementNodes<celltype> evaluate_element_nodes(
      const Core::Elements::Element& ele, const std::vector<double>& disp)
  {
    Discret::ELEMENTS::ElementNodes<celltype> element_nodes;
    for (int i = 0; i < DETAIL::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
      {
        element_nodes.reference_coordinates_(i, d) = ele.Nodes()[i]->X()[d];
        element_nodes.displacements_(i, d) = disp[i * DETAIL::num_dim<celltype> + d];
      }
    }

    return element_nodes;
  }

  /*!
   * @brief Evaluates the nodal coordinates from this iteration
   *
   * @param ele (in) : Reference to the element
   * @param discretization (in) : discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <Core::FE::CellType celltype>
  ElementNodes<celltype> evaluate_element_nodes(const Core::Elements::Element& ele,
      const Core::FE::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& displacements = *discretization.GetState("displacement");

    std::vector<double> mydisp(lm.size());
    Core::FE::ExtractMyValues(displacements, mydisp, lm);

    Discret::ELEMENTS::ElementNodes<celltype> element_nodes =
        evaluate_element_nodes<celltype>(ele, mydisp);

    if constexpr (Core::FE::is_nurbs<celltype>)
    {
      // Obtain the information required for a NURBS element
      bool zero_size = Core::FE::Nurbs::GetMyNurbsKnotsAndWeights(
          discretization, &ele, element_nodes.knots_, element_nodes.weights_);
      if (zero_size)
        FOUR_C_THROW("GetMyNurbsKnotsAndWeights has to return a non zero size NURBS element.");
    }

    return element_nodes;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the Gauss point according the the Gauss rule
   *
   * @tparam celltype : Cell type
   * @param intpoints (in) : Gauss integration points
   * @param gp (in) : id of the Gauss point
   * @return Core::LinAlg::Matrix<num_dim<celltype>, 1> : Coordinates of the Gauss Point in the
   * parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinate(
      const Core::FE::GaussIntegration& intpoints, const int gp)
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = intpoints.Point(gp)[d];

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Hexes
   *
   * Returns xi = [0 0 0].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <Core::FE::CellType celltype,
      std::enable_if_t<Core::FE::is_hex<celltype> | Core::FE::is_nurbs<celltype>, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = 0;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Tets
   *
   * Returns xi = [0.25 0.25 0.25].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <Core::FE::CellType celltype, std::enable_if_t<Core::FE::is_tet<celltype>, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Pyramids
   *
   * Returns xi = [0 0 0.25].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <Core::FE::CellType celltype, std::enable_if_t<Core::FE::is_pyramid<celltype>, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi(true);
    xi(2) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Wedges
   *
   * Returns xi = [1/3 1/3 0].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <Core::FE::CellType celltype, std::enable_if_t<Core::FE::is_wedge<celltype>, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi(true);
    xi(0) = 1.0 / 3.0;
    xi(1) = 1.0 / 3.0;

    return xi;
  }

  /*!
   * @brief Evaluates a point's coordinates in reference configuration
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates_reference (in) : Reference coordinates of the nodes of the element
   * @param shape_functions_point (in) : Shape functions evaluated at the specific point
   * @return Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> : point's reference coordinates
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> EvaluateReferenceCoordinate(
      const Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>&
          nodal_coordinates_reference,
      const Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1>& shape_functions_point)
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> coordinates_reference(true);
    coordinates_reference.MultiplyTN(nodal_coordinates_reference, shape_functions_point);

    return coordinates_reference;
  }

  /*!
   * @brief Evaluates the element centroid's coordinates in reference configuration
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference coordinates of the nodes of the element
   * @return Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> : Element centroid's coordinates in
   * reference configuration
   */
  template <Core::FE::CellType celltype,
      std::enable_if_t<Core::FE::use_lagrange_shapefnct<celltype>, int> = 0>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1> shape_functions_centroid(true);
    Core::FE::shape_function<celltype>(xi_centroid, shape_functions_centroid);

    return EvaluateReferenceCoordinate<celltype>(
        nodal_coordinates.reference_coordinates_, shape_functions_centroid);
  }

  template <Core::FE::CellType celltype, std::enable_if_t<Core::FE::is_nurbs<celltype>, int> = 0>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1> shape_functions_centroid(true);
    Core::FE::Nurbs::nurbs_shape_function_dim(shape_functions_centroid, xi_centroid,
        nodal_coordinates.knots_, nodal_coordinates.weights_, celltype);

    return EvaluateReferenceCoordinate<celltype>(
        nodal_coordinates.reference_coordinates_, shape_functions_centroid);
  }

  /*!
   * @brief Type holding the shape functions and it's first derivatives evaluated at a specific
   * point.
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct ShapeFunctionsAndDerivatives
  {
    Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1> shapefunctions_;
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> derivatives_;
  };

  /*!
   * @brief Evaluates the shape functions and their derivatives at the specified point in the
   * parameter space
   *
   * @tparam celltype : Discretizationt type
   * @param xi (in) : Coordinate in the parameter space
   * @return ShapeFunctionsAndDerivatives<celltype> : An object holding the shape functions and the
   * first derivatives evaluated at the respective point in the parameter space
   */
  template <Core::FE::CellType celltype,
      std::enable_if_t<Core::FE::use_lagrange_shapefnct<celltype>, bool> = true>
  ShapeFunctionsAndDerivatives<celltype> EvaluateShapeFunctionsAndDerivs(
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    Core::FE::shape_function<celltype>(xi, shapefcns.shapefunctions_);
    Core::FE::shape_function_deriv1<celltype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <Core::FE::CellType celltype,
      std::enable_if_t<Core::FE::is_nurbs<celltype>, bool> = true>
  ShapeFunctionsAndDerivatives<celltype> EvaluateShapeFunctionsAndDerivs(
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    Core::FE::Nurbs::nurbs_get_funct_deriv(shapefcns.shapefunctions_, shapefcns.derivatives_, xi,
        nodal_coordinates.knots_, nodal_coordinates.weights_, celltype);

    return shapefcns;
  }

  template <Core::FE::CellType celltype>
  struct JacobianMapping
  {
    /// Determinant of the jacobian
    double determinant_;

    /// Jacobian matrix at a specific point
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> jacobian_;

    /// Inverse jacobian matrix at a specific point
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> inverse_jacobian_;

    /// Derivative of the shape functions w.r.t. the reference coordinates
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> N_XYZ_;
  };

  /*!
   * @brief Evaluates the jacobian mapping of the element
   *
   * @tparam celltype : Cell type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param gp (in) : Id of the Gauss point
   * @return JacobianMapping<celltype> : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ, integration
   * factor)
   */
  template <Core::FE::CellType celltype>
  JacobianMapping<celltype> EvaluateJacobianMapping(
      const ShapeFunctionsAndDerivatives<celltype>& shapefcns,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    JacobianMapping<celltype> jacobian;

    jacobian.jacobian_.Multiply(shapefcns.derivatives_, nodal_coordinates.reference_coordinates_);
    jacobian.inverse_jacobian_ = jacobian.jacobian_;
    jacobian.determinant_ = jacobian.inverse_jacobian_.Invert();
    jacobian.N_XYZ_.Multiply(jacobian.inverse_jacobian_, shapefcns.derivatives_);

    return jacobian;
  }

  /*!
   * @brief Evaluates the jacobian determinant of the element
   *
   * @tparam celltype : Cell type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Jacobian determinant
   */
  template <Core::FE::CellType celltype>
  double EvaluateJacobianDeterminant(const ShapeFunctionsAndDerivatives<celltype>& shapefcns,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> jacobian;
    jacobian.Multiply(shapefcns.derivatives_, nodal_coordinates.reference_coordinates_);

    return jacobian.Determinant();
  }

  /*!
   * @brief Check for negative Jacobian determinants
   *
   * @tparam celltype : Cell type
   * @param element_nodes (in) : Element node information
   */
  template <Core::FE::CellType celltype>
  void ensure_positive_jacobian_determinant_at_element_nodes(
      const ElementNodes<celltype>& element_nodes)
  {
    const Core::LinAlg::SerialDenseMatrix rst =
        Core::FE::getEleNodeNumbering_nodes_paramspace(celltype);

    for (const auto& xi : Core::FE::get_element_nodes_in_parameter_space<celltype>())
    {
      Core::LinAlg::Matrix<3, 1> xi_mat(xi.data(), true);

      ShapeFunctionsAndDerivatives<celltype> shape_functions =
          EvaluateShapeFunctionsAndDerivs(xi_mat, element_nodes);

      JacobianMapping<celltype> jacobian_mapping =
          EvaluateJacobianMapping(shape_functions, element_nodes);

      FOUR_C_THROW_UNLESS(jacobian_mapping.determinant_ > 0,
          "Determinant of jacobian is %f <= 0 at one node of the element.",
          jacobian_mapping.determinant_);
    }
  }

  /*!
   * @brief Evaluate the jacobian mapping at the element centroid
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return JacobianMapping<celltype> : jacobian mapping at the element centroid
   */
  template <Core::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  JacobianMapping<celltype> evaluate_jacobian_mapping_centroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    // shape functions and derivatives evaluated at element centroid
    const ShapeFunctionsAndDerivatives<celltype> shape_functions_centroid =
        EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, nodal_coordinates);

    // jacobian mapping evaluated at centroid
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        EvaluateJacobianMapping(shape_functions_centroid, nodal_coordinates);

    return jacobian_mapping_centroid;
  }

  template <Core::FE::CellType celltype>
  struct SpatialMaterialMapping
  {
    double determinant_deformation_gradient_{};
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>
        deformation_gradient_{};
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>
        inverse_deformation_gradient_{};
  };

  /*!
   * @brief Evaluates the mapping between the spatial and material configuration of the element
   *
   * @tparam celltype : Cell type
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param scale_defgrd (in) : scaling for deformation gradient
   * @param kinematictype (in) : kinematic type of element
   * @return SpatialMaterialMapping<celltype> : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   */
  template <Core::FE::CellType celltype>
  SpatialMaterialMapping<celltype> evaluate_spatial_material_mapping(
      const JacobianMapping<celltype>& jacobian_mapping,
      const ElementNodes<celltype>& nodal_coordinates, const double scale_defgrd = 1.0,
      const Inpar::STR::KinemType& kinematictype = Inpar::STR::KinemType::nonlinearTotLag)
  {
    SpatialMaterialMapping<celltype> spatial_material_mapping;
    spatial_material_mapping.deformation_gradient_ =
        Core::LinAlg::IdentityMatrix<Core::FE::dim<celltype>>();

    if (kinematictype == Inpar::STR::KinemType::nonlinearTotLag)
    {
      spatial_material_mapping.deformation_gradient_.MultiplyTT(
          scale_defgrd, nodal_coordinates.displacements_, jacobian_mapping.N_XYZ_, scale_defgrd);
    }

    spatial_material_mapping.inverse_deformation_gradient_.Invert(
        spatial_material_mapping.deformation_gradient_);
    spatial_material_mapping.determinant_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_.Determinant();

    return spatial_material_mapping;
  }

  /*!
   * @brief Evaluate the determinant of the deformation gradient at the element centroid
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Determinant of the deformation gradient at the centroid
   */
  template <Core::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  double EvaluateDeformationGradientDeterminantCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    // jacobian mapping at centroid of element
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        evaluate_jacobian_mapping_centroid(nodal_coordinates);

    // deformation gradient and strains at centroid of element
    const Discret::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping_centroid =
        evaluate_spatial_material_mapping(jacobian_mapping_centroid, nodal_coordinates);

    return spatial_material_mapping_centroid.determinant_deformation_gradient_;
  }

  /*!
   * @brief Evaluates Green-Lagrange strain from right Cauchy-Green tensor
   *
   * GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
   *
   * @param cauchygreen (in) : Right Cauchy-Green deformation tensor
   * @return Core::LinAlg::Matrix<DETAIL::num_str<celltype>, 1> : Green-Lagrange strain tensor in
   * strain-like Voigt notation
   */
  inline Core::LinAlg::Matrix<6, 1> evaluate_green_lagrange_strain(
      const Core::LinAlg::Matrix<3, 3>& cauchygreen)
  {
    Core::LinAlg::Matrix<6, 1> gl_strain;

    gl_strain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    gl_strain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    gl_strain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    gl_strain(3) = cauchygreen(0, 1);
    gl_strain(4) = cauchygreen(1, 2);
    gl_strain(5) = cauchygreen(2, 0);

    return gl_strain;
  }

  /*!
   * @brief Evaluates right Cauchy-Green deformation tensor
   *
   * @tparam celltype: Cell type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> : Right
   * Cauchy-Green deformation tensor
   */
  template <Core::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> evaluate_cauchy_green(
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> cauchygreen(false);

    cauchygreen.MultiplyTN(spatial_material_mapping.deformation_gradient_,
        spatial_material_mapping.deformation_gradient_);

    return cauchygreen;
  }

  /*!
   * @brief Evaluates the strain gradient (B-Operator) of the specified element
   *
   * @tparam celltype : Cell type
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return Core::LinAlg::Matrix<num_str<celltype>, num_dim<celltype> * num_nodes<celltype>> :
   * B-Operator
   */
  template <Core::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  Core::LinAlg::Matrix<DETAIL::num_str<celltype>,
      DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
  evaluate_strain_gradient(const JacobianMapping<celltype>& jacobian_mapping,
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
  {
    // B-operator
    Core::LinAlg::Matrix<DETAIL::num_str<celltype>,
        DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
        Bop;
    for (int i = 0; i < DETAIL::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
      {
        for (int e = 0; e < DETAIL::num_dim<celltype>; ++e)
        {
          Bop(d, DETAIL::num_dim<celltype> * i + e) =
              spatial_material_mapping.deformation_gradient_(e, d) * jacobian_mapping.N_XYZ_(d, i);
        }
      }

      Bop(3, DETAIL::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, DETAIL::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, DETAIL::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(4, DETAIL::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, DETAIL::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, DETAIL::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(5, DETAIL::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, DETAIL::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, DETAIL::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(2, i);
    }

    return Bop;
  }

  template <Core::FE::CellType celltype>
  struct Stress
  {
    /// Second Piola-Kirchhoff stress tensor in stress-like voigt notation
    Core::LinAlg::Matrix<DETAIL::num_str<celltype>, 1> pk2_;

    /// Linearization of the 2. Piola Kirchhoff stress tensor w.r.t. Green-Lagrange strain tensor in
    /// mixed Voigt notation
    Core::LinAlg::Matrix<DETAIL::num_str<celltype>, DETAIL::num_str<celltype>> cmat_;
  };

  /*!
   * @brief Evaluates the material stress (2. Piola-Kirchhoff stress tensor and the linearization
   * w.r.t. Green-Lagrange strain)
   *
   * @tparam celltype : Cell type
   * @param material (in) : Reference to the material
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @param gl_strain (in) : Green-Lagrange strain
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<celltype> : Object holding the 2. Piola-Kirchhoff stress tensor and the
   * linearization w.r.t. Green-Lagrange strain tensor
   */
  template <Core::FE::CellType celltype>
  Stress<celltype> evaluate_material_stress(Mat::So3Material& material,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>& defgrd,
      const Core::LinAlg::Matrix<DETAIL::num_str<celltype>, 1>& gl_strain,
      Teuchos::ParameterList& params, const int gp, const int eleGID)
  {
    Stress<celltype> stress;

    material.evaluate(&defgrd, &gl_strain, params, &stress.pk2_, &stress.cmat_, gp, eleGID);
    return stress;
  }

  /*!
   * @brief Adds the internal force vector contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_internal_force_vector(const Core::LinAlg::Matrix<DETAIL::num_str<celltype>,
                                     DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>, 1>&
          force_vector)
  {
    force_vector.MultiplyTN(integration_fac, Bop, stress.pk2_, 1.);
  }

  /*!
   * @brief Add elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_elastic_stiffness_matrix(
      const Core::LinAlg::Matrix<DETAIL::num_str<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& stiffness_matrix)
  {
    Core::LinAlg::Matrix<DETAIL::num_str<celltype>,
        DETAIL::num_nodes<celltype> * DETAIL::num_dim<celltype>>
        cb;
    cb.Multiply(stress.cmat_, Bop);
    stiffness_matrix.MultiplyTN(integration_fac, Bop, cb, 1.0);
  }

  /*!
   * @brief Add geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param B_L (in) : B_L operator, i.e. derivatives of the shape functions w.r.t. XYZ
   *                   at the respective Gauss point
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_geometric_stiffness_matrix(
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>>& B_L,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& stiffness_matrix)
  {
    std::array<double, 3> SmB_L;  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < DETAIL::num_nodes<celltype>; ++inod)
    {
      SmB_L[0] = stress.pk2_(0) * B_L(0, inod) + stress.pk2_(3) * B_L(1, inod) +
                 stress.pk2_(5) * B_L(2, inod);
      SmB_L[1] = stress.pk2_(3) * B_L(0, inod) + stress.pk2_(1) * B_L(1, inod) +
                 stress.pk2_(4) * B_L(2, inod);
      SmB_L[2] = stress.pk2_(5) * B_L(0, inod) + stress.pk2_(4) * B_L(1, inod) +
                 stress.pk2_(2) * B_L(2, inod);

      for (int jnod = 0; jnod < DETAIL::num_nodes<celltype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < DETAIL::num_dim<celltype>; ++idim)
          bopstrbop += B_L(idim, jnod) * SmB_L[idim];

        for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
          stiffness_matrix(DETAIL::num_dim<celltype> * inod + d,
              DETAIL::num_dim<celltype> * jnod + d) += integration_fac * bopstrbop;
      }
    }
  }

  /*!
   * @brief Add mass matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param shapefunctions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param integration_factor (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param density (in) : density at the Gauss point
   * @param mass (in/out) : mass matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_mass_matrix(const ShapeFunctionsAndDerivatives<celltype>& shapefunctions,
      const double integration_factor, const double density,
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& mass)
  {
    for (int inod = 0; inod < DETAIL::num_nodes<celltype>; ++inod)
    {
      const double ifactor = shapefunctions.shapefunctions_(inod) * integration_factor * density;
      for (int jnod = 0; jnod < DETAIL::num_nodes<celltype>; ++jnod)
      {
        const double massfactor =
            shapefunctions.shapefunctions_(jnod) * ifactor;  // intermediate factor
        for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
          mass(DETAIL::num_dim<celltype> * inod + d, DETAIL::num_dim<celltype> * jnod + d) +=
              massfactor;
      }
    }
  }

  /*!
   * @brief Calls the @p gp_evaluator for each Gauss point with evaluated jacobian mapping using the
   * integration rule defined by @p integration.
   *
   * @tparam celltype : Cell type known at compile time
   * @tparam GaussPointEvaluator
   * @param nodal_coordinates (in) : The nodal coordinates of the element
   * @param integration (in) : The integration rule to be used.
   * @param gp_evaluator (in) : A callable object (e.g. lambda-function) with signature void(const
   * Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi, const
   * ShapeFunctionsAndDerivatives<celltype>& shape_functions, const JacobianMapping<celltype>&
   * jacobian_mapping, double integration_factor, int gp) that will be called for each integration
   * point.
   */
  template <Core::FE::CellType celltype, typename GaussPointEvaluator>
  inline void ForEachGaussPoint(const ElementNodes<celltype>& nodal_coordinates,
      const Core::FE::GaussIntegration& integration, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < integration.NumPoints(); ++gp)
    {
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi =
          EvaluateParameterCoordinate<celltype>(integration, gp);

      const ShapeFunctionsAndDerivatives<celltype> shape_functions =
          EvaluateShapeFunctionsAndDerivs<celltype>(xi, nodal_coordinates);

      const JacobianMapping<celltype> jacobian_mapping =
          EvaluateJacobianMapping(shape_functions, nodal_coordinates);

      const double integration_factor = jacobian_mapping.determinant_ * integration.Weight(gp);

      gp_evaluator(xi, shape_functions, jacobian_mapping, integration_factor, gp);
    }
  }

  /*!
   * @brief Evaluates the Gauss point coordinates in reference configuration
   * and adds those to the paramater list
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param shape_functions_gp (in) : Shape functions evaluated at the Gauss point
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <Core::FE::CellType celltype>
  void evaluate_gp_coordinates_and_add_to_parameter_list(
      const ElementNodes<celltype>& nodal_coordinates,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions_gp,
      Teuchos::ParameterList& params)
  {
    auto gp_ref_coord = EvaluateReferenceCoordinate<celltype>(
        nodal_coordinates.reference_coordinates_, shape_functions_gp.shapefunctions_);
    params.set("gp_coords_ref", gp_ref_coord);
  }

  /*!
   * @brief Evaluates the element centroid coordinates in reference configuration
   * and adds those to the paramater list
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <Core::FE::CellType celltype>
  void evaluate_centroid_coordinates_and_add_to_parameter_list(
      const ElementNodes<celltype>& nodal_coordinates, Teuchos::ParameterList& params)
  {
    auto element_center = EvaluateReferenceCoordinateCentroid<celltype>(nodal_coordinates);
    params.set("elecenter_coords_ref", element_center);
  }

  /*!
   * @brief For elements with fiber nodes, interpolate fibers to Gauss points and add to the
   * parameter list
   *
   * @tparam celltype : Cell type
   * @param stiffness_matrix_integration (in) : Container holding the integration points
   * @param ele (in) : Reference to the element, possibly having nodal fibers
   * @param params (in/out) : ParameterList the interpolated fibers are added to
   */
  template <Core::FE::CellType celltype>
  inline void InterpolateFibersToGaussPointsAndAddToParameterList(
      const Core::FE::GaussIntegration& stiffness_matrix_integration,
      const Core::Elements::Element& ele, Teuchos::ParameterList& params)
  {
    if (Core::Nodes::HaveNodalFibers<celltype>(ele.Nodes()))
    {
      // This element has fiber nodes.
      // Interpolate fibers to the Gauss points and add them to the parameter list

      // Get shape functions
      const static std::vector<Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1>> shapefcts =
          std::invoke(
              [&]
              {
                std::vector<Core::LinAlg::Matrix<DETAIL::num_nodes<celltype>, 1>> shapefcns(
                    stiffness_matrix_integration.NumPoints());
                for (int gp = 0; gp < stiffness_matrix_integration.NumPoints(); ++gp)
                {
                  Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi(
                      stiffness_matrix_integration.Point(gp), true);
                  Core::FE::shape_function<celltype>(xi, shapefcns[gp]);
                }
                return shapefcns;
              });

      // add fibers to the ParameterList
      Core::Nodes::NodalFiberHolder fiberHolder;

      // Do the interpolation
      Core::Nodes::ProjectFibersToGaussPoints<celltype>(ele.Nodes(), shapefcts, fiberHolder);

      params.set("fiberholder", fiberHolder);
    }
  }

  /*!
   * @brief Linearization of a solid formulation containing the derivatives of the deformation
   * gradient w.r.t. the nodal displacements, xi and the second derivative w.r.t. xi and
   * displacements
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct SolidFormulationLinearization
  {
    /// Derivative of the deformation gradient w.r.t. nodal displacements
    Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>> d_F_dd{};

    /// Derivative of the deformation gradient w.r.t. xi
    Core::LinAlg::Matrix<9, Core::FE::dim<celltype>> d_F_dxi{};

    /// 2. Derivative of the deformation gradient w.r.t. xi and nodal displacements
    Core::LinAlg::Matrix<9,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
        d2_F_dxi_dd{};
  };

}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif