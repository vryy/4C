/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_LIB_HPP
#define BACI_SOLID_ELE_CALC_LIB_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_discretization_fem_general_utils_gauss_point_postprocess.hpp"
#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "baci_fiber_nodal_fiber_holder.hpp"
#include "baci_fiber_utils.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_utils.hpp"
#include "baci_linalg_fixedsizematrix_generators.hpp"
#include "baci_mat_so3_material.hpp"
#include "baci_nurbs_discret_nurbs_utils.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype, typename Enable = void>
  struct ElementNodes;
}  // namespace DRT::ELEMENTS

namespace DRT::ELEMENTS::DETAIL
{
  template <CORE::FE::CellType celltype>
  inline static constexpr int num_nodes = CORE::FE::num_nodes<celltype>;

  template <CORE::FE::CellType celltype>
  inline static constexpr int num_dim = CORE::FE::dim<celltype>;

  template <CORE::FE::CellType celltype>
  inline static constexpr int num_str = num_dim<celltype>*(num_dim<celltype> + 1) / 2;

  template <CORE::FE::CellType celltype>
  inline static constexpr int num_dof_per_ele = num_nodes<celltype>* num_dim<celltype>;
}  // namespace DRT::ELEMENTS::DETAIL

namespace DRT::ELEMENTS
{

  /*!
   * @brief Calculate the lumped mass matrix
   */
  inline void LumpMatrix(CORE::LINALG::SerialDenseMatrix& matrix)
  {
    dsassert(
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
  template <CORE::FE::CellType celltype, typename Enable>
  struct ElementNodes
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>
        reference_coordinates_;

    /*!
     * @brief Displacements of the element nodes
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>> displacements_;
  };

  template <CORE::FE::CellType celltype>
  struct ElementNodes<celltype, std::enable_if_t<CORE::FE::is_nurbs<celltype>>>
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>
        reference_coordinates_;

    /*!
     * @brief Displacements of the element nodes
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>> displacements_;

    /*!
     * @brief Knot span of a NURBS element
     */
    std::vector<CORE::LINALG::SerialDenseVector> knots_;

    /*!
     * @brief Weights of control points
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1, double> weights_;
  };

  /*!
   * @brief Evaluates the nodal coordinates from this iteration
   *
   * @param ele (in) : Reference to the element
   * @param discretization (in) : Discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <CORE::FE::CellType celltype>
  ElementNodes<celltype> EvaluateElementNodes(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& displacements = *discretization.GetState("displacement");

    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(displacements, mydisp, lm);

    DRT::ELEMENTS::ElementNodes<celltype> element_nodes;
    for (int i = 0; i < DETAIL::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
      {
        element_nodes.reference_coordinates_(i, d) = ele.Nodes()[i]->X()[d];
        element_nodes.displacements_(i, d) = mydisp[i * DETAIL::num_dim<celltype> + d];
      }
    }

    if constexpr (CORE::FE::is_nurbs<celltype>)
    {
      // Obtain the information required for a NURBS element
      bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(
          discretization, &ele, element_nodes.knots_, element_nodes.weights_);
      if (zero_size)
        dserror("GetMyNurbsKnotsAndWeights has to return a non zero size NURBS element.");
    }

    return element_nodes;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the Gauss point according the the Gauss rule
   *
   * @tparam celltype : Cell type
   * @param intpoints (in) : Gauss integration points
   * @param gp (in) : id of the Gauss point
   * @return CORE::LINALG::Matrix<num_dim<celltype>, 1> : Coordinates of the Gauss Point in the
   * parameter space
   */
  template <CORE::FE::CellType celltype>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinate(
      const CORE::FE::GaussIntegration& intpoints, const int gp)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = intpoints.Point(gp)[d];

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Hexes
   *
   * Returns xi = [0 0 0].
   *
   * @tparam celltype : Cell type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType celltype,
      std::enable_if_t<CORE::FE::is_hex<celltype> | CORE::FE::is_nurbs<celltype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = 0;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Tets
   *
   * Returns xi = [0.25 0.25 0.25].
   *
   * @tparam celltype : Cell type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::is_tet<celltype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<celltype>; ++d) xi(d) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Pyramids
   *
   * Returns xi = [0 0 0.25].
   *
   * @tparam celltype : Cell type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::is_pyramid<celltype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi(true);
    xi(2) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Wedges
   *
   * Returns xi = [1/3 1/3 0].
   *
   * @tparam celltype : Cell type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::is_wedge<celltype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi(true);
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
   * @return CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> : point's reference coordinates
   */
  template <CORE::FE::CellType celltype>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> EvaluateReferenceCoordinate(
      const CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, DETAIL::num_dim<celltype>>&
          nodal_coordinates_reference,
      const CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1>& shape_functions_point)
  {
    CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> coordinates_reference(true);
    coordinates_reference.MultiplyTN(shape_functions_point, nodal_coordinates_reference);

    return coordinates_reference;
  }

  /*!
   * @brief Evaluates the element centroid's coordinates in reference configuration
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference coordinates of the nodes of the element
   * @return CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> : Element centroid's coordinates in
   * reference configuration
   */
  template <CORE::FE::CellType celltype,
      std::enable_if_t<CORE::FE::use_lagrange_shapefnct<celltype>, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1> shape_functions_centroid(true);
    CORE::FE::shape_function<celltype>(xi_centroid, shape_functions_centroid);

    const CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> centroid_coordinates_reference =
        EvaluateReferenceCoordinate<celltype>(
            nodal_coordinates.reference_coordinates_, shape_functions_centroid);

    return centroid_coordinates_reference;
  }

  template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::is_nurbs<celltype>, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1> shape_functions_centroid(true);
    CORE::FE::NURBS::nurbs_shape_function_dim(shape_functions_centroid, xi_centroid,
        nodal_coordinates.knots_, nodal_coordinates.weights_, celltype);

    const CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype>> centroid_coordinates_reference =
        EvaluateReferenceCoordinate<celltype>(
            nodal_coordinates.reference_coordinates_, shape_functions_centroid);

    return centroid_coordinates_reference;
  }

  /*!
   * @brief Type holding the shape functions and it's first derivatives evaluated at a specific
   * point.
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct ShapeFunctionsAndDerivatives
  {
    CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1> shapefunctions_;
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> derivatives_;
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
  template <CORE::FE::CellType celltype,
      std::enable_if_t<CORE::FE::use_lagrange_shapefnct<celltype>, bool> = true>
  ShapeFunctionsAndDerivatives<celltype> EvaluateShapeFunctionsAndDerivs(
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    CORE::FE::shape_function<celltype>(xi, shapefcns.shapefunctions_);
    CORE::FE::shape_function_deriv1<celltype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <CORE::FE::CellType celltype,
      std::enable_if_t<CORE::FE::is_nurbs<celltype>, bool> = true>
  ShapeFunctionsAndDerivatives<celltype> EvaluateShapeFunctionsAndDerivs(
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    CORE::FE::NURBS::nurbs_get_funct_deriv(shapefcns.shapefunctions_, shapefcns.derivatives_, xi,
        nodal_coordinates.knots_, nodal_coordinates.weights_, celltype);

    return shapefcns;
  }

  template <CORE::FE::CellType celltype>
  struct JacobianMapping
  {
    /// Determinant of the jacobian
    double determinant_;

    /// Jacobian matrix at a specific point
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> jacobian_;

    /// Inverse jacobian matrix at a specific point
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> inverse_jacobian_;

    /// Derivative of the shape functions w.r.t. the reference coordinates
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> N_XYZ_;
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
  template <CORE::FE::CellType celltype>
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
  template <CORE::FE::CellType celltype>
  double EvaluateJacobianDeterminant(const ShapeFunctionsAndDerivatives<celltype>& shapefcns,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> jacobian;
    jacobian.Multiply(shapefcns.derivatives_, nodal_coordinates.reference_coordinates_);

    return jacobian.Determinant();
  }

  /*!
   * @brief Evaluate the jacobian mapping at the element centroid
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return JacobianMapping<celltype> : jacobian mapping at the element centroid
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  JacobianMapping<celltype> EvaluateJacobianMappingCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<celltype>();

    // shape functions and derivatives evaluated at element centroid
    const ShapeFunctionsAndDerivatives<celltype> shape_functions_centroid =
        EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, nodal_coordinates);

    // jacobian mapping evaluated at centroid
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        EvaluateJacobianMapping(shape_functions_centroid, nodal_coordinates);

    return jacobian_mapping_centroid;
  }

  template <CORE::FE::CellType celltype>
  struct SpatialMaterialMapping
  {
    double determinant_deformation_gradient_;
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>
        deformation_gradient_;
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>
        inverse_deformation_gradient_;
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
  template <CORE::FE::CellType celltype>
  SpatialMaterialMapping<celltype> EvaluateSpatialMaterialMapping(
      const JacobianMapping<celltype>& jacobian_mapping,
      const ElementNodes<celltype>& nodal_coordinates, const double scale_defgrd = 1.0,
      const INPAR::STR::KinemType& kinematictype = INPAR::STR::KinemType::nonlinearTotLag)
  {
    SpatialMaterialMapping<celltype> spatial_material_mapping;
    spatial_material_mapping.deformation_gradient_ =
        CORE::LINALG::IdentityMatrix<CORE::FE::dim<celltype>>();

    if (kinematictype == INPAR::STR::KinemType::nonlinearTotLag)
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
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  double EvaluateDeformationGradientDeterminantCentroid(
      const ElementNodes<celltype>& nodal_coordinates)
  {
    // jacobian mapping at centroid of element
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        EvaluateJacobianMappingCentroid(nodal_coordinates);

    // deformation gradient and strains at centroid of element
    const DRT::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping_centroid =
        EvaluateSpatialMaterialMapping(jacobian_mapping_centroid, nodal_coordinates);

    return spatial_material_mapping_centroid.determinant_deformation_gradient_;
  }

  /*!
   * @brief Evaluates Green-Lagrange strain from right Cauchy-Green tensor
   *
   * GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
   *
   * @param cauchygreen (in) : Right Cauchy-Green deformation tensor
   * @return CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> : Green-Lagrange strain tensor in
   * strain-like Voigt notation
   */
  inline CORE::LINALG::Matrix<6, 1> EvaluateGreenLagrangeStrain(
      const CORE::LINALG::Matrix<3, 3>& cauchygreen)
  {
    CORE::LINALG::Matrix<6, 1> gl_strain;

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
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> : Right
   * Cauchy-Green deformation tensor
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> EvaluateCauchyGreen(
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> cauchygreen(false);

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
   * @return CORE::LINALG::Matrix<num_str<celltype>, num_dim<celltype> * num_nodes<celltype>> :
   * B-Operator
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
      DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
  EvaluateStrainGradient(const JacobianMapping<celltype>& jacobian_mapping,
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
  {
    // B-operator
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
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

  template <CORE::FE::CellType celltype>
  struct Stress
  {
    /// Second Piola-Kirchhoff stress tensor in stress-like voigt notation
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> pk2_;

    /// Linearization of the 2. Piola Kirchhoff stress tensor w.r.t. Green-Lagrange strain tensor in
    /// mixed Voigt notation
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>, DETAIL::num_str<celltype>> cmat_;
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
  template <CORE::FE::CellType celltype>
  Stress<celltype> EvaluateMaterialStress(MAT::So3Material& material,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>& defgrd,
      const CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1>& gl_strain,
      Teuchos::ParameterList& params, const int gp, const int eleGID)
  {
    Stress<celltype> stress;

    material.Evaluate(&defgrd, &gl_strain, params, &stress.pk2_, &stress.cmat_, gp, eleGID);
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
  template <CORE::FE::CellType celltype>
  void AddInternalForceVector(const CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
                                  DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>, 1>&
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
  template <CORE::FE::CellType celltype>
  void AddElasticStiffnessMatrix(const CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
                                     DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& stiffness_matrix)
  {
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
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
  template <CORE::FE::CellType celltype>
  void AddGeometricStiffnessMatrix(
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>>& B_L,
      const Stress<celltype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
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
  template <CORE::FE::CellType celltype>
  void AddMassMatrix(const ShapeFunctionsAndDerivatives<celltype>& shapefunctions,
      const double integration_factor, const double density,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
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
   * CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi, const
   * ShapeFunctionsAndDerivatives<celltype>& shape_functions, const JacobianMapping<celltype>&
   * jacobian_mapping, double integration_factor, int gp) that will be called for each integration
   * point.
   */
  template <CORE::FE::CellType celltype, typename GaussPointEvaluator>
  inline void ForEachGaussPoint(const ElementNodes<celltype>& nodal_coordinates,
      const CORE::FE::GaussIntegration& integration, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < integration.NumPoints(); ++gp)
    {
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi =
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
  template <CORE::FE::CellType celltype>
  void EvaluateGPCoordinatesAndAddToParameterList(const ElementNodes<celltype>& nodal_coordinates,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions_gp,
      Teuchos::ParameterList& params)
  {
    auto gp_ref_coord = EvaluateReferenceCoordinate<celltype>(
        nodal_coordinates.reference_coordinates_, shape_functions_gp.shapefunctions_);
    params.set("gprefecoord", gp_ref_coord);
  }

  /*!
   * @brief Evaluates the element centroid coordinates in reference configuration
   * and adds those to the paramater list
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <CORE::FE::CellType celltype>
  void EvaluateCentroidCoordinatesAndAddToParameterList(
      const ElementNodes<celltype>& nodal_coordinates, Teuchos::ParameterList& params)
  {
    auto element_center = EvaluateReferenceCoordinateCentroid<celltype>(nodal_coordinates);
    params.set("elecenter", element_center);
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
  template <CORE::FE::CellType celltype>
  inline void InterpolateFibersToGaussPointsAndAddToParameterList(
      const CORE::FE::GaussIntegration& stiffness_matrix_integration, const DRT::Element& ele,
      Teuchos::ParameterList& params)
  {
    if (DRT::FIBER::UTILS::HaveNodalFibers<celltype>(ele.Nodes()))
    {
      // This element has fiber nodes.
      // Interpolate fibers to the Gauss points and add them to the parameter list

      // Get shape functions
      const static std::vector<CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1>> shapefcts =
          std::invoke(
              [&]
              {
                std::vector<CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1>> shapefcns(
                    stiffness_matrix_integration.NumPoints());
                for (int gp = 0; gp < stiffness_matrix_integration.NumPoints(); ++gp)
                {
                  CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi(
                      stiffness_matrix_integration.Point(gp), true);
                  CORE::FE::shape_function<celltype>(xi, shapefcns[gp]);
                }
                return shapefcns;
              });

      // add fibers to the ParameterList
      DRT::FIBER::NodalFiberHolder fiberHolder;

      // Do the interpolation
      DRT::FIBER::UTILS::ProjectFibersToGaussPoints<celltype>(ele.Nodes(), shapefcts, fiberHolder);

      params.set("fiberholder", fiberHolder);
    }
  }

}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_CALC_LIB_H