/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_LIB_H
#define BACI_SOLID_ELE_CALC_LIB_H

#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"
#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_fiber_nodal_fiber_holder.H"
#include "baci_fiber_utils.H"
#include "baci_lib_discret.H"
#include "baci_lib_element.H"
#include "baci_lib_element_integration_select.H"
#include "baci_lib_utils.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_mat_so3_material.H"
#include "baci_nurbs_discret_nurbs_utils.H"
#include "baci_so3_element_service.H"
#include "baci_solid_ele_utils.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include <numeric>

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType distype, typename Enable = void>
  struct ElementNodes;
}  // namespace DRT::ELEMENTS

namespace DRT::ELEMENTS::DETAIL
{
  template <CORE::FE::CellType distype>
  inline static constexpr int num_nodes = CORE::FE::num_nodes<distype>;

  template <CORE::FE::CellType distype>
  inline static constexpr int num_dim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  template <CORE::FE::CellType distype>
  inline static constexpr int num_str = num_dim<distype>*(num_dim<distype> + 1) / 2;

  template <CORE::FE::CellType distype>
  inline static constexpr int num_dof_per_ele = num_nodes<distype>* num_dim<distype>;
}  // namespace DRT::ELEMENTS::DETAIL

namespace DRT::ELEMENTS::DETAILS
{
  /*!
   * @brief Get the nodal coordinates
   *
   * Extracts the displacement from the previous iteration from the discretization and computes the
   * current and reference position of the nodal coordinates
   *
   * @tparam distype : Discretization type
   * @param ele (in) : Reference to the element
   * @param discretization (in) : Reference to the discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @return ElementNodes<distype> : Nodal coordinates from previous iteration
   */
  template <CORE::FE::CellType distype>
  DRT::ELEMENTS::ElementNodes<distype> GetNodalCoordinates(
      const DRT::Element& ele, const Epetra_Vector& displacements, const std::vector<int>& lm)
  {
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(displacements, mydisp, lm);

    DRT::ELEMENTS::ElementNodes<distype> coordinates;
    for (int i = 0; i < DETAIL::num_nodes<distype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<distype>; ++d)
      {
        coordinates.reference_coordinates_(i, d) = ele.Nodes()[i]->X()[d];
        coordinates.current_coordinates_(i, d) =
            coordinates.reference_coordinates_(i, d) + mydisp[i * DETAIL::num_dim<distype> + d];
      }
    }

    return coordinates;
  }

  /*!
   * @brief Evaluates the nodal coordinates for the desired iteration (state)
   *
   * Extracts the displacement from the discretization and computes the current and reference
   * position of the nodal coordinates
   *
   * @tparam distype : Discretization type
   * @param ele (in) : Reference to the element
   * @param discretization (in) : Reference to the discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @param state_type (in) : State to obtain from discretization
   * @return ElementNodes<distype> : Nodal coordinates
   */

  template <CORE::FE::CellType distype,
      std::enable_if_t<CORE::FE::use_lagrange_shapefnct<distype>, bool> = true>
  DRT::ELEMENTS::ElementNodes<distype> EvaluateElementNodes(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm,
      const std::string& state_type)
  {
    const Epetra_Vector& displacements = *discretization.GetState(state_type);

    return GetNodalCoordinates<distype>(ele, displacements, lm);
  }

  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_nurbs<distype>, bool> = true>
  DRT::ELEMENTS::ElementNodes<distype> EvaluateElementNodes(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm,
      const std::string& state_type)
  {
    if (state_type == "old displacement")
    {
      dserror(
          "The nodal evaluation for the old displacement for the GEMM implementation is not tested "
          "for NURBS elements");
    }

    const Epetra_Vector& displacements = *discretization.GetState(state_type);

    auto coordinates = GetNodalCoordinates<distype>(ele, displacements, lm);

    // Obtain the information required for a NURBS element
    bool zero_size = ::DRT::NURBS::GetMyNurbsKnotsAndWeights(
        discretization, &ele, coordinates.knots_, coordinates.weights_);
    if (zero_size)
      dserror("GetMyNurbsKnotsAndWeights has to return a non zero size NURBS element.");

    return coordinates;
  }
}  // namespace DRT::ELEMENTS::DETAILS

namespace DRT::ELEMENTS
{
  /*!
   * @brief Compare two Gauss integration rules for equality
   */
  inline bool CompareGaussIntegration(const CORE::DRT::UTILS::GaussIntegration& integration_a,
      const CORE::DRT::UTILS::GaussIntegration& integration_b)
  {
    // currently this simple check is sufficient as we only use the same type of gauss integrations.
    return integration_a.NumPoints() == integration_b.NumPoints();
  }

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

  template <typename T>
  inline std::vector<char>& GetStressData(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().StressDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  template <typename T>
  inline std::vector<char>& GetStrainData(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().StrainDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  template <typename T>
  inline INPAR::STR::StressType GetIOStressType(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStressOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress");
    }
  }

  template <typename T>
  inline INPAR::STR::StrainType GetIOStrainType(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStrainOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain");
    }
  }

  /*!
   * @brief Get the default Gauss integration rules for different discretization types.
   *
   * @note It follows the rules defined in DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule,
   * except for the stiffness matrix of tetrahedral elements.
   *
   */
  /// @{

  template <CORE::FE::CellType distype>
  constexpr auto GetGaussRuleMassMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  }

  template <CORE::FE::CellType distype>
  constexpr auto GetGaussRuleStiffnessMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<CORE::FE::CellType::tet10>()
  {
    return CORE::DRT::UTILS::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<CORE::FE::CellType::tet4>()
  {
    return CORE::DRT::UTILS::GaussRule3D::tet_1point;
  }
  /// @}

  /*!
   * @brief Create a Gauss integration interface from a given Gauss rule type
   */
  template <CORE::FE::CellType distype, typename GaussRuleType>
  CORE::DRT::UTILS::GaussIntegration CreateGaussIntegration(GaussRuleType rule)
  {
    constexpr int num_dim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

    // setup default integration
    CORE::DRT::UTILS::IntPointsAndWeights<num_dim> intpoints(rule);

    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<CORE::DRT::UTILS::CollectedGaussPoints> gp =
        Teuchos::rcp(new CORE::DRT::UTILS::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < num_dim; ++d) xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return CORE::DRT::UTILS::GaussIntegration(gp);
  }

  /*!
   * @brief A type holding information of the nodes of an element,
   * such as their reference and current position. Additional information
   * is stored for NURBS elements
   *
   * @tparam distype
   */
  template <CORE::FE::CellType distype, typename Enable>
  struct ElementNodes
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>>
        reference_coordinates_;

    /*!
     * @brief Position of nodes in the current configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>> current_coordinates_;
  };

  template <CORE::FE::CellType distype>
  struct ElementNodes<distype, typename std::enable_if<CORE::FE::is_nurbs<distype>>::type>
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>>
        reference_coordinates_;

    /*!
     * @brief Position of nodes in the current configuration
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>> current_coordinates_;

    /*!
     * @brief Knot span of a NURBS element
     */
    std::vector<CORE::LINALG::SerialDenseVector> knots_;

    /*!
     * @brief Weights of control points
     */
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1, double> weights_;
  };

  /*!
   * @brief Extracts the nodal displacements of the element from the discretization
   *
   * @tparam distype
   * @param discretization (in) : Discretization with displacement state
   * @param lm (in) : Location vector of the element
   * @return CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>> : Nodal
   * displacements of the element
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>> GetNodalDisplacements(
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& displacements = *discretization.GetState("displacement");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(displacements, mydisp, lm);

    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>> nodal_displacements;
    for (int i = 0; i < DETAIL::num_nodes<distype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<distype>; ++d)
      {
        nodal_displacements(i, d) = mydisp[i * DETAIL::num_dim<distype> + d];
      }
    }
    return nodal_displacements;
  }

  /*!
   * @brief Evaluates the nodal coordinates from this iteration
   *
   * @param ele (in) : Reference to the element
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <CORE::FE::CellType distype>
  ElementNodes<distype> EvaluateElementNodes(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    return DETAILS::EvaluateElementNodes<distype>(ele, discretization, lm, "displacement");
  }

  /*!
   * @brief Evaluates the nodal coordinates from the previous iteration (old)
   *
   * @param ele (in) : Reference to the element
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <CORE::FE::CellType distype>
  ElementNodes<distype> EvaluateElementNodesOfPreviousTimestep(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    return DETAILS::EvaluateElementNodes<distype>(ele, discretization, lm, "old displacement");
  }

  /*!
   * @brief Evaluates the parameter coordinate of the Gauss point according the the Gauss rule
   *
   * @tparam distype : Discretization type
   * @param intpoints (in) : Gauss integration points
   * @param gp (in) : id of the Gauss point
   * @return CORE::LINALG::Matrix<num_dim<distype>, 1> : Coordinates of the Gauss Point in the
   * parameter space
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> EvaluateParameterCoordinate(
      const CORE::DRT::UTILS::GaussIntegration& intpoints, const int gp)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<distype>; ++d) xi(d) = intpoints.Point(gp)[d];

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Hexes
   *
   * Returns xi = [0 0 0].
   *
   * @tparam distype : Discretization type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType distype,
      std::enable_if_t<CORE::FE::is_hex<distype> | CORE::FE::is_nurbs<distype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<distype>; ++d) xi(d) = 0;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Tets
   *
   * Returns xi = [0.25 0.25 0.25].
   *
   * @tparam distype : Discretization type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_tet<distype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi;
    for (int d = 0; d < DETAIL::num_dim<distype>; ++d) xi(d) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Pyramids
   *
   * Returns xi = [0 0 0.25].
   *
   * @tparam distype : Discretization type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_pyramid<distype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi(true);
    xi(2) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Wedges
   *
   * Returns xi = [1/3 1/3 0].
   *
   * @tparam distype : Discretization type
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> : Coordinates of the centroid in the
   * parameter space
   */
  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_wedge<distype>, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> EvaluateParameterCoordinateCentroid()
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi(true);
    xi(0) = 1.0 / 3.0;
    xi(1) = 1.0 / 3.0;

    return xi;
  }

  /*!
   * @brief Type holding the shape functions and it's first derivatives evaluated at a specific
   * point.
   *
   * @tparam distype
   */
  template <CORE::FE::CellType distype>
  struct ShapeFunctionsAndDerivatives
  {
    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1> shapefunctions_;
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_nodes<distype>> derivatives_;
  };

  /*!
   * @brief Evaluates the shape functions and their derivatives at the specified point in the
   * parameter space
   *
   * @tparam distype : Discretizationt type
   * @param xi (in) : Coordinate in the parameter space
   * @return ShapeFunctionsAndDerivatives<distype> : An object holding the shape functions and the
   * first derivatives evaluated at the respective point in the parameter space
   */
  template <CORE::FE::CellType distype,
      std::enable_if_t<CORE::FE::use_lagrange_shapefnct<distype>, bool> = true>
  ShapeFunctionsAndDerivatives<distype> EvaluateShapeFunctionsAndDerivs(
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
      const ElementNodes<distype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<distype> shapefcns;
    CORE::DRT::UTILS::shape_function<distype>(xi, shapefcns.shapefunctions_);
    CORE::DRT::UTILS::shape_function_deriv1<distype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_nurbs<distype>, bool> = true>
  ShapeFunctionsAndDerivatives<distype> EvaluateShapeFunctionsAndDerivs(
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
      const ElementNodes<distype>& nodal_coordinates)
  {
    ShapeFunctionsAndDerivatives<distype> shapefcns;
    CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv(shapefcns.shapefunctions_,
        shapefcns.derivatives_, xi, nodal_coordinates.knots_, nodal_coordinates.weights_, distype);

    return shapefcns;
  }

  template <CORE::FE::CellType distype>
  struct JacobianMapping
  {
    /// Determinant of the jacobian
    double determinant_;

    /// Jacobian matrix at a specific point
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> jacobian_;

    /// Inverse jacobian matrix at a specific point
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> inverse_jacobian_;

    /// Derivative of the shape functions w.r.t. the reference coordinates
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_nodes<distype>> N_XYZ_;
  };

  /*!
   * @brief Evaluates the jacobian mapping of the element
   *
   * @tparam distype : Discretization type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param gp (in) : Id of the Gauss point
   * @return JacobianMapping<distype> : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ, integration
   * factor)
   */
  template <CORE::FE::CellType distype>
  JacobianMapping<distype> EvaluateJacobianMapping(
      const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const ElementNodes<distype>& nodal_coordinates)
  {
    JacobianMapping<distype> jacobian;

    jacobian.jacobian_.Multiply(shapefcns.derivatives_, nodal_coordinates.reference_coordinates_);
    jacobian.inverse_jacobian_ = jacobian.jacobian_;
    jacobian.determinant_ = jacobian.inverse_jacobian_.Invert();
    jacobian.N_XYZ_.Multiply(jacobian.inverse_jacobian_, shapefcns.derivatives_);

    return jacobian;
  }

  /*!
   * @brief Evaluates the jacobian determinant of the element
   *
   * @tparam distype : Discretization type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Jacobian determinant
   */
  template <CORE::FE::CellType distype>
  double EvaluateJacobianDeterminant(const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const ElementNodes<distype>& nodal_coordinates)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> jacobian;
    jacobian.Multiply(shapefcns.derivatives_, nodal_coordinates.reference_coordinates_);

    return jacobian.Determinant();
  }

  /*!
   * @brief Evaluate the jacobian mapping at the element centroid
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return JacobianMapping<distype> : jacobian mapping at the element centroid
   */
  template <CORE::FE::CellType distype, std::enable_if_t<DETAIL::num_dim<distype> == 3, int> = 0>
  JacobianMapping<distype> EvaluateJacobianMappingCentroid(
      const ElementNodes<distype>& nodal_coordinates)
  {
    // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<distype>();

    // shape functions and derivatives evaluated at element centroid
    const ShapeFunctionsAndDerivatives<distype> shape_functions_centroid =
        EvaluateShapeFunctionsAndDerivs<distype>(xi_centroid, nodal_coordinates);

    // jacobian mapping evaluated at centroid
    const JacobianMapping<distype> jacobian_mapping_centroid =
        EvaluateJacobianMapping(shape_functions_centroid, nodal_coordinates);

    return jacobian_mapping_centroid;
  }

  /*!
   * @brief Evaluates a point's coordinates in reference configuration
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates_reference (in) : Reference coordinates of the nodes of the element
   * @param shape_functions_point (in) : Shape functions evaluated at the specific point
   * @return CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> : point's reference coordinates
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> EvaluateReferenceCoordinate(
      const CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, DETAIL::num_dim<distype>>&
          nodal_coordinates_reference,
      const CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1>& shape_functions_point)
  {
    CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> coordinates_reference(true);
    coordinates_reference.MultiplyTN(shape_functions_point, nodal_coordinates_reference);

    return coordinates_reference;
  }

  /*!
   * @brief Evaluates the element centroid's coordinates in reference configuration
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference coordinates of the nodes of the element
   * @return CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> : Element centroid's coordinates in
   * reference configuration
   */
  template <CORE::FE::CellType distype,
      std::enable_if_t<CORE::FE::use_lagrange_shapefnct<distype>, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<distype>& nodal_coordinates)
  {
    const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<distype>();

    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1> shape_functions_centroid(true);
    CORE::DRT::UTILS::shape_function<distype>(xi_centroid, shape_functions_centroid);

    const CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> centroid_coordinates_reference =
        EvaluateReferenceCoordinate<distype>(
            nodal_coordinates.reference_coordinates_, shape_functions_centroid);

    return centroid_coordinates_reference;
  }

  template <CORE::FE::CellType distype, std::enable_if_t<CORE::FE::is_nurbs<distype>, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> EvaluateReferenceCoordinateCentroid(
      const ElementNodes<distype>& nodal_coordinates)
  {
    const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<distype>();

    CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1> shape_functions_centroid(true);
    CORE::DRT::NURBS::UTILS::nurbs_shape_function_dim(shape_functions_centroid, xi_centroid,
        nodal_coordinates.knots_, nodal_coordinates.weights_, distype);

    const CORE::LINALG::Matrix<1, DETAIL::num_dim<distype>> centroid_coordinates_reference =
        EvaluateReferenceCoordinate<distype>(
            nodal_coordinates.reference_coordinates_, shape_functions_centroid);

    return centroid_coordinates_reference;
  }

  template <CORE::FE::CellType distype>
  struct SpatialMaterialMapping
  {
    double determinant_deformation_gradient_;
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> deformation_gradient_;
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>
        inverse_deformation_gradient_;
  };

  /*!
   * @brief Evaluates the mapping between the spatial and material configuration of the element
   *
   * @tparam distype : Discretization type
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param scale_defgrd (in) : scaling for deformation gradient
   * @param kinematictype (in) : kinematic type of element
   * @return SpatialMaterialMapping<distype> : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   */
  template <CORE::FE::CellType distype>
  SpatialMaterialMapping<distype> EvaluateSpatialMaterialMapping(
      const JacobianMapping<distype>& jacobian_mapping,
      const ElementNodes<distype>& nodal_coordinates, const double scale_defgrd = 1.0,
      const INPAR::STR::KinemType& kinematictype = INPAR::STR::kinem_nonlinearTotLag)
  {
    SpatialMaterialMapping<distype> spatial_material_mapping;
    spatial_material_mapping.deformation_gradient_.Clear();

    if (kinematictype == INPAR::STR::kinem_linear)
    {
      for (int i = 0; i < DETAIL::num_dim<distype>; i++)
        spatial_material_mapping.deformation_gradient_(i, i) = 1.0;
    }
    else
    {
      spatial_material_mapping.deformation_gradient_.MultiplyTT(
          scale_defgrd, nodal_coordinates.current_coordinates_, jacobian_mapping.N_XYZ_);
    }
    spatial_material_mapping.inverse_deformation_gradient_.Invert(
        spatial_material_mapping.deformation_gradient_);
    spatial_material_mapping.determinant_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_.Determinant();

    return spatial_material_mapping;
  }

  /*!
   * @brief Evaluates Green-Lagrange strain from right Cauchy-Green tensor
   *
   * GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
   *
   * @tparam distype : Discretization type
   * @param cauchygreen (in) : Right Cauchy-Green deformation tensor
   * @return CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> : Green-Lagrange strain tensor in
   * strain-like Voigt notation
   */
  template <CORE::FE::CellType distype, std::enable_if_t<DETAIL::num_dim<distype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> EvaluateGreenLagrangeStrain(
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& cauchygreen)
  {
    CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> gl_strain;

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
   * @tparam distype: Discretization type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> : Right
   * Cauchy-Green deformation tensor
   */
  template <CORE::FE::CellType distype, std::enable_if_t<DETAIL::num_dim<distype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> EvaluateCauchyGreen(
      const SpatialMaterialMapping<distype>& spatial_material_mapping)
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>> cauchygreen(false);

    cauchygreen.MultiplyTN(spatial_material_mapping.deformation_gradient_,
        spatial_material_mapping.deformation_gradient_);

    return cauchygreen;
  }

  /*!
   * @brief Evaluate the determinant of the deformation gradient at the element centroid
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Determinant of the deformation gradient at the centroid
   */
  template <CORE::FE::CellType distype, std::enable_if_t<DETAIL::num_dim<distype> == 3, int> = 0>
  double EvaluateDeformationGradientDeterminantCentroid(
      const ElementNodes<distype>& nodal_coordinates)
  {
    // jacobian mapping at centroid of element
    const JacobianMapping<distype> jacobian_mapping_centroid =
        EvaluateJacobianMappingCentroid(nodal_coordinates);

    // deformation gradient and strains at centroid of element
    const DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping_centroid =
        EvaluateSpatialMaterialMapping(jacobian_mapping_centroid, nodal_coordinates);

    return spatial_material_mapping_centroid.determinant_deformation_gradient_;
  }

  /*!
   * @brief Evaluates the strain gradient (B-Operator) of the specified element
   *
   * @tparam distype : Discretization type
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return CORE::LINALG::Matrix<num_str<distype>, num_dim<distype> * num_nodes<distype>> :
   * B-Operator
   */
  template <CORE::FE::CellType distype, std::enable_if_t<DETAIL::num_dim<distype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_str<distype>,
      DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>
  EvaluateStrainGradient(const JacobianMapping<distype>& jacobian_mapping,
      const SpatialMaterialMapping<distype>& spatial_material_mapping)
  {
    // B-operator
    CORE::LINALG::Matrix<DETAIL::num_str<distype>,
        DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>
        Bop;
    for (int i = 0; i < DETAIL::num_nodes<distype>; ++i)
    {
      for (int d = 0; d < DETAIL::num_dim<distype>; ++d)
      {
        for (int e = 0; e < DETAIL::num_dim<distype>; ++e)
        {
          Bop(d, DETAIL::num_dim<distype> * i + e) =
              spatial_material_mapping.deformation_gradient_(e, d) * jacobian_mapping.N_XYZ_(d, i);
        }
      }

      Bop(3, DETAIL::num_dim<distype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, DETAIL::num_dim<distype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, DETAIL::num_dim<distype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(4, DETAIL::num_dim<distype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, DETAIL::num_dim<distype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, DETAIL::num_dim<distype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(5, DETAIL::num_dim<distype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, DETAIL::num_dim<distype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, DETAIL::num_dim<distype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(2, i);
    }

    return Bop;
  }

  template <CORE::FE::CellType distype>
  struct Stress
  {
    /// Second Piola-Kirchhoff stress tensor in stress-like voigt notation
    CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> pk2_;

    /// Linearization of the 2. Piola Kirchhoff stress tensor w.r.t. Green-Lagrange strain tensor in
    /// mixed Voigt notation
    CORE::LINALG::Matrix<DETAIL::num_str<distype>, DETAIL::num_str<distype>> cmat_;
  };

  /*!
   * @brief Evaluates the material stress (2. Piola-Kirchhoff stress tensor and the linearization
   * w.r.t. Green-Lagrange strain)
   *
   * @tparam distype : Discretization type
   * @param material (in) : Reference to the material
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @param gl_strain (in) : Green-Lagrange strain
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<distype> : Object holding the 2. Piola-Kirchhoff stress tensor and the
   * linearization w.r.t. Green-Lagrange strain tensor
   */
  template <CORE::FE::CellType distype>
  Stress<distype> EvaluateMaterialStress(MAT::So3Material& material,
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& defgrd,
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1>& gl_strain,
      Teuchos::ParameterList& params, const int gp, const int eleGID)
  {
    Stress<distype> stress;

    material.Evaluate(&defgrd, &gl_strain, params, &stress.pk2_, &stress.cmat_, gp, eleGID);
    return stress;
  }

  /*!
   * @brief Evaluates the material stress (2. Piola-Kirchhoff stress tensor and the linearization
   * w.r.t. Green-Lagrange strain)
   *
   * @tparam distype : Discretization type
   * @param material (in) : Reference to the material
   * @param strains (in) : Strain measures of the element
   * @param strains_old (in) : Strain measures of the element old
   * @param strains_m (in) : Strain measures of the element middle
   * @param density (in) : Density
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<distype> : Object holding the 2. Piola-Kirchhoff stress tensor and the
   * linearization w.r.t. Green-Lagrange strain tensor
   */
  template <CORE::FE::CellType distype>
  Stress<distype> EvaluateMaterialStressGEMM(MAT::So3Material& material,
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1>& gl_strains,
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1>& gl_strains_old,
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1>& gl_strains_m,
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& rcg,
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& rcg_old,
      const Teuchos::ParameterList& params, const int gp, const int eleGID)
  {
    Stress<distype> stress;

    material.EvaluateGEMM(&stress.pk2_, &stress.cmat_, &gl_strains_m, &gl_strains, &gl_strains_old,
        &rcg, &rcg_old, gp, eleGID);

    return stress;
  }

  /*!
   * @brief Adds the internal force vector contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddInternalForceVector(const CORE::LINALG::Matrix<DETAIL::num_str<distype>,
                                  DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& Bop,
      const Stress<distype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>, 1>& force_vector)
  {
    force_vector.MultiplyTN(integration_fac, Bop, stress.pk2_, 1.);
  }

  /*!
   * @brief Add elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddElasticStiffnessMatrix(const CORE::LINALG::Matrix<DETAIL::num_str<distype>,
                                     DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& Bop,
      const Stress<distype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& stiffness_matrix)
  {
    CORE::LINALG::Matrix<DETAIL::num_str<distype>,
        DETAIL::num_nodes<distype> * DETAIL::num_dim<distype>>
        cb;
    cb.Multiply(stress.cmat_, Bop);
    stiffness_matrix.MultiplyTN(integration_fac, Bop, cb, 1.0);
  }

  /*!
   * @brief Add GEMM elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param BopM (in) : Mid point strain gradient (B-Operator at mid point)
   * @param stress (in) : Stress measures
   * @param pre_fac (in) : Pre-factor = detJ_w * (1 - alphaf + xi) with the integration factor
   * detJ_w and the GEMM coefficients alphaf and xi
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddElasticStiffnessMatrixGEMM(
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& Bop,
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& BopM,
      const Stress<distype>& stress, const double pre_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& stiffness_matrix)
  {
    // GEMM contribution: (1.0 - gemmalphaf + gemmxi) * detJ_w * (BopM^T . cmat . Bop)
    CORE::LINALG::Matrix<DETAIL::num_str<distype>,
        DETAIL::num_nodes<distype> * DETAIL::num_dim<distype>>
        cb;
    cb.Multiply(stress.cmat_, Bop);
    stiffness_matrix.MultiplyTN(pre_fac, BopM, cb, 1.0);
  }

  /*!
   * @brief Add geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param B_L (in) : B_L operator, i.e. derivatives of the shape functions w.r.t. XYZ
   *                   at the respective Gauss point
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddGeometricStiffnessMatrix(
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_nodes<distype>>& B_L,
      const Stress<distype>& stress, const double integration_fac,
      CORE::LINALG::Matrix<DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& stiffness_matrix)
  {
    std::array<double, 3> SmB_L;  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < DETAIL::num_nodes<distype>; ++inod)
    {
      SmB_L[0] = stress.pk2_(0) * B_L(0, inod) + stress.pk2_(3) * B_L(1, inod) +
                 stress.pk2_(5) * B_L(2, inod);
      SmB_L[1] = stress.pk2_(3) * B_L(0, inod) + stress.pk2_(1) * B_L(1, inod) +
                 stress.pk2_(4) * B_L(2, inod);
      SmB_L[2] = stress.pk2_(5) * B_L(0, inod) + stress.pk2_(4) * B_L(1, inod) +
                 stress.pk2_(2) * B_L(2, inod);

      for (int jnod = 0; jnod < DETAIL::num_nodes<distype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < DETAIL::num_dim<distype>; ++idim)
          bopstrbop += B_L(idim, jnod) * SmB_L[idim];

        for (int d = 0; d < DETAIL::num_dim<distype>; ++d)
          stiffness_matrix(DETAIL::num_dim<distype> * inod + d,
              DETAIL::num_dim<distype> * jnod + d) += integration_fac * bopstrbop;
      }
    }
  }

  /*!
   * @brief Add mass matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param shapefunctions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param integration_factor (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param density (in) : density at the Gauss point
   * @param mass (in/out) : mass matrix where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddMassMatrix(const ShapeFunctionsAndDerivatives<distype>& shapefunctions,
      const double integration_factor, const double density,
      CORE::LINALG::Matrix<DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>,
          DETAIL::num_dim<distype> * DETAIL::num_nodes<distype>>& mass)
  {
    for (int inod = 0; inod < DETAIL::num_nodes<distype>; ++inod)
    {
      const double ifactor = shapefunctions.shapefunctions_(inod) * integration_factor * density;
      for (int jnod = 0; jnod < DETAIL::num_nodes<distype>; ++jnod)
      {
        const double massfactor =
            shapefunctions.shapefunctions_(jnod) * ifactor;  // intermediate factor
        for (int d = 0; d < DETAIL::num_dim<distype>; ++d)
          mass(DETAIL::num_dim<distype> * inod + d, DETAIL::num_dim<distype> * jnod + d) +=
              massfactor;
      }
    }
  }

  /*!
   * @brief Assemble a vector into a matrix row
   *
   * @tparam num_str
   * @param vector (in) : Vector to be assembled into matrix
   * @param data (in/out) : Matrix the vector is assembled into
   * @param row (in) : Matrix row
   */
  template <unsigned num_str>
  void AssembleVectorToMatrixRow(
      CORE::LINALG::Matrix<num_str, 1> vector, CORE::LINALG::SerialDenseMatrix& data, const int row)
  {
    for (unsigned i = 0; i < num_str; ++i) data(row, static_cast<int>(i)) = vector(i);
  }

  /*!
   * @brief Convert Green-Lagrange strains to the desired strain type and assemble to a given matrix
   * row in stress-like Voigt notation
   *
   * @tparam distype : Discretization type
   * @param gl_strain (in) : Green-Lagrange strain
   * @param defgrd (in) : Deformation gradient
   * @param strain_type (in) : Strain type, i.e., Green-Lagrange or Euler-Almansi
   * @param data (in/out) : Matrix the strains are assembled into
   * @param row (in) : Matrix row
   */
  template <CORE::FE::CellType distype>
  void AssembleStrainTypeToMatrixRow(
      const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1>& gl_strain,
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& defgrd,
      const INPAR::STR::StrainType strain_type, CORE::LINALG::SerialDenseMatrix& data,
      const int row)
  {
    switch (strain_type)
    {
      case INPAR::STR::strain_gl:
      {
        CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> gl_strain_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(gl_strain, gl_strain_stress_like);
        AssembleVectorToMatrixRow(gl_strain_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_ea:
      {
        const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> ea =
            STR::UTILS::GreenLagrangeToEulerAlmansi(gl_strain, defgrd);
        CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> ea_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(ea, ea_stress_like);
        AssembleVectorToMatrixRow(ea_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_none:
        return;
      default:
        dserror("strain type not supported");
        break;
    }
  }

  /*!
   * @brief Convert 2nd Piola-Kirchhoff stresses to the desired stress type and assemble to a given
   * matrix row in stress-like Voigt notation
   *
   * @tparam distype : Discretization type
   * @param defgrd (in) : Deformation gradient
   * @param stress (in) : 2nd Piola-Kirchhoff stress
   * @param stress_type (in) : Stress type, i.e., 2nd Piola-Kirchhoff or Cauchy
   * @param data (in/out) : Matrix the stresses are assembled into
   * @param row (in) : Matrix row
   */
  template <CORE::FE::CellType distype>
  void AssembleStressTypeToMatrixRow(
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>& defgrd,
      const Stress<distype>& stress, const INPAR::STR::StressType stress_type,
      CORE::LINALG::SerialDenseMatrix& data, const int row)
  {
    switch (stress_type)
    {
      case INPAR::STR::stress_2pk:
      {
        AssembleVectorToMatrixRow(stress.pk2_, data, row);
        return;
      }
      case INPAR::STR::stress_cauchy:
      {
        CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> cauchy;
        STR::UTILS::Pk2ToCauchy(stress.pk2_, defgrd, cauchy);
        AssembleVectorToMatrixRow(cauchy, data, row);
        return;
      }
      case INPAR::STR::stress_none:

        return;
      default:
        dserror("stress type not supported");
        break;
    }
  }

  /*!
   * @brief Serialize a matrix by conversion to a vector representation

   * @param matrix (in) : Matrix
   * @param serialized_matrix (in/out) : Serialized matrix
   */
  inline void Serialize(
      const CORE::LINALG::SerialDenseMatrix& matrix, std::vector<char>& serialized_matrix)
  {
    DRT::PackBuffer packBuffer;
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    packBuffer.StartPacking();
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    std::copy(packBuffer().begin(), packBuffer().end(), std::back_inserter(serialized_matrix));
  }

  /*!
   * @brief Calls the @p gp_evaluator for each Gauss point with evaluated jacobian mapping using the
   * integration rule defined by @p integration.
   *
   * @tparam distype : Discretization type known at compile time
   * @tparam GaussPointEvaluator
   * @param nodal_coordinates (in) : The nodal coordinates of the element
   * @param integration (in) : The integration rule to be used.
   * @param gp_evaluator (in) : A callable object (e.g. lambda-function) with signature void(const
   * CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi, const
   * ShapeFunctionsAndDerivatives<distype>& shape_functions, const JacobianMapping<distype>&
   * jacobian_mapping, double integration_factor, int gp) that will be called for each integration
   * point.
   */
  template <CORE::FE::CellType distype, typename GaussPointEvaluator>
  inline void ForEachGaussPoint(const ElementNodes<distype>& nodal_coordinates,
      const CORE::DRT::UTILS::GaussIntegration& integration, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < integration.NumPoints(); ++gp)
    {
      const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi =
          EvaluateParameterCoordinate<distype>(integration, gp);

      const ShapeFunctionsAndDerivatives<distype> shape_functions =
          EvaluateShapeFunctionsAndDerivs<distype>(xi, nodal_coordinates);

      const JacobianMapping<distype> jacobian_mapping =
          EvaluateJacobianMapping(shape_functions, nodal_coordinates);

      const double integration_factor = jacobian_mapping.determinant_ * integration.Weight(gp);

      gp_evaluator(xi, shape_functions, jacobian_mapping, integration_factor, gp);
    }
  }

  /*!
   * @brief Evaluates the Gauss point coordinates in reference configuration
   * and adds those to the paramater list
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param shape_functions_gp (in) : Shape functions evaluated at the Gauss point
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <CORE::FE::CellType distype>
  void EvaluateGPCoordinatesAndAddToParameterList(const ElementNodes<distype>& nodal_coordinates,
      const ShapeFunctionsAndDerivatives<distype>& shape_functions_gp,
      Teuchos::ParameterList& params)
  {
    auto gp_ref_coord = EvaluateReferenceCoordinate<distype>(
        nodal_coordinates.reference_coordinates_, shape_functions_gp.shapefunctions_);
    params.set("gprefecoord", gp_ref_coord);
  }

  /*!
   * @brief Evaluates the element centroid coordinates in reference configuration
   * and adds those to the paramater list
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <CORE::FE::CellType distype>
  void EvaluateCentroidCoordinatesAndAddToParameterList(
      const ElementNodes<distype>& nodal_coordinates, Teuchos::ParameterList& params)
  {
    auto element_center = EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates);
    params.set("elecenter", element_center);
  }

  /*!
   * @brief Asks the material for the Gauss Point output quantities and adds the information to
   * the Gauss point output data manager
   *
   * @param num_gp (in) : Number of Gauss Points of the element
   * @param solid_material (in) : Solid material of the element
   * @param gp_data_output_manager (in/out) : Gauss point data output manager
   *                                          (only for new structure time integration)
   */
  inline void AskAndAddQuantitiesToGaussPointDataOutput(const int num_gp,
      const MAT::So3Material& solid_material,
      STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Save number of Gauss Points of the element for gauss point data output
    gp_data_output_manager.AddElementNumberOfGaussPoints(num_gp);

    // holder for output quantity names and their size
    std::unordered_map<std::string, int> quantities_map{};

    // Ask material for the output quantity names and sizes
    solid_material.RegisterOutputDataNames(quantities_map);

    // Add quantities to the Gauss point output data manager (if they do not already exist)
    gp_data_output_manager.MergeQuantities(quantities_map);
  }

  /*!
   * @brief Collect Gauss Point output data from material and assemble/interpolate depending on
   * output type to element center, Gauss Points, or nodes
   *
   * @tparam distype : Discretization type
   * @param stiffness_matrix_integration (in) : Container holding the integration points
   * @param solid_material (in) : Solid material of the element
   * @param ele (in) : Reference to the element
   * @param gp_data_output_manager (in/out) : Gauss point data output manager
   *                                          (only for new structure time integration)
   */
  template <CORE::FE::CellType distype>
  inline void CollectAndAssembleGaussPointDataOutput(
      const CORE::DRT::UTILS::GaussIntegration& stiffness_matrix_integration,
      const MAT::So3Material& solid_material, const DRT::Element& ele,
      STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Collection and assembly of gauss point data
    for (const auto& quantity : gp_data_output_manager.GetQuantities())
    {
      const std::string& quantity_name = quantity.first;
      const int quantity_size = quantity.second;

      // Step 1: Collect the data for each Gauss point for the material
      CORE::LINALG::SerialDenseMatrix gp_data(
          stiffness_matrix_integration.NumPoints(), quantity_size, true);
      bool data_available = solid_material.EvaluateOutputData(quantity_name, gp_data);

      // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
      // point)
      if (data_available)
      {
        switch (gp_data_output_manager.GetOutputType())
        {
          case INPAR::STR::GaussPointDataOutputType::element_center:
          {
            // compute average of the quantities
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.GetElementCenterData().at(quantity_name);
            CORE::DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::nodes:
          {
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.GetNodalData().at(quantity_name);

            Epetra_IntVector& global_nodal_element_count =
                *gp_data_output_manager.GetNodalDataCount().at(quantity_name);

            CORE::DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<distype>(
                ele, gp_data, *global_data, false, stiffness_matrix_integration);
            DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::gauss_points:
          {
            std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                gp_data_output_manager.GetGaussPointData().at(quantity_name);
            DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::none:
            dserror(
                "You specified a Gauss point data output type of none, so you should not end up "
                "here.");
          default:
            dserror("Unknown Gauss point data output type.");
        }
      }
    }
  }

  /*!
   * @brief For elements with fiber nodes, interpolate fibers to Gauss points and add to the
   * parameter list
   *
   * @tparam distype : Discretization type
   * @param stiffness_matrix_integration (in) : Container holding the integration points
   * @param ele (in) : Reference to the element, possibly having nodal fibers
   * @param params (in/out) : ParameterList the interpolated fibers are added to
   */
  template <CORE::FE::CellType distype>
  inline void InterpolateFibersToGaussPointsAndAddToParameterList(
      const CORE::DRT::UTILS::GaussIntegration& stiffness_matrix_integration,
      const DRT::Element& ele, Teuchos::ParameterList& params)
  {
    if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
    {
      // This element has fiber nodes.
      // Interpolate fibers to the Gauss points and add them to the parameter list

      // Get shape functions
      const static std::vector<CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1>> shapefcts =
          std::invoke(
              [&]
              {
                std::vector<CORE::LINALG::Matrix<DETAIL::num_nodes<distype>, 1>> shapefcns(
                    stiffness_matrix_integration.NumPoints());
                for (int gp = 0; gp < stiffness_matrix_integration.NumPoints(); ++gp)
                {
                  CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1> xi(
                      stiffness_matrix_integration.Point(gp), true);
                  CORE::DRT::UTILS::shape_function<distype>(xi, shapefcns[gp]);
                }
                return shapefcns;
              });

      // add fibers to the ParameterList
      DRT::FIBER::NodalFiberHolder fiberHolder;

      // Do the interpolation
      DRT::FIBER::UTILS::ProjectFibersToGaussPoints<distype>(ele.Nodes(), shapefcts, fiberHolder);

      params.set("fiberholder", fiberHolder);
    }
  }

}  // namespace DRT::ELEMENTS

#endif  // SOLID_ELE_CALC_LIB_H