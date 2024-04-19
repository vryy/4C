/*! \file

\brief A library of free functions for the shell7p element calculation

\level 3
*/

#ifndef FOUR_C_SHELL7P_ELE_CALC_LIB_HPP
#define FOUR_C_SHELL7P_ELE_CALC_LIB_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_lib_element_integration_select.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_fixedsizematrix_tensor_transformation.hpp"
#include "baci_mat_so3_material.hpp"
#include "baci_shell7p_ele.hpp"

#include <Teuchos_ParameterList.hpp>

#include <numeric>

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS::SHELL
{
  /*!
   * @brief An object holding the nodal coordinates in reference and current configuration
   *
   * @tparam distype : The discretization type known at compile time
   */
  template <CORE::FE::CellType distype>
  struct NodalCoordinates
  {
    CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim> x_refe_;
    CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim> x_curr_;
    CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim> a3_refe_;
    CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim> a3_curr_;
  };

  /*!
   * \brief Create matrix with spatial configuration
   *
   * @tparam distype :  The discretization type known at compile time
   * @param x       (out)  : Nodal coords in spatial frame
   * @param x_refe  (out)  : Nodal coords in reference frame
   * @param disp    (int)  : Displacements
   * @param index   (int)  : Integer to consider either displacements or director displacmenets
   */
  template <CORE::FE::CellType distype>
  void SpatialConfiguration(CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim>& x,
      const CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim>& x_refe,
      const std::vector<double>& disp, const int index)
  {
    const int nodedof = DRT::ELEMENTS::SHELL::DETAIL::node_dof;
    for (int i = 0; i < DETAIL::num_node<distype>; ++i)
    {
      x(i, 0) = x_refe(i, 0) + disp[i * nodedof + 0 + index];
      x(i, 1) = x_refe(i, 1) + disp[i * nodedof + 1 + index];
      x(i, 2) = x_refe(i, 2) + disp[i * nodedof + 2 + index];
    }
  }
  /*!
   * \brief Evaluates coordinates of element nodes at the shell surface and approximated nodal
   * points at the shell thickness. The thickness is represented by nodal directors a3_reference.
   *
   * @tparam distype :  The discretization type known at compile time
   * @param nodes (in)  : Nodal coords in reference frame
   * @param disp (in)  : Displacements
   * @param thicknesss (int)  : Nodal thickness values in refernce frame
   * @param a3_reference (int)  : Nodal directors in reference frame
   * @param factor (int)  : Scaling factor due to SDC
   */
  template <CORE::FE::CellType distype>
  DRT::ELEMENTS::SHELL::NodalCoordinates<distype> EvaluateNodalCoordinates(DRT::Node** nodes,
      std::vector<double>& disp, const double& thickness,
      const CORE::LINALG::SerialDenseMatrix& a3_reference, const double factor)
  {
    DRT::ELEMENTS::SHELL::NodalCoordinates<distype> coordinates;
    for (auto i = 0; i < DRT::ELEMENTS::SHELL::DETAIL::num_node<distype>; ++i)
    {
      const double h2 = thickness * factor * 0.5;

      const auto& x = nodes[i]->X();
      coordinates.x_refe_(i, 0) = x[0];
      coordinates.x_refe_(i, 1) = x[1];
      coordinates.x_refe_(i, 2) = x[2];

      coordinates.a3_refe_(i, 0) = a3_reference(i, 0) * h2;
      coordinates.a3_refe_(i, 1) = a3_reference(i, 1) * h2;
      coordinates.a3_refe_(i, 2) = a3_reference(i, 2) * h2;

      const int nodedof = DETAIL::node_dof;
      coordinates.x_curr_(i, 0) = coordinates.x_refe_(i, 0) + disp[i * nodedof + 0];
      coordinates.x_curr_(i, 1) = coordinates.x_refe_(i, 1) + disp[i * nodedof + 1];
      coordinates.x_curr_(i, 2) = coordinates.x_refe_(i, 2) + disp[i * nodedof + 2];

      coordinates.a3_curr_(i, 0) = coordinates.a3_refe_(i, 0) + disp[i * nodedof + 3];
      coordinates.a3_curr_(i, 1) = coordinates.a3_refe_(i, 1) + disp[i * nodedof + 4];
      coordinates.a3_curr_(i, 2) = coordinates.a3_refe_(i, 2) + disp[i * nodedof + 5];
    }
    return coordinates;
  }

  /*!
   * @brief A struct holding the the strains measures (deformation gradient and Green-lagrange
   * strains)
   */
  struct Strains
  {
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> defgrd_;
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> gl_strain_;
  };

  /*!
   * @brief A struct holding the the stress measures (2-Piola-Kirchhoff stresses and elasticity
   * matrix)
   */
  template <int numint>
  struct Stress
  {
    CORE::LINALG::Matrix<numint, 1> pk2_;
    CORE::LINALG::Matrix<numint, numint> cmat_;
  };

  /*!
   * @brief A struct holding the enhanced stress measures (2-Piola-Kirchhoff stresses and
   * enhanced material matrix)
   */
  struct StressEnhanced
  {
    CORE::LINALG::SerialDenseVector stress_;
    CORE::LINALG::SerialDenseMatrix dmat_;
  };


  /*!
   * @brief Evaluates the material stress in a global cartesian coordinate system (2. Piola
   * Kirchhoff stress tensor and the linearization w.r.t Green-Lagrange strain)
   *
   * @tparam dim : Dimension
   * @param material (in) : Reference to the material
   * @param strains (in) : Strain measures of the element
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<MAT::NUM_STRESS_3D> : Object holding the 2. Piola Kirchhoff stress tensor and
   * the linearization w.r.t. Green Lagrange strain tensor
   */
  template <int dim>
  Stress<MAT::NUM_STRESS_3D> EvaluateMaterialStressCartesianSystem(MAT::So3Material& material,
      const Strains& strains, Teuchos::ParameterList& params, int gp, int eleGID)
  {
    if (dim != 3) FOUR_C_THROW("stop: this currently only works for 3D");
    DRT::ELEMENTS::SHELL::Stress<MAT::NUM_STRESS_3D> stress;

    material.Evaluate(
        &strains.defgrd_, &strains.gl_strain_, params, &stress.pk2_, &stress.cmat_, gp, eleGID);

    return stress;
  }

  /*!
   * @brief Returns the optimal gauss integration rule based on the element discretization type
   */
  template <CORE::FE::CellType distype>
  constexpr auto GetGaussRule()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  }

  /*!
   * @brief Returns the gauss integration based on the optimal gauss integration rule
   *
   * @tparam distype : The discretization type known at compile time
   * @tparam GaussRuleType (in) : Optimal gauss rule integration type
   * @return CORE::FE::IntegrationPoints2D : Integration points
   */
  template <CORE::FE::CellType distype, typename GaussRuleType>
  CORE::FE::IntegrationPoints2D CreateGaussIntegrationPoints(GaussRuleType gaussrule)
  {
    return CORE::FE::IntegrationPoints2D(gaussrule);
  }

  /*!
   * \brief Calculate the deformation gradient that is consistent to the EAS enhanced strains.
   *
   * @tparam dim: dimension
   *
   * @param defgrd_disp [in]: displacement based deformation gradient
   * @param glstrain_enh [in]: enhanced strains
   * @param defgrd_enh [in/out]: enhanced (consistent) deformation gradient
   */
  template <unsigned dim>
  void CalcConsistentDefgrd(const CORE::LINALG::Matrix<dim, dim>& defgrd_disp,
      const CORE::LINALG::Matrix<dim*(dim + 1) / 2, 1>& glstrain_enh,
      CORE::LINALG::Matrix<dim, dim>& defgrd_enh)
  {
    CORE::LINALG::Matrix<dim, dim> R;       // rotation tensor
    CORE::LINALG::Matrix<dim, dim> U_enh;   // modified right stretch tensor
    CORE::LINALG::Matrix<dim, dim> U_disp;  // displacement-based right stretch tensor
    CORE::LINALG::Matrix<dim, dim> EW;      // temporarily store eigenvalues
    CORE::LINALG::Matrix<dim, dim> tmp;     // temporary matrix for matrix matrix matrix products
    CORE::LINALG::Matrix<dim, dim> tmp2;    // temporary matrix for matrix matrix matrix products

    // We calculate the "enhanced" deformation gradient from the enhanced GL strains with the help
    // of two polar decompositions
    if (dim != 3) FOUR_C_THROW("stop: this currently only works for 3D");
    // First step: calculate enhanced right stretch tensor  U_enh from C_enh=U_enh^T*U_enh
    // -> get C_enh from enhanced GL strains
    for (unsigned i = 0; i < dim; i++) U_enh(i, i) = 2. * glstrain_enh(i) + 1.;
    U_enh(0, 1) = glstrain_enh(dim);
    U_enh(1, 0) = glstrain_enh(dim);
    U_enh(1, 2) = glstrain_enh(4);
    U_enh(2, 1) = glstrain_enh(4);
    U_enh(0, 2) = glstrain_enh(5);
    U_enh(2, 0) = glstrain_enh(5);

    CORE::LINALG::SYEV(U_enh, EW, U_enh);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_enh, EW);
    tmp2.MultiplyNT(tmp, U_enh);
    U_enh.Update(tmp2);

    // Second step: calculate displacement-based right stretch tensor
    U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

    CORE::LINALG::SYEV(U_disp, EW, U_disp);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_disp, EW);
    tmp2.MultiplyNT(tmp, U_disp);
    U_disp.Update(tmp2);

    // Third step: compose consistent deformation gradient
    U_disp.Invert();
    R.Multiply(defgrd_disp, U_disp);
    defgrd_enh.Multiply(R, U_enh);

    const double det_defgrd_enh = defgrd_enh.Determinant();
    if (det_defgrd_enh <= 0.0)
    {
      FOUR_C_THROW("Negative jacobian determinant of modified deformation gradient");
    }
  }

  /*!
   * @brief A struct holding the mass matrix integration factors
   */
  struct MassMatrixVariables
  {
    double factor_v_ = 0.;
    double factor_w_ = 0.;
    double factor_vw_ = 0.;
  };

  /*!
   * @brief An object holding the shape functions and the first derivatives evaluated at the
   * respective point in the parameter space
   *
   * @tparam distype :  The discretization type known at compile time
   */
  template <CORE::FE::CellType distype>
  struct ShapefunctionsAndDerivatives
  {
    CORE::LINALG::Matrix<DETAIL::num_node<distype>, 1> shapefunctions_;
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_node<distype>> derivatives_;
  };

  /*!
   * @brief Evaluates the shape functions and the first derivatives w.r.t spatial coordinates at the
   * specified point in the parameter space
   *
   * @tparam distype :  The discretization type known at compile time
   * @param xi_gp (in) : Coordinate of the integration point in the parameter space
   * @return ShapeFunctionsAndDerivatives<distype> : An object holding the shape functions and the
   * first derivatives evaluated at the respective point in the parameter space
   */
  template <CORE::FE::CellType distype>
  ShapefunctionsAndDerivatives<distype> EvaluateShapefunctionsAndDerivs(
      const std::array<double, 2>& xi_gp)
  {
    ShapefunctionsAndDerivatives<distype> shapefunctions;
    CORE::FE::shape_function_2D(shapefunctions.shapefunctions_, xi_gp[0], xi_gp[1], distype);
    CORE::FE::shape_function_2D_deriv1(shapefunctions.derivatives_, xi_gp[0], xi_gp[1], distype);
    return shapefunctions;
  }

  /*!
   * @brief An object holding the basis vectors and metric tensors
   *
   * @tparam distype :  The discretization type known at compile time
   */
  template <CORE::FE::CellType distype>
  struct BasisVectorsAndMetrics
  {
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> kovariant_;
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> kontravariant_;
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> metric_kovariant_;
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> metric_kontravariant_;
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> partial_derivative_;
    double detJ_;
  };

  /*!
   * @brief Evaluates the strain gradient (B-Operator) of the specified element
   *
   * @tparam distype :  The discretization type known at compile time
   * @param akov (in) : kovariant basis vectors
   * @param da3kov (in) : Partial derivatives of the kovariant basis vectors w.r.t spatial
   * coordinates
   * @param shapefunctions_derivatives (in) : An object holding the shape functions and the
   * first derivatives w.r.t spatial coordinates evaluated at the respective point in the parameter
   * space
   * @return CORE::LINALG::SerialDenseMatrix: B-operator Matrix
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::SerialDenseMatrix CalcBOperator(
      const CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim>& akov,
      const CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim>& da3kov,
      const DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>& shapefunctions_derivatives)
  {
    const CORE::LINALG::Matrix<DETAIL::num_node<distype>, 1> shapefunctions =
        shapefunctions_derivatives.shapefunctions_;
    const CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_node<distype>>& derivs =
        shapefunctions_derivatives.derivatives_;
    CORE::LINALG::SerialDenseMatrix Bop(
        DETAIL::num_internal_variables, DETAIL::numdofperelement<distype>);
    const int nodedof = DETAIL::node_dof;
    for (int i = 0; i < DETAIL::num_node<distype>; ++i)
    {
      Bop(0, nodedof * i + 0) = derivs(0, i) * akov(0, 0);
      Bop(0, nodedof * i + 1) = derivs(0, i) * akov(0, 1);
      Bop(0, nodedof * i + 2) = derivs(0, i) * akov(0, 2);
      Bop(0, nodedof * i + 3) = 0.0;
      Bop(0, nodedof * i + 4) = 0.0;
      Bop(0, nodedof * i + 5) = 0.0;

      Bop(1, nodedof * i + 0) = derivs(1, i) * akov(0, 0) + derivs(0, i) * akov(1, 0);
      Bop(1, nodedof * i + 1) = derivs(1, i) * akov(0, 1) + derivs(0, i) * akov(1, 1);
      Bop(1, nodedof * i + 2) = derivs(1, i) * akov(0, 2) + derivs(0, i) * akov(1, 2);
      Bop(1, nodedof * i + 3) = 0.0;
      Bop(1, nodedof * i + 4) = 0.0;
      Bop(1, nodedof * i + 5) = 0.0;

      Bop(2, nodedof * i + 0) = derivs(0, i) * akov(2, 0);
      Bop(2, nodedof * i + 1) = derivs(0, i) * akov(2, 1);
      Bop(2, nodedof * i + 2) = derivs(0, i) * akov(2, 2);
      Bop(2, nodedof * i + 3) = shapefunctions(i) * akov(0, 0);
      Bop(2, nodedof * i + 4) = shapefunctions(i) * akov(0, 1);
      Bop(2, nodedof * i + 5) = shapefunctions(i) * akov(0, 2);

      Bop(3, nodedof * i + 0) = derivs(1, i) * akov(1, 0);
      Bop(3, nodedof * i + 1) = derivs(1, i) * akov(1, 1);
      Bop(3, nodedof * i + 2) = derivs(1, i) * akov(1, 2);
      Bop(3, nodedof * i + 3) = 0.0;
      Bop(3, nodedof * i + 4) = 0.0;
      Bop(3, nodedof * i + 5) = 0.0;

      Bop(4, nodedof * i + 0) = derivs(1, i) * akov(2, 0);
      Bop(4, nodedof * i + 1) = derivs(1, i) * akov(2, 1);
      Bop(4, nodedof * i + 2) = derivs(1, i) * akov(2, 2);
      Bop(4, nodedof * i + 3) = shapefunctions(i) * akov(1, 0);
      Bop(4, nodedof * i + 4) = shapefunctions(i) * akov(1, 1);
      Bop(4, nodedof * i + 5) = shapefunctions(i) * akov(1, 2);

      Bop(5, nodedof * i + 0) = 0.0;
      Bop(5, nodedof * i + 1) = 0.0;
      Bop(5, nodedof * i + 2) = 0.0;
      Bop(5, nodedof * i + 3) = shapefunctions(i) * akov(2, 0);
      Bop(5, nodedof * i + 4) = shapefunctions(i) * akov(2, 1);
      Bop(5, nodedof * i + 5) = shapefunctions(i) * akov(2, 2);

      Bop(6, nodedof * i + 0) = derivs(0, i) * da3kov(0, 0);
      Bop(6, nodedof * i + 1) = derivs(0, i) * da3kov(0, 1);
      Bop(6, nodedof * i + 2) = derivs(0, i) * da3kov(0, 2);
      Bop(6, nodedof * i + 3) = derivs(0, i) * akov(0, 0);
      Bop(6, nodedof * i + 4) = derivs(0, i) * akov(0, 1);
      Bop(6, nodedof * i + 5) = derivs(0, i) * akov(0, 2);

      Bop(7, nodedof * i + 0) = derivs(0, i) * da3kov(1, 0) + derivs(1, i) * da3kov(0, 0);
      Bop(7, nodedof * i + 1) = derivs(0, i) * da3kov(1, 1) + derivs(1, i) * da3kov(0, 1);
      Bop(7, nodedof * i + 2) = derivs(0, i) * da3kov(1, 2) + derivs(1, i) * da3kov(0, 2);
      Bop(7, nodedof * i + 3) = derivs(0, i) * akov(1, 0) + derivs(1, i) * akov(0, 0);
      Bop(7, nodedof * i + 4) = derivs(0, i) * akov(1, 1) + derivs(1, i) * akov(0, 1);
      Bop(7, nodedof * i + 5) = derivs(0, i) * akov(1, 2) + derivs(1, i) * akov(0, 2);

      Bop(8, nodedof * i + 0) = 0.0;
      Bop(8, nodedof * i + 1) = 0.0;
      Bop(8, nodedof * i + 2) = 0.0;
      Bop(8, nodedof * i + 3) = shapefunctions(i) * da3kov(0, 0) + derivs(0, i) * akov(2, 0);
      Bop(8, nodedof * i + 4) = shapefunctions(i) * da3kov(0, 1) + derivs(0, i) * akov(2, 1);
      Bop(8, nodedof * i + 5) = shapefunctions(i) * da3kov(0, 2) + derivs(0, i) * akov(2, 2);

      Bop(9, nodedof * i + 0) = derivs(1, i) * da3kov(1, 0);
      Bop(9, nodedof * i + 1) = derivs(1, i) * da3kov(1, 1);
      Bop(9, nodedof * i + 2) = derivs(1, i) * da3kov(1, 2);
      Bop(9, nodedof * i + 3) = derivs(1, i) * akov(1, 0);
      Bop(9, nodedof * i + 4) = derivs(1, i) * akov(1, 1);
      Bop(9, nodedof * i + 5) = derivs(1, i) * akov(1, 2);

      Bop(10, nodedof * i + 0) = 0.0;
      Bop(10, nodedof * i + 1) = 0.0;
      Bop(10, nodedof * i + 2) = 0.0;
      Bop(10, nodedof * i + 3) = shapefunctions(i) * da3kov(1, 0) + derivs(1, i) * akov(2, 0);
      Bop(10, nodedof * i + 4) = shapefunctions(i) * da3kov(1, 1) + derivs(1, i) * akov(2, 1);
      Bop(10, nodedof * i + 5) = shapefunctions(i) * da3kov(1, 2) + derivs(1, i) * akov(2, 2);

      Bop(11, nodedof * i + 0) = 0.0;
      Bop(11, nodedof * i + 1) = 0.0;
      Bop(11, nodedof * i + 2) = 0.0;
      Bop(11, nodedof * i + 3) = 0.0;
      Bop(11, nodedof * i + 4) = 0.0;
      Bop(11, nodedof * i + 5) = 0.0;
    }
    return Bop;
  }

  /*!
   * @brief Modifiy the strain gradient (B-Operator) of the specified element due to transverse
   * shear strain ANS  (B-Bar method)
   *
   * @tparam distype : The discretization type known at compile time
   * @param Bop (in) : B-Operator matrix
   * @param shapefunctions_ans (in) : Vector holding the ANS shapefunctions and its derivatves w.r.t
   * spatial coordinates
   * @param shapefunctions_q (in) : Vector holding the shapefunctions and the first derivatives
   * w.r.t spatial coordinates evaluated at collocation points
   * @param metric_currq (in) : Metric tensor in material configuration evaluated at collocation
   * points
   * @param numans (in) : Number of ANS collocation points
   */
  template <CORE::FE::CellType distype>
  void ModifyBOperatorAns(CORE::LINALG::SerialDenseMatrix& Bop,
      const std::vector<double>& shapefunctions_ans,
      const std::vector<ShapefunctionsAndDerivatives<distype>>& shapefunctions_q,
      const std::vector<BasisVectorsAndMetrics<distype>>& metric_currq, const int& numans)
  {
    const int nodedof = DETAIL::node_dof;
    for (int i = 0; i < DETAIL::num_node<distype>; ++i)
    {
      Bop(2, nodedof * i + 0) = 0.0;
      Bop(2, nodedof * i + 1) = 0.0;
      Bop(2, nodedof * i + 2) = 0.0;
      Bop(2, nodedof * i + 3) = 0.0;
      Bop(2, nodedof * i + 4) = 0.0;
      Bop(2, nodedof * i + 5) = 0.0;
      //
      Bop(4, nodedof * i + 0) = 0.0;
      Bop(4, nodedof * i + 1) = 0.0;
      Bop(4, nodedof * i + 2) = 0.0;
      Bop(4, nodedof * i + 3) = 0.0;
      Bop(4, nodedof * i + 4) = 0.0;
      Bop(4, nodedof * i + 5) = 0.0;

      for (int j = 0; j < numans; ++j)
      {
        const double a1x1 = metric_currq[j].kovariant_(0, 0);
        const double a1y1 = metric_currq[j].kovariant_(0, 1);
        const double a1z1 = metric_currq[j].kovariant_(0, 2);

        const double a3x1 = metric_currq[j].kovariant_(2, 0);
        const double a3y1 = metric_currq[j].kovariant_(2, 1);
        const double a3z1 = metric_currq[j].kovariant_(2, 2);

        const double a2x2 = metric_currq[j + numans].kovariant_(1, 0);
        const double a2y2 = metric_currq[j + numans].kovariant_(1, 1);
        const double a2z2 = metric_currq[j + numans].kovariant_(1, 2);

        const double a3x2 = metric_currq[j + numans].kovariant_(2, 0);
        const double a3y2 = metric_currq[j + numans].kovariant_(2, 1);
        const double a3z2 = metric_currq[j + numans].kovariant_(2, 2);

        const double N1 = shapefunctions_q[j].shapefunctions_(i);
        const double N2 = shapefunctions_q[j + numans].shapefunctions_(i);

        const double dN1d1 = shapefunctions_q[j].derivatives_(0, i);
        const double dN2d2 = shapefunctions_q[j + numans].derivatives_(1, i);

        const double Nans1 = shapefunctions_ans[j];
        const double Nans2 = shapefunctions_ans[j + numans];

        // E_13 const remedy of transverse shear locking
        Bop(2, nodedof * i + 0) += dN1d1 * a3x1 * Nans1;
        Bop(2, nodedof * i + 1) += dN1d1 * a3y1 * Nans1;
        Bop(2, nodedof * i + 2) += dN1d1 * a3z1 * Nans1;
        Bop(2, nodedof * i + 3) += N1 * a1x1 * Nans1;
        Bop(2, nodedof * i + 4) += N1 * a1y1 * Nans1;
        Bop(2, nodedof * i + 5) += N1 * a1z1 * Nans1;

        // E_23 const remedy of transverse shear locking
        Bop(4, nodedof * i + 0) += dN2d2 * a3x2 * Nans2;
        Bop(4, nodedof * i + 1) += dN2d2 * a3y2 * Nans2;
        Bop(4, nodedof * i + 2) += dN2d2 * a3z2 * Nans2;
        Bop(4, nodedof * i + 3) += N2 * a2x2 * Nans2;
        Bop(4, nodedof * i + 4) += N2 * a2y2 * Nans2;
        Bop(4, nodedof * i + 5) += N2 * a2z2 * Nans2;
      }
    }
  }

  /*!
   * @brief Updates the current thickness of an element at a gaussian point
   *
   * @tparam distype : The discretization type known at compile time
   * qparam a3current (in) : Nodal directors in spatial configuraton
   * @param shape_functions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   */
  template <CORE::FE::CellType distype>
  double UpdateGaussPointThickness(
      const CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim>& a3current,
      const CORE::LINALG::Matrix<DETAIL::num_node<distype>, 1>& shapefunctions)
  {
    double current_thickness = 0;
    for (int i = 0; i < DETAIL::num_node<distype>; ++i)
    {
      double nodal_thickness = 0;
      for (int d = 0; d < DETAIL::num_dim; ++d)
        nodal_thickness += a3current(i, d) * a3current(i, d);
      nodal_thickness = 2 * std::sqrt(nodal_thickness);
      current_thickness += shapefunctions(i) * nodal_thickness;
    }
    return std::abs(current_thickness);
  }

  /*!
   * @brief Evaluates the kovariant basis vectors and metric tensors of the element
   *
   * @tparam distype :  The discretization type known at compile time
   * @param shape_functions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param basis_and_metrics_reference (in/out) : An object holding the basis vectors and metric
   * tensors of the element in reference configuration
   * @param basis_and_metrics_current (in/out) : An object holding the basis vectors and metric
   * tensors of the element in current configuration
   * @param nodal_coordinates (in) : Coordinates of the nodes of the element
   * @param zeta (in) : Thickness coordinate of gaussian point (scaled via SDC)
   */
  template <CORE::FE::CellType distype>
  void EvaluateMetrics(
      const DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& basis_and_metrics_reference,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& basis_and_metrics_current,
      const DRT::ELEMENTS::SHELL::NodalCoordinates<distype>& nodal_coordinates, const double zeta)
  {
    EvaluateKovariantVectorsAndMetrics(shape_functions, basis_and_metrics_reference,
        nodal_coordinates.x_refe_, nodal_coordinates.a3_refe_, zeta);
    EvaluateKontravariantVectorsAndMetrics(basis_and_metrics_reference);
    EvaluateKovariantVectorsAndMetrics(shape_functions, basis_and_metrics_current,
        nodal_coordinates.x_curr_, nodal_coordinates.a3_curr_, zeta);
    EvaluateKontravariantVectorsAndMetrics(basis_and_metrics_current);
  }

  /*!
   * @brief Evaluates the kovariant basis vectors and metric tensors of the element
   *
   * @tparam distype :  The discretization type known at compile time
   * @param shape_functions(in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param basis_and_metrics (in/out) : An object holding the basis vectors and metric tensors of
   * the element
   * @param x (in) : Coordinates of the nodes of the element
   * @param a3 (in) : Coordinates of the nodal directors of the element
   * @param zeta (in) : Thickness coordinate of gaussian point (scaled via SDC)
   */
  template <CORE::FE::CellType distype>
  void EvaluateKovariantVectorsAndMetrics(
      const DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& basis_and_metrics,
      const CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim>& x,
      const CORE::LINALG::Matrix<DETAIL::num_node<distype>, DETAIL::num_dim>& a3, const double zeta)
  {
    // interpolation of kovariant a1,a2
    basis_and_metrics.kovariant_.MultiplyNN(1.0, shape_functions.derivatives_, x, 0.0);

    // displacement field in the third dimension is approximated by assuming a linear variation of
    // the displacements across the thickness
    if (zeta)
    {
      CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> tmp(true);
      tmp.MultiplyNN(zeta, shape_functions.derivatives_, a3, 0.0);
      basis_and_metrics.kovariant_.Update(1.0, tmp, 1.0);
    }

    // interpolation of a3
    for (int idim = 0; idim < DETAIL::num_dim; idim++)
    {
      basis_and_metrics.kovariant_(2, idim) = 0.0;
      for (int inode = 0; inode < DETAIL::num_node<distype>; inode++)
        basis_and_metrics.kovariant_(2, idim) +=
            shape_functions.shapefunctions_(inode) * a3(inode, idim);
    }
    // get kovariant metric tensor
    basis_and_metrics.metric_kovariant_.MultiplyNT(
        1.0, basis_and_metrics.kovariant_, basis_and_metrics.kovariant_, 0.0);

    // get partial derivatives of the basis vector in thickness direction
    basis_and_metrics.partial_derivative_.MultiplyNN(1.0, shape_functions.derivatives_, a3, 0.0);
  }

  /*!
   * @brief Evaluates the kontravariant basis vectors and metric tensors of the element
   *
   * @tparam distype : The discretization type known at compile time
   *
   * @param metrics (in/out) : An object holding the basis vectors and metric tensors of the element
   */
  template <CORE::FE::CellType distype>
  void EvaluateKontravariantVectorsAndMetrics(
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& metrics)
  {
    // get kontravariant basis vectors g1,g2,g3 (inverse transpose of kov)
    metrics.kontravariant_ = metrics.kovariant_;
    metrics.detJ_ = metrics.kontravariant_.Invert();

    // get the transpose
    for (int i = 0; i < DETAIL::num_dim; ++i)
    {
      for (int j = i + 1; j < DETAIL::num_dim; ++j)
      {
        const double tmp = metrics.kontravariant_(j, i);
        metrics.kontravariant_(j, i) = metrics.kontravariant_(i, j);
        metrics.kontravariant_(i, j) = tmp;
      }
    }
    metrics.metric_kontravariant_ = metrics.metric_kovariant_;
    // get  kontravariant metrictensor
    metrics.metric_kontravariant_.Invert();
  }

  /*!
   * @brief Modify the kovariant metric tensors of the element to neglect the quadratic terms of the
   * normal strain in thickness direction.
   *
   * @Note : After the modification for the current metric tensor holds g_ij != g_i*g_j
   *
   * @tparam distype :  The discretization type known at compile time
   * @param g_reference (in/out) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in/out) : An object holding the current basis vectors and metric
   * tensors of the shell body
   *  @param a_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the midsurface
   * @param a_current (in) : An object holding the current basis vectors and metric tensors
   * tensors of the midsurface
   * @param zeta (in) : Thickness coordinate of gaussian point (scaled via SDC)
   */
  template <CORE::FE::CellType distype>
  void ModifyKovariantMetricsStandart(
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_reference,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_current,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_current, const double& zeta)
  {
    double b11c = 0.0;
    double b12c = 0.0;
    double b21c = 0.0;
    double b22c = 0.0;
    double b31c = 0.0;
    double b32c = 0.0;
    double b11r = 0.0;
    double b12r = 0.0;
    double b21r = 0.0;
    double b22r = 0.0;
    double b31r = 0.0;
    double b32r = 0.0;

    for (int i = 0; i < DETAIL::num_dim; ++i)
    {
      b11c += a_current.kovariant_(0, i) * a_current.partial_derivative_(0, i);
      b12c += a_current.kovariant_(0, i) * a_current.partial_derivative_(1, i);
      b21c += a_current.kovariant_(1, i) * a_current.partial_derivative_(0, i);
      b22c += a_current.kovariant_(1, i) * a_current.partial_derivative_(1, i);
      b31c += a_current.kovariant_(2, i) * a_current.partial_derivative_(0, i);
      b32c += a_current.kovariant_(2, i) * a_current.partial_derivative_(1, i);

      b11r += a_reference.kovariant_(0, i) * a_reference.partial_derivative_(0, i);
      b12r += a_reference.kovariant_(0, i) * a_reference.partial_derivative_(1, i);
      b21r += a_reference.kovariant_(1, i) * a_reference.partial_derivative_(0, i);
      b22r += a_reference.kovariant_(1, i) * a_reference.partial_derivative_(1, i);
      b31r += a_reference.kovariant_(2, i) * a_reference.partial_derivative_(0, i);
      b32r += a_reference.kovariant_(2, i) * a_reference.partial_derivative_(1, i);
    }
    g_current.metric_kovariant_(0, 0) =
        g_reference.metric_kovariant_(0, 0) +
        (a_current.metric_kovariant_(0, 0) - a_reference.metric_kovariant_(0, 0)) +
        zeta * 2.0 * (b11c - b11r);
    g_current.metric_kovariant_(1, 1) =
        g_reference.metric_kovariant_(1, 1) +
        (a_current.metric_kovariant_(1, 1) - a_reference.metric_kovariant_(1, 1)) +
        zeta * 2.0 * (b22c - b22r);
    g_current.metric_kovariant_(2, 2) =
        g_reference.metric_kovariant_(2, 2) +
        (a_current.metric_kovariant_(2, 2) - a_reference.metric_kovariant_(2, 2));

    g_current.metric_kovariant_(0, 1) =
        g_reference.metric_kovariant_(0, 1) +
        (a_current.metric_kovariant_(0, 1) - a_reference.metric_kovariant_(0, 1)) +
        zeta * (b21c + b12c - b21r - b12r);
    g_current.metric_kovariant_(0, 2) =
        g_reference.metric_kovariant_(0, 2) +
        (a_current.metric_kovariant_(0, 2) - a_reference.metric_kovariant_(0, 2)) +
        zeta * (b31c - b31r);
    g_current.metric_kovariant_(1, 2) =
        g_reference.metric_kovariant_(1, 2) +
        (a_current.metric_kovariant_(1, 2) - a_reference.metric_kovariant_(1, 2)) +
        zeta * (b32c - b32r);

    g_current.metric_kovariant_(2, 0) = g_current.metric_kovariant_(0, 2);
    g_current.metric_kovariant_(2, 1) = g_current.metric_kovariant_(1, 2);
    g_current.metric_kovariant_(1, 0) = g_current.metric_kovariant_(0, 1);

    g_current.metric_kontravariant_.Update(g_current.metric_kovariant_);
    double detJ = g_current.metric_kontravariant_.Invert();
    g_current.detJ_ = std::sqrt(detJ);
    g_reference.metric_kontravariant_.Update(g_reference.metric_kovariant_);
    detJ = g_reference.metric_kontravariant_.Invert();
    g_reference.detJ_ = std::sqrt(detJ);
  }

  /*!
   * @brief Modify the kovariant metric tensors of the element due to transverse shear strain ANS
   *
   * @Note : After the modification: g_ij != g_i*g_j
   *
   * @tparam distype :  The discretization type known at compile time
   * @param g_reference (in/out) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in/out) : An object holding the current basis vectors and metric
   * tensors of the shell body
   * @param a_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the midsurface
   * @param a_current (in) : An object holding the current basis vectors and metric tensors
   * tensors of the midsurface
   * @param zeta (in) : thickness coordinate of gaussian point (scaled via SDC)
   * @param shapefunctions_ans (in) : Vector holding the ANS shapefunctions and derivatives
   * @param metrics_collocation_reference (in) : vector holding the reference basis vectors and
   * metric tensors of the midsurface at collocation points
   * @param metrics_collocation_current (in) : Vector holding the current basis vectors and metric
   * tensors of the midsurface at collocation points
   * @param numansq (in) : Number of ANS collocation points
   */
  template <CORE::FE::CellType distype>
  void ModifyKovariantMetricsAns(DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_reference,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_current,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_current, const double& zeta,
      const std::vector<double>& shapefunctions_ans,
      const std::vector<DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>>&
          metrics_collocation_reference,
      const std::vector<DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>>&
          metrics_collocation_current,
      const int& numansq)
  {
    double b11c = 0.0;
    double b12c = 0.0;
    double b21c = 0.0;
    double b22c = 0.0;
    double b31c = 0.0;
    double b32c = 0.0;
    double b11r = 0.0;
    double b12r = 0.0;
    double b21r = 0.0;
    double b22r = 0.0;
    double b31r = 0.0;
    double b32r = 0.0;

    for (int i = 0; i < DETAIL::num_dim; ++i)
    {
      b11c += a_current.kovariant_(0, i) * a_current.partial_derivative_(0, i);
      b12c += a_current.kovariant_(0, i) * a_current.partial_derivative_(1, i);
      b21c += a_current.kovariant_(1, i) * a_current.partial_derivative_(0, i);
      b22c += a_current.kovariant_(1, i) * a_current.partial_derivative_(1, i);
      b31c += a_current.kovariant_(2, i) * a_current.partial_derivative_(0, i);
      b32c += a_current.kovariant_(2, i) * a_current.partial_derivative_(1, i);

      b11r += a_reference.kovariant_(0, i) * a_reference.partial_derivative_(0, i);
      b12r += a_reference.kovariant_(0, i) * a_reference.partial_derivative_(1, i);
      b21r += a_reference.kovariant_(1, i) * a_reference.partial_derivative_(0, i);
      b22r += a_reference.kovariant_(1, i) * a_reference.partial_derivative_(1, i);
      b31r += a_reference.kovariant_(2, i) * a_reference.partial_derivative_(0, i);
      b32r += a_reference.kovariant_(2, i) * a_reference.partial_derivative_(1, i);
    }
    g_current.metric_kovariant_(0, 0) =
        g_reference.metric_kovariant_(0, 0) +
        (a_current.metric_kovariant_(0, 0) - a_reference.metric_kovariant_(0, 0)) +
        zeta * 2.0 * (b11c - b11r);
    g_current.metric_kovariant_(1, 1) =
        g_reference.metric_kovariant_(1, 1) +
        (a_current.metric_kovariant_(1, 1) - a_reference.metric_kovariant_(1, 1)) +
        zeta * 2.0 * (b22c - b22r);
    g_current.metric_kovariant_(2, 2) =
        g_reference.metric_kovariant_(2, 2) +
        (a_current.metric_kovariant_(2, 2) - a_reference.metric_kovariant_(2, 2));

    g_current.metric_kovariant_(0, 1) =
        g_reference.metric_kovariant_(0, 1) +
        (a_current.metric_kovariant_(0, 1) - a_reference.metric_kovariant_(0, 1)) +
        zeta * (b21c + b12c - b21r - b12r);
    g_current.metric_kovariant_(0, 2) = g_reference.metric_kovariant_(0, 2) + zeta * (b31c - b31r);
    g_current.metric_kovariant_(1, 2) = g_reference.metric_kovariant_(1, 2) + zeta * (b32c - b32r);
    for (int qp = 0; qp < numansq; ++qp)
    {
      g_current.metric_kovariant_(0, 2) +=
          (metrics_collocation_current[qp].metric_kovariant_(0, 2) -
              metrics_collocation_reference[qp].metric_kovariant_(0, 2)) *
          shapefunctions_ans[qp];
      g_current.metric_kovariant_(1, 2) +=
          (metrics_collocation_current[qp + numansq].metric_kovariant_(1, 2) -
              metrics_collocation_reference[qp + numansq].metric_kovariant_(1, 2)) *
          shapefunctions_ans[qp + numansq];
    }
    g_current.metric_kovariant_(2, 0) = g_current.metric_kovariant_(0, 2);
    g_current.metric_kovariant_(2, 1) = g_current.metric_kovariant_(1, 2);
    g_current.metric_kovariant_(1, 0) = g_current.metric_kovariant_(0, 1);

    g_current.metric_kontravariant_.Update(g_current.metric_kovariant_);
    double detJ = g_current.metric_kontravariant_.Invert();
    g_current.detJ_ = std::sqrt(detJ);
    g_reference.metric_kontravariant_.Update(g_reference.metric_kovariant_);
    detJ = g_reference.metric_kontravariant_.Invert();
    g_reference.detJ_ = std::sqrt(detJ);
  }


  /*!
   * @brief Modify the kovariant metric tensors of the element
   *
   *
   * @tparam distype :  The discretization type known at compile time
   * @param g_reference (in/out) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in/out) : An object holding the current basis vectors and metric
   * tensors of the shell body
   * @param a_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the midsurface
   * @param a_current (in) : An object holding the current basis vectors and metric tensors
   * tensors of the midsurface
   * @param zeta (in) : thickness coordinate of gaussian point (scaled via SDC)
   * @param shapefunctions_ans (in) : Vector holding the ANS shapefunctions and derivatives
   * @param metrics_collocation_reference (in) : vector holding the reference basis vectors and
   * metric tensors of the midsurface at collocation points
   * @param metrics_collocation_current (in) : Vector holding the current basis vectors and metric
   * tensors of the midsurface at collocation points
   * @param numansq (in) : Number of ANS collocation points
   */
  template <CORE::FE::CellType distype>
  void ModifyKovariantMetrics(DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_reference,
      DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_current,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_current, const double& zeta,
      const std::vector<double>& shapefunctions_ans,
      const std::vector<DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>>&
          metrics_collocation_reference,
      const std::vector<DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>>&
          metrics_collocation_current,
      const int& numansq)
  {
    std::invoke(
        [&]()
        {
          if (numansq == 0)
          {
            ModifyKovariantMetricsStandart(g_reference, g_current, a_reference, a_current, zeta);
          }
          else
          {
            // modify the current kovariant metric tensor due to transverse shear strain ANS
            ModifyKovariantMetricsAns(g_reference, g_current, a_reference, a_current, zeta,
                shapefunctions_ans, metrics_collocation_reference, metrics_collocation_current,
                numansq);
          }
        });
  }

  /*!
   * @brief Evaluates displacement based Deformation gradient tensor from metric tensors
   *
   * The deformamtion gradient tensor is defined as :
   * \f[
   *    \mathbf{F} = \mathbf{g}_i \otimes \mathbf{G}^i
   * \f]
   * This can be transformed to a cartesian system by \f$ F_{ij}=\mathbf{E}_i\cdot \mathbf{F}
   * \mathbf{E}_j$\f. As we used a global cartesian basis with the basis vectors  \f$ E_1=[1, 0,
   * 0]^T, E_2=[0, 1, 0]^T, E_3=[0, 0, 1]^T$\f, no transformation is necessary
   *
   * @tparam distype : The discretization type known at compile time
   * @params strains (in/out) :  Strain measures of the element (deformation gradient,
   * Green-Lagrange strain tensor)
   * @param g_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in) : An object holding the current basis vectors and metric
   * tensors of the shell body
   */
  template <CORE::FE::CellType distype>
  void EvaluateDeformationGradient(DRT::ELEMENTS::SHELL::Strains& strains,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_current)
  {
    strains.defgrd_.MultiplyNT(1.0, g_current.kovariant_, g_reference.kontravariant_);
  }

  /*!
   * @brief Evaluates Green-Lagrange strains from metrics
   *
   *  The reen-Lagrange strains are defined as :
   * \f[
   *    \mathbf{E} = E_{ij} \mathbf{G}^i \otimes \mathbf{G}^j
   * \f]
   * GL strain vector glstrain = [E_{11},E_{22}, E_{33}, 2*E_{12}, 2*E_{23}, 2*E_{31}]^T
   * @tparam distype : The discretization type known at compile time
   * @param strains (in/out) : Strain measures (deformation graident tensor, Green-Lagrange strain
   * tensor)
   * @param g_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in) : An object holding the current basis vectors and metric
   * tensors of the shell body
   */
  template <CORE::FE::CellType distype>
  void EvaluateGreenLagrangeStrain(Strains& strains,
      const BasisVectorsAndMetrics<distype>& g_reference,
      const BasisVectorsAndMetrics<distype>& g_current)
  {
    //  evaluate strain tensor in curvilinear coordinate system E_ij = 0.5 (g_ij-G_ij)
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> gl_strain_tensor(true);
    for (int i = 0; i < DETAIL::num_dim; ++i)
    {
      for (int j = 0; j < DETAIL::num_dim; ++j)
        gl_strain_tensor(i, j) =
            0.5 * (g_current.metric_kovariant_(i, j) - g_reference.metric_kovariant_(i, j));
    }
    // map gl strains from curvilinear system to global cartesian system
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> gl_strain_tensor_cartesian(true);
    CORE::LINALG::TENSOR::InverseTensorRotation<DETAIL::num_dim>(
        g_reference.kontravariant_, gl_strain_tensor, gl_strain_tensor_cartesian);
    // GL strain vector glstrain for solid material E={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::VOIGT::Strains::MatrixToVector(gl_strain_tensor_cartesian, strains.gl_strain_);
  }

  /*!
   * @brief Maps the stress quantities (2. Piola-Kirchoff stress tensor and the linearization
   * w.r.t Green-Lagrange strain) from the global cartesioan system
   * back to the curvilinear system and re-arranged the order for the specific shell formulation
   *
   * @Note: PK Stress are ordered here like
   * \f[
   *    \mathbf{S} = [S_{11},S_{12}, S_{13}, S_{22}, S_{23}, S_{33}]^T
   * \f]
   *
   * @tparam distype :  The discretization type known at compile time
   * @params stress (in/out) : Stress quantities of the element (2. Piola-Kirchoff stress tensor
   * and the linearization w.r.t Green-Lagrange strain)
   * @param g_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   */
  template <CORE::FE::CellType distype>
  void MapMaterialStressToCurvilinearSystem(
      Stress<MAT::NUM_STRESS_3D>& stress, const BasisVectorsAndMetrics<distype>& g_reference)
  {
    // transform Piola-Kirchhoff-stresses from global cartesian coordinate system back to local
    // curvilinear coordinate system
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> stress_tensor(true);
    CORE::LINALG::VOIGT::Stresses::VectorToMatrix(stress.pk2_, stress_tensor);
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> tmp(true);
    CORE::LINALG::TENSOR::TensorRotation(g_reference.kontravariant_, stress_tensor, tmp);

    // re-arrange indices for shell element formulation:
    // PK Stress=[S_{11},S_{12}, S_{13}, S_{22}, S_{23}, S_{33}]^T
    stress.pk2_(0) = tmp(0, 0);
    stress.pk2_(1) = tmp(0, 1);
    stress.pk2_(2) = tmp(0, 2);
    stress.pk2_(3) = tmp(1, 1);
    stress.pk2_(4) = tmp(1, 2);
    stress.pk2_(5) = tmp(2, 2);

    // transform elasticity matrix from global cartesian coordinate system back to local
    // curvilinear coordinate system
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> Cmat(true);
    CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> g_metrics_trans(true);
    g_metrics_trans.UpdateT(g_reference.kontravariant_);
    CORE::LINALG::TENSOR::InverseFourthTensorRotation(g_metrics_trans, stress.cmat_, Cmat);

    // re-arrange indices for shell element formulation
    static constexpr std::array voigt_inconsistent_ = {0, 3, 5, 1, 4, 2};
    for (unsigned int i = 0; i < 6; ++i)
      for (unsigned int j = 0; j < 6; ++j)
        stress.cmat_(i, j) = Cmat(voigt_inconsistent_[i], voigt_inconsistent_[j]);
  }

  /*!
   * @brief Get the shapefunctions evaluated at the gaussian points to remedy transverse shear
   * strains within the ANS method
   *
   * @tparam distype : The discretization type known at compile time
   * @params xi_gp (in) : Coordinate of the integration point in the parameter space
   * @return shapefunctions for ANS
   */
  template <CORE::FE::CellType distype>
  static std::vector<double> SetShapefunctionsForAns(const std::array<double, 2>& xi_gp)
  {
    std::vector<double> shapefunctions;
    double r = xi_gp[0];
    double s = xi_gp[1];

    if (distype == CORE::FE::CellType::quad4)
    {
      shapefunctions.resize(4);
      shapefunctions[0] = 0.5 * (1.0 - s);
      shapefunctions[1] = 0.5 * (1.0 + s);
      shapefunctions[2] = 0.5 * (1.0 - r);
      shapefunctions[3] = 0.5 * (1.0 + r);
    }
    else if (distype == CORE::FE::CellType::quad9)
    {
      shapefunctions.resize(12);
      const double rthreei = 1.0 / (sqrt(3.0));
      std::array<double, 3> pr;
      std::array<double, 3> ps;
      std::array<double, 2> qr;
      std::array<double, 2> qs;

      pr[0] = -0.5 * s * (1.0 - s);
      pr[1] = (1.0 - s) * (1.0 + s);
      pr[2] = 0.5 * s * (1.0 + s);

      qr[0] = 0.5 * (1.0 - r / rthreei);
      qr[1] = 0.5 * (1.0 + r / rthreei);

      ps[0] = -0.5 * r * (1.0 - r);
      ps[1] = (1.0 - r) * (1.0 + r);
      ps[2] = 0.5 * r * (1.0 + r);

      qs[0] = 0.5 * (1.0 - s / rthreei);
      qs[1] = 0.5 * (1.0 + s / rthreei);

      shapefunctions[0] = pr[0] * qr[0];
      shapefunctions[1] = pr[1] * qr[0];
      shapefunctions[2] = pr[2] * qr[0];
      shapefunctions[3] = pr[0] * qr[1];
      shapefunctions[4] = pr[1] * qr[1];
      shapefunctions[5] = pr[2] * qr[1];

      shapefunctions[6] = ps[0] * qs[0];
      shapefunctions[7] = ps[1] * qs[0];
      shapefunctions[8] = ps[2] * qs[0];
      shapefunctions[9] = ps[0] * qs[1];
      shapefunctions[10] = ps[1] * qs[1];
      shapefunctions[11] = ps[2] * qs[1];
    }
    return shapefunctions;
  }

  /*!
   * @brief Get the shapefunctions evaluated at the gaussian points to remedy transverse shear
   * strains within the ANS method
   *
   * @tparam distype : The discretization type known at compile time
   * @params xi_gp (in) : Coordinate of the integration point in the parameter space
   * @params num_ans (in) : Number of ANS collocation points
   * @return shapefunctions for ANS
   */
  template <CORE::FE::CellType distype>
  std::vector<double> GetShapefunctionsForAns(
      const std::array<double, 2>& xi_gp, const int& num_ans)
  {
    const std::vector<double> shape_functions_ans = std::invoke(
        [&]()
        {
          if (num_ans == 0)
          {
            return std::vector<double>{};
          }
          else
          {
            return SHELL::SetShapefunctionsForAns<distype>(xi_gp);
          }
        });
    return shape_functions_ans;
  }
  /*!
   * @brief Transform the Green-Lagrange strains to Euler-Almansi strains
   *
   * @tparam distype : The discretization type known at compile time
   * @params gl (in) :  Green-Lagrange strains
   * @param defgrd (in) : Deformations gradient tensor
   * @param ea (in/out) : Euler-Almansi strains
   */
  template <CORE::FE::CellType distype>
  void GreenLagrangeToEulerAlmansi(const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& gl,
      const CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim>& defgrd,
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& ea)
  {
    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> invdefgrd(defgrd);
    invdefgrd.Invert();

    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> E_matrix;
    CORE::LINALG::VOIGT::Strains::VectorToMatrix(gl, E_matrix);

    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> invFTE;
    invFTE.MultiplyTN(invdefgrd, E_matrix);

    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> ea_matrix;
    ea_matrix.MultiplyNN(invFTE, invdefgrd);

    CORE::LINALG::VOIGT::Strains::MatrixToVector(ea_matrix, ea);
  }

  /*!
   * @brief Transform the 2. Piola-Kirchoff stresses to Cauchy stresses
   *
   * @tparam distype : The discretization type known at compile time
   * @params pk2 (in) :  2. Piola-Kirchoff stresses
   * @param defgrd (in) : Deformations gradient tensor
   * @param cauchy (in/out) : Cauchy stresses
   */
  template <CORE::FE::CellType distype>
  void Pk2ToCauchy(const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& pk2,
      const CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim>& defgrd,
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& cauchy)
  {
    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> S_matrix;
    CORE::LINALG::VOIGT::Stresses::VectorToMatrix(pk2, S_matrix);

    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> FS;
    FS.MultiplyNN(defgrd, S_matrix);

    CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> cauchy_matrix;
    cauchy_matrix.MultiplyNT(1.0 / defgrd.Determinant(), FS, defgrd, 0.0);

    CORE::LINALG::VOIGT::Stresses::MatrixToVector(cauchy_matrix, cauchy);
  }

  /*!
   * @brief Assembles data in voigt vector notation to matrix
   *
   * @tparam num_internal_variables : number of internal parameters
   * @params vector (in) :  Vector in voigt notation
   * @param data (in/out) : Data to which the vector will be assembled
   * @param thickness_weight (in) : Weighting factor to consider thickness integration
   */
  template <unsigned num_internal_variables>
  void AssembleVectorToMatrixRow(CORE::LINALG::Matrix<num_internal_variables, 1> vector,
      CORE::LINALG::SerialDenseMatrix& data, int row, double thickness_weight)
  {
    for (unsigned int i = 0; i < num_internal_variables; ++i)
      data(row, i) += thickness_weight * vector(i);
  }

  /*!
   * @brief Assembles strain data in voigt vector notation to matrix
   *
   * @tparam distype : The discretization type known at compile time
   * @params strains (in) :  Strain vector
   * @params strains_type (in) :  Strain type
   * @param data (in/out) : Data to which the vector will be assembled
   * @param row (in) : Row number
   * @param thickness_weight (in) : Weighting factor to consider thickness integration
   */
  template <CORE::FE::CellType distype>
  void AssembleStrainTypeToMatrixRow(const Strains& strains, INPAR::STR::StrainType strain_type,
      CORE::LINALG::SerialDenseMatrix& data, int row, const double thickness_weight)
  {
    switch (strain_type)
    {
      case INPAR::STR::strain_gl:
      {
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> gl_strain_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(strains.gl_strain_, gl_strain_stress_like);
        AssembleVectorToMatrixRow(gl_strain_stress_like, data, row, thickness_weight);
        return;
      }
      case INPAR::STR::strain_ea:
      {
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> ea;
        GreenLagrangeToEulerAlmansi<distype>(strains.gl_strain_, strains.defgrd_, ea);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> ea_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(ea, ea_stress_like);
        AssembleVectorToMatrixRow(ea_stress_like, data, row, thickness_weight);
        return;
      }
      case INPAR::STR::strain_none:
        return;
      default:
        FOUR_C_THROW("strain type not supported");
    }
  }
  /*!
   * @brief Assembles stress data in voigt vector notation to matrix
   *
   * @tparam distype : The discretization type known at compile time
   * @params stress (in) :  Stress vector
   * @params stress_type (in) :  Stress type
   * @param data (in/out) : Data to which the vector will be assembled
   * @param row (in) : Row number
   * @param thickness_weight (in) : Weighting factor to consider thickness integration
   */
  template <CORE::FE::CellType distype>
  void AssembleStressTypeToMatrixRow(const Strains& strains,
      const Stress<MAT::NUM_STRESS_3D> stress, INPAR::STR::StressType stress_type,
      CORE::LINALG::SerialDenseMatrix& data, int row, const double thickness_weight)
  {
    switch (stress_type)
    {
      case INPAR::STR::stress_2pk:
      {
        AssembleVectorToMatrixRow(stress.pk2_, data, row, thickness_weight);
        return;
      }
      case INPAR::STR::stress_cauchy:
      {
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cauchy;
        Pk2ToCauchy<distype>(stress.pk2_, strains.defgrd_, cauchy);
        AssembleVectorToMatrixRow(cauchy, data, row, thickness_weight);
        return;
      }
      case INPAR::STR::stress_none:
        return;
      default:
        FOUR_C_THROW("stress type not supported");
    }
  }

  inline void Serialize(const CORE::LINALG::SerialDenseMatrix& matrix, std::vector<char>& data)
  {
    CORE::COMM::PackBuffer packBuffer;
    CORE::COMM::ParObject::AddtoPack(packBuffer, matrix);
    packBuffer.StartPacking();
    CORE::COMM::ParObject::AddtoPack(packBuffer, matrix);
    std::copy(packBuffer().begin(), packBuffer().end(), std::back_inserter(data));
  }


  /*!
   * @brief Returns the local coordinates of the collocation points for ANS.
   *
   * This functions returns the local coordinates of the collocation points for ANS to avoid
   * transverse shear locking. The choice of coordinates depends on the discretization types of the
   * element. Here, they are chosen such that the points lie on the middle of each edge.
   *
   * @tparam distype : Discretization type
   * @params qp (in) :  Index of collocation point
   */
  template <CORE::FE::CellType distype>
  void GetCoordinatesOfAnsCollocationPoints(std::vector<std::array<double, 2>>& collocation_points)
  {
    if (distype == CORE::FE::CellType::quad4)
    {
      collocation_points.resize(4);
      // for a_13
      collocation_points[0] = {0.0, -1.0};
      collocation_points[1] = {0.0, 1.0};
      // for a_23
      collocation_points[2] = {-1.0, 0.0};
      collocation_points[3] = {1.0, 0.0};
    }
    else if (distype == CORE::FE::CellType::quad9)
    {
      const double sqrt3inv = 1.0 / (sqrt(3.0));
      collocation_points.resize(DETAIL::num_internal_variables);
      // for a_13
      collocation_points[0] = {-sqrt3inv, -1.0};
      collocation_points[1] = {-sqrt3inv, 0.0};
      collocation_points[2] = {-sqrt3inv, 1.0};
      collocation_points[3] = {sqrt3inv, -1.0};
      collocation_points[4] = {sqrt3inv, 0.0};
      collocation_points[5] = {sqrt3inv, 1.0};
      // for a_23
      collocation_points[6] = {-1.0, -sqrt3inv};
      collocation_points[7] = {0.0, -sqrt3inv};
      collocation_points[8] = {1.0, -sqrt3inv};
      collocation_points[9] = {-1.0, sqrt3inv};
      collocation_points[10] = {0.0, sqrt3inv};
      collocation_points[11] = {1.0, sqrt3inv};
    }
  }

  /*!
   * @brief Setup shapefunction, the first derivatives w.r.t spatial coordinates and the metrics
   * evaluated at collocation points for ANS
   *
   * @tparam distype : The discretization type known at compile time
   * @params shapefunctions_collocation (in/out) :  An object holding the shape functions and the
   * first derivatives w.r.t spatial coordinates
   * @params metrics_collocation_reference (in/out) :  An object holding the reference basis vectors
   * and metric tensors of the shell body
   * @param metrics_collocation_current (in/out) : An object holding the current basis vectors and
   * metric tensors of the shell body
   * @param nodal_coordinates (in) : Nodal coordinates in reference and current frame
   * @param num_collocation_points (in) : Number of collocation points
   */
  template <CORE::FE::CellType distype>
  void SetupANS(
      std::vector<SHELL::ShapefunctionsAndDerivatives<distype>>& shapefunctions_collocation,
      std::vector<SHELL::BasisVectorsAndMetrics<distype>>& metrics_collocation_reference,
      std::vector<SHELL::BasisVectorsAndMetrics<distype>>& metrics_collocation_current,
      const NodalCoordinates<distype> nodal_coordinates, const int& num_collocation_points)
  {
    // get coordinates of collocations points
    std::vector<std::array<double, 2>> collocation_points;
    SHELL::GetCoordinatesOfAnsCollocationPoints<distype>(collocation_points);
    for (int qp = 0; qp < num_collocation_points; ++qp)
    {
      // get shape functions and derivatives for ANS at each collocation point
      shapefunctions_collocation[qp] =
          SHELL::EvaluateShapefunctionsAndDerivs<distype>(collocation_points[qp]);
      //  get basis vectors and metric at collocation points
      SHELL::EvaluateMetrics(shapefunctions_collocation[qp], metrics_collocation_reference[qp],
          metrics_collocation_current[qp], nodal_coordinates, 0.0);
    }
  }

  /*!
   * @brief This function perfoms the preintegration of the the 2. Piola-Kirchhoff stresses and
   * the material tensor through the shell thickness.
   *
   * @tparam distype : The discretization type known at compile time
   * @param stress_enh (in/out) : An object holding the enhanced stress resultants
   * @param stress (in) :  An object holding the stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param zeta (in) : Thickness coordinate of gaussian point (scaled via SDC)
   * to
   */
  template <CORE::FE::CellType distype>
  void ThicknessIntegration(DRT::ELEMENTS::SHELL::StressEnhanced& stress_enh,
      const DRT::ELEMENTS::SHELL::Stress<MAT::NUM_STRESS_3D>& stress,
      const double& integration_factor, const double& zeta)
  {
    for (int i = 0; i < DETAIL::node_dof; ++i)
    {
      const double stress_fact = stress.pk2_(i) * integration_factor;
      stress_enh.stress_(i) += stress_fact;
      stress_enh.stress_(i + DETAIL::node_dof) += stress_fact * zeta;

      for (int j = 0; j < DETAIL::node_dof; ++j)
      {
        const double C_fact = stress.cmat_(i, j) * integration_factor;
        stress_enh.dmat_(i, j) += C_fact;
        stress_enh.dmat_(i + DETAIL::node_dof, j) += C_fact * zeta;
        stress_enh.dmat_(i + DETAIL::node_dof, j + DETAIL::node_dof) += C_fact * zeta * zeta;
      }
    }
    // symmetrize enhanced material tensor
    for (int i = 0; i < DETAIL::num_internal_variables; i++)
      for (int j = i + 1; j < DETAIL::num_internal_variables; j++)
        stress_enh.dmat_(i, j) = stress_enh.dmat_(j, i);
  }

  /*!
   * @brief Evaluates strain measures
   *
   * @tparam distype : The discretization type known at compile time
   * @param g_reference (in) : An object holding the reference basis vectors and metric
   * tensors of the shell body
   * @param g_current (in) : An object holding the current basis vectors and metric
   * tensors of the shell body
   * @return Strains<distype> : Strain measures of the element (deformation gradient, Green-Lagrange
   * strain tensor)
   */
  template <CORE::FE::CellType distype>
  DRT::ELEMENTS::SHELL::Strains EvaluateStrains(
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_metric_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g_metric_current)
  {
    Strains strains;
    EvaluateDeformationGradient(strains, g_metric_reference, g_metric_current);
    EvaluateGreenLagrangeStrain(strains, g_metric_reference, g_metric_current);
    return strains;
  }

  /*!
   * @brief Adds elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : The discretization type known at compile time
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param D (in) : Material tensor (equals integration of the material tensor
   * through the shell thickness)
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian times scaled director factor)
   * @param stiffness_matrix (in/out) : Stiffness matrix where the local contribution is added
   * to
   */
  template <CORE::FE::CellType distype>
  void AddElasticStiffnessMatrix(const CORE::LINALG::SerialDenseMatrix& Bop,
      const CORE::LINALG::SerialDenseMatrix& Dmat, const double& integration_fac,
      CORE::LINALG::SerialDenseMatrix& stiffness_matrix)
  {
    // calculate Ke = integration_factor * B^TDB
    CORE::LINALG::SerialDenseMatrix DB(
        DETAIL::num_internal_variables, DETAIL::numdofperelement<distype>);
    DB.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Dmat, Bop, 0.0);
    stiffness_matrix.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, integration_fac, Bop, DB, 1.0);
  }

  /*!
   * @brief Adds geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : The discretization type known at compile time
   *  @param shapefunctions_collocation (in) : Shapefunctions and the first derivatives w.r.t
   * spatial coordinates evaluated at collocation points
   * @param shapefunctions_ans (in) : Shapefunctions for transverse shear strain ANS
   * @param shapefunctions (in) : Shapefunctions and the first derivatives w.r.t spatial coordinates
   * (midsurface)
   * @param stress_enh (in) :  An object holding the enhanced stress resultants
   * @param numansq (in) : Number of ANS collocation points
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param stiffness_matrix (in/out) : Stiffness matrix where the local contribution is added
   * to
   */
  template <CORE::FE::CellType distype>
  void AddGeometricStiffnessMatrix(
      const std::vector<DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>>&
          shapefunctions_q,
      const std::vector<double>& shapefunctions_ans,
      const DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>& shapefunctions,
      const CORE::LINALG::SerialDenseVector& stress_enh, const int& numans,
      const double& integration_fac, CORE::LINALG::SerialDenseMatrix& stiffness_matrix)
  {
    const int nodedof = DETAIL::node_dof;
    const int num_dim = DETAIL::num_dim;
    for (int inod = 0; inod < DETAIL::num_node<distype>; ++inod)
    {
      for (int jnod = 0; jnod <= inod; ++jnod)
      {
        const double Ni = shapefunctions.shapefunctions_(inod);
        const double Nj = shapefunctions.shapefunctions_(jnod);
        const double dN11 =
            shapefunctions.derivatives_(0, inod) * shapefunctions.derivatives_(0, jnod);
        const double dN12 =
            shapefunctions.derivatives_(0, inod) * shapefunctions.derivatives_(1, jnod);
        const double dN21 =
            shapefunctions.derivatives_(1, inod) * shapefunctions.derivatives_(0, jnod);
        const double dN22 =
            shapefunctions.derivatives_(1, inod) * shapefunctions.derivatives_(1, jnod);

        const double tmp1 =
            (dN11 * stress_enh(0) + (dN12 + dN21) * stress_enh(1) + dN22 * stress_enh(3)) *
            integration_fac;
        const double tmp2 =
            (dN11 * stress_enh(6) + (dN12 + dN21) * stress_enh(7) + dN22 * stress_enh(9)) *
            integration_fac;

        double tmp3 = 0.0;
        double tmp4 = 0.0;
        const double dN1dij = shapefunctions.derivatives_(0, inod) * Nj;
        const double dN1dji = shapefunctions.derivatives_(0, jnod) * Ni;
        const double dN2dij = shapefunctions.derivatives_(1, inod) * Nj;
        const double dN2dji = shapefunctions.derivatives_(1, jnod) * Ni;

        if (!numans)
        {
          tmp3 = (dN1dji * stress_enh(2) + dN2dji * stress_enh(4)) * integration_fac;
          tmp4 = (dN1dij * stress_enh(2) + dN2dij * stress_enh(4)) * integration_fac;
        }
        // modification due to transverse shear strain ANS
        else
        {
          for (int q = 0; q < numans; ++q)
          {
            double dNq1dij = shapefunctions_q[q].derivatives_(0, inod) *
                             shapefunctions_q[q].shapefunctions_(jnod) * shapefunctions_ans[q];
            double dNq1dji = shapefunctions_q[q].derivatives_(0, jnod) *
                             shapefunctions_q[q].shapefunctions_(inod) * shapefunctions_ans[q];
            int index_q23 = q + numans;
            double dNq2dij = shapefunctions_q[index_q23].derivatives_(1, inod) *
                             shapefunctions_q[index_q23].shapefunctions_(jnod) *
                             shapefunctions_ans[index_q23];
            double dNq2dji = shapefunctions_q[index_q23].derivatives_(1, jnod) *
                             shapefunctions_q[index_q23].shapefunctions_(inod) *
                             shapefunctions_ans[index_q23];
            tmp3 += (dNq1dji * stress_enh(2) + dNq2dji * stress_enh(4)) * integration_fac;
            tmp4 += (dNq1dij * stress_enh(2) + dNq2dij * stress_enh(4)) * integration_fac;
          }
        }

        const double tmp5 =
            ((dN1dij + dN1dji) * stress_enh(8) + (dN2dij + dN2dji) * stress_enh(10)) *
            integration_fac;
        const double tmp6 = (Ni * Nj * stress_enh(5)) * integration_fac;

        for (int d = 0; d < num_dim; ++d)
        {
          stiffness_matrix(inod * nodedof + d, jnod * nodedof + d) += tmp1;
          stiffness_matrix(inod * nodedof + d + num_dim, jnod * nodedof + d) += (tmp2 + tmp3);
          stiffness_matrix(inod * nodedof + d, jnod * nodedof + d + num_dim) += (tmp2 + tmp4);
          stiffness_matrix(inod * nodedof + d + num_dim, jnod * nodedof + d + num_dim) +=
              (tmp5 + tmp6);
        }
        if (inod != jnod)
        {
          for (int d = 0; d < num_dim; ++d)
          {
            stiffness_matrix(jnod * nodedof + d, inod * nodedof + d) += tmp1;
            stiffness_matrix(jnod * nodedof + d, inod * nodedof + d + num_dim) += (tmp2 + tmp3);
            stiffness_matrix(jnod * nodedof + d + num_dim, inod * nodedof + d) += (tmp2 + tmp4);
            stiffness_matrix(jnod * nodedof + d + num_dim, inod * nodedof + d + num_dim) +=
                (tmp5 + tmp6);
          }
        }
      }
    }
  }

  /*!
   * @brief Adds the internal force vector contribution of one Gauss point
   *
   * @tparam distype : The discretization type known at compile time
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress_enh (in) :  An object holding the enhanced stress resultants
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddInternalForceVector(const CORE::LINALG::SerialDenseMatrix& Bop,
      const CORE::LINALG::SerialDenseVector& stress_enh, const double integration_fac,
      CORE::LINALG::SerialDenseVector& force_vector)
  {
    force_vector.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, integration_fac, Bop, stress_enh, 1.);
  }
  /*!
   * @brief Adds mass matrix contribution of one Gauss point
   *
   * @tparam distype : The discretization type known at compile time
   * @param shapefunctions (in) : Shape functions and the first derivatives w.r.t spatial
   * coordinates evaluated at the respective point in the parameter space
   * @param mass_matrix_variables (in) : Thickness Integration factor for mass matrix
   * @param mass_matrix (in/out) : Mass matrix where the local contribution is added to
   */
  template <CORE::FE::CellType distype>
  void AddMassMatrix(
      const DRT::ELEMENTS::SHELL::ShapefunctionsAndDerivatives<distype>& shapefunctions,
      DRT::ELEMENTS::SHELL::MassMatrixVariables& mass_matrix_variables, const double& thickness,
      CORE::LINALG::SerialDenseMatrix& massmatrix)
  {
    // half element thickness at gaussian point
    double half_thickness = 0.0;
    for (int i = 0; i < DETAIL::num_node<distype>; ++i)
      half_thickness += thickness * shapefunctions.shapefunctions_(i);
    half_thickness *= 0.5;

    const double squared_half_thickness = half_thickness * half_thickness;

    // integrate consistent mass matrix
    double massfactor;
    const int nodedof = DETAIL::node_dof;
    for (int inod = 0; inod < DETAIL::num_node<distype>; ++inod)
    {
      for (int jnod = 0; jnod < DETAIL::num_node<distype>; ++jnod)
      {
        massfactor = shapefunctions.shapefunctions_(inod) * shapefunctions.shapefunctions_(jnod);
        double massfactor_v = massfactor * mass_matrix_variables.factor_v_;

        massmatrix(nodedof * jnod + 0, nodedof * inod + 0) += massfactor_v;
        massmatrix(nodedof * jnod + 1, nodedof * inod + 1) += massfactor_v;
        massmatrix(nodedof * jnod + 2, nodedof * inod + 2) += massfactor_v;

        double massfactor_w = massfactor * mass_matrix_variables.factor_w_ * squared_half_thickness;
        massmatrix(nodedof * jnod + 3, nodedof * inod + 3) += massfactor_w;
        massmatrix(nodedof * jnod + 4, nodedof * inod + 4) += massfactor_w;
        massmatrix(nodedof * jnod + 5, nodedof * inod + 5) += massfactor_w;

        if (std::abs(mass_matrix_variables.factor_vw_) > 1.0e-14)
        {
          double massfactor_vw = massfactor * mass_matrix_variables.factor_vw_ * half_thickness;
          massmatrix(nodedof * jnod + 0, nodedof * inod + 3) += massfactor_vw;
          massmatrix(nodedof * jnod + 1, nodedof * inod + 4) += massfactor_vw;
          massmatrix(nodedof * jnod + 2, nodedof * inod + 5) += massfactor_vw;
          massmatrix(nodedof * jnod + 3, nodedof * inod + 0) += massfactor_vw;
          massmatrix(nodedof * jnod + 4, nodedof * inod + 1) += massfactor_vw;
          massmatrix(nodedof * jnod + 5, nodedof * inod + 2) += massfactor_vw;
        }
      }
    }
  }

  /*!
   * @brief Calls the @p gp_evaluator to evaluate the jacobian mapping at each integration point of
   * @p intpoints_ of the mid-surface of the shell.
   *
   * @tparam distype  : The discretization type known at compile time
   * @param nodal_coordinates (in) : The nodal coordinates of the element
   * @param intpoints_ (in) : Integratiopn points .
   * @param gp_evaluator (in) : A callable object (e.g. lambda-function) with signature void(const
   * std::vector<double>& xi_gp, const SHELL::ShapefunctionsAndDerivatives<distype>&
   * shape_functions, SHELL::BasisVectorsAndMetrics<distype>& a_current,
   * SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp) that
   * will be called for each integration point.
   */
  template <CORE::FE::CellType distype, typename GaussPointEvaluator>
  inline void ForEachGaussPoint(DRT::ELEMENTS::SHELL::NodalCoordinates<distype>& nodal_coordinates,
      const CORE::FE::IntegrationPoints2D& intpoints_, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < intpoints_.NumPoints(); ++gp)
    {
      // get gauss points from integration rule
      std::array<double, 2> xi_gp;
      xi_gp[0] = intpoints_.qxg[gp][0];
      xi_gp[1] = intpoints_.qxg[gp][1];

      // get gauss weight at current gp
      double gpweight = intpoints_.qwgt[gp];

      // get shape functions and derivatives at gaussian points
      ShapefunctionsAndDerivatives<distype> shapefunctions =
          EvaluateShapefunctionsAndDerivs<distype>(xi_gp);

      // basis vectors and metric on mid-surface
      BasisVectorsAndMetrics<distype> a_reference;
      BasisVectorsAndMetrics<distype> a_current;

      EvaluateMetrics(shapefunctions, a_reference, a_current, nodal_coordinates, 0.0);

      // make h as cross product in ref configuration to get area da on shell mid-surface
      CORE::LINALG::Matrix<DETAIL::num_dim, 1> h(true);
      {
        CORE::LINALG::Matrix<DETAIL::num_dim, DETAIL::num_dim> akovrefe = a_reference.kovariant_;
        h(0) = akovrefe(0, 1) * akovrefe(1, 2) - akovrefe(0, 2) * akovrefe(1, 1);
        h(1) = akovrefe(0, 2) * akovrefe(1, 0) - akovrefe(0, 0) * akovrefe(1, 2);
        h(2) = akovrefe(0, 0) * akovrefe(1, 1) - akovrefe(0, 1) * akovrefe(1, 0);
      }
      // make director unit length and get mid-surface area da from it
      double da = h.Norm2();

      gp_evaluator(xi_gp, shapefunctions, a_current, a_reference, gpweight, da, gp);
    }
  }

}  // namespace DRT::ELEMENTS::SHELL

FOUR_C_NAMESPACE_CLOSE

#endif
