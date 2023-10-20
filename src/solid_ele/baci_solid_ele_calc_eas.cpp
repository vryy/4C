/*! \file

\brief Implementation of routines for calculation of solid element with EAS element technology

\level 1
*/

#include "baci_solid_ele_calc_eas.H"

#include "baci_lib_utils.H"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_mat_so3_material.H"
#include "baci_solid_ele.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_utils.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>


namespace
{
  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int num_nodes =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int num_dim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int num_str = num_dim<distype>*(num_dim<distype> + 1) / 2;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int num_dof_per_ele = num_nodes<distype>* num_dim<distype>;

  /*!
   * @brief Solve for the inverse of a matrix and throw errors if unsuccessful
   *
   * @tparam dim : matrix dimensions
   * @param matrix(in/out) : matrix to be inverted
   */
  template <int dim>
  void SolveForInverse(CORE::LINALG::Matrix<dim, dim>& matrix)
  {
    CORE::LINALG::FixedSizeSerialDenseSolver<dim, dim, 1> solve_for_inverse;
    solve_for_inverse.SetMatrix(matrix);

    int err_fac = solve_for_inverse.Factor();
    if (err_fac != 0) dserror("Factorization of matrix during inversion failed");

    int err_inv = solve_for_inverse.Invert();
    if (err_inv != 0) dserror("Inversion of matrix failed");
  }

  template <DRT::Element::DiscretizationType distype>
  struct CentroidTransformation
  {
    // transformation matrix T0^{-T}, which maps the matrix M from parameter space to the material
    // configuration see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805
    CORE::LINALG::Matrix<num_str<distype>, num_str<distype>> T0invT_;

    // Jacobi determinant evaluated at the element centroid
    double detJ0_;
  };

  /*!
   * @brief Evaluates and returns the transformation matrix T0^{-T} which maps the matrix M from
   * parameter space to the material configuration
   *
   * For details, see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805.
   *
   * @tparam distype : Discretization type
   * @param jacobian_centroid(in) : Jacobian mapping evaluated at the element centroid
   * @return double : transformation matrix
   */
  template <DRT::Element::DiscretizationType distype>
  CORE::LINALG::Matrix<num_str<distype>, num_str<distype>> EvaluateT0invT(
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_centroid)
  {
    // build T0^T (based on strain-like Voigt notation: xx,yy,zz,xy,yz,xz)
    // currently only works in 3D
    CORE::LINALG::Matrix<num_str<distype>, num_str<distype>> T0invT(false);
    T0invT(0, 0) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(0, 0);
    T0invT(1, 0) = jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(1, 0);
    T0invT(2, 0) = jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(2, 0);
    T0invT(3, 0) = 2 * jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(1, 0);
    T0invT(4, 0) = 2 * jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(2, 0);
    T0invT(5, 0) = 2 * jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(2, 0);

    T0invT(0, 1) = jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(0, 1);
    T0invT(1, 1) = jacobian_centroid.jacobian_(1, 1) * jacobian_centroid.jacobian_(1, 1);
    T0invT(2, 1) = jacobian_centroid.jacobian_(2, 1) * jacobian_centroid.jacobian_(2, 1);
    T0invT(3, 1) = 2 * jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(1, 1);
    T0invT(4, 1) = 2 * jacobian_centroid.jacobian_(1, 1) * jacobian_centroid.jacobian_(2, 1);
    T0invT(5, 1) = 2 * jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(2, 1);

    T0invT(0, 2) = jacobian_centroid.jacobian_(0, 2) * jacobian_centroid.jacobian_(0, 2);
    T0invT(1, 2) = jacobian_centroid.jacobian_(1, 2) * jacobian_centroid.jacobian_(1, 2);
    T0invT(2, 2) = jacobian_centroid.jacobian_(2, 2) * jacobian_centroid.jacobian_(2, 2);
    T0invT(3, 2) = 2 * jacobian_centroid.jacobian_(0, 2) * jacobian_centroid.jacobian_(1, 2);
    T0invT(4, 2) = 2 * jacobian_centroid.jacobian_(1, 2) * jacobian_centroid.jacobian_(2, 2);
    T0invT(5, 2) = 2 * jacobian_centroid.jacobian_(0, 2) * jacobian_centroid.jacobian_(2, 2);

    T0invT(0, 3) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(0, 1);
    T0invT(1, 3) = jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(1, 1);
    T0invT(2, 3) = jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(2, 1);
    T0invT(3, 3) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(1, 1) +
                   jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(0, 1);
    T0invT(4, 3) = jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(2, 1) +
                   jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(1, 1);
    T0invT(5, 3) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(2, 1) +
                   jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(0, 1);

    T0invT(0, 4) = jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(0, 2);
    T0invT(1, 4) = jacobian_centroid.jacobian_(1, 1) * jacobian_centroid.jacobian_(1, 2);
    T0invT(2, 4) = jacobian_centroid.jacobian_(2, 1) * jacobian_centroid.jacobian_(2, 2);
    T0invT(3, 4) = jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(1, 2) +
                   jacobian_centroid.jacobian_(1, 1) * jacobian_centroid.jacobian_(0, 2);
    T0invT(4, 4) = jacobian_centroid.jacobian_(1, 1) * jacobian_centroid.jacobian_(2, 2) +
                   jacobian_centroid.jacobian_(2, 1) * jacobian_centroid.jacobian_(1, 2);
    T0invT(5, 4) = jacobian_centroid.jacobian_(0, 1) * jacobian_centroid.jacobian_(2, 2) +
                   jacobian_centroid.jacobian_(2, 1) * jacobian_centroid.jacobian_(0, 2);

    T0invT(0, 5) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(0, 2);
    T0invT(1, 5) = jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(1, 2);
    T0invT(2, 5) = jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(2, 2);
    T0invT(3, 5) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(1, 2) +
                   jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(0, 2);
    T0invT(4, 5) = jacobian_centroid.jacobian_(1, 0) * jacobian_centroid.jacobian_(2, 2) +
                   jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(1, 2);
    T0invT(5, 5) = jacobian_centroid.jacobian_(0, 0) * jacobian_centroid.jacobian_(2, 2) +
                   jacobian_centroid.jacobian_(2, 0) * jacobian_centroid.jacobian_(0, 2);

    // evaluate the inverse T0^{-T} with solver
    SolveForInverse<num_str<distype>>(T0invT);

    return T0invT;
  }

  /*!
   * @brief Evaluates and returns the centroid transformation quantities, i.e., the jacobi
   * determinant at the element centroid and the transformation matrix T0^{-T}
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates(in) : reference and current coordinates of the nodes of the element
   * @return CentroidTransformation<distype> : Jacobi determinant at the element centroid and
   * transformation matrix T0^{-T}
   */
  template <DRT::Element::DiscretizationType distype>
  CentroidTransformation<distype> EvaluateCentroidTransformation(
      const DRT::ELEMENTS::NodalCoordinates<distype>& nodal_coordinates)
  {
    CentroidTransformation<distype> centroid_transformation;

    // 1) compute jacobian at element centroid
    const DRT::ELEMENTS::JacobianMapping<distype> jacobian_mapping_centroid =
        DRT::ELEMENTS::EvaluateJacobianMappingCentroid(nodal_coordinates);

    centroid_transformation.detJ0_ = jacobian_mapping_centroid.determinant_;

    // 2) compute matrix T0^{-T}: T0^{-T} maps the matrix M from local to global coordinates
    centroid_transformation.T0invT_ = EvaluateT0invT(jacobian_mapping_centroid);

    return centroid_transformation;
  }

  /*!
   * @brief Extracts and returns the residual displacement
   *
   * @tparam distype : Discretization type
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @return double : residual displacement or displacement increment
   */
  template <DRT::Element::DiscretizationType distype>
  CORE::LINALG::Matrix<num_dof_per_ele<distype>, 1> GetDisplacementIncrement(
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    auto residual_from_dis = discretization.GetState("residual displacement");
    std::vector<double> residual(lm.size());
    DRT::UTILS::ExtractMyValues(*residual_from_dis, residual, lm);
    CORE::LINALG::Matrix<num_dof_per_ele<distype>, 1> displ_inc(false);
    for (int i = 0; i < num_dof_per_ele<distype>; ++i) displ_inc(i) = residual[i];

    return displ_inc;
  }

  /*!
   * @brief Evaluates and returns the enhanced strains scalar increment
   *
   * @tparam distype, eastype
   * @param displ_inc(in) : displacement increment delta_D_{i+1}
   * @param eas_iteration_data(in) : EAS matrices and vectors from iteration i
   * @return double : enhanced strains scalar increment delta_alpha_{i+1}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1> EvaluateAlphaIncrement(
      const CORE::LINALG::Matrix<num_dof_per_ele<distype>, 1>& displ_inc,
      const DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data)
  {
    // the enhanced strains scalar increment is computed to:
    // delta_alpha_{i+1} = - invKaa_{i} (s_{i} + Kad_{i} delta_D_{i+1})
    CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1> alpha_inc(true);

    // init as enhancement vector s_{i}
    CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1> tmp(
        eas_iteration_data.s_);

    // addition of Kad_{i} delta_D_{i+1}
    tmp.MultiplyTN(1.0, eas_iteration_data.Kda_, displ_inc, 1.0);

    // multiplication with (- invKaa_{i})
    alpha_inc.Multiply(-1.0, eas_iteration_data.invKaa_, tmp);

    return alpha_inc;
  }

  /*!
   * @brief Evaluates the enhanced strain scalars and updates eas_iteration_data.alpha_
   * accordingly
   *
   * @tparam distype, eastype
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void EvaluateAlpha(DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    // residual displacement at the previous step
    CORE::LINALG::Matrix<num_dof_per_ele<distype>, 1> displ_inc(false);
    displ_inc = GetDisplacementIncrement<distype>(discretization, lm);

    // compute the enhanced strain scalar increment delta_alpha
    CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1> alpha_inc =
        EvaluateAlphaIncrement<distype, eastype>(displ_inc, eas_iteration_data);

    // update alpha_i with the increment delta_alpha such that alpha_{i+1} = alpha_{i} + delta_alpha
    eas_iteration_data.alpha_.Update(1.0, alpha_inc, 1.0);
  }

  /*!
   * @brief Compute the matrix M which is the element-wise matrix of the shape functions for the
   * enhanced strains in the parameter space
   *
   * @tparam distype, eastype
   * @param xi(in) : coordinate in the parameter space
   * @return CORE::LINALG::Matrix<num_str, num_eas> : enhanced strains shape function matrix in
   * parameter space
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
  EvaluateEASShapeFunctionsParameterSpace(const CORE::LINALG::Matrix<num_dim<distype>, 1>& xi)
  {
    CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas> M(
        true);

    switch (eastype)
    {
      /* easmild is the EAS interpolation of 9 modes, based on
      **            r 0 0   0 0 0 0 0 0
      **            0 s 0   0 0 0 0 0 0
      **    M =     0 0 t   0 0 0 0 0 0
      **            0 0 0   r s 0 0 0 0
      **            0 0 0   0 0 s t 0 0
      **            0 0 0   0 0 0 0 r t
      */
      case STR::ELEMENTS::EasType::eastype_h8_9:
      {
        M(0, 0) = xi(0);
        M(1, 1) = xi(1);
        M(2, 2) = xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);
        M(4, 5) = xi(1);
        M(4, 6) = xi(2);
        M(5, 7) = xi(0);
        M(5, 8) = xi(2);

        break;
      }
      /* easfull is the EAS interpolation of 21 modes, based on
      **            r 0 0   0 0 0 0 0 0   0  0  0  0  0  0   rs rt 0  0  0  0
      **            0 s 0   0 0 0 0 0 0   0  0  0  0  0  0   0  0  rs st 0  0
      **    M =     0 0 t   0 0 0 0 0 0   0  0  0  0  0  0   0  0  0  0  rt st
      **            0 0 0   r s 0 0 0 0   rt st 0  0  0  0   0  0  0  0  0  0
      **            0 0 0   0 0 s t 0 0   0  0  rs rt 0  0   0  0  0  0  0  0
      **            0 0 0   0 0 0 0 r t   0  0  0  0  rs st  0  0  0  0  0  0
      */
      case STR::ELEMENTS::EasType::eastype_h8_21:
      {
        M(0, 0) = xi(0);
        M(0, 15) = xi(0) * xi(1);
        M(0, 16) = xi(0) * xi(2);
        M(1, 1) = xi(1);
        M(1, 17) = xi(0) * xi(1);
        M(1, 18) = xi(1) * xi(2);
        M(2, 2) = xi(2);
        M(2, 19) = xi(0) * xi(2);
        M(2, 20) = xi(1) * xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);
        M(3, 9) = xi(0) * xi(2);
        M(3, 10) = xi(1) * xi(2);
        M(4, 5) = xi(1);
        M(4, 6) = xi(2);
        M(4, 11) = xi(0) * xi(1);
        M(4, 12) = xi(0) * xi(2);
        M(5, 7) = xi(0);
        M(5, 8) = xi(2);
        M(5, 13) = xi(0) * xi(1);
        M(5, 14) = xi(1) * xi(2);

        break;
      }
      /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
      ** consisting of 7 modes, based on
      **            r 0 0   0 0 0  0
      **            0 s 0   0 0 0  0
      **    M =     0 0 t   0 0 rt st
      **            0 0 0   r s 0  0
      **            0 0 0   0 0 0  0
      **            0 0 0   0 0 0  0
      */
      case STR::ELEMENTS::EasType::eastype_sh8_7:
      {
        /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
        ** consisting of 7 modes, based on
        **            r 0 0   0 0 0  0
        **            0 s 0   0 0 0  0
        **    M =     0 0 t   0 0 rt st
        **            0 0 0   r s 0  0
        **            0 0 0   0 0 0  0
        **            0 0 0   0 0 0  0
        */
        M(0, 0) = xi(0);
        M(1, 1) = xi(1);
        M(2, 2) = xi(2);
        M(2, 5) = xi(0) * xi(2);
        M(2, 6) = xi(1) * xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);

        break;
      }

      default:
        dserror("unknown EAS type");
        break;
    }
    return M;
  }

  /*!
   * @brief Map the matrix M in the parameter space to Mtilde in the material configuration and
   * return Mtilde
   *
   * @tparam distype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant at
   * element centroid
   * @param M(in) : matrix M in the parameter space
   * @return CORE::LINALG::Matrix<num_str, num_eas> : matrix Mtilde in the material configuration
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
  MapEASShapeFunctionsToMaterialConfig(const double detJ,
      const CentroidTransformation<distype>& centroid_transformation,
      const CORE::LINALG::Matrix<num_str<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& M)
  {
    CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas> Mtilde;

    // Mtilde = detJ0/detJ T0^{-T} M
    Mtilde.Multiply(centroid_transformation.detJ0_ / detJ, centroid_transformation.T0invT_, M);

    return Mtilde;
  }

  /*!
   * @brief Evaluate the element-wise matrix of the shape functions for the enhanced strains in the
   * parameter space Mtilde. Therefore set up M (in the material configuration) and map M to Mtilde
   * via T0^{-T}.
   *
   * @tparam distype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant at
   * element centroid
   * @param xi(in) : coordinate in the parameter space
   * @return CORE::LINALG::Matrix<num_str, num_eas> : matrix Mtilde in the material configuration
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
  EvaluateEASShapeFunctionsMaterialConfig(const double detJ,
      const CentroidTransformation<distype>& centroid_transformation,
      const CORE::LINALG::Matrix<num_dim<distype>, 1>& xi)
  {
    CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas> M(
        EvaluateEASShapeFunctionsParameterSpace<distype, eastype>(xi));
    CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
        Mtilde = MapEASShapeFunctionsToMaterialConfig<distype, eastype>(
            detJ, centroid_transformation, M);
    return Mtilde;
  }

  /*!
   * @brief Add the enhanced assumed Green-Lagrange strains E^{enh} = Mtilde alpha to the
   * conventional Green-Lagrange strains E^{u}
   *
   * Background: Choose deformation gradient F as sum of displacement-based F^{u} and enhanced
   * gradient F^{enh}. Considering F_0 the deformation gradient evaluated at the element centroid,
   * F^{enh} is computed to F^{enh} = F_0^{u} Mtilde alpha.
   *
   * @tparam distype, eastype
   * @param gl_strain(in) : Green-Lagrange strains E^{u}
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return CORE::LINALG::Matrix<num_str, 1>  : enhanced Green-Lagrange strains E^{enh}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<num_str<distype>, 1> EvaluateEnhancedAssumedGLStrains(
      const CORE::LINALG::Matrix<num_str<distype>, 1>& gl_strain,
      const CORE::LINALG::Matrix<num_str<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    CORE::LINALG::Matrix<num_str<distype>, 1> enhanced_gl_strain(gl_strain);
    enhanced_gl_strain.Multiply(1.0, Mtilde, alpha, 1.0);
    return enhanced_gl_strain;
  }

  /*!
   * @brief Evaluate the enhanced assumed Green-Lagrange strains E^{enh}

   * @tparam distype, eastype
   * @param displacement_based_mapping(in) : displacement-based spatial mapping
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return CORE::LINALG::Matrix<num_str, 1> : Enhanced Green-Lagrange strains E^{enh}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  CORE::LINALG::Matrix<num_str<distype>, 1> EvaluateEnhancedAssumedGLStrains(
      const DRT::ELEMENTS::SpatialMaterialMapping<distype>& displacement_based_mapping,
      const CORE::LINALG::Matrix<num_str<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    const CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> displacement_based_cauchygreen =
        DRT::ELEMENTS::EvaluateCauchyGreen<distype>(displacement_based_mapping);

    const CORE::LINALG::Matrix<num_str<distype>, 1> gl_strain =
        DRT::ELEMENTS::EvaluateGreenLagrangeStrain<distype>(displacement_based_cauchygreen);

    return EvaluateEnhancedAssumedGLStrains<distype, eastype>(gl_strain, Mtilde, alpha);
  }

  /*!
   * @brief Compute the enhanced deformation gradient F^{enh}

   * @tparam dim
   * @param defgrd_disp(in) : displacement-based deformation gradient F^{u}
   * @param enhanced_gl_strain(in) : enhanced Green-Lagrange strains E^{enh}
   * @return CORE::LINALG::Matrix<dim, dim> : enhanced deformation gradient F^{enh}
   */
  template <unsigned dim>
  CORE::LINALG::Matrix<dim, dim> EvaluateConsistentDefgrd(
      const CORE::LINALG::Matrix<dim, dim>& defgrd_disp,
      const CORE::LINALG::Matrix<dim*(dim + 1) / 2, 1>& enhanced_gl_strain)
  {
    CORE::LINALG::Matrix<dim, dim> R;       // rotation tensor
    CORE::LINALG::Matrix<dim, dim> U_enh;   // enhanced right stretch tensor
    CORE::LINALG::Matrix<dim, dim> U_disp;  // displacement-based right stretch tensor
    CORE::LINALG::Matrix<dim, dim> EW;      // temporarily store eigenvalues
    CORE::LINALG::Matrix<dim, dim> tmp;     // temporary matrix for matrix matrix matrix products
    CORE::LINALG::Matrix<dim, dim> tmp2;    // temporary matrix for matrix matrix matrix products

    // calculate modified right stretch tensor
    if (dim != 3) dserror("stop: this currently only works for 3D");
    for (unsigned i = 0; i < dim; i++) U_enh(i, i) = 2. * enhanced_gl_strain(i) + 1.;
    U_enh(0, 1) = enhanced_gl_strain(dim);
    U_enh(1, 0) = enhanced_gl_strain(dim);
    U_enh(1, 2) = enhanced_gl_strain(4);
    U_enh(2, 1) = enhanced_gl_strain(4);
    U_enh(0, 2) = enhanced_gl_strain(5);
    U_enh(2, 0) = enhanced_gl_strain(5);

    CORE::LINALG::SYEV(U_enh, EW, U_enh);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_enh, EW);
    tmp2.MultiplyNT(tmp, U_enh);
    U_enh.Update(tmp2);

    // calculate displacement-based right stretch tensor
    U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

    CORE::LINALG::SYEV(U_disp, EW, U_disp);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_disp, EW);
    tmp2.MultiplyNT(tmp, U_disp);
    U_disp.Update(tmp2);

    // compose consistent deformation gradient
    U_disp.Invert();
    R.Multiply(defgrd_disp, U_disp);

    CORE::LINALG::Matrix<dim, dim> defgrd_enh;
    defgrd_enh.Multiply(R, U_enh);
    return defgrd_enh;
  }

  /*!
   * @brief Integrate the EAS stiffness matrices
   *
   * @tparam distype, eastype
   * @param stress(in) : 2. Piola Kirchhoff stress tensor and material tangent
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param Bop(in) : B-operator
   * @param integration_factor(in) : integration factor (Gauss point weight times Jacobi
   * determinant)
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void IntegrateEAS(const DRT::ELEMENTS::Stress<distype>& stress,
      const CORE::LINALG::Matrix<num_str<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const CORE::LINALG::Matrix<num_str<distype>, num_dof_per_ele<distype>>& Bop,
      const double integration_factor,
      DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data)
  {
    // integrate Kaa: Kaa += (Mtilde^T . cmat . Mtilde) * detJ * w(gp)
    // IMPORTANT: We save this in invKaa_ here since after the loop over all Gauss points, we
    // invert the matrix. At this point, this is still Kaa and NOT invKaa.
    CORE::LINALG::Matrix<num_str<distype>, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas> cmatM(
        true);
    cmatM.Multiply(stress.cmat_, Mtilde);
    eas_iteration_data.invKaa_.MultiplyTN(integration_factor, Mtilde, cmatM, 1.);

    // integrate Kda: Kda += (B^T . cmat . Mtilde) * detJ * w(gp)
    eas_iteration_data.Kda_.MultiplyTN(integration_factor, Bop, cmatM, 1.);

    // integrate s: s += (Mtilde^T . S) * detJ * w(gp)
    eas_iteration_data.s_.MultiplyTN(integration_factor, Mtilde, stress.pk2_, 1.);
  }

  /*!
   * @brief Add EAS internal force contribution of one Gauss point
   *
   * The EAS internal force contribution is $f_{eas} = - K_{da} K_{aa}^{-1} s$.
   *
   * @tparam distype : Discretization type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param s(in) : enhancement vector s
   * @param force(in/out) : internal force vector where the contribution is added to
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void AddEASInternalForce(const CORE::LINALG::Matrix<num_dof_per_ele<distype>,
                               STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const CORE::LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1>& s,
      CORE::LINALG::Matrix<num_dof_per_ele<distype>, 1>& force_vector)
  {
    force_vector.MultiplyNN(1.0, minusKdainvKaa, s, 1.0);
  }

  /*!
   * @brief Add EAS stiffness matrix contribution of one Gauss point
   *
   * The EAS stiffness matrix contribution is $- K_{da} K_{aa}^{-1} K_{da}^T$.
   *
   * @tparam distype : Discretization type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param Kda(in) : EAS stiffness matrix part K_{da}
   * @param stiffness_matrix(in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void AddEASStiffnessMatrix(const CORE::LINALG::Matrix<num_dof_per_ele<distype>,
                                 STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const CORE::LINALG::Matrix<num_dof_per_ele<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>& Kda,
      CORE::LINALG::Matrix<num_dof_per_ele<distype>, num_dof_per_ele<distype>>& stiffness_matrix)
  {
    stiffness_matrix.MultiplyNT(1.0, minusKdainvKaa, Kda, 1.0);
  }
}  // namespace

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::SolidEleCalcEas()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Pack(DRT::PackBuffer& data) const
{
  constexpr int num_dof_per_element =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement *
      CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1>(
      data, eas_iteration_data_.alpha_);
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas, 1>(
      data, eas_iteration_data_.s_);
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>(data, eas_iteration_data_.invKaa_);
  DRT::ELEMENTS::Solid::AddtoPack<num_dof_per_element,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>(data, eas_iteration_data_.Kda_);
};

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.alpha_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.s_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.invKaa_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.Kda_);
};

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff = {};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass = {};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);

  double mean_density = 0.0;

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  // clear for integration
  eas_iteration_data_.invKaa_.Clear();
  eas_iteration_data_.Kda_.Clear();
  eas_iteration_data_.s_.Clear();

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, displacement_based_spatial_material_mapping);

        const CORE::LINALG::Matrix<num_str_, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
            Mtilde = EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);

        const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>
            consistent_defgrd = EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        const Stress<distype> stress = EvaluateMaterialStress<distype>(
            solid_material, consistent_defgrd, enhanced_gl_strain, params, gp, ele.Id());

        IntegrateEAS<distype, eastype>(
            stress, Mtilde, Bop, integration_factor, eas_iteration_data_);

        if (force.has_value())
        {
          AddInternalForceVector(Bop, stress, integration_factor, *force);
        }

        if (stiff.has_value())
        {
          AddElasticStiffnessMatrix(Bop, stress, integration_factor, *stiff);
          AddGeometricStiffnessMatrix(jacobian_mapping.N_XYZ_, stress, integration_factor, *stiff);
        }

        if (mass.has_value())
        {
          if (equal_integration_mass_stiffness)
          {
            AddMassMatrix(shape_functions, integration_factor, solid_material.Density(gp), *mass);
          }
          else
          {
            mean_density += solid_material.Density(gp) / stiffness_matrix_integration_.NumPoints();
          }
        }
      });

  // invert Kaa with solver. eas_iteration_data_.invKaa_ then is Kaa^{-1}
  SolveForInverse<STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>(eas_iteration_data_.invKaa_);

  // compute the product (- Kda Kaa^{-1}) which is later needed for force and stiffness update
  CORE::LINALG::Matrix<num_dof_per_ele_, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
      minusKdainvKaa(true);
  minusKdainvKaa.MultiplyNN(-1.0, eas_iteration_data_.Kda_, eas_iteration_data_.invKaa_);

  if (force.has_value())
  {
    AddEASInternalForce<distype, eastype>(minusKdainvKaa, eas_iteration_data_.s_, *force);
  }

  if (stiff.has_value())
  {
    AddEASStiffnessMatrix<distype, eastype>(minusKdainvKaa, eas_iteration_data_.Kda_, *stiff);
  }

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    ForEachGaussPoint<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateNonlinearForceStiffnessMassGEMM(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  dserror("GEMM is not implemented for EAS elements.");
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  dserror("Recovering is not yet implemented for EAS elements");
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);
  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  // No need to update alpha here. Update is called to copy states from t_{n+1} to t_{n} after the
  // time step and output. Hence, there are no more Newton iterations that would require an update
  // of alpha

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const CORE::LINALG::Matrix<num_str_, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
            Mtilde = EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);

        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>
            consistent_defgrd = EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        solid_material.Update(consistent_defgrd, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::CalculateStress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const StressIO& stressIO, const StrainIO& strainIO,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const CORE::LINALG::Matrix<num_str_, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
            Mtilde = EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);


        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, DETAIL::num_dim<distype>>
            consistent_defgrd = EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        const Stress<distype> stress = EvaluateMaterialStress<distype>(
            solid_material, consistent_defgrd, enhanced_gl_strain, params, gp, ele.Id());

        AssembleStrainTypeToMatrixRow<distype>(
            enhanced_gl_strain, consistent_defgrd, strainIO.type, strain_data, gp);
        AssembleStressTypeToMatrixRow(consistent_defgrd, stress, stressIO.type, stress_data, gp);
      });

  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
double DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::CalculateInternalEnergy(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const CORE::LINALG::Matrix<num_str_, STR::ELEMENTS::EasTypeToNumEas<eastype>::num_eas>
            Mtilde = EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);

        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        double psi = 0.0;
        solid_material.StrainEnergy(enhanced_gl_strain, psi, gp, ele.Id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<distype>(
      stiffness_matrix_integration_, ele, params);

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::InitializeGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<distype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8,
    STR::ELEMENTS::EasType::eastype_h8_9>;
template class DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8,
    STR::ELEMENTS::EasType::eastype_h8_21>;

static_assert(
    DRT::ELEMENTS::IsPackable<
        DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, STR::ELEMENTS::EasType::eastype_h8_9>*>,
    "EAS needs to implement the method Pack(DRT::PackBuffer&) to be able to store history data!");
static_assert(
    DRT::ELEMENTS::IsUnpackable<
        DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, STR::ELEMENTS::EasType::eastype_h8_9>*>,
    "EAS needs to implement the method Unpack(std::size_t, std::vector<char>&) to be able to store "
    "history data!");