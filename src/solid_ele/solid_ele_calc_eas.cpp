/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element with EAS element technology
\level 1
*/
/*----------------------------------------------------------------------*/

#include "linalg_utils_densematrix_eigen.H"
#include "solid_ele_calc_lib.H"
#include "solid_ele_calc_eas.H"
#include <Teuchos_ParameterList.hpp>
#include <memory>
#include <optional>
#include "utils_exceptions.H"
#include "lib_utils.H"
#include "solid_ele.H"
#include "lib_discret.H"
#include "mat_so3_material.H"
#include "solid_ele_utils.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "fiber_utils.H"
#include "fiber_nodal_fiber_holder.H"
#include "discretization_fem_general_utils_gauss_point_postprocess.H"
#include "discretization_fem_general_utils_gauss_point_extrapolation.H"

#include "structure_new_gauss_point_data_output_manager.H"
#include "so3_element_service.H"


namespace
{

  template <DRT::Element::DiscretizationType distype>
  struct CentroidTransformation
  {
    // transformation matrix T0^{-T}, which maps the matrix M from parameter space to the material
    // configuration see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, DRT::ELEMENTS::DETAIL::numstr<distype>>
        T0invT_;

    // Jacobi determinant evaluated at the element centroid
    double detJ0_;
  };

  /*!
   * @brief Evaluates and returns the transformation matrix T0^{-T} which maps the matrix M from
   * parameter space to the material configuration
   *
   * For details, see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805.
   *
   * @tparam distype
   * @param jacobian_centroid(in) : Jacobian mapping evaluated at the element centroid
   * @return T0invT : transformation matrix
   */
  template <DRT::Element::DiscretizationType distype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, DRT::ELEMENTS::DETAIL::numstr<distype>>
  EvaluateT0invT(const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_centroid)
  {
    // build T0^T (based on strain-like Voigt notation: xx,yy,zz,xy,yz,xz)
    // currently only works in 3D
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, DRT::ELEMENTS::DETAIL::numstr<distype>>
        T0invT(false);
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
    LINALG::FixedSizeSerialDenseSolver<DRT::ELEMENTS::DETAIL::numstr<distype>,
        DRT::ELEMENTS::DETAIL::numstr<distype>, 1>
        solve_for_inverseT0;
    solve_for_inverseT0.SetMatrix(T0invT);
    int err2 = solve_for_inverseT0.Factor();
    int err = solve_for_inverseT0.Invert();
    if ((err != 0) || (err2 != 0)) dserror("Inversion of T0inv (Jacobian0) failed");

    return T0invT;
  }

  /*!
   * @brief Evaluates and returns the centroid transformation quantities, i.e., the jacobi
   * determinant at the element centroid and the transformation matrix T0^{-T}
   *
   * @tparam distype
   * @param nodal_coordinates(in) : reference and current coordinates of the nodes of the element
   * @return centroid_transformation : Jacobi determinant at the element centroid and
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
   * @tparam distype
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : location vector of this element
   * @return displ_inc : residual displacement or displacement increment
   */
  template <DRT::Element::DiscretizationType distype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> GetDisplacementIncrement(
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    auto residual_from_dis = discretization.GetState("residual displacement");
    std::vector<double> residual(lm.size());
    DRT::UTILS::ExtractMyValues(*residual_from_dis, residual, lm);
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> displ_inc(false);
    for (int i = 0; i < DRT::ELEMENTS::DETAIL::numdofperelement<distype>; ++i)
      displ_inc(i) = residual[i];

    return displ_inc;
  }

  /*!
   * @brief Evaluates and returns the enhanced strains scalar increment
   *
   * @tparam distype, eastype
   * @param displ_inc(in) : displacement increment delta_D_{i+1}
   * @param eas_iteration_data(in) : EAS matrices and vectors from iteration i
   * @return alpha_inc : enhanced strains scalar increment delta_alpha_{i+1}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1> EvaluateAlphaIncrement(
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1>& displ_inc,
      const DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data)
  {
    // the enhanced strains scalar increment is computed to:
    // delta_alpha_{i+1} = - invKaa_{i} (s_{i} + Kad_{i} delta_D_{i+1})
    LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1> alpha_inc(true);

    // init as enhancement vector s_{i} (EAS portion of internal forces)
    LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1> tmp(eas_iteration_data.feas_);

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
   * @param lm(in) : location vector of this element
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void EvaluateAlpha(DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    // residual displacement at the previous step
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> displ_inc(false);
    displ_inc = GetDisplacementIncrement<distype>(discretization, lm);

    // compute the enhanced strain scalar increment delta_alpha
    LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1> alpha_inc =
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
   * @return M : enhanced strains shape function matrix in parameter space
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
  EvaluateEASShapeFunctionsParameterSpace(
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, 1>& xi)
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
        STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
        M(true);

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
   * @return Mtilde : matrix Mtilde in the material configuration
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
  MapEASShapeFunctionsToMaterialConfig(const double detJ,
      const CentroidTransformation<distype>& centroid_transformation,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& M)
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
        STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
        Mtilde;

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
   * @return Mtilde : matrix Mtilde in the material configuration
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
  EvaluateEASShapeFunctionsMaterialConfig(const double detJ,
      const CentroidTransformation<distype>& centroid_transformation,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, 1>& xi)
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
        STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
        M(EvaluateEASShapeFunctionsParameterSpace<distype, eastype>(xi));
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
        STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
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
   * @return enhanced_gl_strain : enhanced Green-Lagrange strains E^{enh}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> EvaluateEnhancedAssumedGLStrains(
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1>& gl_strain,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& Mtilde,
      const LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>& alpha)
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> enhanced_gl_strain(gl_strain);
    enhanced_gl_strain.Multiply(1.0, Mtilde, alpha, 1.0);
    return enhanced_gl_strain;
  }

  /*!
   * @brief Evaluate the enhanced assumed Green-Lagrange strains E^{enh}

   * @tparam distype, eastype
   * @param displacement_based_mapping(in) : displacement-based spatial mapping
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return Enhanced Green-Lagrange strains E^{enh}
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> EvaluateEnhancedAssumedGLStrains(
      const DRT::ELEMENTS::SpatialMaterialMapping<distype>& displacement_based_mapping,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& Mtilde,
      const LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>& alpha)
  {
    const DRT::ELEMENTS::CauchyGreen<distype> displacement_based_cauchygreen =
        EvaluateCauchyGreen(displacement_based_mapping);

    const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> gl_strain =
        EvaluateGreenLagrangeStrain(displacement_based_cauchygreen);

    return EvaluateEnhancedAssumedGLStrains<distype, eastype>(gl_strain, Mtilde, alpha);
  }

  /*!
   * @brief Compute the enhanced deformation gradient F^{enh}

   * @tparam dim
   * @param defgrd_disp(in) : displacement-based deformation gradient F^{u}
   * @param enhanced_gl_strain(in) : enhanced Green-Lagrange strains E^{enh}
   * @return defgrd_enh : enhanced deformation gradient F^{enh}
   */
  template <unsigned dim>
  LINALG::Matrix<dim, dim> EvaluateConsistentDefgrd(const LINALG::Matrix<dim, dim>& defgrd_disp,
      const LINALG::Matrix<dim*(dim + 1) / 2, 1>& enhanced_gl_strain)
  {
    LINALG::Matrix<dim, dim> R;       // rotation tensor
    LINALG::Matrix<dim, dim> U_enh;   // enhanced right stretch tensor
    LINALG::Matrix<dim, dim> U_disp;  // displacement-based right stretch tensor
    LINALG::Matrix<dim, dim> EW;      // temporarily store eigenvalues
    LINALG::Matrix<dim, dim> tmp;     // temporary matrix for matrix matrix matrix products
    LINALG::Matrix<dim, dim> tmp2;    // temporary matrix for matrix matrix matrix products

    // calculate modified right stretch tensor
    if (dim != 3) dserror("stop: this currently only works for 3D");
    for (unsigned i = 0; i < dim; i++) U_enh(i, i) = 2. * enhanced_gl_strain(i) + 1.;
    U_enh(0, 1) = enhanced_gl_strain(dim);
    U_enh(1, 0) = enhanced_gl_strain(dim);
    U_enh(1, 2) = enhanced_gl_strain(4);
    U_enh(2, 1) = enhanced_gl_strain(4);
    U_enh(0, 2) = enhanced_gl_strain(5);
    U_enh(2, 0) = enhanced_gl_strain(5);

    LINALG::SYEV(U_enh, EW, U_enh);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_enh, EW);
    tmp2.MultiplyNT(tmp, U_enh);
    U_enh.Update(tmp2);

    // calculate displacement-based right stretch tensor
    U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

    LINALG::SYEV(U_disp, EW, U_disp);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_disp, EW);
    tmp2.MultiplyNT(tmp, U_disp);
    U_disp.Update(tmp2);

    // compose consistent deformation gradient
    U_disp.Invert();
    R.Multiply(defgrd_disp, U_disp);

    LINALG::Matrix<dim, dim> defgrd_enh;
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
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& Mtilde,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
          DRT::ELEMENTS::DETAIL::numdofperelement<distype>>& Bop,
      const double integration_factor,
      DRT::ELEMENTS::EasIterationData<distype, eastype>& eas_iteration_data)
  {
    // integrate Kaa: Kaa += (Mtilde^T . cmat . Mtilde) * detJ * w(gp)
    // IMPORTANT: We save this in invKaa_ here since after the loop over all Gauss points, we
    // invert the matrix. At this point, this is still Kaa and NOT invKaa.
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
        STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>
        cmatM(true);
    cmatM.Multiply(stress.cmat_, Mtilde);
    eas_iteration_data.invKaa_.MultiplyTN(integration_factor, Mtilde, cmatM, 1.);

    // integrate Kda: Kda += (B^T . cmat . Mtilde) * detJ * w(gp)
    eas_iteration_data.Kda_.MultiplyTN(integration_factor, Bop, cmatM, 1.);

    // integrate feas: feas += (Mtilde^T . sigma) * detJ * w(gp)
    eas_iteration_data.feas_.MultiplyTN(integration_factor, Mtilde, stress.pk2_, 1.);
  }

  /*!
   * @brief Add EAS internal force contribution of one Gauss point
   *
   * The EAS internal force contribution is $- K_{da} K_{aa}^{-1} feas$.
   *
   * @tparam distype
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param feas(in) : EAS portion of internal forces, also called enhancement vector
   * @param force(in/out) : internal force vector where the contribution is added to
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void AddEASInternalForce(const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>,
                               STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& minusKdainvKaa,
      const LINALG::Matrix<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>& feas,
      LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1>& force_vector)
  {
    force_vector.MultiplyNN(1.0, minusKdainvKaa, feas, 1.0);
  }

  /*!
   * @brief Add EAS stiffness matrix contribution of one Gauss point
   *
   * The EAS stiffness matrix contribution is $- K_{da} K_{aa}^{-1} K_{da}^T$.
   *
   * @tparam distype
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param Kda(in) : EAS stiffness matrix part K_{da}
   * @param stiffness_matrix(in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
  void AddEASStiffnessMatrix(const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>,
                                 STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& minusKdainvKaa,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>,
          STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>& Kda,
      LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>,
          DRT::ELEMENTS::DETAIL::numdofperelement<distype>>& stiffness_matrix)
  {
    stiffness_matrix.MultiplyNT(1.0, minusKdainvKaa, Kda, 1.0);
  }
}  // namespace

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::SolidEleCalcEas()
    : DRT::ELEMENTS::SolidEleCalcInterface::SolidEleCalcInterface(),
      stiffness_matrix_integration_(
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
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>(
      data, eas_iteration_data_.alpha_);
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>(
      data, eas_iteration_data_.feas_);
  DRT::ELEMENTS::Solid::AddtoPack<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>(data, eas_iteration_data_.invKaa_);
  DRT::ELEMENTS::Solid::AddtoPack<num_dof_per_element,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas>(data, eas_iteration_data_.Kda_);
};

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.alpha_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.feas_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.invKaa_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.Kda_);
};

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector* force_vector,
    Epetra_SerialDenseMatrix* stiffness_matrix, Epetra_SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<LINALG::Matrix<numdofperelement_, numdofperelement_>> stiff = {};
  std::optional<LINALG::Matrix<numdofperelement_, numdofperelement_>> mass = {};
  std::optional<LINALG::Matrix<numdofperelement_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // TODO: This is a quite unsafe check, whether the same integrations are used
  bool equal_integration_mass_stiffness =
      mass_matrix_integration_.NumPoints() == stiffness_matrix_integration_.NumPoints();

  double mean_density = 0.0;

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  // clear for integration
  eas_iteration_data_.invKaa_.Clear();
  eas_iteration_data_.Kda_.Clear();
  eas_iteration_data_.feas_.Clear();

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        LINALG::Matrix<numstr_, numdofperelement_> Bop =
            EvaluateStrainGradient(jacobian_mapping, displacement_based_spatial_material_mapping);

        const LINALG::Matrix<numstr_, STR::ELEMENTS::EasTypeToNumEas<eastype>::neas> Mtilde =
            EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);

        const LINALG::Matrix<DETAIL::numstr<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> consistent_defgrd =
            EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

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
  LINALG::FixedSizeSerialDenseSolver<STR::ELEMENTS::EasTypeToNumEas<eastype>::neas,
      STR::ELEMENTS::EasTypeToNumEas<eastype>::neas, 1>
      solve_for_invKaa;
  solve_for_invKaa.SetMatrix(eas_iteration_data_.invKaa_);
  int err2 = solve_for_invKaa.Factor();
  int err = solve_for_invKaa.Invert();
  if ((err != 0) || (err2 != 0)) dserror("Inversion of Kaa failed");

  // compute the product (- Kda Kaa^{-1}) which is later needed for force and stiffness update
  LINALG::Matrix<numdofperelement_, STR::ELEMENTS::EasTypeToNumEas<eastype>::neas> minusKdainvKaa(
      true);
  minusKdainvKaa.MultiplyNN(-1.0, eas_iteration_data_.Kda_, eas_iteration_data_.invKaa_);

  if (force.has_value())
  {
    AddEASInternalForce<distype, eastype>(minusKdainvKaa, eas_iteration_data_.feas_, *force);
  }

  if (stiff.has_value())
  {
    AddEASStiffnessMatrix<distype, eastype>(minusKdainvKaa, eas_iteration_data_.Kda_, *stiff);
  }

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {  // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateNonlinearForceStiffnessMassGEMM(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector* force_vector,
    Epetra_SerialDenseMatrix* stiffness_matrix, Epetra_SerialDenseMatrix* mass_matrix)
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

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const LINALG::Matrix<numstr_, STR::ELEMENTS::EasTypeToNumEas<eastype>::neas> Mtilde =
            EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);


        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const LINALG::Matrix<DETAIL::numstr<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> consistent_defgrd =
            EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

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
  // TODO: If we get rid of post_drt_*, we don't need this here anymore. We could directly use
  // InitializeGaussPointDataOutput and EvaluateGaussPointDataOutput and write the stresses there.
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  Epetra_SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), numstr_);
  Epetra_SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), numstr_);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const LINALG::Matrix<numstr_, STR::ELEMENTS::EasTypeToNumEas<eastype>::neas> Mtilde =
            EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);


        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const LINALG::Matrix<DETAIL::numstr<distype>, 1> enhanced_gl_strain =
            EvaluateEnhancedAssumedGLStrains<distype, eastype>(
                displacement_based_spatial_material_mapping, Mtilde, eas_iteration_data_.alpha_);

        const LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> consistent_defgrd =
            EvaluateConsistentDefgrd(
                displacement_based_spatial_material_mapping.deformation_gradient_,
                enhanced_gl_strain);

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
  // need update
  double intenergy = 0.0;
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  CentroidTransformation<distype> centroid_transformation =
      EvaluateCentroidTransformation<distype>(nodal_coordinates);

  EvaluateAlpha<distype, eastype>(eas_iteration_data_, discretization, lm);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const LINALG::Matrix<numstr_, STR::ELEMENTS::EasTypeToNumEas<eastype>::neas> Mtilde =
            EvaluateEASShapeFunctionsMaterialConfig<distype, eastype>(
                jacobian_mapping.determinant_, centroid_transformation, xi);

        const SpatialMaterialMapping<distype> displacement_based_spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const LINALG::Matrix<DETAIL::numstr<distype>, 1> enhanced_gl_strain =
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
  if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const static std::vector<LINALG::Matrix<nen_, 1>> shapefcts = std::invoke(
        [&]
        {
          std::vector<LINALG::Matrix<nen_, 1>> shapefcns(stiffness_matrix_integration_.NumPoints());
          for (int gp = 0; gp < stiffness_matrix_integration_.NumPoints(); ++gp)
          {
            LINALG::Matrix<nsd_, 1> xi(stiffness_matrix_integration_.Point(gp), true);
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

  // Save number of Gauss of the element for gauss point data output
  gp_data_output_manager.AddElementNumberOfGaussPoints(stiffness_matrix_integration_.NumPoints());

  // holder for output quantity names and their size
  std::unordered_map<std::string, int> quantities_map{};

  // Ask material for the output quantity names and sizes
  solid_material.RegisterVtkOutputDataNames(quantities_map);

  // Add quantities to the Gauss point output data manager (if they do not already exist)
  gp_data_output_manager.MergeQuantities(quantities_map);
}

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::EvaluateGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  // Collection and assembly of gauss point data
  for (const auto& quantity : gp_data_output_manager.GetQuantities())
  {
    const std::string& quantity_name = quantity.first;
    const int quantity_size = quantity.second;

    // Step 1: Collect the data for each Gauss point for the material
    LINALG::SerialDenseMatrix gp_data(
        stiffness_matrix_integration_.NumPoints(), quantity_size, true);
    bool data_available = solid_material.EvaluateVtkOutputData(quantity_name, gp_data);

    // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
    // point)
    if (data_available)
    {
      switch (gp_data_output_manager.GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          // compute average of the quantities
          Teuchos::RCP<Epetra_MultiVector> global_data =
              gp_data_output_manager.GetMutableElementCenterData().at(quantity_name);
          CORE::DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> global_data =
              gp_data_output_manager.GetMutableNodalData().at(quantity_name);

          Epetra_IntVector& global_nodal_element_count =
              *gp_data_output_manager.GetMutableNodalDataCount().at(quantity_name);

          CORE::DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<distype>(
              ele, gp_data, *global_data, false, stiffness_matrix_integration_);
          DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
              gp_data_output_manager.GetMutableGaussPointData().at(quantity_name);
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

template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
void DRT::ELEMENTS::SolidEleCalcEas<distype, eastype>::ResetAll(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetAll(stiffness_matrix_integration_.NumPoints());

  eas_iteration_data_ = {};
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
