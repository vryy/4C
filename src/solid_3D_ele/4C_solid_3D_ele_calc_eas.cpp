// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_3D_ele_calc_eas.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"
#include "4C_utils_exceptions.hpp"
#include "internal/4C_solid_3D_ele_calc_eas_helpers.hpp"

#include <Teuchos_dyn_cast.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief Solve for the inverse of a matrix and ignore any errors
   *
   * @tparam dim : matrix dimensions
   * @param matrix(in/out) : matrix to be inverted
   */
  template <unsigned int dim>
  void solve_for_inverse_ignoring_errors(Core::LinAlg::Matrix<dim, dim>& matrix)
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<dim, dim, 1> solve_for_inverse;
    solve_for_inverse.set_matrix(matrix);

    solve_for_inverse.invert();
  }

  /*!
   * @brief Solve for the inverse of a matrix and throw errors if unsuccessful
   *
   * @tparam dim : matrix dimensions
   * @param matrix(in/out) : matrix to be inverted
   */
  template <unsigned int dim>
  void solve_for_inverse(Core::LinAlg::Matrix<dim, dim>& matrix)
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<dim, dim, 1> solve_for_inverse;
    solve_for_inverse.set_matrix(matrix);

    int err_inv = solve_for_inverse.invert();
    if (err_inv != 0) FOUR_C_THROW("Inversion of matrix failed with LAPACK error code %d", err_inv);
  }

  /*!
   * @brief Evaluates and returns the centroid transformation quantities, i.e., the jacobi
   * determinant at the element centroid and the transformation matrix T0^{-T}
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates(in) : reference and current coordinates of the nodes of the element
   * @return CentroidTransformation<celltype> : Jacobi determinant at the element centroid and
   * transformation matrix T0^{-T}
   */
  template <Core::FE::CellType celltype>
  Discret::Elements::CentroidTransformation<celltype> evaluate_centroid_transformation(
      const Discret::Elements::ElementNodes<celltype>& nodal_coordinates)
  {
    Discret::Elements::CentroidTransformation<celltype> centroid_transformation;

    // 1) compute jacobian at element centroid
    const Discret::Elements::JacobianMapping<celltype> jacobian_mapping_centroid =
        Discret::Elements::evaluate_jacobian_mapping_centroid(nodal_coordinates);

    centroid_transformation.detJ0_ = jacobian_mapping_centroid.determinant_;

    // 2) compute matrix T0^{-T}: T0^{-T} maps the matrix M from local to global coordinates
    centroid_transformation.T0invT_ =
        evaluate_voigt_transformation_matrix(jacobian_mapping_centroid);

    return centroid_transformation;
  }

  /*!
   * @brief Extracts and returns the residual displacement
   *
   * @tparam celltype : Cell type
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @return double : residual displacement or displacement increment
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1> get_displacement_increment(
      const Core::FE::Discretization& discretization, const std::vector<int>& lm)
  {
    auto residual_from_dis = discretization.get_state("residual displacement");
    std::vector<double> residual(lm.size());
    Core::FE::extract_my_values(*residual_from_dis, residual, lm);
    Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1> displ_inc(false);
    for (int i = 0; i < Discret::Elements::num_dof_per_ele<celltype>; ++i)
      displ_inc(i) = residual[i];

    return displ_inc;
  }

  /*!
   * @brief Updates the enhanced strains scalar increment
   *
   * @tparam celltype, eastype
   * @param displ_inc(in) : displacement increment delta_D_{i+1}
   * @param eas_iteration_data(in) : EAS matrices and vectors from iteration i
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void update_alpha_increment(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1>& displ_inc,
      Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    // the enhanced strains scalar increment is computed to:
    // delta_alpha_{i+1} = - invKaa_{i} (s_{i} + Kad_{i} delta_D_{i+1})

    // init as enhancement vector s_{i}
    Core::LinAlg::Matrix<Solid::Elements::EasTypeToNumEas<eastype>::num_eas, 1> tmp(
        eas_iteration_data.s_);

    // addition of Kad_{i} delta_D_{i+1}
    tmp.multiply_tn(1.0, eas_iteration_data.Kda_, displ_inc, 1.0);

    // multiplication with (- invKaa_{i})
    eas_iteration_data.alpha_inc_.multiply(-1.0, eas_iteration_data.invKaa_, tmp);
  }

  /*!
   * @brief Updates the enhanced strain scalars, alpha_inc and alpha in the iteration data
   * accordingly
   *
   * @tparam celltype, eastype
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void update_alpha(Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data,
      const Core::FE::Discretization& discretization, const std::vector<int>& lm,
      const double step_length = 1.0)
  {
    // residual displacement at the previous step
    Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1> displ_inc(false);
    displ_inc = get_displacement_increment<celltype>(discretization, lm);

    // compute the enhanced strain scalar increment delta_alpha
    update_alpha_increment<celltype, eastype>(displ_inc, eas_iteration_data);

    // update alpha_i with the increment delta_alpha such that alpha_{i+1} = alpha_{i} + delta_alpha
    eas_iteration_data.alpha_.update(step_length, eas_iteration_data.alpha_inc_, 1.0);
  }

  /*!
   * @brief Correct alpha in a line search step by adapting the step length
   *
   * @tparam celltype
   * @tparam eastype
   * @param eas_iteration_data (in/out) : EAS iteration data
   * @param new_step_length (in) : new step length
   * @param old_step_length (in) : old step length
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void correct_alpha(Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data,
      const double new_step_length, const double old_step_length)
  {
    eas_iteration_data.alpha_.update(
        new_step_length - old_step_length, eas_iteration_data.alpha_inc_, 1.0);
  }

  /*!
   * @brief Compute the matrix M which is the element-wise matrix of the shape functions for the
   * enhanced strains in the parameter space
   *
   * @tparam celltype, eastype
   * @param xi(in) : coordinate in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : enhanced strains shape function matrix in
   * parameter space
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
  evaluate_eas_shape_functions_parameter_space(
      const Core::LinAlg::Matrix<Discret::Elements::num_dim<celltype>, 1>& xi)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
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
      case Solid::Elements::EasType::eastype_h8_9:
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
      case Solid::Elements::EasType::eastype_h8_21:
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
      case Solid::Elements::EasType::eastype_sh8_7:
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
        FOUR_C_THROW("unknown EAS type");
        break;
    }
    return M;
  }

  /*!
   * @brief Map the matrix M in the parameter space to Mtilde in the material configuration and
   * return Mtilde
   *
   * @tparam celltype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant at
   * element centroid
   * @param M(in) : matrix M in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : matrix Mtilde in the material configuration
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
  map_eas_shape_functions_to_material_config(const double detJ,
      const Discret::Elements::CentroidTransformation<celltype>& centroid_transformation,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& M)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
        Mtilde;

    // Mtilde = detJ0/detJ T0^{-T} M
    Mtilde.multiply(centroid_transformation.detJ0_ / detJ, centroid_transformation.T0invT_, M);

    return Mtilde;
  }

  /*!
   * @brief Evaluate the element-wise matrix of the shape functions for the enhanced strains in the
   * parameter space Mtilde. Therefore set up M (in the material configuration) and map M to Mtilde
   * via T0^{-T}.
   *
   * @tparam celltype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant at
   * element centroid
   * @param xi(in) : coordinate in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : matrix Mtilde in the material configuration
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
  evaluate_eas_shape_functions_material_config(const double detJ,
      const Discret::Elements::CentroidTransformation<celltype>& centroid_transformation,
      const Core::LinAlg::Matrix<Discret::Elements::num_dim<celltype>, 1>& xi)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
        M(evaluate_eas_shape_functions_parameter_space<celltype, eastype>(xi));
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
        Mtilde = map_eas_shape_functions_to_material_config<celltype, eastype>(
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
   * @tparam celltype, eastype
   * @param gl_strain(in) : Green-Lagrange strains E^{u}
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return Core::LinAlg::Matrix<num_str, 1>  : enhanced Green-Lagrange strains E^{enh}
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1>
  evaluate_enhanced_assumed_gl_strains(
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1>& gl_strain,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Solid::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1> enhanced_gl_strain(gl_strain);
    enhanced_gl_strain.multiply(1.0, Mtilde, alpha, 1.0);
    return enhanced_gl_strain;
  }

  /*!
   * @brief Evaluate the enhanced assumed Green-Lagrange strains E^{enh}

   * @tparam celltype, eastype
   * @param displacement_based_mapping(in) : displacement-based spatial mapping
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return Core::LinAlg::Matrix<num_str, 1> : Enhanced Green-Lagrange strains E^{enh}
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1>
  evaluate_enhanced_assumed_gl_strains(
      const Discret::Elements::SpatialMaterialMapping<celltype>& displacement_based_mapping,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Solid::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    const Core::LinAlg::Matrix<Discret::Elements::num_dim<celltype>,
        Discret::Elements::num_dim<celltype>>
        displacement_based_cauchygreen =
            Discret::Elements::evaluate_cauchy_green<celltype>(displacement_based_mapping);

    const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1> gl_strain =
        Discret::Elements::evaluate_green_lagrange_strain(displacement_based_cauchygreen);

    return evaluate_enhanced_assumed_gl_strains<celltype, eastype>(gl_strain, Mtilde, alpha);
  }

  /*!
   * @brief Integrate the EAS stiffness matrices
   *
   * @tparam celltype, eastype
   * @param stress(in) : 2. Piola Kirchhoff stress tensor and material tangent
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param Bop(in) : B-operator
   * @param integration_factor(in) : integration factor (Gauss point weight times Jacobi
   * determinant)
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void integrate_eas(const Discret::Elements::Stress<celltype>& stress,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::num_dof_per_ele<celltype>>& Bop,
      const double integration_factor,
      Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    // integrate Kaa: Kaa += (Mtilde^T . cmat . Mtilde) * detJ * w(gp)
    // IMPORTANT: We save this in invKaa_ here since after the loop over all Gauss points, we
    // invert the matrix. At this point, this is still Kaa and NOT invKaa.
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
        cmatM(true);
    cmatM.multiply(stress.cmat_, Mtilde);
    eas_iteration_data.invKaa_.multiply_tn(integration_factor, Mtilde, cmatM, 1.);

    // integrate Kda: Kda += (B^T . cmat . Mtilde) * detJ * w(gp)
    eas_iteration_data.Kda_.multiply_tn(integration_factor, Bop, cmatM, 1.);

    // integrate s: s += (Mtilde^T . S) * detJ * w(gp)
    eas_iteration_data.s_.multiply_tn(integration_factor, Mtilde, stress.pk2_, 1.);
  }

  /*!
   * @brief Add EAS internal force contribution of one Gauss point
   *
   * The EAS internal force contribution is $f_{eas} = - K_{da} K_{aa}^{-1} s$.
   *
   * @tparam celltype : Cell type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param s(in) : enhancement vector s
   * @param force(in/out) : internal force vector where the contribution is added to
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void add_eas_internal_force(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const Core::LinAlg::Matrix<Solid::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& s,
      Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1>& force_vector)
  {
    force_vector.multiply_nn(1.0, minusKdainvKaa, s, 1.0);
  }

  /*!
   * @brief Add EAS stiffness matrix contribution of one Gauss point
   *
   * The EAS stiffness matrix contribution is $- K_{da} K_{aa}^{-1} K_{da}^T$.
   *
   * @tparam celltype : Cell type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param Kda(in) : EAS stiffness matrix part K_{da}
   * @param stiffness_matrix(in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  void add_eas_stiffness_matrix(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Solid::Elements::EasTypeToNumEas<eastype>::num_eas>& Kda,
      Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Discret::Elements::num_dof_per_ele<celltype>>& stiffness_matrix)
  {
    stiffness_matrix.multiply_nt(1.0, minusKdainvKaa, Kda, 1.0);
  }

  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype>
  struct EASKinematics
  {
    static constexpr int num_str = Core::FE::dim<celltype> * (Core::FE::dim<celltype> + 1) / 2;
    static constexpr int num_dof_per_ele = Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>;
    Core::LinAlg::Matrix<num_str, num_dof_per_ele> b_op{};
    Core::LinAlg::Matrix<num_str, 1> enhanced_gl{};
    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        enhanced_deformation_gradient{};
    Core::LinAlg::Matrix<num_str, Solid::Elements::EasTypeToNumEas<eastype>::num_eas> m_tilde{};
  };

  template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
      Inpar::Solid::KinemType kinematic_type>
  EASKinematics<celltype, eastype> evaluate_eas_kinematics(
      const Discret::Elements::ElementNodes<celltype> nodal_coordinates,
      const Discret::Elements::CentroidTransformation<celltype>& centeroid_transformation,
      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1>& xi,
      const Discret::Elements::JacobianMapping<celltype>& jacobian_mapping,
      const Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    EASKinematics<celltype, eastype> eas_kinematics{};

    if constexpr (kinematic_type == Inpar::Solid::KinemType::nonlinearTotLag)
    {
      const Discret::Elements::SpatialMaterialMapping<celltype>
          displacement_based_spatial_material_mapping =
              evaluate_spatial_material_mapping(jacobian_mapping, nodal_coordinates);

      eas_kinematics.b_op = Discret::Elements::evaluate_strain_gradient(
          jacobian_mapping, displacement_based_spatial_material_mapping);

      eas_kinematics.m_tilde = evaluate_eas_shape_functions_material_config<celltype, eastype>(
          jacobian_mapping.determinant_, centeroid_transformation, xi);

      eas_kinematics.enhanced_gl = evaluate_enhanced_assumed_gl_strains<celltype, eastype>(
          displacement_based_spatial_material_mapping, eas_kinematics.m_tilde,
          eas_iteration_data.alpha_);

      eas_kinematics.enhanced_deformation_gradient =
          Discret::Elements::compute_deformation_gradient_from_gl_strains(
              displacement_based_spatial_material_mapping.deformation_gradient_,
              eas_kinematics.enhanced_gl);
    }
    else if constexpr (kinematic_type == Inpar::Solid::KinemType::linear)
    {
      eas_kinematics.b_op = Discret::Elements::evaluate_linear_strain_gradient(jacobian_mapping);

      eas_kinematics.m_tilde = evaluate_eas_shape_functions_material_config<celltype, eastype>(
          jacobian_mapping.determinant_, centeroid_transformation, xi);

      Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1> gl_strain_displacement_based =
          evaluate_linear_gl_strain(nodal_coordinates, eas_kinematics.b_op);

      eas_kinematics.enhanced_gl = evaluate_enhanced_assumed_gl_strains<celltype, eastype>(
          gl_strain_displacement_based, eas_kinematics.m_tilde, eas_iteration_data.alpha_);

      eas_kinematics.enhanced_deformation_gradient =
          Core::LinAlg::identity_matrix<Core::FE::dim<celltype>>();
    }
    else
    {
      FOUR_C_THROW("Unknown kinematic type!");
    }

    return eas_kinematics;
  }
}  // namespace

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::SolidEleCalcEas()
    : stiffness_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_stiffness_matrix<celltype>())),
      mass_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_mass_matrix<celltype>()))
{
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::pack(
    Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, eas_iteration_data_.alpha_inc_);
  add_to_pack(data, eas_iteration_data_.alpha_);
  add_to_pack(data, eas_iteration_data_.s_);
  add_to_pack(data, eas_iteration_data_.invKaa_);
  add_to_pack(data, eas_iteration_data_.Kda_);
  add_to_pack(data, old_step_length_);
};

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, eas_iteration_data_.alpha_inc_);
  extract_from_pack(buffer, eas_iteration_data_.alpha_);
  extract_from_pack(buffer, eas_iteration_data_.s_);
  extract_from_pack(buffer, eas_iteration_data_.invKaa_);
  extract_from_pack(buffer, eas_iteration_data_.Kda_);
  extract_from_pack(buffer, old_step_length_);
};

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
    Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseVector* force_vector,
    Core::LinAlg::SerialDenseMatrix* stiffness_matrix, Core::LinAlg::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff = {};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass = {};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      compare_gauss_integration(mass_matrix_integration_, stiffness_matrix_integration_);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  if (!ele.is_params_interface())
  {
    // Update alpha only in old time integration scheme
    update_alpha<celltype, eastype>(eas_iteration_data_, discretization, lm);
  }

  // clear for integration
  eas_iteration_data_.invKaa_.clear();
  eas_iteration_data_.Kda_.clear();
  eas_iteration_data_.s_.clear();

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  double element_mass = 0.0;
  double element_volume = 0.0;
  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        const Stress<celltype> stress = evaluate_material_stress<celltype>(solid_material,
            kinematic_quantitites.enhanced_deformation_gradient, kinematic_quantitites.enhanced_gl,
            params, gp, ele.id());

        integrate_eas<celltype, eastype>(stress, kinematic_quantitites.m_tilde,
            kinematic_quantitites.b_op, integration_factor, eas_iteration_data_);

        if (force.has_value())
        {
          add_internal_force_vector(kinematic_quantitites.b_op, stress, integration_factor, *force);
        }

        if (stiff.has_value())
        {
          add_elastic_stiffness_matrix(
              kinematic_quantitites.b_op, stress, integration_factor, *stiff);
          add_geometric_stiffness_matrix(
              jacobian_mapping.N_XYZ_, stress, integration_factor, *stiff);
        }

        if (mass.has_value())
        {
          if (equal_integration_mass_stiffness)
          {
            add_mass_matrix(shape_functions, integration_factor, solid_material.density(gp), *mass);
          }
          else
          {
            element_mass += solid_material.density(gp) * integration_factor;
            element_volume += integration_factor;
          }
        }
      });

  // invert Kaa with solver. eas_iteration_data_.invKaa_ then is Kaa^{-1}
  solve_for_inverse_ignoring_errors(eas_iteration_data_.invKaa_);

  // compute the product (- Kda Kaa^{-1}) which is later needed for force and stiffness update
  Core::LinAlg::Matrix<num_dof_per_ele_, FourC::Solid::Elements::EasTypeToNumEas<eastype>::num_eas>
      minusKdainvKaa(true);
  minusKdainvKaa.multiply_nn(-1.0, eas_iteration_data_.Kda_, eas_iteration_data_.invKaa_);

  if (force.has_value())
  {
    add_eas_internal_force<celltype, eastype>(minusKdainvKaa, eas_iteration_data_.s_, *force);
  }

  if (stiff.has_value())
  {
    add_eas_stiffness_matrix<celltype, eastype>(minusKdainvKaa, eas_iteration_data_.Kda_, *stiff);
  }

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    FOUR_C_ASSERT(element_mass > 0, "It looks like the element mass is 0.0");
    Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates, mass_matrix_integration_,
        [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp) {
          add_mass_matrix(
              shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::recover(
    Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  FourC::Solid::Elements::ParamsInterface& params_interface =
      *Teuchos::rcp_dynamic_cast<FourC::Solid::Elements::ParamsInterface>(
          ele.params_interface_ptr());

  const double step_length = params_interface.get_step_length();

  if (params_interface.is_default_step())
  {
    params_interface.sum_into_my_previous_sol_norm(NOX::Nln::StatusTest::quantity_eas,
        FourC::Solid::Elements::EasTypeToNumEas<eastype>::num_eas,
        &eas_iteration_data_.alpha_(0, 0), ele.owner());

    // Update alpha
    update_alpha(eas_iteration_data_, discretization, lm, step_length);
  }
  else
  {
    correct_alpha(eas_iteration_data_, step_length, old_step_length_);
  }

  // store old step length
  old_step_length_ = step_length;

  params_interface.sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas,
      FourC::Solid::Elements::EasTypeToNumEas<eastype>::num_eas,
      &eas_iteration_data_.alpha_inc_(0, 0), &eas_iteration_data_.alpha_(0, 0), step_length,
      ele.owner());
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::update(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);
  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  // No need to update alpha here. Update is called to copy states from t_{n+1} to t_{n} after the
  // time step and output. Hence, there are no more Newton iterations that would require an update
  // of alpha

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        solid_material.update(
            kinematic_quantitites.enhanced_deformation_gradient, gp, params, ele.id());
      });

  solid_material.update();
}


template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::calculate_stress(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material, const StressIO& stressIO,
    const StrainIO& strainIO, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  Core::LinAlg::SerialDenseMatrix stress_data(stiffness_matrix_integration_.num_points(), num_str_);
  Core::LinAlg::SerialDenseMatrix strain_data(stiffness_matrix_integration_.num_points(), num_str_);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        const Stress<celltype> stress = evaluate_material_stress<celltype>(solid_material,
            kinematic_quantitites.enhanced_deformation_gradient, kinematic_quantitites.enhanced_gl,
            params, gp, ele.id());

        assemble_strain_type_to_matrix_row<celltype>(kinematic_quantitites.enhanced_gl,
            kinematic_quantitites.enhanced_deformation_gradient, strainIO.type, strain_data, gp);
        assemble_stress_type_to_matrix_row(kinematic_quantitites.enhanced_deformation_gradient,
            stress, stressIO.type, stress_data, gp);
      });

  serialize(stress_data, serialized_stress_data);
  serialize(strain_data, serialized_strain_data);
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
double
Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::calculate_internal_energy(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        double psi = 0.0;
        solid_material.strain_energy(kinematic_quantitites.enhanced_gl, psi, gp, ele.id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::setup(
    Mat::So3Material& solid_material, const Core::IO::InputParameterContainer& container)
{
  solid_material.setup(stiffness_matrix_integration_.num_points(), container);
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::material_post_setup(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  interpolate_fibers_to_gauss_points_and_add_to_parameter_list<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call post_setup of material
  solid_material.post_setup(params, ele.id());
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::initialize_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  ask_and_add_quantities_to_gauss_point_data_output(
      stiffness_matrix_integration_.num_points(), solid_material, gp_data_output_manager);
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  collect_and_assemble_gauss_point_data_output<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::reset_to_last_converged(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <Core::FE::CellType celltype, Solid::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::for_each_gauss_point(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    const std::function<void(Mat::So3Material&, double, int)>& integrator) const
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  Discret::Elements::for_each_gauss_point(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      { integrator(solid_material, integration_factor, gp); });
}

// template classes
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Solid::Elements::EasType::eastype_h8_21, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Solid::Elements::EasType::eastype_sh8_7, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Solid::Elements::EasType::eastype_h8_21, Inpar::Solid::KinemType::linear>;

static_assert(
    Core::Communication::is_packable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>>,
    "EAS needs to implement the method pack(Core::Communication::PackBuffer&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_unpackable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>>,
    "EAS needs to implement the method unpack(std::size_t, std::vector<char>&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_packable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>>,
    "EAS needs to implement the method pack(Core::Communication::PackBuffer&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_unpackable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>>,
    "EAS needs to implement the method unpack(std::size_t, std::vector<char>&) to be able to store "
    "history data!");
FOUR_C_NAMESPACE_CLOSE
