// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_CALC_SHELL_ANS_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_SHELL_ANS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_eas_helpers.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  template <Core::FE::CellType celltype>
  struct SamplingPointData
  {
    static constexpr int num_ans = 3;
    ShapeFunctionsAndDerivatives<celltype> shape_functions{};
    Core::LinAlg::Matrix<3, 3> reference_jacobian{};
    Core::LinAlg::Matrix<3, 3> current_jacobian{};

    // modified B-operator at local parametric space
    Core::LinAlg::Matrix<num_ans, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
        bop_ans_local{};
  };

  template <Core::FE::CellType celltype>
  struct ShellANSPreparationData
  {
    static constexpr int num_sampling_points = 8;
    std::array<SamplingPointData<celltype>, num_sampling_points> sampling_point_data{};
  };

  namespace Internal
  {
    template <Core::FE::CellType celltype>
    Core::LinAlg::Matrix<num_str<celltype>, num_dof_per_ele<celltype>> evaluate_local_b_operator(
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const ShellANSPreparationData<celltype>& preparation_data)
    {
      Core::LinAlg::Matrix<num_str<celltype>, num_dof_per_ele<celltype>> bop_loc{};
      Core::LinAlg::Matrix<3, 3> current_jacobian(jacobian_mapping.jacobian_);
      current_jacobian.multiply(
          1.0, shape_functions.derivatives_, element_nodes.displacements, 1.0);

      for (int inode = 0; inode < Core::FE::num_nodes<celltype>; ++inode)
      {
        for (int dim = 0; dim < Core::FE::dim<celltype>; ++dim)
        {
          // rr
          bop_loc(0, inode * 3 + dim) =
              shape_functions.derivatives_(0, inode) * current_jacobian(0, dim);

          // ss
          bop_loc(1, inode * 3 + dim) =
              shape_functions.derivatives_(1, inode) * current_jacobian(1, dim);

          // rs
          bop_loc(3, inode * 3 + dim) =
              shape_functions.derivatives_(0, inode) * current_jacobian(1, dim) +
              shape_functions.derivatives_(1, inode) * current_jacobian(0, dim);

          // interpolate along (r x s) of bop_ans_local(tt) (sampling points 4, 5, 6, 7)
          bop_loc(2, inode * 3 + dim) =
              0.25 * (1 - xi(0)) * (1 - xi(1)) *
                  preparation_data.sampling_point_data[4].bop_ans_local(0, inode * 3 + dim) +
              0.25 * (1 + xi(0)) * (1 - xi(1)) *
                  preparation_data.sampling_point_data[5].bop_ans_local(0, inode * 3 + dim) +
              0.25 * (1 + xi(0)) * (1 + xi(1)) *
                  preparation_data.sampling_point_data[6].bop_ans_local(0, inode * 3 + dim) +
              0.25 * (1 - xi(0)) * (1 + xi(1)) *
                  preparation_data.sampling_point_data[7].bop_ans_local(0, inode * 3 + dim);

          // interpolate along (r x s) of bop_ans_local(st) (sampling points 1, 3)
          bop_loc(4, inode * 3 + dim) =
              0.5 * (1 + xi(0)) *
                  preparation_data.sampling_point_data[1].bop_ans_local(1, inode * 3 + dim) +
              0.5 * (1 - xi(0)) *
                  preparation_data.sampling_point_data[3].bop_ans_local(1, inode * 3 + dim);

          // interpolate along (r x s) of bop_ans_local(rt) (sampling points 0, 2)
          bop_loc(5, inode * 3 + dim) =
              0.5 * (1 - xi(1)) *
                  preparation_data.sampling_point_data[0].bop_ans_local(2, inode * 3 + dim) +
              0.5 * (1 + xi(1)) *
                  preparation_data.sampling_point_data[2].bop_ans_local(2, inode * 3 + dim);
          ;
        }
      }

      return bop_loc;
    }

    template <Core::FE::CellType celltype>
    Core::LinAlg::Matrix<num_str<celltype>, 1> evaluate_local_glstrain(
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const ShellANSPreparationData<celltype>& preparation_data)
    {
      Core::LinAlg::Matrix<3, 3> current_jacobian(jacobian_mapping.jacobian_);
      current_jacobian.multiply(
          1.0, shape_functions.derivatives_, element_nodes.displacements, 1.0);

      Core::LinAlg::Matrix<num_str<celltype>, 1> glstrain;
      // evaluate glstrains in local(parameter) coords
      // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
      glstrain(0) =
          0.5 * (+(current_jacobian(0, 0) * current_jacobian(0, 0) +
                     current_jacobian(0, 1) * current_jacobian(0, 1) +
                     current_jacobian(0, 2) * current_jacobian(0, 2)) -
                    (jacobian_mapping.jacobian_(0, 0) * jacobian_mapping.jacobian_(0, 0) +
                        jacobian_mapping.jacobian_(0, 1) * jacobian_mapping.jacobian_(0, 1) +
                        jacobian_mapping.jacobian_(0, 2) * jacobian_mapping.jacobian_(0, 2)));
      // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
      glstrain(1) =
          0.5 * (+(current_jacobian(1, 0) * current_jacobian(1, 0) +
                     current_jacobian(1, 1) * current_jacobian(1, 1) +
                     current_jacobian(1, 2) * current_jacobian(1, 2)) -
                    (jacobian_mapping.jacobian_(1, 0) * jacobian_mapping.jacobian_(1, 0) +
                        jacobian_mapping.jacobian_(1, 1) * jacobian_mapping.jacobian_(1, 1) +
                        jacobian_mapping.jacobian_(1, 2) * jacobian_mapping.jacobian_(1, 2)));
      // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
      glstrain(3) = (+(current_jacobian(0, 0) * current_jacobian(1, 0) +
                         current_jacobian(0, 1) * current_jacobian(1, 1) +
                         current_jacobian(0, 2) * current_jacobian(1, 2)) -
                     (jacobian_mapping.jacobian_(0, 0) * jacobian_mapping.jacobian_(1, 0) +
                         jacobian_mapping.jacobian_(0, 1) * jacobian_mapping.jacobian_(1, 1) +
                         jacobian_mapping.jacobian_(0, 2) * jacobian_mapping.jacobian_(1, 2)));

      // ANS modification of strains ************************************** ANS
      double dxdt_A = 0.0;
      double dXdt_A = 0.0;
      double dydt_B = 0.0;
      double dYdt_B = 0.0;
      double dxdt_C = 0.0;
      double dXdt_C = 0.0;
      double dydt_D = 0.0;
      double dYdt_D = 0.0;

      double dzdt_E = 0.0;
      double dZdt_E = 0.0;
      double dzdt_F = 0.0;
      double dZdt_F = 0.0;
      double dzdt_G = 0.0;
      double dZdt_G = 0.0;
      double dzdt_H = 0.0;
      double dZdt_H = 0.0;

      // vector product of rows of jacobians at corresponding sampling point
      for (int dim = 0; dim < Core::FE::dim<celltype>; ++dim)
      {
        dxdt_A += preparation_data.sampling_point_data[0].current_jacobian(0, dim) *
                  preparation_data.sampling_point_data[0].current_jacobian(2, dim);  // g_13^A
        dXdt_A += preparation_data.sampling_point_data[0].reference_jacobian(0, dim) *
                  preparation_data.sampling_point_data[0].reference_jacobian(2, dim);  // G_13^A
        dydt_B += preparation_data.sampling_point_data[1].current_jacobian(1, dim) *
                  preparation_data.sampling_point_data[1].current_jacobian(2, dim);  // g_23^B
        dYdt_B += preparation_data.sampling_point_data[1].reference_jacobian(1, dim) *
                  preparation_data.sampling_point_data[1].reference_jacobian(2, dim);  // G_23^B
        dxdt_C += preparation_data.sampling_point_data[2].current_jacobian(0, dim) *
                  preparation_data.sampling_point_data[2].current_jacobian(2, dim);  // g_13^C
        dXdt_C += preparation_data.sampling_point_data[2].reference_jacobian(0, dim) *
                  preparation_data.sampling_point_data[2].reference_jacobian(2, dim);  // G_13^C
        dydt_D += preparation_data.sampling_point_data[3].current_jacobian(1, dim) *
                  preparation_data.sampling_point_data[3].current_jacobian(2, dim);  // g_23^D
        dYdt_D += preparation_data.sampling_point_data[3].reference_jacobian(1, dim) *
                  preparation_data.sampling_point_data[3].reference_jacobian(2, dim);  // G_23^D

        dzdt_E += preparation_data.sampling_point_data[4].current_jacobian(2, dim) *
                  preparation_data.sampling_point_data[4].current_jacobian(2, dim);
        dZdt_E += preparation_data.sampling_point_data[4].reference_jacobian(2, dim) *
                  preparation_data.sampling_point_data[4].reference_jacobian(2, dim);
        dzdt_F += preparation_data.sampling_point_data[5].current_jacobian(2, dim) *
                  preparation_data.sampling_point_data[5].current_jacobian(2, dim);
        dZdt_F += preparation_data.sampling_point_data[5].reference_jacobian(2, dim) *
                  preparation_data.sampling_point_data[5].reference_jacobian(2, dim);
        dzdt_G += preparation_data.sampling_point_data[6].current_jacobian(2, dim) *
                  preparation_data.sampling_point_data[6].current_jacobian(2, dim);
        dZdt_G += preparation_data.sampling_point_data[6].reference_jacobian(2, dim) *
                  preparation_data.sampling_point_data[6].reference_jacobian(2, dim);
        dzdt_H += preparation_data.sampling_point_data[7].current_jacobian(2, dim) *
                  preparation_data.sampling_point_data[7].current_jacobian(2, dim);
        dZdt_H += preparation_data.sampling_point_data[7].reference_jacobian(2, dim) *
                  preparation_data.sampling_point_data[7].reference_jacobian(2, dim);
      }

      // E33: remedy of curvature thickness locking
      // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
      glstrain(2) = 0.5 * (0.25 * (1 - xi(0)) * (1 - xi(1)) * (dzdt_E - dZdt_E) +
                              0.25 * (1 + xi(0)) * (1 - xi(1)) * (dzdt_F - dZdt_F) +
                              0.25 * (1 + xi(0)) * (1 + xi(1)) * (dzdt_G - dZdt_G) +
                              0.25 * (1 - xi(0)) * (1 + xi(1)) * (dzdt_H - dZdt_H));
      // E23: remedy of transverse shear locking
      // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
      glstrain(4) = 0.5 * (1 + xi(0)) * (dydt_B - dYdt_B) + 0.5 * (1 - xi(0)) * (dydt_D - dYdt_D);
      // E13: remedy of transverse shear locking
      // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
      glstrain(5) = 0.5 * (1 - xi(1)) * (dxdt_A - dXdt_A) + 0.5 * (1 + xi(1)) * (dxdt_C - dXdt_C);

      return glstrain;
    }

    template <Core::FE::CellType celltype>
    void add_ans_geometric_stiffness(const Core::LinAlg::Matrix<3, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const Stress<celltype>& stress, const double integration_factor,
        const Core::LinAlg::Matrix<Internal::num_str<celltype>, Internal::num_str<celltype>>& TinvT,
        const ShellANSPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Internal::num_dim<celltype> * Internal::num_nodes<celltype>,
            Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& stiffness_matrix)
    {
      for (int inod = 0; inod < Core::FE::num_nodes<celltype>; ++inod)
      {
        for (int jnod = 0; jnod < Core::FE::num_nodes<celltype>; ++jnod)
        {
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> G_ij;
          G_ij(0) = shape_functions.derivatives_(0, inod) *
                    shape_functions.derivatives_(0, jnod);  // rr-dir
          G_ij(1) = shape_functions.derivatives_(1, inod) *
                    shape_functions.derivatives_(1, jnod);  // ss-dir
          G_ij(3) = shape_functions.derivatives_(0, inod) * shape_functions.derivatives_(1, jnod) +
                    shape_functions.derivatives_(1, inod) *
                        shape_functions.derivatives_(0, jnod);  // rs-dir


          // ANS modification in tt-dir
          G_ij(2) =
              0.25 * (1 - xi(0)) * (1 - xi(1)) *
                  preparation_data.sampling_point_data[4].shape_functions.derivatives_(2, inod) *
                  preparation_data.sampling_point_data[4].shape_functions.derivatives_(2, jnod) +
              0.25 * (1 + xi(0)) * (1 - xi(1)) *
                  preparation_data.sampling_point_data[5].shape_functions.derivatives_(2, inod) *
                  preparation_data.sampling_point_data[5].shape_functions.derivatives_(2, jnod) +
              0.25 * (1 + xi(0)) * (1 + xi(1)) *
                  preparation_data.sampling_point_data[6].shape_functions.derivatives_(2, inod) *
                  preparation_data.sampling_point_data[6].shape_functions.derivatives_(2, jnod) +
              0.25 * (1 - xi(0)) * (1 + xi(1)) *
                  preparation_data.sampling_point_data[7].shape_functions.derivatives_(2, inod) *
                  preparation_data.sampling_point_data[7].shape_functions.derivatives_(2, jnod);
          // ANS modification in st-dir
          G_ij(4) =
              0.5 *
              ((1 + xi(0)) *
                      (preparation_data.sampling_point_data[1].shape_functions.derivatives_(
                           1, inod) *
                              preparation_data.sampling_point_data[1].shape_functions.derivatives_(
                                  2, jnod) +
                          preparation_data.sampling_point_data[1].shape_functions.derivatives_(
                              2, inod) *
                              preparation_data.sampling_point_data[1].shape_functions.derivatives_(
                                  1, jnod)) +
                  (1 - xi(0)) *
                      (preparation_data.sampling_point_data[3].shape_functions.derivatives_(
                           1, inod) *
                              preparation_data.sampling_point_data[3].shape_functions.derivatives_(
                                  2, jnod) +
                          preparation_data.sampling_point_data[3].shape_functions.derivatives_(
                              2, inod) *
                              preparation_data.sampling_point_data[3].shape_functions.derivatives_(
                                  1, jnod)));
          // ANS modification in rt-dir
          G_ij(5) =
              0.5 *
              ((1 - xi(1)) *
                      (preparation_data.sampling_point_data[0].shape_functions.derivatives_(
                           0, inod) *
                              preparation_data.sampling_point_data[0].shape_functions.derivatives_(
                                  2, jnod) +
                          preparation_data.sampling_point_data[0].shape_functions.derivatives_(
                              2, inod) *
                              preparation_data.sampling_point_data[0].shape_functions.derivatives_(
                                  0, jnod)) +
                  (1 + xi(1)) *
                      (preparation_data.sampling_point_data[2].shape_functions.derivatives_(
                           0, inod) *
                              preparation_data.sampling_point_data[2].shape_functions.derivatives_(
                                  2, jnod) +
                          preparation_data.sampling_point_data[2].shape_functions.derivatives_(
                              2, inod) *
                              preparation_data.sampling_point_data[2].shape_functions.derivatives_(
                                  0, jnod)));


          // transformation of local(parameter) space 'back' to global(material) space
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> G_ij_glob;
          G_ij_glob.multiply(TinvT, G_ij);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          const double Gij = integration_factor * stress.pk2_.dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          stiffness_matrix(
              Core::FE::dim<celltype> * inod + 0, Core::FE::dim<celltype> * jnod + 0) += Gij;
          stiffness_matrix(
              Core::FE::dim<celltype> * inod + 1, Core::FE::dim<celltype> * jnod + 1) += Gij;
          stiffness_matrix(
              Core::FE::dim<celltype> * inod + 2, Core::FE::dim<celltype> * jnod + 2) += Gij;
        }
      }
    }
  }  // namespace Internal


  template <Core::FE::CellType celltype>
  struct ShellANSLinearizationContainer
  {
    Core::LinAlg::Matrix<Internal::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
        Bop_{};

    Core::LinAlg::Matrix<Internal::num_str<celltype>, Internal::num_str<celltype>> TinvT{};
  };

  /*!
   * @brief A solid-shell element formulation with ANS
   *
   * @note The ordering of the nodes within the element determine the direction of the thickness of
   * the shell. The thickness-direction is assumed to be the t-direction of the element's parameter
   * space.
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct ShellANSFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = true;
    static constexpr bool has_condensed_contribution = false;

    using LinearizationContainer = ShellANSLinearizationContainer<celltype>;
    using PreparationData = ShellANSPreparationData<celltype>;

    static ShellANSPreparationData<celltype> prepare(
        const Core::Elements::Element& ele, const ElementNodes<celltype>& nodal_coordinates)
    {
      constexpr std::array<std::array<double, 3>, 8> ans_sampling_points{
          {{{0.0, -1.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}},
              {{-1.0, -1.0, 0.0}}, {{1.0, -1.0, 0.0}}, {{1.0, 1.0, 0.0}}, {{-1.0, 1.0, 0.0}}}};

      ShellANSPreparationData<celltype> shell_ans_data{};
      std::size_t i = 0;
      for (const std::array<double, 3>& sampling_point : ans_sampling_points)
      {
        Core::LinAlg::Matrix<3, 1> xi_view(sampling_point.data(), true);
        // evaluate derivative of the shape functions
        shell_ans_data.sampling_point_data[i].shape_functions =
            evaluate_shape_functions_and_derivs(xi_view, nodal_coordinates);

        const ShapeFunctionsAndDerivatives<celltype>& shape_functions =
            shell_ans_data.sampling_point_data[i].shape_functions;

        shell_ans_data.sampling_point_data[i].reference_jacobian.multiply(
            shape_functions.derivatives_, nodal_coordinates.reference_coordinates);

        shell_ans_data.sampling_point_data[i].current_jacobian =
            shell_ans_data.sampling_point_data[i].reference_jacobian;
        shell_ans_data.sampling_point_data[i].current_jacobian.multiply(
            1.0, shape_functions.derivatives_, nodal_coordinates.displacements, 1.0);


        const Core::LinAlg::Matrix<3, 3>& current_jacobian =
            shell_ans_data.sampling_point_data[i].current_jacobian;
        // build local ans b-operator
        for (int inode = 0; inode < Core::FE::num_nodes<celltype>; ++inode)
        {
          for (int dim = 0; dim < Core::FE::dim<celltype>; ++dim)
          {
            // modified b-operator in tt
            shell_ans_data.sampling_point_data[i].bop_ans_local(0, inode * 3 + dim) =
                shape_functions.derivatives_(2, inode) * current_jacobian(2, dim);

            // modified b-operator in st
            shell_ans_data.sampling_point_data[i].bop_ans_local(1, inode * 3 + dim) =
                shape_functions.derivatives_(1, inode) * current_jacobian(2, dim) +
                shape_functions.derivatives_(2, inode) * current_jacobian(1, dim);

            // modified b-operator in rt
            shell_ans_data.sampling_point_data[i].bop_ans_local(2, inode * 3 + dim) =
                shape_functions.derivatives_(0, inode) * current_jacobian(2, dim) +
                shape_functions.derivatives_(2, inode) * current_jacobian(0, dim);
          }
        }


        i += 1;
      }

      return shell_ans_data;
    }

    template <typename Evaluator>
    static inline auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const ShellANSPreparationData<celltype>& preparation_data, Evaluator evaluator)
    {
      // evaluate local b-op
      Core::LinAlg::Matrix<num_str<celltype>, num_dof_per_ele<celltype>> bop_local =
          Internal::evaluate_local_b_operator(
              element_nodes, xi, shape_functions, jacobian_mapping, preparation_data);

      const ShellANSLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            ShellANSLinearizationContainer<celltype> linearization{};
            linearization.TinvT = evaluate_voigt_transformation_matrix(jacobian_mapping);
            linearization.Bop_.multiply(linearization.TinvT, bop_local);
            return linearization;
          });

      Core::LinAlg::Matrix<6, 1> gl_strain_local = Internal::evaluate_local_glstrain(
          element_nodes, xi, shape_functions, jacobian_mapping, preparation_data);

      Core::LinAlg::Matrix<6, 1> gl_strain;
      gl_strain.multiply(linearization.TinvT, gl_strain_local);


      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          evaluate_spatial_material_mapping(jacobian_mapping, element_nodes);

      Core::LinAlg::Matrix<3, 3> consistent_defgrd = compute_deformation_gradient_from_gl_strains(
          spatial_material_mapping.deformation_gradient_, gl_strain);

      return evaluator(consistent_defgrd, gl_strain, linearization);
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>&
            deformation_gradient,
        const ShellANSPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This derivative of the deformation gradient w.r.t. the displacements is not "
          "implemented");
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>&
            deformation_gradient,
        const ShellANSPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW("This derivative of the deformation gradient w.r.t. xi is not implemented");
    }

    static inline Core::LinAlg::Matrix<9,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>&
            deformation_gradient,
        const ShellANSPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This second derivative of the deformation gradient w.r.t. the displacements and xi "
          "is "
          "not implemented");
    }

    static Core::LinAlg::Matrix<Internal::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    get_linear_b_operator(const ShellANSLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop_;
    }

    static void add_internal_force_vector(
        const ShellANSLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const ShellANSPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::Elements::add_internal_force_vector(
          linearization.Bop_, stress, integration_factor, force_vector);
    }

    static void add_stiffness_matrix(const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const ShellANSLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const ShellANSPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
            Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::Elements::add_elastic_stiffness_matrix(
          linearization.Bop_, stress, integration_factor, stiffness_matrix);


      Internal::add_ans_geometric_stiffness(xi, shape_functions, stress, integration_factor,
          linearization.TinvT, preparation_data, stiffness_matrix);
    }
  };

  template <Core::FE::CellType celltype>
  using ANSSolidShellIntegrator = SolidEleCalc<celltype, ShellANSFormulation<celltype>>;


}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE
#endif
