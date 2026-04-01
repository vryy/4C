// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_solid_3D_ele_calc_lib_plane.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_element_dof_matrix.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_fiber_node_utils.hpp"
#include "4C_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_nurbs_discretization_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_tensor_symmetric_einstein.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_local_newton.hpp"

#include <Teuchos_ParameterList.hpp>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

namespace
{
  Mat::EvaluationContext<3> translate_context(const Mat::EvaluationContext<2>& context,
      Core::LinAlg::Tensor<double, 3>& xi, Core::LinAlg::Tensor<double, 3>& ref_coords)
  {
    if (context.xi)
    {
      xi(0) = (*context.xi)(0);
      xi(1) = (*context.xi)(1);
      xi(2) = 0.0;
    }
    if (context.ref_coords)
    {
      ref_coords(0) = (*context.ref_coords)(0);
      ref_coords(1) = (*context.ref_coords)(1);
      ref_coords(2) = 0.0;
    }
    return {
        .total_time = context.total_time,
        .time_step_size = context.time_step_size,
        .xi = context.xi ? &xi : nullptr,
        .ref_coords = context.ref_coords ? &ref_coords : nullptr,
    };
  }

  Core::LinAlg::SymmetricTensor<double, 2, 2> extract_2d_part(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& tensor_3d)
  {
    return Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{
        {{tensor_3d(0, 0), tensor_3d(0, 1)}, {tensor_3d(1, 0), tensor_3d(1, 1)}}});
  }

  Core::LinAlg::SymmetricTensor<double, 2, 2, 2, 2> extract_2d_part(
      const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& tensor_3d)
  {
    Core::LinAlg::SymmetricTensor<double, 2, 2, 2, 2> result{};

    for (int i = 0; i < 2; ++i)
      for (int j = i; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          for (int l = k; l < 2; ++l) result(i, j, k, l) = tensor_3d(i, j, k, l);

    return result;
  }

  Core::LinAlg::Tensor<double, 3, 3> make_3d_tensor(
      const Core::LinAlg::Tensor<double, 2, 2>& tensor_2d, double diagonal_value = 0.0)
  {
    return Core::LinAlg::Tensor<double, 3, 3>{{
        {tensor_2d(0, 0), tensor_2d(0, 1), 0},
        {tensor_2d(1, 0), tensor_2d(1, 1), 0},
        {0, 0, diagonal_value},
    }};
  }

  Core::LinAlg::SymmetricTensor<double, 3, 3> make_3d_tensor(
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& tensor_2d, double diagonal_value = 0.0)
  {
    return Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
        {tensor_2d(0, 0), tensor_2d(0, 1), 0},
        {tensor_2d(1, 0), tensor_2d(1, 1), 0},
        {0, 0, diagonal_value},
    }});
  }


  template <Core::FE::CellType celltype>
  Discret::Elements::Stress<celltype> evaluate_material_plane_strain(Mat::So3Material& material,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID)
  {
    // make 3D tensors out of the 2D ones
    const Core::LinAlg::Tensor<double, 3, 3> defgrd_3d = make_3d_tensor(defgrd, 1.0);
    const Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d = make_3d_tensor(gl_strain, 0.0);

    // making 3D context out of 2D context
    Core::LinAlg::Tensor<double, 3> xi, gp_ref_coord;
    const Mat::EvaluationContext<3> context_3d = translate_context(context, xi, gp_ref_coord);

    Core::LinAlg::SymmetricTensor<double, 3, 3> pk2_3d{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_3d{};
    material.evaluate(&defgrd_3d, gl_strain_3d, params, context_3d, pk2_3d, cmat_3d, gp, eleGID);

    return {.pk2_ = extract_2d_part(pk2_3d), .cmat_ = extract_2d_part(cmat_3d)};
  }

  template <Core::FE::CellType celltype>
  struct PlaneStressQuantities
  {
    Discret::Elements::Stress<celltype> stress;
    Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d;
  };

  template <Core::FE::CellType celltype>
  PlaneStressQuantities<celltype> evaluate_material_plane_stress(Mat::So3Material& material,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID)
  {
    // making 3D context out of 2D context
    Core::LinAlg::Tensor<double, 3> xi, gp_ref_coord;
    const Mat::EvaluationContext<3> context_3d = translate_context(context, xi, gp_ref_coord);

    Core::LinAlg::SymmetricTensor<double, 3, 3> pk2_3d{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_3d{};

    auto residuum_and_jacobian = [&](const Core::LinAlg::Tensor<double, 3>& gl_strain_components_3)
        -> std::tuple<Core::LinAlg::Tensor<double, 3>, Core::LinAlg::Tensor<double, 3, 3>>
    {
      Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d =
          Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
              {gl_strain(0, 0), gl_strain(0, 1), gl_strain_components_3(0)},
              {gl_strain(1, 0), gl_strain(1, 1), gl_strain_components_3(1)},
              {gl_strain_components_3(0), gl_strain_components_3(1), gl_strain_components_3(2)},
          }});

      pk2_3d = {};
      cmat_3d = {};
      material.evaluate(nullptr, gl_strain_3d, params, context_3d, pk2_3d, cmat_3d, gp, eleGID);

      Core::LinAlg::Tensor<double, 3> residuum = {{pk2_3d(0, 2), pk2_3d(1, 2), pk2_3d(2, 2)}};
      Core::LinAlg::Tensor<double, 3, 3> jacobian = {{
          {cmat_3d(0, 2, 0, 2), cmat_3d(0, 2, 1, 2), cmat_3d(0, 2, 2, 2)},
          {cmat_3d(1, 2, 0, 2), cmat_3d(1, 2, 1, 2), cmat_3d(1, 2, 2, 2)},
          {cmat_3d(2, 2, 0, 2), cmat_3d(2, 2, 1, 2), cmat_3d(2, 2, 2, 2)},
      }};
      return {residuum, jacobian};
    };

    auto [gl_strain_components_3, jacobian] = Core::Utils::solve_local_newton_and_return_jacobian(
        residuum_and_jacobian, Core::LinAlg::Tensor<double, 3>{{0.0, 0.0, 0.0}}, 1e-9);

    Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
            {gl_strain(0, 0), gl_strain(0, 1), gl_strain_components_3(0)},
            {gl_strain(1, 0), gl_strain(1, 1), gl_strain_components_3(1)},
            {gl_strain_components_3(0), gl_strain_components_3(1), gl_strain_components_3(2)},
        }});

    Core::LinAlg::Tensor<double, 2, 2, 3> cfr{};
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 3; ++k) cfr(i, j, k) = cmat_3d(i, j, k, 2);

    Core::LinAlg::SymmetricTensor<double, 2, 2, 2, 2> plane_stress_linearization =
        Core::LinAlg::assume_symmetry(
            cfr * Core::LinAlg::inv(jacobian) * Core::LinAlg::einsum<"bca">(cfr));

    return {.stress = {.pk2_ = extract_2d_part(pk2_3d),
                .cmat_ = extract_2d_part(cmat_3d) - plane_stress_linearization},
        .gl_strain_3d = gl_strain_3d};
  }
}  // namespace

template <Core::FE::CellType celltype>
  requires(Core::FE::dim<celltype> == 2)
Discret::Elements::Stress<celltype> Discret::Elements::evaluate_material_stress(
    Mat::So3Material& material, const ElementProperties<celltype>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID)
{
  // Note: This is a legacy evaluation for 2D materials. 4C currently does not have pure 2D
  // materials, so we use normal 3D materials instead for now.
  Core::Mat::Material* base_material = &material;
  // if we have a poro-material, we are only interested in the solid contribution
  if (material.material_type() == Core::Materials::m_structporo or
      material.material_type() == Core::Materials::m_structpororeaction or
      material.material_type() == Core::Materials::m_structpororeactionECM)
  {
    base_material = dynamic_cast<Mat::StructPoro&>(material).get_material().get();
  }


  switch (base_material->material_type())
  {
    case Core::Materials::m_stvenant:
    {
      // The St. Venant-Kirchhoff material can be simplified a lot for both plane stress and plane
      // strain assumptions
      const Mat::StVenantKirchhoff& actmat =
          dynamic_cast<const Mat::StVenantKirchhoff&>(*base_material);
      const double ym = actmat.youngs();
      const double nu = actmat.poisson_ratio();
      const double lamb = element_properties.plane_assumption == PlaneAssumption::plane_stress
                              ? ym * nu / (1 - nu * nu)
                              : ym * nu / ((1 + nu) * (1 - 2 * nu));
      const double mue = ym / (2 * (1 + nu));

      Stress<celltype> stress{};
      stress.pk2_ = lamb * Core::LinAlg::trace(gl_strain) *
                        Core::LinAlg::TensorGenerators::identity<double, 2, 2> +
                    2 * mue * gl_strain;
      stress.cmat_ =
          2 * mue * Core::LinAlg::TensorGenerators::symmetric_identity<double, 2, 2, 2, 2> +
          lamb * Core::LinAlg::dyadic(Core::LinAlg::TensorGenerators::identity<double, 2, 2>,
                     Core::LinAlg::TensorGenerators::identity<double, 2, 2>);
      return stress;
    }
    default:
    {
      switch (element_properties.plane_assumption)
      {
        case PlaneAssumption::plane_stress:
          return evaluate_material_plane_stress<celltype>(
              material, gl_strain, params, context, gp, eleGID)
              .stress;
        case PlaneAssumption::plane_strain:
          return evaluate_material_plane_strain<celltype>(
              material, defgrd, gl_strain, params, context, gp, eleGID);
        default:
          FOUR_C_THROW("Unknown plane assumption for 2D solid element.");
      }
    }
  }
}

template <Core::FE::CellType celltype>
  requires(Core::FE::dim<celltype> == 2)
void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<celltype>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID)
{
  Core::LinAlg::Tensor<double, 3, 3> F_consistent{};

  switch (element_properties.plane_assumption)
  {
    case PlaneAssumption::plane_stress:
    {
      // compute 2d strains
      Core::LinAlg::SymmetricTensor<double, 2, 2> gl_strain_2d =
          0.5 * (Core::LinAlg::assume_symmetry(Core::LinAlg::transpose(defgrd) * defgrd) -
                    Core::LinAlg::TensorGenerators::identity<double, 2, 2>);

      // solve local system to obtain the consistent 3d strains for plane stress
      Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d =
          evaluate_material_plane_stress<celltype>(
              material, gl_strain_2d, params, context, gp, eleGID)
              .gl_strain_3d;

      // compute consistent deformation gradient
      F_consistent = Discret::Elements::compute_deformation_gradient_from_gl_strains(
          make_3d_tensor(defgrd, 1.0), gl_strain_3d);
      break;
    }
    case PlaneAssumption::plane_strain:
    {
      F_consistent = make_3d_tensor(defgrd, 1.0);
      break;
    }
    default:
      FOUR_C_THROW("Unknown plane assumption for 2D solid element.");
  }

  // making 3D context out of 2D context
  Core::LinAlg::Tensor<double, 3> xi, gp_ref_coord;
  const Mat::EvaluationContext<3> context_3d = translate_context(context, xi, gp_ref_coord);

  // call consistent update
  material.update(F_consistent, gp, params, context_3d, eleGID);
}

template <Core::FE::CellType celltype>
  requires(Core::FE::dim<celltype> == 2)
double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<celltype>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID)
{
  Core::LinAlg::SymmetricTensor<double, 3, 3> gl_strain_3d;
  switch (element_properties.plane_assumption)
  {
    case PlaneAssumption::plane_stress:
    {
      // solve local system to obtain the consistent 3d strains for plane stress
      gl_strain_3d =
          evaluate_material_plane_stress<celltype>(material, gl_strain, params, context, gp, eleGID)
              .gl_strain_3d;
      break;
    }
    case PlaneAssumption::plane_strain:
    {
      gl_strain_3d = make_3d_tensor(gl_strain, 0.0);
      break;
    }
    default:
      FOUR_C_THROW("Unknown plane assumption for 2D solid element.");
  }

  // making 3D context out of 2D context
  Core::LinAlg::Tensor<double, 3> xi, gp_ref_coord;
  const Mat::EvaluationContext<3> context_3d = translate_context(context, xi, gp_ref_coord);

  return material.strain_energy(gl_strain_3d, context_3d, gp, eleGID);
}

template Discret::Elements::Stress<Core::FE::CellType::quad4>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::quad4>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad4>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template Discret::Elements::Stress<Core::FE::CellType::quad8>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::quad8>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad8>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template Discret::Elements::Stress<Core::FE::CellType::quad9>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::quad9>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad9>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template Discret::Elements::Stress<Core::FE::CellType::nurbs9>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::nurbs9>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::nurbs9>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template Discret::Elements::Stress<Core::FE::CellType::tri3>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::tri3>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri3>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template Discret::Elements::Stress<Core::FE::CellType::tri6>
Discret::Elements::evaluate_material_stress<Core::FE::CellType::tri6>(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri6>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad4>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad8>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad9>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::nurbs9>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri3>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template void Discret::Elements::update_material(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri6>& element_properties,
    const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad4>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad8>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::quad9>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::nurbs9>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri3>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);
template double Discret::Elements::evaluate_material_strain_energy(Mat::So3Material& material,
    const ElementProperties<Core::FE::CellType::tri6>& element_properties,
    const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
    const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

FOUR_C_NAMESPACE_CLOSE