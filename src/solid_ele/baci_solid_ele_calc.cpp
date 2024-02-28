/*! \file

\brief Implementation of routines for calculation of solid element with templated element
formulation

\level 1
*/

#include "baci_solid_ele_calc.hpp"

#include "baci_discretization_fem_general_cell_type.hpp"
#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_mat_so3_material.hpp"
#include "baci_solid_ele_calc_displacement_based.hpp"
#include "baci_solid_ele_calc_displacement_based_linear_kinematics.hpp"
#include "baci_solid_ele_calc_fbar.hpp"
#include "baci_solid_ele_calc_lib.hpp"
#include "baci_solid_ele_calc_lib_integration.hpp"
#include "baci_solid_ele_calc_lib_io.hpp"
#include "baci_solid_ele_calc_lib_nitsche.hpp"
#include "baci_solid_ele_interface_serializable.hpp"

#include <Teuchos_ParameterList.hpp>

#include <optional>

BACI_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType celltype>
  DRT::ELEMENTS::CauchyNDirAndLinearization<3> EvaluateCauchyNDirAndDerivatives(
      MAT::So3Material& mat,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
          deformation_gradient,
      const DRT::ELEMENTS::SolidFormulationLinearization<celltype>& solid_formulation_linearization,
      const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir, int eleGID)
  {
    CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF(true);
    CORE::LINALG::Matrix<9, 9> d2_cauchyndir_dF2(true);
    CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>> d2_cauchyndir_dF_dn(true);
    CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>> d2_cauchyndir_dF_ddir(true);

    DRT::ELEMENTS::CauchyNDirAndLinearization<3> cauchy_n_dir_with_linearization{};

    mat.EvaluateCauchyNDirAndDerivatives(deformation_gradient, n, dir,
        cauchy_n_dir_with_linearization.cauchy_n_dir,
        &cauchy_n_dir_with_linearization.d_cauchyndir_dn,
        &cauchy_n_dir_with_linearization.d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2,
        &d2_cauchyndir_dF_dn, &d2_cauchyndir_dF_ddir, -1, eleGID, nullptr, nullptr, nullptr,
        nullptr);

    // evaluate first derivative w.r.t. displacements
    DRT::ELEMENTS::EvaluateDCauchyNDirDDisplacements<celltype>(
        solid_formulation_linearization.d_F_dd, d_cauchyndir_dF,
        cauchy_n_dir_with_linearization.d_cauchyndir_dd);


    // evaluate second derivative w.r.t. displacements, normal
    DRT::ELEMENTS::EvaluateD2CauchyNDirDDisplacementsDNormal<celltype>(
        solid_formulation_linearization.d_F_dd, d2_cauchyndir_dF_dn,
        cauchy_n_dir_with_linearization.d2_cauchyndir_dd_dn);

    // evaluate second derivative w.r.t. displacements, direction
    DRT::ELEMENTS::EvaluateD2CauchyNDirDDisplacementsDDir<celltype>(
        solid_formulation_linearization.d_F_dd, d2_cauchyndir_dF_ddir,
        cauchy_n_dir_with_linearization.d2_cauchyndir_dd_ddir);

    // evaluate second derivative w.r.t. displacements, displacements
    DRT::ELEMENTS::EvaluateD2CauchyNDirDDisplacements2<celltype>(
        solid_formulation_linearization.d_F_dd, d2_cauchyndir_dF2,
        cauchy_n_dir_with_linearization.d2_cauchyndir_dd2);

    // evaluate first derivative w.r.t. xi
    DRT::ELEMENTS::EvaluateDCauchyNDirDXi<celltype>(solid_formulation_linearization.d_F_dxi,
        d_cauchyndir_dF, cauchy_n_dir_with_linearization.d_cauchyndir_dxi);

    // evaluate second derivative w.r.t. displacements, xi
    DRT::ELEMENTS::EvaluateD2CauchyNDirDDisplacementsDXi<celltype>(
        solid_formulation_linearization.d2_F_dxi_dd, d_cauchyndir_dF,
        cauchy_n_dir_with_linearization.d2_cauchyndir_dd_dxi);

    return cauchy_n_dir_with_linearization;
  }
}  // namespace

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::SolidEleCalc()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<celltype>(GetGaussRuleStiffnessMatrix<celltype>())),
      mass_matrix_integration_(CreateGaussIntegration<celltype>(GetGaussRuleMassMatrix<celltype>()))
{
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData, HistoryData>::Pack(
    CORE::COMM::PackBuffer& data) const
{
  ElementFormulation::Pack(history_data_, data);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
{
  ElementFormulation::Unpack(position, data, history_data_);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::EvaluateNonlinearForceStiffnessMass(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff{};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass{};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, 1>> force{};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);


  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, nodal_coordinates, history_data_);

  double element_mass = 0.0;
  double element_volume = 0.0;
  ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        ElementFormulation::Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping,
            preparation_data, history_data_,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              const Stress<celltype> stress = EvaluateMaterialStress<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              if (force.has_value())
              {
                ElementFormulation::AddInternalForceVector(linearization, stress,
                    integration_factor, preparation_data, history_data_, *force);
              }

              if (stiff.has_value())
              {
                ElementFormulation::AddStiffnessMatrix(linearization, jacobian_mapping, stress,
                    integration_factor, preparation_data, history_data_, *stiff);
              }

              if (mass.has_value())
              {
                if (equal_integration_mass_stiffness)
                {
                  AddMassMatrix(
                      shape_functions, integration_factor, solid_material.Density(gp), *mass);
                }
                else
                {
                  element_mass += solid_material.Density(gp) * integration_factor;
                  element_volume += integration_factor;
                }
              }
            });
      });

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(element_mass > 0, "It looks like the element mass is 0.0");
    ForEachGaussPoint<celltype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp) {
          AddMassMatrix(shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::Update(const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, nodal_coordinates, history_data_);

  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        ElementFormulation::Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping,
            preparation_data, history_data_,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            { solid_material.Update(deformation_gradient, gp, params, ele.Id()); });
      });

  solid_material.Update();
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
double DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::CalculateInternalEnergy(const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, nodal_coordinates, history_data_);

  double intenergy = 0;
  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        ElementFormulation::Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping,
            preparation_data, history_data_,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              double psi = 0.0;
              solid_material.StrainEnergy(gl_strain, psi, gp, ele.Id());
              intenergy += psi * integration_factor;
            });
      });

  return intenergy;
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
    const StressIO& stressIO, const StrainIO& strainIO, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, nodal_coordinates, history_data_);

  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        ElementFormulation::Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping,
            preparation_data, history_data_,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              const Stress<celltype> stress = EvaluateMaterialStress<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              AssembleStrainTypeToMatrixRow<celltype>(
                  gl_strain, deformation_gradient, strainIO.type, strain_data, gp);
              AssembleStressTypeToMatrixRow(
                  deformation_gradient, stress, stressIO.type, stress_data, gp);
            });
      });

  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData, HistoryData>::Setup(
    MAT::So3Material& solid_material, INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::MaterialPostSetup(const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::InitializeGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
void DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::ResetToLastConverged(const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

template <CORE::FE::CellType celltype, typename ElementFormulation, typename PreparationData,
    typename HistoryData>
DRT::ELEMENTS::CauchyNDirAndLinearization<3>
DRT::ELEMENTS::SolidEleCalc<celltype, ElementFormulation, PreparationData,
    HistoryData>::GetCauchyNDirAndDerivativesAtXi(const DRT::Element& ele,
    MAT::So3Material& solid_material, const std::vector<double>& disp,
    const CORE::LINALG::Matrix<3, 1>& xi, const CORE::LINALG::Matrix<3, 1>& n,
    const CORE::LINALG::Matrix<3, 1>& dir)
{
  ElementNodes<celltype> element_nodes = EvaluateElementNodes<celltype>(ele, disp);

  const ShapeFunctionsAndDerivatives<celltype> shape_functions =
      EvaluateShapeFunctionsAndDerivs<celltype>(xi, element_nodes);

  const JacobianMapping<celltype> jacobian_mapping =
      EvaluateJacobianMapping(shape_functions, element_nodes);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, element_nodes, history_data_);

  return ElementFormulation::Evaluate(ele, element_nodes, xi, shape_functions, jacobian_mapping,
      preparation_data, history_data_,
      [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
              deformation_gradient,
          const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
      {
        SolidFormulationLinearization<celltype> solid_linearization =
            ElementFormulation::EvaluateFullLinearization(ele, element_nodes, xi, shape_functions,
                jacobian_mapping, deformation_gradient, preparation_data, history_data_);

        return EvaluateCauchyNDirAndDerivatives<celltype>(
            solid_material, deformation_gradient, solid_linearization, n, dir, ele.Id());
      });
}

template <CORE::FE::CellType... celltypes>
struct VerifyPackable
{
  static constexpr bool are_all_packable =
      (DRT::ELEMENTS::IsPackable<DRT::ELEMENTS::SolidEleCalc<celltypes,
              DRT::ELEMENTS::DisplacementBasedFormulation<celltypes>,
              DRT::ELEMENTS::DisplacementBasedPreparationData,
              DRT::ELEMENTS::DisplacementBasedHistoryData>*> &&
          ...);

  static constexpr bool are_all_unpackable =
      (DRT::ELEMENTS::IsUnpackable<DRT::ELEMENTS::SolidEleCalc<celltypes,
              DRT::ELEMENTS::DisplacementBasedFormulation<celltypes>,
              DRT::ELEMENTS::DisplacementBasedPreparationData,
              DRT::ELEMENTS::DisplacementBasedHistoryData>*> &&
          ...);

  void StaticAsserts() const
  {
    static_assert(are_all_packable);
    static_assert(are_all_unpackable);
  }
};

template struct VerifyPackable<CORE::FE::CellType::hex8, CORE::FE::CellType::hex18,
    CORE::FE::CellType::hex20, CORE::FE::CellType::hex27, CORE::FE::CellType::nurbs27,
    CORE::FE::CellType::tet4, CORE::FE::CellType::tet10, CORE::FE::CellType::pyramid5,
    CORE::FE::CellType::wedge6, CORE::FE::CellType::hex8, CORE::FE::CellType::hex8,
    CORE::FE::CellType::hex8, CORE::FE::CellType::hex8, CORE::FE::CellType::hex8>;

// explicit instantiations of template classes
// for displacement based formulation
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex8>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex18,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex18>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex20,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex20>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex27,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex27>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::nurbs27,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::nurbs27>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet4,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::tet4>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet10,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::tet10>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::pyramid5,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::pyramid5>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::wedge6,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::wedge6>,
    DRT::ELEMENTS::DisplacementBasedPreparationData, DRT::ELEMENTS::DisplacementBasedHistoryData>;

// for displacement based formulation with linear kinematics
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::hex8>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex18,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::hex18>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex20,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::hex20>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex27,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::hex27>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::nurbs27,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::nurbs27>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet4,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::tet4>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet10,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::tet10>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::pyramid5,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::pyramid5>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::wedge6,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsFormulation<CORE::FE::CellType::wedge6>,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsPreparationData,
    DRT::ELEMENTS::DisplacementBasedLinearKinematicsHistoryData>;

// Fbar element technology
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::FBarFormulation<CORE::FE::CellType::hex8>,
    DRT::ELEMENTS::FBarPreparationData<CORE::FE::CellType::hex8>, DRT::ELEMENTS::FBarHistoryData>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::pyramid5,
    DRT::ELEMENTS::FBarFormulation<CORE::FE::CellType::pyramid5>,
    DRT::ELEMENTS::FBarPreparationData<CORE::FE::CellType::pyramid5>,
    DRT::ELEMENTS::FBarHistoryData>;

BACI_NAMESPACE_CLOSE