/*! \file

\brief Implementation of routines for calculation of solid element with templated element
formulation

\level 1
*/

#include "baci_solid_ele_calc.H"

#include "baci_discretization_fem_general_cell_type.H"
#include "baci_mat_so3_material.H"
#include "baci_solid_ele_calc_displacement_based.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_calc_lib_integration.H"
#include "baci_solid_ele_calc_lib_io.H"
#include "baci_solid_ele_interface_serializable.H"

#include <Teuchos_ParameterList.hpp>

#include <optional>

BACI_NAMESPACE_OPEN


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

  double mean_density = 0.0;

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData preparation_data =
      ElementFormulation::Prepare(ele, nodal_coordinates, history_data_);

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
                  mean_density +=
                      solid_material.Density(gp) / stiffness_matrix_integration_.NumPoints();
                }
              }
            });
      });

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    ForEachGaussPoint<celltype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
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
  if (discretization.Comm().MyPID() != ele.Owner()) return;

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
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
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

BACI_NAMESPACE_CLOSE