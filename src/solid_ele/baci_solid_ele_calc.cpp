/*! \file

\brief Implementation of routines for calculation of solid element simple displacement based

\level 1
*/

#include "baci_solid_ele_calc.H"

#include "baci_mat_so3_material.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_calc_lib_integration.H"
#include "baci_solid_ele_calc_lib_io.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>


template <CORE::FE::CellType distype>
DRT::ELEMENTS::SolidEleCalc<distype>::SolidEleCalc()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNonlinearForceStiffnessMass(
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

  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);

  double mean_density = 0.0;

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<num_dim_, num_dim_> cauchygreen =
            EvaluateCauchyGreen<distype>(spatial_material_mapping);

        const CORE::LINALG::Matrix<num_str_, 1> gl_strain =
            EvaluateGreenLagrangeStrain<distype>(cauchygreen);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        const Stress<distype> stress = EvaluateMaterialStress<distype>(solid_material,
            spatial_material_mapping.deformation_gradient_, gl_strain, params, gp, ele.Id());

        if (force.has_value()) AddInternalForceVector(Bop, stress, integration_factor, *force);

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


  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    ForEachGaussPoint<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNonlinearForceStiffnessMassGEMM(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  const double gemmalphaf = params.get<double>("alpha f");
  const double gemmxi = params.get<double>("xi");

  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff = {};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass = {};
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  const ElementNodes<distype> nodal_coordinates_old =
      EvaluateElementNodesOfPreviousTimestep<distype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);

  double mean_density = 0.0;

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        SpatialMaterialMapping<distype> spatial_material_mapping_old =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates_old);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop_old =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping_old);

        CORE::LINALG::Matrix<num_dim_, num_dim_> cauchy_green =
            EvaluateCauchyGreen(spatial_material_mapping);
        CORE::LINALG::Matrix<num_dim_, num_dim_> cauchy_green_old =
            EvaluateCauchyGreen(spatial_material_mapping_old);

        CORE::LINALG::Matrix<num_str_, 1> gl_strain =
            EvaluateGreenLagrangeStrain<distype>(cauchy_green);
        CORE::LINALG::Matrix<num_str_, 1> gl_strain_old =
            EvaluateGreenLagrangeStrain<distype>(cauchy_green_old);

        // computed averaged mid-point quantities
        // 1. non-linear mid-B-operator from Bop and Bop_old
        // B_{m} = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> BopM(true);
        BopM.Update(1.0 - gemmalphaf, Bop, gemmalphaf, Bop_old);

        // 2. mid-strain GL-vector from gl_strains and gl_strains_old
        // E_{m} = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_{n}
        CORE::LINALG::Matrix<num_str_, 1> gl_strains_m(true);
        gl_strains_m.Update(
            1.0 - gemmalphaf + gemmxi, gl_strain, gemmalphaf - gemmxi, gl_strain_old);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        const Stress<distype> stress =
            EvaluateMaterialStressGEMM<distype>(solid_material, gl_strain, gl_strain_old,
                gl_strains_m, cauchy_green, cauchy_green_old, params, gp, ele.Id());

        if (force.has_value()) AddInternalForceVector(BopM, stress, integration_factor, *force);

        if (stiff.has_value())
        {
          AddElasticStiffnessMatrixGEMM(
              Bop, BopM, stress, integration_factor * (1.0 - gemmalphaf + gemmxi), *stiff);
          AddGeometricStiffnessMatrix(
              jacobian_mapping.N_XYZ_, stress, integration_factor * (1.0 - gemmalphaf), *stiff);
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


  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    ForEachGaussPoint<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        solid_material.Update(spatial_material_mapping.deformation_gradient_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <CORE::FE::CellType distype>
double DRT::ELEMENTS::SolidEleCalc<distype>::CalculateInternalEnergy(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<num_dim_, num_dim_> cauchygreen =
            EvaluateCauchyGreen<distype>(spatial_material_mapping);

        const CORE::LINALG::Matrix<num_str_, 1> gl_strain =
            EvaluateGreenLagrangeStrain<distype>(cauchygreen);

        double psi = 0.0;
        solid_material.StrainEnergy(gl_strain, psi, gp, ele.Id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::CalculateStress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const StressIO& stressIO, const StrainIO& strainIO,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CORE::LINALG::Matrix<num_dim_, num_dim_> cauchygreen =
            EvaluateCauchyGreen<distype>(spatial_material_mapping);

        const CORE::LINALG::Matrix<num_str_, 1> gl_strain =
            EvaluateGreenLagrangeStrain<distype>(cauchygreen);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        const Stress<distype> stress = EvaluateMaterialStress<distype>(solid_material,
            spatial_material_mapping.deformation_gradient_, gl_strain, params, gp, ele.Id());

        AssembleStrainTypeToMatrixRow<distype>(gl_strain,
            spatial_material_mapping.deformation_gradient_, strainIO.type, strain_data, gp);
        AssembleStressTypeToMatrixRow(
            spatial_material_mapping.deformation_gradient_, stress, stressIO.type, stress_data, gp);
      });


  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<distype>(
      stiffness_matrix_integration_, ele, params);

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<distype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex18>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::nurbs27>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::SolidEleCalc<CORE::FE::CellType::wedge6>;
