/*! \file

\brief Implementation of routines for calculation of solid element displacement
based with MULF prestressing

\level 1
*/

#include "baci_solid_ele_calc_mulf.H"

#include "baci_lib_parobject.H"
#include "baci_linalg_fixedsizematrix_generators.H"
#include "baci_mat_so3_material.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>

namespace
{
  template <CORE::FE::CellType distype>
  inline static constexpr int num_dim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  template <CORE::FE::CellType distype>
  inline static constexpr int num_nodes = CORE::FE::num_nodes<distype>;

  template <CORE::FE::CellType distype>
  CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> EvaluateMulfDeformationGradientUpdate(
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping,
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<distype>& shape_functions,
      const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>>& nodal_displacements,
      const DRT::ELEMENTS::MulfHistoryData<distype>& mulf_history_data,
      const INPAR::STR::KinemType& kinematictype = INPAR::STR::kinem_nonlinearTotLag)
  {
    if (kinematictype == INPAR::STR::kinem_linear)
    {
      dserror(
          "Linear kinematics have not yet been implemented for MULF "
          "prestressing");
    }

    CORE::LINALG::Matrix<num_dim<distype>, num_nodes<distype>> N_xyz;

    N_xyz.Multiply(mulf_history_data.inverse_jacobian, shape_functions.derivatives_);

    CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> defgrd =
        CORE::LINALG::IdentityMatrix<num_dim<distype>>();

    defgrd.MultiplyTT(1.0, nodal_displacements, N_xyz, 1.0);

    return defgrd;
  }

  template <CORE::FE::CellType distype>
  DRT::ELEMENTS::SpatialMaterialMapping<distype> EvaluateMulfSpatialMaterialMapping(
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping,
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<distype>& shape_functions,
      const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>>& nodal_displacements,
      const DRT::ELEMENTS::MulfHistoryData<distype>& mulf_history_data,
      const INPAR::STR::KinemType& kinematictype = INPAR::STR::kinem_nonlinearTotLag)
  {
    DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping;

    CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> defgrd =
        EvaluateMulfDeformationGradientUpdate(jacobian_mapping, shape_functions,
            nodal_displacements, mulf_history_data, kinematictype);

    spatial_material_mapping.deformation_gradient_.Multiply(
        defgrd, mulf_history_data.deformation_gradient);

    spatial_material_mapping.inverse_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_;
    spatial_material_mapping.inverse_deformation_gradient_.Invert();

    return spatial_material_mapping;
  }

}  // namespace

template <CORE::FE::CellType distype>
DRT::ELEMENTS::SolidEleCalcMulf<distype>::SolidEleCalcMulf()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::Pack(DRT::PackBuffer& data) const
{
  data.AddtoPack(history_data_.size());
  for (const auto& item : history_data_)
  {
    DRT::ParObject::AddtoPack(data, item.inverse_jacobian);
    DRT::ParObject::AddtoPack(data, item.deformation_gradient);
    DRT::ParObject::AddtoPack(data, item.is_setup);
  }
};

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  std::size_t size;
  DRT::ParObject::ExtractfromPack(position, data, size);
  history_data_.resize(size);

  for (auto& item : history_data_)
  {
    DRT::ParObject::ExtractfromPack(position, data, item.inverse_jacobian);
    DRT::ParObject::ExtractfromPack(position, data, item.deformation_gradient);
    int is_setup_int;
    DRT::ParObject::ExtractfromPack(position, data, is_setup_int);
    item.is_setup = static_cast<bool>(is_setup_int);
  }
};

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff;
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass;
  std::optional<CORE::LINALG::Matrix<num_dof_per_ele_, 1>> force;
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>> nodal_displacements =
      GetNodalDisplacements<distype>(discretization, lm);

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
        if (!history_data_[gp].is_setup)
        {
          history_data_[gp].inverse_jacobian = jacobian_mapping.inverse_jacobian_;
          history_data_[gp].is_setup = true;
        }
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateMulfSpatialMaterialMapping(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);

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
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::EvaluateNonlinearForceStiffnessMassGEMM(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  dserror("GEMM not implemented for MULF prestressing");
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>> nodal_displacements =
      GetNodalDisplacements<distype>(discretization, lm);
  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateMulfSpatialMaterialMapping(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);

        EvaluateGPCoordinatesAndAddToParameterList<distype>(
            nodal_coordinates, shape_functions, params);

        solid_material.Update(spatial_material_mapping.deformation_gradient_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::UpdatePrestress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>> nodal_displacements =
      GetNodalDisplacements<distype>(discretization, lm);
  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);


  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateMulfSpatialMaterialMapping(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);


        CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> delta_defgrd =
            EvaluateMulfDeformationGradientUpdate(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);

        // update mulf history data only if prestress is active
        CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> inv_delta_defgrd(delta_defgrd);
        inv_delta_defgrd.Invert();

        CORE::LINALG::Matrix<num_dim<distype>, num_dim<distype>> invJ_new;

        invJ_new.MultiplyTN(inv_delta_defgrd, history_data_[gp].inverse_jacobian);

        history_data_[gp].deformation_gradient = spatial_material_mapping.deformation_gradient_;
        history_data_[gp].inverse_jacobian = std::move(invJ_new);


        solid_material.Update(spatial_material_mapping.deformation_gradient_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <CORE::FE::CellType distype>
double DRT::ELEMENTS::SolidEleCalcMulf<distype>::CalculateInternalEnergy(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  double intenergy = 0.0;

  const CORE::LINALG::Matrix<num_nodes_, num_dim_> nodal_displacements =
      GetNodalDisplacements<distype>(discretization, lm);

  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateMulfSpatialMaterialMapping(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);

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
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::CalculateStress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const StressIO& stressIO, const StrainIO& strainIO,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const CORE::LINALG::Matrix<num_nodes<distype>, num_dim<distype>> nodal_displacements =
      GetNodalDisplacements<distype>(discretization, lm);

  const ElementNodes<distype> nodal_coordinates =
      EvaluateElementNodes<distype>(ele, discretization, lm);

  EvaluateCentroidCoordinatesAndAddToParameterList<distype>(nodal_coordinates, params);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateMulfSpatialMaterialMapping(
                jacobian_mapping, shape_functions, nodal_displacements, history_data_[gp]);

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
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  history_data_.resize(stiffness_matrix_integration_.NumPoints(), {});
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points
  // and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<distype>(
      stiffness_matrix_integration_, ele, params);

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::InitializeGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called "
      "from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called "
      "from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<distype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::SolidEleCalcMulf<distype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex18>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::nurbs27>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::wedge6>;

static_assert(DRT::ELEMENTS::IsPrestressUpdateable<
                  DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex8>*>,
    "MULF calculation interface needs to be prestress updatable. Please carefully check the "
    "signature of UpdatePrestress(...)");
static_assert(DRT::ELEMENTS::IsPackable<DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex8>*>,
    "MULF calculation interface needs to be packable. Please carefully check the "
    "signature of Pack(...)");
static_assert(
    DRT::ELEMENTS::IsUnpackable<DRT::ELEMENTS::SolidEleCalcMulf<CORE::FE::CellType::hex8>*>,
    "MULF calculation interface needs to be unpackable. Please carefully check the "
    "signature of Unpack(...)");