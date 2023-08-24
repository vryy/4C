/*! \file

\brief Implementation of routines for calculation of solid element simple displacement based

\level 1
*/

#include "baci_solid_ele_calc.H"

#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_fiber_nodal_fiber_holder.H"
#include "baci_fiber_node.H"
#include "baci_fiber_utils.H"
#include "baci_lib_discret.H"
#include "baci_lib_utils.H"
#include "baci_lib_voigt_notation.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_element_service.H"
#include "baci_solid_ele.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_utils.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>


template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalc<distype>::SolidEleCalc()
    : DRT::ELEMENTS::SolidEleCalcInterface::SolidEleCalcInterface(),
      stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>> stiff = {};
  std::optional<CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>> mass = {};
  std::optional<CORE::LINALG::Matrix<numdofperelement_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // TODO: This is a quite unsafe check, whether the same integrations are used
  bool equal_integration_mass_stiffness =
      mass_matrix_integration_.NumPoints() == stiffness_matrix_integration_.NumPoints();

  double mean_density = 0.0;

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        const CORE::LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain =
            EvaluateGreenLagrangeStrain(cauchygreen);

        CORE::LINALG::Matrix<numstr_, numdofperelement_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

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
  {  // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNonlinearForceStiffnessMassGEMM(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  const double gemmalphaf = params.get<double>("alpha f");
  const double gemmxi = params.get<double>("xi");

  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>> stiff = {};
  std::optional<CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>> mass = {};
  std::optional<CORE::LINALG::Matrix<numdofperelement_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  const NodalCoordinates<distype> nodal_coordinates_old =
      EvaluateNodalCoordinatesOld<distype>(ele, discretization, lm);

  // TODO: This is a quite unsafe check, whether the same integrations are used
  bool equal_integration_mass_stiffness =
      mass_matrix_integration_.NumPoints() == stiffness_matrix_integration_.NumPoints();

  double mean_density = 0.0;

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        SpatialMaterialMapping<distype> spatial_material_mapping_old =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates_old);

        CORE::LINALG::Matrix<numstr_, numdofperelement_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        CORE::LINALG::Matrix<numstr_, numdofperelement_> Bop_old =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping_old);

        CauchyGreen<distype> cauchy_green = EvaluateCauchyGreen(spatial_material_mapping);
        CauchyGreen<distype> cauchy_green_old = EvaluateCauchyGreen(spatial_material_mapping_old);

        CORE::LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain =
            EvaluateGreenLagrangeStrain(cauchy_green);
        CORE::LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain_old =
            EvaluateGreenLagrangeStrain(cauchy_green_old);

        // computed averaged mid-point quantities
        // 1. non-linear mid-B-operator from Bop and Bop_old
        // B_{m} = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
        CORE::LINALG::Matrix<numstr_, numdofperelement_> BopM(true);
        BopM.Update(1.0 - gemmalphaf, Bop, gemmalphaf, Bop_old);

        // 2. mid-strain GL-vector from gl_strains and gl_strains_old
        // E_{m} = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_{n}
        CORE::LINALG::Matrix<numstr_, 1> gl_strains_m(true);
        gl_strains_m.Update(
            1.0 - gemmalphaf + gemmxi, gl_strain, gemmalphaf - gemmxi, gl_strain_old);

        const Stress<distype> stress = EvaluateMaterialStressGEMM(solid_material, gl_strain,
            gl_strain_old, gl_strains_m, cauchy_green, cauchy_green_old, params, gp, ele.Id());

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
  {  // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // nothing to do for a standard element
  // except...
  // TODO: We need to recover history information of materials!
  // which was also not implemented in the old elements
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        solid_material.Update(spatial_material_mapping.deformation_gradient_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::SolidEleCalc<distype>::CalculateInternalEnergy(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        const CORE::LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain =
            EvaluateGreenLagrangeStrain(cauchygreen);

        double psi = 0.0;
        solid_material.StrainEnergy(gl_strain, psi, gp, ele.Id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::CalculateStress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const StressIO& stressIO, const StrainIO& strainIO,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // TODO: If we get rid of post_drt_*, we don't need this here anymore. We could directly use
  // InitializeGaussPointDataOutput and EvaluateGaussPointDataOutput and write the stresses there.
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), numstr_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), numstr_);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        const CORE::LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain =
            EvaluateGreenLagrangeStrain(cauchygreen);

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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};
  if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const static std::vector<CORE::LINALG::Matrix<nen_, 1>> shapefcts = std::invoke(
        [&]
        {
          std::vector<CORE::LINALG::Matrix<nen_, 1>> shapefcns(
              stiffness_matrix_integration_.NumPoints());
          for (int gp = 0; gp < stiffness_matrix_integration_.NumPoints(); ++gp)
          {
            CORE::LINALG::Matrix<nsd_, 1> xi(stiffness_matrix_integration_.Point(gp), true);
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
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
    CORE::LINALG::SerialDenseMatrix gp_data(
        stiffness_matrix_integration_.NumPoints(), quantity_size, true);
    bool data_available = solid_material.EvaluateVtkOutputData(quantity_name, gp_data);

    // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::ResetAll(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetAll(stiffness_matrix_integration_.NumPoints());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex18>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::wedge6>;
