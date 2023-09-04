/*! \file

\brief Implementation of routines for calculation of solid element with fbar element technology

\level 1
*/

#include "baci_solid_ele_calc_fbar.H"

#include "baci_mat_so3_material.H"
#include "baci_solid_ele_calc_lib.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>

namespace
{
  /*!
   * @brief Evaluate the fbar factor \f[ \frac{\mathbf{F}_{\mathrm{centroid}}}{\mathbf{F}}^{1/3} \f]
   *
   * @param defgrd_centroid (in) : Deformation gradient evaluated at the element centroid
   * @param defgrd_gp (in) : Deformation gradient evaluated at the Gauss point
   * @return double : Fbar factor
   */
  double EvaluateFbarFactor(const double& defgrd_centroid, const double& defgrd_gp)
  {
    const double fbar_factor = std::pow(defgrd_centroid / defgrd_gp, 1.0 / 3.0);
    return fbar_factor;
  }

  /*!
   * @brief Evaluates the H-Operator used in F-bar of the specified element
   *
   * @tparam distype : Discretization type
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping evaluated at the Gauss point
   * @param jacobian_mapping_centroid (in) : Quantities of the jacobian mapping evaluated at the
   * element centroid
   * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient) evaluated at the Gauss point
   * @param spatial_material_mapping_centroid (in) : An object holding quantities of the spatial
   * material mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient) evaluated at the element centroid
   * @return CORE::LINALG::Matrix<num_dof_per_ele, 1> : H-Operator
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::ELEMENTS::DETAIL::num_dim<distype> == 3, int> = 0>
  const CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>, 1>
  EvaluateFbarHOperator(const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping,
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping_centroid,
      const DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping,
      const DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping_centroid)
  {
    // inverse deformation gradient at centroid
    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        invdefgrd_centroid;
    invdefgrd_centroid.Invert(spatial_material_mapping_centroid.deformation_gradient_);

    // inverse deformation gradient at gp
    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        invdefgrd;
    invdefgrd.Invert(spatial_material_mapping.deformation_gradient_);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>, 1> Hop(true);
    for (int idof = 0; idof < DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>; idof++)
    {
      for (int idim = 0; idim < DRT::ELEMENTS::DETAIL::num_dim<distype>; idim++)
      {
        Hop(idof) +=
            invdefgrd_centroid(idim, idof % DRT::ELEMENTS::DETAIL::num_dim<distype>) *
            jacobian_mapping_centroid.N_XYZ_(idim, idof / DRT::ELEMENTS::DETAIL::num_dim<distype>);
        Hop(idof) -= invdefgrd(idim, idof % DRT::ELEMENTS::DETAIL::num_dim<distype>) *
                     jacobian_mapping.N_XYZ_(idim, idof / DRT::ELEMENTS::DETAIL::num_dim<distype>);
      }
    }

    return Hop;
  }

  /*!
   * @brief Add fbar stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param Hop (in) : H-Operator
   * @param f_bar_factor (in) : f_bar_factor
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param cauchyGreen (in) : An object holding the right Cauchy-Green deformation tensor and
   * its inverse
   * @param stress_bar (in) : Deviatoric part of stress measures
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddFbarStiffnessMatrix(const CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>,
                                  DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>>& Bop,
      const CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>, 1>& Hop,
      const double f_bar_factor, const double integration_fac,
      const DRT::ELEMENTS::CauchyGreen<distype> cauchyGreen,
      const DRT::ELEMENTS::Stress<distype> stress_bar,
      CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>,
          DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>>& stiffness_matrix)
  {
    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>, 1> rcg_bar_voigt;
    ::UTILS::VOIGT::Strains::MatrixToVector(cauchyGreen.right_cauchy_green_, rcg_bar_voigt);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>, 1> ccg;
    ccg.MultiplyNN(stress_bar.cmat_, rcg_bar_voigt);

    // auxiliary integrated stress_bar
    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>, 1> bopccg(false);
    bopccg.MultiplyTN(integration_fac * f_bar_factor / 3.0, Bop, ccg);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>, 1> bops(false);
    bops.MultiplyTN(-integration_fac / f_bar_factor / 3.0, Bop, stress_bar.pk2_);

    for (int idof = 0; idof < DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>; idof++)
    {
      for (int jdof = 0; jdof < DRT::ELEMENTS::DETAIL::num_dof_per_ele<distype>; jdof++)
      {
        stiffness_matrix(idof, jdof) += Hop(jdof) * (bops(idof, 0) + bopccg(idof, 0));
      }
    }
  }

  template <DRT::Element::DiscretizationType distype>
  CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>, 1> EvaluateStrainsBar(
      const DRT::ELEMENTS::NodalCoordinates<distype>& nodal_coordinates,
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping, double detF_centroid)
  {
    const DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping =
        EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

    // factor (detF0/detF)^1/3
    const double fbar_factor = EvaluateFbarFactor(
        detF_centroid, spatial_material_mapping.determinant_deformation_gradient_);

    // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
    const DRT::ELEMENTS::SpatialMaterialMapping<distype> spatial_material_mapping_fbar_factor =
        EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

    const DRT::ELEMENTS::CauchyGreen<distype> cauchygreen_fbar_factor =
        EvaluateCauchyGreen(spatial_material_mapping_fbar_factor);

    return EvaluateGreenLagrangeStrain(cauchygreen_fbar_factor);
  }
}  // namespace

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalcFbar<distype>::SolidEleCalcFbar()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::EvaluateNonlinearForceStiffnessMass(
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

  // get current nodal coordinates
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);

  // initialize element density
  double mean_density = 0.0;

  // jacobian mapping evaluated at element centroid
  const JacobianMapping<distype> jacobian_mapping_centroid =
      EvaluateJacobianMappingCentroid(nodal_coordinates);

  // deformation gradient at element centroid
  const SpatialMaterialMapping<distype> spatial_material_mapping_centroid =
      EvaluateSpatialMaterialMapping(jacobian_mapping_centroid, nodal_coordinates);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        // factor (detF0/detF)^1/3
        const double fbar_factor =
            EvaluateFbarFactor(spatial_material_mapping_centroid.determinant_deformation_gradient_,
                spatial_material_mapping.determinant_deformation_gradient_);

        const CORE::LINALG::Matrix<num_dof_per_ele_, 1> Hop =
            EvaluateFbarHOperator(jacobian_mapping, jacobian_mapping_centroid,
                spatial_material_mapping, spatial_material_mapping_centroid);

        // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
        const SpatialMaterialMapping<distype> spatial_material_mapping_bar =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

        const CauchyGreen<distype> cauchygreen_bar =
            EvaluateCauchyGreen(spatial_material_mapping_bar);

        CORE::LINALG::Matrix<DETAIL::num_str<distype>, 1> gl_strain_bar =
            EvaluateGreenLagrangeStrain(cauchygreen_bar);

        // stress stress_bar evaluated using strains_bar
        const Stress<distype> stress_bar = EvaluateMaterialStress<distype>(solid_material,
            spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, params, gp,
            ele.Id());

        // evaluate internal force vector using only the deviatoric component
        if (force.has_value())
          AddInternalForceVector(Bop, stress_bar, integration_factor / fbar_factor, *force);

        if (stiff.has_value())
        {
          // evaluate elastic stiffness matrix using only the deviatoric component Fbar
          AddElasticStiffnessMatrix(Bop, stress_bar, integration_factor * fbar_factor, *stiff);

          // evaluate geometric stiffness matrix using only the deviatoric component Fbar
          AddGeometricStiffnessMatrix(
              jacobian_mapping.N_XYZ_, stress_bar, integration_factor / fbar_factor, *stiff);

          // additional stiffness matrix needed for fbar method
          AddFbarStiffnessMatrix(
              Bop, Hop, fbar_factor, integration_factor, cauchygreen, stress_bar, *stiff);
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
        [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  // Loop over all Gauss points
  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const double fbar_factor = EvaluateFbarFactor(
            detF_centroid, spatial_material_mapping.determinant_deformation_gradient_);

        const SpatialMaterialMapping<distype> spatial_material_mapping_bar =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

        solid_material.Update(
            spatial_material_mapping_bar.deformation_gradient_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::CalculateStress(const DRT::Element& ele,
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

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  // Loop over all Gauss points
  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>, 1> gl_strains_bar =
            EvaluateStrainsBar(nodal_coordinates, jacobian_mapping, detF_centroid);

        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

        const double fbar_factor = EvaluateFbarFactor(
            detF_centroid, spatial_material_mapping.determinant_deformation_gradient_);

        const SpatialMaterialMapping<distype> spatial_material_mapping_bar =
            EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

        const Stress<distype> stress_bar = EvaluateMaterialStress<distype>(solid_material,
            spatial_material_mapping_bar.deformation_gradient_, gl_strains_bar, params, gp,
            ele.Id());

        AssembleStrainTypeToMatrixRow<distype>(gl_strains_bar,
            spatial_material_mapping_bar.deformation_gradient_, strainIO.type, strain_data, gp);
        AssembleStressTypeToMatrixRow(spatial_material_mapping_bar.deformation_gradient_,
            stress_bar, stressIO.type, stress_data, gp);
      });

  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::SolidEleCalcFbar<distype>::CalculateInternalEnergy(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  ForEachGaussPoint<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_str<distype>, 1> gl_strains_bar =
            EvaluateStrainsBar(nodal_coordinates, jacobian_mapping, detF_centroid);

        double psi = 0.0;
        solid_material.StrainEnergy(gl_strains_bar, psi, gp, ele.Id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<distype>(
      stiffness_matrix_integration_, ele, params);

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::InitializeGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<distype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::ResetAll(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetAll(stiffness_matrix_integration_.NumPoints());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalcFbar<DRT::Element::hex8>;
template class DRT::ELEMENTS::SolidEleCalcFbar<DRT::Element::pyramid5>;