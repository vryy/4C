/*! \file

\brief Implementation of routines for calculation of a coupled solid-scatra element with templated
solid formulation

\level 1
*/

#include "4C_solid_scatra_3D_ele_calc.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType celltype>
  inline static constexpr int num_str = CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2;

  template <CORE::FE::CellType celltype>
  CORE::LINALG::Matrix<num_str<celltype>, 1> EvaluateDMaterialStressDScalar(
      MAT::So3Material& solid_material,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
          deformation_gradient,
      const CORE::LINALG::Matrix<num_str<celltype>, 1>& gl_strain, Teuchos::ParameterList& params,
      const int gp, const int eleGID)
  {
    CORE::LINALG::Matrix<num_str<celltype>, 1> dStressDScalar(true);

    // The derivative of the solid stress w.r.t. the scalar is implemented in the normal material
    // Evaluate call by not passing the linearization matrix.
    solid_material.Evaluate(
        &deformation_gradient, &gl_strain, params, &dStressDScalar, nullptr, gp, eleGID);

    return dStressDScalar;
  }

  template <CORE::FE::CellType celltype>
  auto ProjectQuantityToGaussPoint(
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const std::vector<CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, 1>>& nodal_quantities)
  {
    std::vector<double> quantities_at_gp(nodal_quantities.size(), 0.0);

    for (std::size_t k = 0; k < nodal_quantities.size(); ++k)
    {
      quantities_at_gp[k] = shape_functions.shapefunctions_.Dot(nodal_quantities[k]);
    }
    return quantities_at_gp;
  }

  template <CORE::FE::CellType celltype>
  auto ProjectQuantityToGaussPoint(
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, 1>& nodal_quantity)
  {
    return shape_functions.shapefunctions_.Dot(nodal_quantity);
  }


  template <bool is_scalar, CORE::FE::CellType celltype>
  auto GetElementQuantities(const int num_scalars, const std::vector<double>& quantities_at_dofs)
  {
    if constexpr (is_scalar)
    {
      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, 1> nodal_quantities(num_scalars);
      for (int i = 0; i < CORE::FE::num_nodes<celltype>; ++i)
        nodal_quantities(i, 0) = quantities_at_dofs.at(i);

      return nodal_quantities;
    }
    else
    {
      std::vector<CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, 1>> nodal_quantities(
          num_scalars);

      for (int k = 0; k < num_scalars; ++k)
        for (int i = 0; i < CORE::FE::num_nodes<celltype>; ++i)
          (nodal_quantities[k])(i, 0) = quantities_at_dofs.at(num_scalars * i + k);

      return nodal_quantities;
    }
  }

  template <bool is_scalar, CORE::FE::CellType celltype>
  void PrepareScatraQuantityInParameterList(const DRT::Discretization& discretization,
      const DRT::Element::LocationArray& la,
      const DRT::ELEMENTS::ElementNodes<celltype>& element_nodes, const std::string& field_name,
      const int field_index, const int num_scalars,
      const CORE::FE::GaussIntegration& gauss_integration, Teuchos::ParameterList& params,
      const std::string& target_name)
  {
    FOUR_C_ASSERT(discretization.HasState(field_index, field_name),
        "Could not find the requested field in the discretization.");

    FOUR_C_ASSERT(
        !is_scalar || num_scalars == 1, "numscalars must be 1 if result type is not a vector!");

    // get quantitiy from discretization
    Teuchos::RCP<const Epetra_Vector> quantitites_np =
        discretization.GetState(field_index, field_name);
    if (quantitites_np == Teuchos::null)
      FOUR_C_THROW("Cannot get state vector '%s' ", field_name.c_str());

    // extract my values
    auto my_quantities = std::vector<double>(la[field_index].lm_.size(), 0.0);
    CORE::FE::ExtractMyValues(*quantitites_np, my_quantities, la[field_index].lm_);

    // get nodal quantities for the scalars
    auto nodal_quantities = GetElementQuantities<is_scalar, celltype>(num_scalars, my_quantities);

    // Determine quantity type at Gauss point (can be double (for scalar) or a std::vector of
    // scalars)
    using gp_quantity_type = decltype(ProjectQuantityToGaussPoint(
        std::declval<const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>&>(),
        nodal_quantities));

    // the material expects the gp-quantities at the Gauss points in a rcp-std::vector
    auto quantity_at_gp =
        Teuchos::rcp(new std::vector<gp_quantity_type>(gauss_integration.NumPoints()));

    DRT::ELEMENTS::ForEachGaussPoint(element_nodes, gauss_integration,
        [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, 1>& xi,
            const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
            double integration_factor, int gp)
        {
          // Project to Gauss point
          (*quantity_at_gp)[gp] = ProjectQuantityToGaussPoint(shape_functions, nodal_quantities);
        });

    params.set(target_name, quantity_at_gp);
  };

  template <CORE::FE::CellType celltype>
  void PrepareScatraQuantitiesInParameterList(const DRT::Element& element,
      const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
      const DRT::ELEMENTS::ElementNodes<celltype>& element_nodes,
      const CORE::FE::GaussIntegration& gauss_integration, Teuchos::ParameterList& params)
  {
    if (la.Size() > 1)
    {
      // prepare data from the scatra-field
      if (discretization.HasState(1, "scalarfield"))
      {
        const int num_scalars = discretization.NumDof(1, element.Nodes()[0]);
        constexpr bool is_scalar = false;
        PrepareScatraQuantityInParameterList<is_scalar>(discretization, la, element_nodes,
            "scalarfield", 1, num_scalars, gauss_integration, params, "gp_conc");
      }

      // additionally prepare temperature-filed if available
      if (discretization.NumDofSets() == 3 && discretization.HasState(2, "tempfield"))
      {
        constexpr bool is_scalar = true;
        PrepareScatraQuantityInParameterList<is_scalar>(discretization, la, element_nodes,
            "tempfield", 2, 1, gauss_integration, params, "gp_temp");
      }
    }
  }
}  // namespace

template <CORE::FE::CellType celltype, typename SolidFormulation>
DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::SolidScatraEleCalc()
    : stiffness_matrix_integration_(
          CreateGaussIntegration<celltype>(GetGaussRuleStiffnessMatrix<celltype>())),
      mass_matrix_integration_(CreateGaussIntegration<celltype>(GetGaussRuleMassMatrix<celltype>()))
{
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Pack(
    CORE::COMM::PackBuffer& data) const
{
  DRT::ELEMENTS::Pack(data, history_data_);
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  DRT::ELEMENTS::Unpack(position, data, history_data_);
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_nonlinear_force_stiffness_mass(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const DRT::Element::LocationArray& la, Teuchos::ParameterList& params,
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
      EvaluateElementNodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  PrepareScatraQuantitiesInParameterList(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  bool equal_integration_mass_stiffness =
      CompareGaussIntegration(mass_matrix_integration_, stiffness_matrix_integration_);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  double element_mass = 0.0;
  double element_volume = 0.0;
  ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              const Stress<celltype> stress = EvaluateMaterialStress<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              if (force.has_value())
              {
                add_internal_force_vector(linearization, stress, integration_factor,
                    preparation_data, history_data_, gp, *force);
              }

              if (stiff.has_value())
              {
                AddStiffnessMatrix(linearization, jacobian_mapping, stress, integration_factor,
                    preparation_data, history_data_, gp, *stiff);
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
    FOUR_C_ASSERT(element_mass > 0, "It looks like the element mass is 0.0");
    ForEachGaussPoint<celltype>(nodal_coordinates, mass_matrix_integration_,
        [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp) {
          AddMassMatrix(shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::evaluate_d_stress_d_scalar(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& stiffness_matrix_dScalar)
{
  const int scatra_column_stride = std::invoke(
      [&]()
      {
        if (params.isParameter("numscatradofspernode"))
        {
          return params.get<int>("numscatradofspernode");
        }
        return 1;
      });


  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  PrepareScatraQuantitiesInParameterList(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              CORE::LINALG::Matrix<6, 1> dSdc = EvaluateDMaterialStressDScalar<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              // linear B-opeartor
              const CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
                  CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
                  bop = SolidFormulation::GetLinearBOperator(linearization);

              constexpr int num_dof_per_ele =
                  CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>;

              // Assemble matrix
              // k_dS = B^T . dS/dc * detJ * N * w(gp)
              CORE::LINALG::Matrix<num_dof_per_ele, 1> BdSdc(true);
              BdSdc.MultiplyTN(integration_factor, bop, dSdc);

              // loop over rows
              for (int rowi = 0; rowi < num_dof_per_ele; ++rowi)
              {
                const double BdSdc_rowi = BdSdc(rowi, 0);
                // loop over columns
                for (int coli = 0; coli < CORE::FE::num_nodes<celltype>; ++coli)
                {
                  stiffness_matrix_dScalar(rowi, coli * scatra_column_stride) +=
                      BdSdc_rowi * shape_functions.shapefunctions_(coli, 0);
                }
              }
            });
      });
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
    Teuchos::ParameterList& params)
{
  // nothing needs to be done for simple displacement based elements
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const DRT::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  PrepareScatraQuantitiesInParameterList(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>&
                    deformation_gradient,
                const CORE::LINALG::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            { solid_material.Update(deformation_gradient, gp, params, ele.Id()); });
      });

  solid_material.Update();
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
double DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::calculate_internal_energy(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  PrepareScatraQuantitiesInParameterList(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  double intenergy = 0;
  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
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

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::CalculateStress(
    const DRT::Element& ele, MAT::So3Material& solid_material, const StressIO& stressIO,
    const StrainIO& strainIO, const DRT::Discretization& discretization,
    const DRT::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  CORE::LINALG::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  PrepareScatraQuantitiesInParameterList(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  EvaluateCentroidCoordinatesAndAddToParameterList(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  DRT::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        EvaluateGPCoordinatesAndAddToParameterList(nodal_coordinates, shape_functions, params);
        Evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
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

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Setup(
    MAT::So3Material& solid_material, INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::material_post_setup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call post_setup of material
  solid_material.post_setup(params, ele.Id());
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::initialize_gauss_point_data_output(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  AskAndAddQuantitiesToGaussPointDataOutput(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_gauss_point_data_output(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  CollectAndAssembleGaussPointDataOutput<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <CORE::FE::CellType celltype, typename SolidFormulation>
void DRT::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::reset_to_last_converged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <CORE::FE::CellType... celltypes>
struct VerifyPackable
{
  static constexpr bool are_all_packable =
      (DRT::ELEMENTS::IsPackable<DRT::ELEMENTS::SolidScatraEleCalc<celltypes,
              DRT::ELEMENTS::DisplacementBasedFormulation<celltypes>>*> &&
          ...);

  static constexpr bool are_all_unpackable =
      (DRT::ELEMENTS::IsUnpackable<DRT::ELEMENTS::SolidScatraEleCalc<celltypes,
              DRT::ELEMENTS::DisplacementBasedFormulation<celltypes>>*> &&
          ...);

  void StaticAsserts() const
  {
    static_assert(are_all_packable);
    static_assert(are_all_unpackable);
  }
};

template struct VerifyPackable<CORE::FE::CellType::hex8, CORE::FE::CellType::hex27,
    CORE::FE::CellType::tet4, CORE::FE::CellType::tet10>;

// explicit instantiations of template classes
// for displacement based formulation
template class DRT::ELEMENTS::SolidScatraEleCalc<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex8>>;
template class DRT::ELEMENTS::SolidScatraEleCalc<CORE::FE::CellType::hex27,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::hex27>>;
template class DRT::ELEMENTS::SolidScatraEleCalc<CORE::FE::CellType::tet4,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::tet4>>;
template class DRT::ELEMENTS::SolidScatraEleCalc<CORE::FE::CellType::tet10,
    DRT::ELEMENTS::DisplacementBasedFormulation<CORE::FE::CellType::tet10>>;

FOUR_C_NAMESPACE_CLOSE