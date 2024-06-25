/*! \file

\brief Implementation of routines for calculation of a coupled solid-scatra element with templated
solid formulation

\level 1
*/

#include "4C_solid_scatra_3D_ele_calc.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_formulation.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename T>
  T* get_ptr(std::optional<T>& opt)
  {
    return opt.has_value() ? &opt.value() : nullptr;
  }
  template <typename T>
  const T* get_data(const std::optional<std::vector<T>>& opt)
  {
    return opt.has_value() ? opt.value().data() : nullptr;
  }

  template <Core::FE::CellType celltype>
  inline static constexpr int num_str = Core::FE::dim<celltype>*(Core::FE::dim<celltype> + 1) / 2;

  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<num_str<celltype>, 1> EvaluateDMaterialStressDScalar(
      Mat::So3Material& solid_material,
      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          deformation_gradient,
      const Core::LinAlg::Matrix<num_str<celltype>, 1>& gl_strain, Teuchos::ParameterList& params,
      const int gp, const int eleGID)
  {
    Core::LinAlg::Matrix<num_str<celltype>, 1> dStressDScalar(true);

    // The derivative of the solid stress w.r.t. the scalar is implemented in the normal material
    // Evaluate call by not passing the linearization matrix.
    solid_material.evaluate(
        &deformation_gradient, &gl_strain, params, &dStressDScalar, nullptr, gp, eleGID);

    return dStressDScalar;
  }

  template <Core::FE::CellType celltype>
  auto ProjectQuantityToGaussPoint(
      const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, 1>>& nodal_quantities)
  {
    std::vector<double> quantities_at_gp(nodal_quantities.size(), 0.0);

    for (std::size_t k = 0; k < nodal_quantities.size(); ++k)
    {
      quantities_at_gp[k] = shape_functions.shapefunctions_.dot(nodal_quantities[k]);
    }
    return quantities_at_gp;
  }

  template <Core::FE::CellType celltype>
  auto ProjectQuantityToGaussPoint(
      const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, 1>& nodal_quantity)
  {
    return shape_functions.shapefunctions_.dot(nodal_quantity);
  }


  template <bool is_scalar, Core::FE::CellType celltype>
  auto GetElementQuantities(const int num_scalars, const std::vector<double>& quantities_at_dofs)
  {
    if constexpr (is_scalar)
    {
      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, 1> nodal_quantities(num_scalars);
      for (int i = 0; i < Core::FE::num_nodes<celltype>; ++i)
        nodal_quantities(i, 0) = quantities_at_dofs.at(i);

      return nodal_quantities;
    }
    else
    {
      std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, 1>> nodal_quantities(
          num_scalars);

      for (int k = 0; k < num_scalars; ++k)
        for (int i = 0; i < Core::FE::num_nodes<celltype>; ++i)
          (nodal_quantities[k])(i, 0) = quantities_at_dofs.at(num_scalars * i + k);

      return nodal_quantities;
    }
  }

  template <bool is_scalar, Core::FE::CellType celltype>
  void PrepareScatraQuantityInParameterList(const Core::FE::Discretization& discretization,
      const Core::Elements::Element::LocationArray& la,
      const Discret::ELEMENTS::ElementNodes<celltype>& element_nodes, const std::string& field_name,
      const int field_index, const int num_scalars,
      const Core::FE::GaussIntegration& gauss_integration, Teuchos::ParameterList& params,
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
    Core::FE::ExtractMyValues(*quantitites_np, my_quantities, la[field_index].lm_);

    // get nodal quantities for the scalars
    auto nodal_quantities = GetElementQuantities<is_scalar, celltype>(num_scalars, my_quantities);

    // Determine quantity type at Gauss point (can be double (for scalar) or a std::vector of
    // scalars)
    using gp_quantity_type = decltype(ProjectQuantityToGaussPoint(
        std::declval<const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>&>(),
        nodal_quantities));

    // the material expects the gp-quantities at the Gauss points in a rcp-std::vector
    auto quantity_at_gp =
        Teuchos::rcp(new std::vector<gp_quantity_type>(gauss_integration.NumPoints()));

    Discret::ELEMENTS::ForEachGaussPoint(element_nodes, gauss_integration,
        [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1>& xi,
            const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const Discret::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
            double integration_factor, int gp)
        {
          // Project to Gauss point
          (*quantity_at_gp)[gp] = ProjectQuantityToGaussPoint(shape_functions, nodal_quantities);
        });

    params.set(target_name, quantity_at_gp);
  };

  template <Core::FE::CellType celltype>
  void prepare_scatra_quantities_in_parameter_list(const Core::Elements::Element& element,
      const Core::FE::Discretization& discretization,
      const Core::Elements::Element::LocationArray& la,
      const Discret::ELEMENTS::ElementNodes<celltype>& element_nodes,
      const Core::FE::GaussIntegration& gauss_integration, Teuchos::ParameterList& params)
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



  template <Core::FE::CellType celltype, typename SolidFormulation>
  double EvaluateCauchyNDirAtXi(Mat::So3Material& mat,
      Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype> shape_functions,
      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          deformation_gradient,
      const std::optional<std::vector<double>>& scalars_at_xi, const Core::LinAlg::Matrix<3, 1>& n,
      const Core::LinAlg::Matrix<3, 1>& dir, int eleGID,
      const Discret::ELEMENTS::ElementFormulationDerivativeEvaluator<celltype, SolidFormulation>&
          evaluator,
      Discret::ELEMENTS::SolidScatraCauchyNDirLinearizations<3>& linearizations)
  {
    Discret::ELEMENTS::CauchyNDirLinearizationDependencies<celltype> linearization_dependencies =
        Discret::ELEMENTS::get_initialized_cauchy_n_dir_linearization_dependencies(
            evaluator, linearizations);

    double cauchy_n_dir = 0;
    mat.evaluate_cauchy_n_dir_and_derivatives(deformation_gradient, n, dir, cauchy_n_dir,
        linearizations.solid.d_cauchyndir_dn, linearizations.solid.d_cauchyndir_ddir,
        get_ptr(linearization_dependencies.d_cauchyndir_dF),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF2),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF_dn),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF_ddir), -1, eleGID,
        get_data(scalars_at_xi), nullptr, nullptr, nullptr);

    // Evaluate pure solid linearizations
    Discret::ELEMENTS::evaluate_cauchy_n_dir_linearizations<celltype>(
        linearization_dependencies, linearizations.solid);

    // Evaluate ssi-linearizations
    if (linearizations.d_cauchyndir_ds)
    {
      FOUR_C_ASSERT(linearization_dependencies.d_cauchyndir_dF, "Not all tensors are computed!");
      FOUR_C_ASSERT(
          scalars_at_xi.has_value(), "Scalar needs to have a value if the derivatives are needed!");
      linearizations.d_cauchyndir_ds->shape(Core::FE::num_nodes<celltype>, 1);

      static Core::LinAlg::Matrix<9, 1> d_F_dc(true);
      mat.evaluate_linearization_od(deformation_gradient, (*scalars_at_xi)[0], &d_F_dc);

      double d_cauchyndir_ds_gp = (*linearization_dependencies.d_cauchyndir_dF).dot(d_F_dc);

      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, 1>(
          linearizations.d_cauchyndir_ds->values(), true)
          .update(d_cauchyndir_ds_gp, shape_functions.shapefunctions_, 1.0);
    }
    return cauchy_n_dir;
  }
}  // namespace

template <Core::FE::CellType celltype, typename SolidFormulation>
Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::SolidScatraEleCalc()
    : stiffness_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_stiffness_matrix<celltype>())),
      mass_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_mass_matrix<celltype>()))
{
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::pack(
    Core::Communication::PackBuffer& data) const
{
  Discret::ELEMENTS::pack(data, history_data_);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  Discret::ELEMENTS::unpack(position, data, history_data_);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
    Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseVector* force_vector,
    Core::LinAlg::SerialDenseMatrix* stiffness_matrix, Core::LinAlg::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff{};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass{};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, 1>> force{};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  prepare_scatra_quantities_in_parameter_list(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  bool equal_integration_mass_stiffness =
      compare_gauss_integration(mass_matrix_integration_, stiffness_matrix_integration_);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  double element_mass = 0.0;
  double element_volume = 0.0;
  ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        evaluate_gp_coordinates_and_add_to_parameter_list(
            nodal_coordinates, shape_functions, params);
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                    deformation_gradient,
                const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              const Stress<celltype> stress = evaluate_material_stress<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              if (force.has_value())
              {
                add_internal_force_vector(linearization, stress, integration_factor,
                    preparation_data, history_data_, gp, *force);
              }

              if (stiff.has_value())
              {
                add_stiffness_matrix(linearization, jacobian_mapping, stress, integration_factor,
                    preparation_data, history_data_, gp, *stiff);
              }

              if (mass.has_value())
              {
                if (equal_integration_mass_stiffness)
                {
                  add_mass_matrix(
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
        [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp) {
          add_mass_matrix(
              shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::evaluate_d_stress_d_scalar(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& stiffness_matrix_dScalar)
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
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  prepare_scatra_quantities_in_parameter_list(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  // Check for negative Jacobian determinants
  ensure_positive_jacobian_determinant_at_element_nodes(nodal_coordinates);

  ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        evaluate_gp_coordinates_and_add_to_parameter_list(
            nodal_coordinates, shape_functions, params);
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                    deformation_gradient,
                const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              Core::LinAlg::Matrix<6, 1> dSdc = EvaluateDMaterialStressDScalar<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              // linear B-opeartor
              const Core::LinAlg::Matrix<Details::num_str<celltype>,
                  Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
                  bop = SolidFormulation::get_linear_b_operator(linearization);

              constexpr int num_dof_per_ele =
                  Core::FE::dim<celltype> * Core::FE::num_nodes<celltype>;

              // Assemble matrix
              // k_dS = B^T . dS/dc * detJ * N * w(gp)
              Core::LinAlg::Matrix<num_dof_per_ele, 1> BdSdc(true);
              BdSdc.multiply_tn(integration_factor, bop, dSdc);

              // loop over rows
              for (int rowi = 0; rowi < num_dof_per_ele; ++rowi)
              {
                const double BdSdc_rowi = BdSdc(rowi, 0);
                // loop over columns
                for (int coli = 0; coli < Core::FE::num_nodes<celltype>; ++coli)
                {
                  stiffness_matrix_dScalar(rowi, coli * scatra_column_stride) +=
                      BdSdc_rowi * shape_functions.shapefunctions_(coli, 0);
                }
              }
            });
      });
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Recover(
    const Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  // nothing needs to be done for simple displacement based elements
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::Update(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  prepare_scatra_quantities_in_parameter_list(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  Discret::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        evaluate_gp_coordinates_and_add_to_parameter_list(
            nodal_coordinates, shape_functions, params);
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                    deformation_gradient,
                const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            { solid_material.update(deformation_gradient, gp, params, ele.Id()); });
      });

  solid_material.update();
}

template <Core::FE::CellType celltype, typename SolidFormulation>
double Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::calculate_internal_energy(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  prepare_scatra_quantities_in_parameter_list(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  double intenergy = 0;
  Discret::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        evaluate_gp_coordinates_and_add_to_parameter_list(
            nodal_coordinates, shape_functions, params);
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                    deformation_gradient,
                const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              double psi = 0.0;
              solid_material.StrainEnergy(gl_strain, psi, gp, ele.Id());
              intenergy += psi * integration_factor;
            });
      });

  return intenergy;
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::calculate_stress(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material, const StressIO& stressIO,
    const StrainIO& strainIO, const Core::FE::Discretization& discretization,
    const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  Core::LinAlg::SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), num_str_);
  Core::LinAlg::SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), num_str_);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // prepare scatra data in the parameter list
  prepare_scatra_quantities_in_parameter_list(
      ele, discretization, la, nodal_coordinates, stiffness_matrix_integration_, params);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      Prepare(ele, nodal_coordinates, history_data_);

  Discret::ELEMENTS::ForEachGaussPoint(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        evaluate_gp_coordinates_and_add_to_parameter_list(
            nodal_coordinates, shape_functions, params);
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                    deformation_gradient,
                const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
            {
              const Stress<celltype> stress = evaluate_material_stress<celltype>(
                  solid_material, deformation_gradient, gl_strain, params, gp, ele.Id());

              assemble_strain_type_to_matrix_row<celltype>(
                  gl_strain, deformation_gradient, strainIO.type, strain_data, gp);
              assemble_stress_type_to_matrix_row(
                  deformation_gradient, stress, stressIO.type, stress_data, gp);
            });
      });

  serialize(stress_data, serialized_stress_data);
  serialize(strain_data, serialized_strain_data);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
double
Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::get_normal_cauchy_stress_at_xi(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const std::vector<double>& disp, const std::optional<std::vector<double>>& scalars,
    const Core::LinAlg::Matrix<3, 1>& xi, const Core::LinAlg::Matrix<3, 1>& n,
    const Core::LinAlg::Matrix<3, 1>& dir, SolidScatraCauchyNDirLinearizations<3>& linearizations)
{
  if constexpr (has_gauss_point_history<SolidFormulation>)
  {
    FOUR_C_THROW(
        "Cannot evaluate the Cauchy stress at xi with an element formulation with Gauss point "
        "history. The element formulation is %s.",
        Core::UTILS::TryDemangle(typeid(SolidFormulation).name()).c_str());
  }
  else
  {
    // project scalar values to xi
    const auto scalar_values_at_xi = std::invoke(
        [&]() -> std::optional<std::vector<double>>
        {
          if (!scalars.has_value()) return std::nullopt;

          return Discret::ELEMENTS::ProjectNodalQuantityToXi<celltype>(xi, *scalars);
        });

    ElementNodes<celltype> element_nodes = evaluate_element_nodes<celltype>(ele, disp);

    const ShapeFunctionsAndDerivatives<celltype> shape_functions =
        EvaluateShapeFunctionsAndDerivs<celltype>(xi, element_nodes);

    const JacobianMapping<celltype> jacobian_mapping =
        EvaluateJacobianMapping(shape_functions, element_nodes);

    const PreparationData<SolidFormulation> preparation_data =
        Prepare(ele, element_nodes, history_data_);

    return evaluate(ele, element_nodes, xi, shape_functions, jacobian_mapping, preparation_data,
        history_data_,
        [&](const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                deformation_gradient,
            const Core::LinAlg::Matrix<num_str_, 1>& gl_strain, const auto& linearization)
        {
          const ElementFormulationDerivativeEvaluator<celltype, SolidFormulation> evaluator(ele,
              element_nodes, xi, shape_functions, jacobian_mapping, deformation_gradient,
              preparation_data, history_data_);

          return EvaluateCauchyNDirAtXi<celltype>(solid_material, shape_functions,
              deformation_gradient, scalar_values_at_xi, n, dir, ele.Id(), evaluator,
              linearizations);
        });
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::setup(
    Mat::So3Material& solid_material, Input::LineDefinition* linedef)
{
  solid_material.setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::material_post_setup(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  InterpolateFibersToGaussPointsAndAddToParameterList<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call post_setup of material
  solid_material.post_setup(params, ele.Id());
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::initialize_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  ask_and_add_quantities_to_gauss_point_data_output(
      stiffness_matrix_integration_.NumPoints(), solid_material, gp_data_output_manager);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  collect_and_assemble_gauss_point_data_output<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::ELEMENTS::SolidScatraEleCalc<celltype, SolidFormulation>::reset_to_last_converged(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <Core::FE::CellType... celltypes>
struct VerifyPackable
{
  static constexpr bool are_all_packable =
      (Discret::ELEMENTS::IsPackable<Discret::ELEMENTS::SolidScatraEleCalc<celltypes,
              Discret::ELEMENTS::DisplacementBasedFormulation<celltypes>>*> &&
          ...);

  static constexpr bool are_all_unpackable =
      (Discret::ELEMENTS::IsUnpackable<Discret::ELEMENTS::SolidScatraEleCalc<celltypes,
              Discret::ELEMENTS::DisplacementBasedFormulation<celltypes>>*> &&
          ...);

  void StaticAsserts() const
  {
    static_assert(are_all_packable);
    static_assert(are_all_unpackable);
  }
};

template struct VerifyPackable<Core::FE::CellType::hex8, Core::FE::CellType::hex27,
    Core::FE::CellType::tet4, Core::FE::CellType::tet10>;

// explicit instantiations of template classes
// for displacement based formulation
template class Discret::ELEMENTS::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::ELEMENTS::DisplacementBasedFormulation<Core::FE::CellType::hex8>>;
template class Discret::ELEMENTS::SolidScatraEleCalc<Core::FE::CellType::hex27,
    Discret::ELEMENTS::DisplacementBasedFormulation<Core::FE::CellType::hex27>>;
template class Discret::ELEMENTS::SolidScatraEleCalc<Core::FE::CellType::tet4,
    Discret::ELEMENTS::DisplacementBasedFormulation<Core::FE::CellType::tet4>>;
template class Discret::ELEMENTS::SolidScatraEleCalc<Core::FE::CellType::tet10,
    Discret::ELEMENTS::DisplacementBasedFormulation<Core::FE::CellType::tet10>>;

// FBar based formulation
template class Discret::ELEMENTS::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::ELEMENTS::FBarFormulation<Core::FE::CellType::hex8>>;

FOUR_C_NAMESPACE_CLOSE