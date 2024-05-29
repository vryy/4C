/*! \file

\brief Implementation of routines for calculation of shell element simple displacement based

\level 3
*/

#include "4C_shell7p_ele_calc.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_gauss_point_extrapolation.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_shell7p_ele.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function_of_time.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


template <CORE::FE::CellType distype>
DRT::ELEMENTS::Shell7pEleCalc<distype>::Shell7pEleCalc()
    : DRT::ELEMENTS::Shell7pEleCalcInterface::Shell7pEleCalcInterface(),
      intpoints_midsurface_(
          SHELL::CreateGaussIntegrationPoints<distype>(SHELL::get_gauss_rule<distype>()))
{
  cur_thickness_.resize(intpoints_midsurface_.NumPoints(), shell_data_.thickness);
}


template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::Setup(CORE::Elements::Element& ele,
    MAT::So3Material& solid_material, INPUT::LineDefinition* linedef,
    const STR::ELEMENTS::ShellLockingTypes& locking_types,
    const STR::ELEMENTS::ShellData& shell_data)
{
  shell_data_ = shell_data;
  // initialize current thickness at all gp
  cur_thickness_.resize(intpoints_midsurface_.NumPoints(), shell_data_.thickness);
  //  set up of materials with GP data (e.g., history variables)
  solid_material.Setup(intpoints_midsurface_.NumPoints(), linedef);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.sdc);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.thickness);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.num_ans);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, cur_thickness_);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  CORE::COMM::ParObject::ExtractfromPack(position, data, shell_data_.sdc);
  CORE::COMM::ParObject::ExtractfromPack(position, data, shell_data_.thickness);
  CORE::COMM::ParObject::ExtractfromPack(position, data, shell_data_.num_ans);
  CORE::COMM::ParObject::ExtractfromPack(position, data, cur_thickness_);
}


template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::material_post_setup(
    CORE::Elements::Element& ele, MAT::So3Material& solid_material)
{
  // element/nodal wise defined data
  Teuchos::ParameterList params{};

  // Call post_setup of material
  solid_material.post_setup(params, ele.Id());
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::reset_to_last_converged(
    CORE::Elements::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <CORE::FE::CellType distype>
double DRT::ELEMENTS::Shell7pEleCalc<distype>::calculate_internal_energy(
    CORE::Elements::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  // need update
  double intenergy = 0.0;

  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  if (disp == Teuchos::null || res == Teuchos::null)
    FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
  std::vector<double> displacement(dof_index_array.size());
  CORE::FE::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  CORE::FE::ExtractMyValues(*res, residual, dof_index_array);

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;

  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  std::vector<double> shape_functions_ans(true);

  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;

  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  // init enhanced strains for shell
  CORE::LINALG::SerialDenseVector strain_enh(SHELL::DETAIL::num_internal_variables);


  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        double integration_factor = gpweight * da;

        const std::vector<double> shape_functions_ans =
            SHELL::GetShapefunctionsForAns<distype>(xi_gp, shell_data_.num_ans);

        // integration loop in thickness direction, here we prescribe 2 integration points
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;

          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // modify the current kovariant metric tensor to neglect the quadratic terms in thickness
          // directions
          SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta,
              shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
              shell_data_.num_ans);

          // evaluate Green-Lagrange strains and deformation gradient in cartesian coordinate system
          auto strains = EvaluateStrains(g_reference, g_current);

          // call material for evaluation of strain energy function
          double psi = 0.0;
          solid_material.StrainEnergy(strains.gl_strain_, psi, gp, ele.Id());

          double thickness = 0.0;
          for (int i = 0; i < SHELL::DETAIL::num_node<distype>; ++i)
            thickness += thickness * shape_functions.shapefunctions_(i);

          intenergy += psi * integration_factor * 0.5 * thickness;
        }
      });

  return intenergy;
}


template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::calculate_stresses_strains(
    CORE::Elements::Element& ele, MAT::So3Material& solid_material, const ShellStressIO& stressIO,
    const ShellStrainIO& strainIO, const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(
      intpoints_midsurface_.NumPoints(), MAT::NUM_STRESS_3D);
  CORE::LINALG::SerialDenseMatrix strain_data(
      intpoints_midsurface_.NumPoints(), MAT::NUM_STRESS_3D);


  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  if (disp == Teuchos::null || res == Teuchos::null)
    FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
  std::vector<double> displacement(dof_index_array.size());
  CORE::FE::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  CORE::FE::ExtractMyValues(*res, residual, dof_index_array);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;

  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // metric of element at centroid point (for EAS)
  std::array<double, 2> centroid_point = {0.0, 0.0};
  SHELL::ShapefunctionsAndDerivatives<distype> shapefunctions_centroid =
      SHELL::EvaluateShapefunctionsAndDerivs<distype>(centroid_point);
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_reference;
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_current;

  SHELL::EvaluateMetrics(shapefunctions_centroid, metrics_centroid_reference,
      metrics_centroid_current, nodal_coordinates, 0.0);

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;
  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  // init enhanced strains for shell:
  CORE::LINALG::SerialDenseVector strain_enh(SHELL::DETAIL::num_internal_variables);

  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        const std::vector<double> shape_functions_ans =
            SHELL::GetShapefunctionsForAns<distype>(xi_gp, shell_data_.num_ans);

        // integration loop in thickness direction, here we prescribe 2 integration points
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;
          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // modify the current kovariant metric tensor to neglect the quadratic terms in thickness
          // directions
          SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta,
              shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
              shell_data_.num_ans);

          // evaluate Green-Lagrange strains and deformation gradient in global cartesian coordinate
          // system
          auto strains = SHELL::EvaluateStrains(g_reference, g_current);

          // update the deformation gradient
          CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> defgrd_enh(false);
          SHELL::calc_consistent_defgrd<SHELL::DETAIL::num_dim>(
              strains.defgrd_, strains.gl_strain_, defgrd_enh);
          strains.defgrd_ = defgrd_enh;

          // evaluate stress in global cartesian system
          auto stress = SHELL::EvaluateMaterialStressCartesianSystem<SHELL::DETAIL::num_dim>(
              solid_material, strains, params, gp, ele.Id());
          SHELL::AssembleStrainTypeToMatrixRow<distype>(
              strains, strainIO.type, strain_data, gp, 0.5);
          SHELL::AssembleStressTypeToMatrixRow<distype>(
              strains, stress, stressIO.type, stress_data, gp, 0.5);
        }
      });
  SHELL::Serialize(stress_data, serialized_stress_data);
  SHELL::Serialize(strain_data, serialized_strain_data);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::evaluate_nonlinear_force_stiffness_mass(
    CORE::Elements::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  if (disp == Teuchos::null || res == Teuchos::null)
    FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
  std::vector<double> displacement(dof_index_array.size());
  CORE::FE::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  CORE::FE::ExtractMyValues(*res, residual, dof_index_array);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;


  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;

  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  // init enhanced strain for shell
  constexpr auto num_internal_variables = SHELL::DETAIL::num_internal_variables;
  CORE::LINALG::SerialDenseVector strain_enh(num_internal_variables);
  SHELL::StressEnhanced stress_enh;

  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        double integration_factor = gpweight * da;

        // update current thickness at gauss point
        cur_thickness_[gp] = SHELL::UpdateGaussPointThickness<distype>(
            nodal_coordinates.a3_curr_, shape_functions.shapefunctions_);

        // reset mid-surface material tensor and stress resultants to zero
        stress_enh.dmat_.shape(num_internal_variables, num_internal_variables);
        stress_enh.stress_.size(num_internal_variables);

        // init mass matrix variables
        SHELL::MassMatrixVariables mass_matrix_variables;

        // calculate B-operator for compatible strains (displacement)
        CORE::LINALG::SerialDenseMatrix Bop = SHELL::CalcBOperator<distype>(
            a_current.kovariant_, a_current.partial_derivative_, shape_functions);

        const std::vector<double> shape_functions_ans =
            SHELL::GetShapefunctionsForAns<distype>(xi_gp, shell_data_.num_ans);

        // modifications due to ANS with B-bar method (Hughes (1980))
        std::invoke(
            [&]()
            {
              if (shell_data_.num_ans > 0)
              {
                ModifyBOperatorAns(Bop, shape_functions_ans, shapefunctions_collocation,
                    metrics_collocation_current, shell_data_.num_ans);
              }
            });


        // integration loop in thickness direction, here we prescribe 2 integration points
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;
          double factor = intpoints_thickness_.qwgt[gpt];

          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // modify the current kovariant metric tensor to neglect the quadratic terms in thickness
          // directions
          SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta,
              shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
              shell_data_.num_ans);

          // calc shell shifter and put it in the integration factor
          const double shifter = (1.0 / condfac) * (g_reference.detJ_ / da);
          factor *= shifter;

          // evaluate Green-Lagrange strains and deformation gradient in cartesian coordinate system
          auto strains = EvaluateStrains(g_reference, g_current);

          // update the deformation gradient (if needed?)
          if (solid_material.needs_defgrd())
          {
            CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> defgrd_enh(false);
            SHELL::calc_consistent_defgrd<SHELL::DETAIL::num_dim>(
                strains.defgrd_, strains.gl_strain_, defgrd_enh);
            strains.defgrd_ = defgrd_enh;
          }

          auto stress = SHELL::EvaluateMaterialStressCartesianSystem<SHELL::DETAIL::num_dim>(
              solid_material, strains, params, gp, ele.Id());
          SHELL::MapMaterialStressToCurvilinearSystem(stress, g_reference);
          SHELL::ThicknessIntegration<distype>(stress_enh, stress, factor, zeta);

          // thickness integration of mass matrix variables
          if (mass_matrix != nullptr)
          {
            double tmp_integration_factor = intpoints_thickness_.qwgt[gpt] * g_reference.detJ_;
            mass_matrix_variables.factor_v_ += tmp_integration_factor;
            mass_matrix_variables.factor_w_ += tmp_integration_factor *
                                               intpoints_thickness_.qxg[gpt][0] *
                                               intpoints_thickness_.qxg[gpt][0];
            mass_matrix_variables.factor_vw_ +=
                tmp_integration_factor * intpoints_thickness_.qxg[gpt][0];
          }
        }
        // add stiffness matrix
        if (stiffness_matrix != nullptr)
        {
          // elastic stiffness matrix Ke
          SHELL::AddElasticStiffnessMatrix<distype>(
              Bop, stress_enh.dmat_, integration_factor, *stiffness_matrix);
          // geometric stiffness matrix Kg
          SHELL::AddGeometricStiffnessMatrix<distype>(shapefunctions_collocation,
              shape_functions_ans, shape_functions, stress_enh.stress_, shell_data_.num_ans,
              integration_factor, *stiffness_matrix);
          // make stiffness matrix absolute symmetric
          for (int i = 0; i < SHELL::DETAIL::numdofperelement<distype>; ++i)
          {
            for (int j = i + 1; j < SHELL::DETAIL::numdofperelement<distype>; ++j)
            {
              const double average = 0.5 * ((*stiffness_matrix)(i, j) + (*stiffness_matrix)(j, i));
              (*stiffness_matrix)(i, j) = average;
              (*stiffness_matrix)(j, i) = average;
            }
          }
        }
        // add internal force vector
        if (force_vector != nullptr)
        {
          SHELL::add_internal_force_vector<distype>(
              Bop, stress_enh.stress_, integration_factor, *force_vector);
        }
        // add internal mass_matrix
        if (mass_matrix != nullptr)
        {
          double density = solid_material.Density(gp);
          mass_matrix_variables.factor_v_ *= gpweight * density;
          mass_matrix_variables.factor_w_ *= gpweight * density;
          mass_matrix_variables.factor_vw_ *= gpweight * density;
          SHELL::AddMassMatrix(
              shape_functions, mass_matrix_variables, shell_data_.thickness, *mass_matrix);
        }
      });
}


template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::Recover(CORE::Elements::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params, STR::ELEMENTS::ParamsInterface& str_interface)
{
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::Update(CORE::Elements::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement' ");
  std::vector<double> displacement(dof_index_array.size());
  CORE::FE::ExtractMyValues(*disp, displacement, dof_index_array);

  // calculate and update inelastic deformation gradient if needed
  if (solid_material.UsesExtendedUpdate())
  {
    const double condfac = shell_data_.sdc;

    // init gauss point in thickness direction that will be modified via SDC
    double zeta = 0.0;

    // get nodal coordinates
    SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
        ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

    // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
    // for a_13 and a_23 each
    const int total_ansq = 2 * shell_data_.num_ans;
    std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(
        total_ansq);
    std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
    std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

    if (shell_data_.num_ans > 0)
    {
      SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
          metrics_collocation_current, nodal_coordinates, total_ansq);
    }
    // init metric tensor and basis vectors of mid-surface
    SHELL::BasisVectorsAndMetrics<distype> a_reference;
    SHELL::BasisVectorsAndMetrics<distype> a_current;

    // init metric tensor and basis vectors of shell body
    SHELL::BasisVectorsAndMetrics<distype> g_reference;
    SHELL::BasisVectorsAndMetrics<distype> g_current;

    // enhanced strains
    CORE::LINALG::SerialDenseVector strain_enh(SHELL::DETAIL::num_internal_variables);

    SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
        [&](const std::array<double, 2>& xi_gp,
            const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
            SHELL::BasisVectorsAndMetrics<distype>& a_current,
            SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
        {
          const std::vector<double> shape_functions_ans =
              SHELL::GetShapefunctionsForAns<distype>(xi_gp, shell_data_.num_ans);

          // integration loop in thickness direction, here we prescribe 2 integration points
          for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
          {
            zeta = intpoints_thickness_.qxg[gpt][0] / condfac;

            SHELL::EvaluateMetrics(
                shape_functions, g_reference, g_current, nodal_coordinates, zeta);

            // modify the current kovariant metric tensor to neglect the quadratic terms in
            // thickness directions
            SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta,
                shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
                shell_data_.num_ans);

            auto strains = EvaluateStrains(g_reference, g_current);

            // calculate deformation gradient consistent with modified GL strain tensor
            if (solid_material.needs_defgrd())
            {
              // update the deformation gradient (if needed?)
              CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> defgrd_enh(
                  false);
              SHELL::calc_consistent_defgrd<SHELL::DETAIL::num_dim>(
                  strains.defgrd_, strains.gl_strain_, defgrd_enh);
              strains.defgrd_ = defgrd_enh;
            }
            solid_material.Update(strains.defgrd_, gp, params, ele.Id());
          }
        });
  }
  solid_material.Update();
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Shell7pEleCalc<distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  if (name == "thickness")
  {
    if (data.size() != 1) FOUR_C_THROW("size mismatch");
    for (auto& thickness_data : cur_thickness_)
    {
      data[0] += thickness_data;
    }
    data[0] = data[0] / intpoints_midsurface_.NumPoints();
  }

}  // VisData()

// template classes
template class DRT::ELEMENTS::Shell7pEleCalc<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::Shell7pEleCalc<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::Shell7pEleCalc<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::Shell7pEleCalc<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::Shell7pEleCalc<CORE::FE::CellType::tri6>;

FOUR_C_NAMESPACE_CLOSE
