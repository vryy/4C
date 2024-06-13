/*! \file

\brief Implementation of routines for calculation of solidporo element with pressure based
implementation

\level 1
*/
#include "4C_solid_poro_3D_ele_calc_pressure_based.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_poro_3D_ele_calc_lib.hpp"

#include <optional>

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType celltype>
Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>::SolidPoroPressureBasedEleCalc()
    : gauss_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_stiffness_matrix_poro<celltype>()))
{
}

template <Core::FE::CellType celltype>
void Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>::poro_setup(
    Mat::StructPoro& porostructmat, Input::LineDefinition* linedef)
{
  porostructmat.poro_setup(gauss_integration_.NumPoints(), linedef);
}

template <Core::FE::CellType celltype>
void Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>::evaluate_nonlinear_force_stiffness(
    const Core::Elements::Element& ele, Mat::StructPoro& porostructmat,
    Mat::FluidPoroMultiPhase& porofluidmat, const Inpar::STR::KinemType& kinematictype,
    const Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector* force_vector,
    Core::LinAlg::SerialDenseMatrix* stiffness_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<Core::LinAlg::Matrix<num_dim_ * num_nodes_, num_dim_* num_nodes_>> stiff = {};
  std::optional<Core::LinAlg::Matrix<num_dim_ * num_nodes_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  // get primary variables of multiphase porous medium flow
  std::vector<double> fluidmultiphase_ephi(la[1].Size());
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
  Core::FE::ExtractMyValues(*matrix_state, fluidmultiphase_ephi, la[1].lm_);

  // Initialize variables of multiphase porous medium flow
  const int nummultifluiddofpernode = porofluidmat.NumMat();
  const int numfluidphases = porofluidmat.NumFluidPhases();
  const int numvolfrac = porofluidmat.NumVolFrac();
  const bool hasvolfracs = (nummultifluiddofpernode > numfluidphases);

  // get nodal coordinates current and reference
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // Loop over all Gauss points
  ForEachGaussPoint(nodal_coordinates, gauss_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const SpatialMaterialMapping<celltype> spatial_material_mapping =
            evaluate_spatial_material_mapping<celltype>(
                jacobian_mapping, nodal_coordinates, 1.0, kinematictype);

        const CauchyGreenAndInverse<celltype> cauchygreen =
            evaluate_cauchy_green_and_inverse(spatial_material_mapping);

        Core::LinAlg::Matrix<num_str_, num_dof_per_ele_> Bop =
            evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

        Core::LinAlg::Matrix<num_str_, num_dof_per_ele_> dInverseRightCauchyGreen_dDisp =
            EvaluateInverseCauchyGreenLinearization(
                cauchygreen, jacobian_mapping, spatial_material_mapping);

        const double volchange = ComputeVolumeChange<celltype>(spatial_material_mapping,
            jacobian_mapping, ele, discretization, la[0].lm_, kinematictype);

        Core::LinAlg::Matrix<1, num_dof_per_ele_> dDetDefGrad_dDisp =
            ComputeLinearizationOfDetDefGradWrtDisp<celltype>(
                spatial_material_mapping, jacobian_mapping, kinematictype);

        const Core::LinAlg::Matrix<1, num_dof_per_ele_> dVolchange_dDisp =
            ComputeLinearizationOfVolchangeWrtDisp<celltype>(
                dDetDefGrad_dDisp, jacobian_mapping, kinematictype);

        std::vector<double> fluidmultiphase_phiAtGP =
            ComputeFluidMultiPhasePrimaryVariablesAtGP<celltype>(
                fluidmultiphase_ephi, nummultifluiddofpernode, shape_functions);

        double solidpressure = compute_sol_pressure_at_gp<celltype>(
            nummultifluiddofpernode, numfluidphases, fluidmultiphase_phiAtGP, porofluidmat);
        // derivative of press w.r.t. displacements (only in case of volfracs)
        Core::LinAlg::Matrix<1, num_dof_per_ele_> dSolidpressure_dDisp(true);

        if (hasvolfracs)
        {
          Core::LinAlg::Matrix<1, num_dof_per_ele_> dPorosity_dDisp;
          double porosity = 0.0;

          compute_porosity_and_linearization<celltype>(porostructmat, params, solidpressure, gp,
              volchange, porosity, dDetDefGrad_dDisp, dPorosity_dDisp);

          // save the pressure coming from the fluid S_i*p_i (old solidpressure, without accounting
          // for volfracs)
          const double fluidpress = solidpressure;

          solidpressure = recalculate_sol_pressure_at_gp(fluidpress, porosity,
              nummultifluiddofpernode, numfluidphases, numvolfrac, fluidmultiphase_phiAtGP);

          RecalculateLinearizationOfSolPressWrtDisp<celltype>(fluidpress, porosity,
              nummultifluiddofpernode, numfluidphases, numvolfrac, fluidmultiphase_phiAtGP,
              dPorosity_dDisp, dSolidpressure_dDisp);
        }

        // inverse Right Cauchy-Green tensor as vector in voigt notation
        Core::LinAlg::Matrix<num_str_, 1> C_inv_vec(false);
        Core::LinAlg::Voigt::Stresses::matrix_to_vector(
            cauchygreen.inverse_right_cauchy_green_, C_inv_vec);

        // B^T . C^-1
        Core::LinAlg::Matrix<num_dof_per_ele_, 1> BopCinv(true);
        BopCinv.MultiplyTN(Bop, C_inv_vec);

        // update internal force vector
        if (force.has_value())
        {
          UpdateInternalForceVectorMultiPhasePressureBased<celltype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, BopCinv,
              *force);
        }

        // update stiffness matrix
        if (stiff.has_value())
        {
          UpdateElasticStiffnessMatrixMultiPhasePressureBased<celltype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, BopCinv,
              Bop, dDetDefGrad_dDisp, dSolidpressure_dDisp, dInverseRightCauchyGreen_dDisp, *stiff);

          UpdateGeometricStiffnessMatrixMultiPhasePressureBased<celltype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, C_inv_vec,
              jacobian_mapping.N_XYZ_, *stiff);
        }
      });
}

template <Core::FE::CellType celltype>
void Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>::coupling_poroelast(
    const Core::Elements::Element& ele, Mat::StructPoro& porostructmat,
    Mat::FluidPoroMultiPhase& porofluidmat, const Inpar::STR::KinemType& kinematictype,
    const Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Teuchos::ParameterList& params, Core::LinAlg::SerialDenseMatrix& stiffness_matrix)
{
  // get primary variables of multiphase porous medium flow
  std::vector<double> fluidmultiphase_ephi(la[1].Size());
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
  Core::FE::ExtractMyValues(*matrix_state, fluidmultiphase_ephi, la[1].lm_);

  // Initialize variables of multiphase porous medium flow
  const int nummultifluiddofpernode = porofluidmat.NumMat();
  const int numfluidphases = porofluidmat.NumFluidPhases();
  const int numvolfrac = porofluidmat.NumVolFrac();
  const bool hasvolfracs = (nummultifluiddofpernode > numfluidphases);

  // get nodal coordinates current and reference
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  // Loop over all Gauss points
  ForEachGaussPoint(nodal_coordinates, gauss_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp

      )
      {
        const SpatialMaterialMapping<celltype> spatial_material_mapping =
            evaluate_spatial_material_mapping<celltype>(
                jacobian_mapping, nodal_coordinates, 1.0, kinematictype);

        const CauchyGreenAndInverse<celltype> cauchygreen =
            evaluate_cauchy_green_and_inverse(spatial_material_mapping);

        Core::LinAlg::Matrix<num_str_, num_dof_per_ele_> Bop =
            evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

        // volume change (used for porosity law). Same as J in nonlinear theory.
        const double volchange = ComputeVolumeChange<celltype>(spatial_material_mapping,
            jacobian_mapping, ele, discretization, la[0].lm_, kinematictype);

        std::vector<double> fluidmultiphase_phiAtGP =
            ComputeFluidMultiPhasePrimaryVariablesAtGP<celltype>(
                fluidmultiphase_ephi, nummultifluiddofpernode, shape_functions);

        std::vector<double> solidpressurederiv = ComputeSolidPressureDeriv<celltype>(
            porofluidmat, fluidmultiphase_phiAtGP, numfluidphases);

        if (hasvolfracs)
        {
          double solidpressure = compute_sol_pressure_at_gp<celltype>(
              nummultifluiddofpernode, numfluidphases, fluidmultiphase_phiAtGP, porofluidmat);

          double porosity =
              compute_porosity<celltype>(porostructmat, params, solidpressure, volchange, gp);

          recalculate_sol_pressure_deriv(fluidmultiphase_phiAtGP, nummultifluiddofpernode,
              numfluidphases, numvolfrac, solidpressure, porosity, solidpressurederiv);
        }

        const double detJ_w = jacobian_mapping.determinant_ * gauss_integration_.Weight(gp);

        // inverse Right Cauchy-Green tensor as vector in voigt notation
        Core::LinAlg::Matrix<num_str_, 1> C_inv_vec(false);
        Core::LinAlg::Voigt::Stresses::matrix_to_vector(
            cauchygreen.inverse_right_cauchy_green_, C_inv_vec);

        // B^T . C^-1
        Core::LinAlg::Matrix<num_dof_per_ele_, 1> BopCinv(true);
        BopCinv.MultiplyTN(Bop, C_inv_vec);

        UpdateStiffnessMatrixCouplingMultiPhasePressureBased<celltype>(detJ_w, solidpressurederiv,
            BopCinv, shape_functions, spatial_material_mapping.determinant_deformation_gradient_,
            nummultifluiddofpernode, stiffness_matrix);
      });
}

template <Core::FE::CellType celltype>
void Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>::coupling_stress(
    const Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  auto iocouplingstress = Core::UTILS::GetAsEnum<Inpar::STR::StressType>(
      params, "iocouplstress", Inpar::STR::stress_none);

  // check for output of coupling stress
  if (iocouplingstress == Inpar::STR::stress_none)
  {
    // nothing to do for calculation of effective stress
    return;
  }
  else
  {
    FOUR_C_THROW("coupling stress poroelast not yet implemented for pressure-based variant");
  }
}

// template classes
template class Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<Core::FE::CellType::tet10>;

FOUR_C_NAMESPACE_CLOSE
