/*! \file

\brief Implementation of routines for calculation of solidporo element with pressure based
implementation

\level 1
*/
#include "baci_solid_ele_poro_calc_pressure_based.H"

#include "baci_lib_discret.H"
#include "baci_lib_utils.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_poro_calc_lib.H"

#include <optional>

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>::SolidPoroPressureBasedEleCalc()
    : gauss_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrixPoro<distype>()))
{
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>::PoroSetup(
    MAT::StructPoro& porostructmat, DRT::INPUT::LineDefinition* linedef)
{
  porostructmat.PoroSetup(gauss_integration_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>::EvaluateNonlinearForceStiffness(
    const DRT::Element& ele, MAT::StructPoro& porostructmat, MAT::FluidPoroMultiPhase& porofluidmat,
    const INPAR::STR::KinemType& kinematictype, const DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<CORE::LINALG::Matrix<num_dim_ * num_nodes_, num_dim_* num_nodes_>> stiff = {};
  std::optional<CORE::LINALG::Matrix<num_dim_ * num_nodes_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  // get primary variables of multiphase porous medium flow
  std::vector<double> fluidmultiphase_ephi(la[1].Size());
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
  DRT::UTILS::ExtractMyValues(*matrix_state, fluidmultiphase_ephi, la[1].lm_);

  // Initialize variables of multiphase porous medium flow
  const int nummultifluiddofpernode = porofluidmat.NumMat();
  const int numfluidphases = porofluidmat.NumFluidPhases();
  const int numvolfrac = porofluidmat.NumVolFrac();
  const bool hasvolfracs = (nummultifluiddofpernode > numfluidphases);

  // get nodal coordinates current and reference
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, la[0].lm_);

  // Loop over all Gauss points
  ForEachGaussPoint(nodal_coordinates, gauss_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp

      )
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping<distype>(
                jacobian_mapping, nodal_coordinates, 1.0, kinematictype);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> dInverseRightCauchyGreen_dDisp =
            EvaluateInverseCauchyGreenLinearization(
                cauchygreen, jacobian_mapping, spatial_material_mapping);

        const double volchange = ComputeVolumeChange<distype>(spatial_material_mapping,
            jacobian_mapping, ele, discretization, la[0].lm_, kinematictype);

        CORE::LINALG::Matrix<1, num_dof_per_ele_> dDetDefGrad_dDisp =
            ComputeLinearizationOfDetDefGradWrtDisp<distype>(
                spatial_material_mapping, jacobian_mapping, kinematictype);

        const CORE::LINALG::Matrix<1, num_dof_per_ele_> dVolchange_dDisp =
            ComputeLinearizationOfVolchangeWrtDisp<distype>(
                dDetDefGrad_dDisp, jacobian_mapping, kinematictype);

        std::vector<double> fluidmultiphase_phiAtGP =
            ComputeFluidMultiPhasePrimaryVariablesAtGP<distype>(
                fluidmultiphase_ephi, nummultifluiddofpernode, shape_functions);

        double solidpressure = ComputeSolPressureAtGP<distype>(
            nummultifluiddofpernode, numfluidphases, fluidmultiphase_phiAtGP, porofluidmat);
        // derivative of press w.r.t. displacements (only in case of volfracs)
        CORE::LINALG::Matrix<1, num_dof_per_ele_> dSolidpressure_dDisp(true);

        if (hasvolfracs)
        {
          CORE::LINALG::Matrix<1, num_dof_per_ele_> dPorosity_dDisp;
          double porosity = 0.0;

          ComputePorosityAndLinearization<distype>(porostructmat, params, solidpressure, gp,
              volchange, porosity, dDetDefGrad_dDisp, dPorosity_dDisp);

          // save the pressure coming from the fluid S_i*p_i (old solidpressure, without accounting
          // for volfracs)
          const double fluidpress = solidpressure;

          solidpressure = RecalculateSolPressureAtGP(fluidpress, porosity, nummultifluiddofpernode,
              numfluidphases, numvolfrac, fluidmultiphase_phiAtGP);

          RecalculateLinearizationOfSolPressWrtDisp<distype>(fluidpress, porosity,
              nummultifluiddofpernode, numfluidphases, numvolfrac, fluidmultiphase_phiAtGP,
              dPorosity_dDisp, dSolidpressure_dDisp);
        }

        // inverse Right Cauchy-Green tensor as vector in voigt notation
        CORE::LINALG::Matrix<num_str_, 1> C_inv_vec(false);
        CORE::LINALG::VOIGT::Stresses::MatrixToVector(
            cauchygreen.inverse_right_cauchy_green, C_inv_vec);

        // B^T . C^-1
        CORE::LINALG::Matrix<num_dof_per_ele_, 1> BopCinv(true);
        BopCinv.MultiplyTN(Bop, C_inv_vec);

        // update internal force vector
        if (force.has_value())
        {
          UpdateInternalForceVectorMultiPhasePressureBased<distype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, BopCinv,
              *force);
        }

        // update stiffness matrix
        if (stiff.has_value())
        {
          UpdateElasticStiffnessMatrixMultiPhasePressureBased<distype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, BopCinv,
              Bop, dDetDefGrad_dDisp, dSolidpressure_dDisp, dInverseRightCauchyGreen_dDisp, *stiff);

          UpdateGeometricStiffnessMatrixMultiPhasePressureBased<distype>(integration_factor,
              solidpressure, spatial_material_mapping.determinant_deformation_gradient_, C_inv_vec,
              jacobian_mapping.N_XYZ_, *stiff);
        }
      });
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>::CouplingPoroelast(
    const DRT::Element& ele, MAT::StructPoro& porostructmat, MAT::FluidPoroMultiPhase& porofluidmat,
    const INPAR::STR::KinemType& kinematictype, const DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseMatrix& stiffness_matrix)
{
  // get primary variables of multiphase porous medium flow
  std::vector<double> fluidmultiphase_ephi(la[1].Size());
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
  DRT::UTILS::ExtractMyValues(*matrix_state, fluidmultiphase_ephi, la[1].lm_);

  // Initialize variables of multiphase porous medium flow
  const int nummultifluiddofpernode = porofluidmat.NumMat();
  const int numfluidphases = porofluidmat.NumFluidPhases();
  const int numvolfrac = porofluidmat.NumVolFrac();
  const bool hasvolfracs = (nummultifluiddofpernode > numfluidphases);

  // get nodal coordinates current and reference
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, la[0].lm_);

  // Loop over all Gauss points
  ForEachGaussPoint(nodal_coordinates, gauss_integration_,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp

      )
      {
        const SpatialMaterialMapping<distype> spatial_material_mapping =
            EvaluateSpatialMaterialMapping<distype>(
                jacobian_mapping, nodal_coordinates, 1.0, kinematictype);

        const CauchyGreen<distype> cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

        CORE::LINALG::Matrix<num_str_, num_dof_per_ele_> Bop =
            EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

        // volume change (used for porosity law). Same as J in nonlinear theory.
        const double volchange = ComputeVolumeChange<distype>(spatial_material_mapping,
            jacobian_mapping, ele, discretization, la[0].lm_, kinematictype);

        std::vector<double> fluidmultiphase_phiAtGP =
            ComputeFluidMultiPhasePrimaryVariablesAtGP<distype>(
                fluidmultiphase_ephi, nummultifluiddofpernode, shape_functions);

        std::vector<double> solidpressurederiv = ComputeSolidPressureDeriv<distype>(
            porofluidmat, fluidmultiphase_phiAtGP, numfluidphases);

        if (hasvolfracs)
        {
          double solidpressure = ComputeSolPressureAtGP<distype>(
              nummultifluiddofpernode, numfluidphases, fluidmultiphase_phiAtGP, porofluidmat);

          double porosity =
              ComputePorosity<distype>(porostructmat, params, solidpressure, volchange, gp);

          RecalculateSolPressureDeriv(fluidmultiphase_phiAtGP, nummultifluiddofpernode,
              numfluidphases, numvolfrac, solidpressure, porosity, solidpressurederiv);
        }

        const double detJ_w = jacobian_mapping.determinant_ * gauss_integration_.Weight(gp);

        // inverse Right Cauchy-Green tensor as vector in voigt notation
        CORE::LINALG::Matrix<num_str_, 1> C_inv_vec(false);
        CORE::LINALG::VOIGT::Stresses::MatrixToVector(
            cauchygreen.inverse_right_cauchy_green, C_inv_vec);

        // B^T . C^-1
        CORE::LINALG::Matrix<num_dof_per_ele_, 1> BopCinv(true);
        BopCinv.MultiplyTN(Bop, C_inv_vec);

        UpdateStiffnessMatrixCouplingMultiPhasePressureBased<distype>(detJ_w, solidpressurederiv,
            BopCinv, shape_functions, spatial_material_mapping.determinant_deformation_gradient_,
            nummultifluiddofpernode, stiffness_matrix);
      });
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>::CouplingStress(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  dserror("coupling stress poroelast not yet implemented for pressure-based variant");
}

// template classes
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<
    DRT::Element::DiscretizationType::hex18>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<
    DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<
    DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<
    DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<
    DRT::Element::DiscretizationType::pyramid5>;