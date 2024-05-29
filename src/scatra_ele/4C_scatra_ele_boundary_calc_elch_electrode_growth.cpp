/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer
growth, e.g., lithium plating

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_boundary_calc_elch_electrode_growth.hpp"

#include "4C_discretization_fem_general_utils_boundary_integration.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode_growth_utils.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode_utils.hpp"
#include "4C_scatra_ele_parameter_boundary.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](int numdofpernode, int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::ScaTraEleBoundaryCalcElchElectrodeGrowth(const int numdofpernode, const int numscal,
    const std::string& disname)
    : myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode, numscal, disname),
      egrowth_(true)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::evaluate_min_max_overpotential(const CORE::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->parent_element()->Material());
  if (matelectrode == Teuchos::null)
    FOUR_C_THROW("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  extract_node_values(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::extract_node_values(emasterphinp, discretization, la, "imasterphinp");

  if (my::scatraparamsboundary_->ConditionType() != CORE::Conditions::S2IKineticsGrowth)
    FOUR_C_THROW("Received illegal condition type!");

  // access input parameters associated with condition
  const double faraday = myelch::elchparams_->Faraday();
  const double resistivity = my::scatraparamsboundary_->Resistivity();
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    FOUR_C_THROW(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    my::eval_shape_func_and_int_fac(intpoints, gpid);

    // evaluate factor F/RT
    const double frt = myelch::elchparams_->FRT();

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint * resistivity;

    double eta = 0.0;

    switch (kineticmodel)
    {
      case INPAR::S2I::growth_kinetics_butlervolmer:
      {
        const double alphaa = my::scatraparamsboundary_->AlphaA();
        const double kr = my::scatraparamsboundary_->charge_transfer_constant();
        const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
        const double epd = 0.0;  // equilibrium potential is 0 for the plating reaction

        // compute exchange mass flux density
        const double j0 = kr * std::pow(emasterphiint, alphaa);

        // compute mass flux density of growth kinetics via Newton-Raphson iteration
        const double j = CalculateGrowthMassFluxDensity(j0, frt, eslavepotint, emasterpotint, epd,
            eslaveresistanceint, eslavegrowthint, faraday, my::scatraparams_,
            my::scatraparamsboundary_);

        // calculate electrode-electrolyte overpotential at integration point
        eta = eslavepotint - emasterpotint - j * faraday * eslaveresistanceint;

        break;
      }
      default:
      {
        FOUR_C_THROW("Model for scatra-scatra interface growth kinetics not implemented!");
      }
    }

    // check for minimality and update result if applicable
    auto& etagrowthmin = params.get<double>("etagrowthmin");
    if (eta < etagrowthmin) etagrowthmin = eta;

    // check for maximality and update result if applicable
    auto& etagrowthmax = params.get<double>("etagrowthmax");
    if (eta > etagrowthmax) etagrowthmax = eta;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::evaluate_s2_i_coupling(const CORE::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix,
    CORE::LINALG::SerialDenseMatrix& emastermatrix, CORE::LINALG::SerialDenseVector& eslaveresidual)
{
  // safety checks
  if (my::numscal_ != 1) FOUR_C_THROW("Invalid number of transported scalars!");
  if (my::numdofpernode_ != 2) FOUR_C_THROW("Invalid number of degrees of freedom per node!");
  if (myelch::elchparams_->EquPot() != INPAR::ELCH::equpot_divi)
    FOUR_C_THROW("Invalid closing equation for electric potential!");

  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->parent_element()->Material());
  if (matelectrode == Teuchos::null)
    FOUR_C_THROW("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  extract_node_values(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::extract_node_values(emasterphinp, discretization, la, "imasterphinp");

  // extract condition type
  const CORE::Conditions::ConditionType& s2iconditiontype =
      my::scatraparamsboundary_->ConditionType();
  if (s2iconditiontype != CORE::Conditions::S2IKinetics and
      s2iconditiontype != CORE::Conditions::S2IKineticsGrowth)
    FOUR_C_THROW("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->charge_transfer_constant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_int_fac(intpoints, gpid);
    const double detF = my::calculate_det_f_of_parent_element(ele, intpoints.Point(gpid));

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) FOUR_C_THROW("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelch::elchparams_->FRT();

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint * resistivity;

    switch (s2iconditiontype)
    {
      case CORE::Conditions::S2IKinetics:
      {
        // equilibrium electric potential difference and its derivative w.r.t. concentration at
        // electrode surface
        const double epd =
            matelectrode->compute_open_circuit_potential(eslavephiint, faraday, frt, detF);
        const double epdderiv = matelectrode->compute_d_open_circuit_potential_d_concentration(
            eslavephiint, faraday, frt, detF);

        // compute exchange mass flux density
        const double j0 = CalculateButlerVolmerExchangeMassFluxDensity(
            kr, alphaa, alphac, cmax, eslavephiint, emasterphiint, kineticmodel, s2iconditiontype);

        // compute Butler-Volmer mass flux density via Newton-Raphson iteration
        const double j = CalculateModifiedButlerVolmerMassFluxDensity(j0, alphaa, alphac, frt,
            eslavepotint, emasterpotint, epd, eslaveresistanceint,
            my::scatraparams_->int_layer_growth_ite_max(),
            my::scatraparams_->int_layer_growth_conv_tol(), faraday);

        // continue with evaluation of linearizations and residual contributions only in case of
        // non-zero Butler-Volmer mass flux density to avoid unnecessary effort and to consistently
        // enforce the lithium plating condition
        if (std::abs(j) > 1.0e-16)
        {
          // calculate electrode-electrolyte overpotential at integration point and regularization
          // factor
          const double eta = eslavepotint - emasterpotint - epd - j * faraday * eslaveresistanceint;
          // exponential Butler-Volmer terms
          const double expterm1 = std::exp(alphaa * frt * eta);
          const double expterm2 = std::exp(-alphac * frt * eta);

          double dj_dc_slave(0.0), dj_dc_master(0.0), dj_dpot_slave(0.0), dj_dpot_master(0.0);
          CalculateButlerVolmerElchLinearizations(kineticmodel, j0, frt, epdderiv, alphaa, alphac,
              eslaveresistanceint, expterm1, expterm2, kr, faraday, emasterphiint, eslavephiint,
              cmax, eta, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

          calculate_rhs_and_linearization(numelectrons, timefacfac, timefacrhsfac, j, dj_dc_slave,
              dj_dc_master, dj_dpot_slave, dj_dpot_master, eslavematrix, emastermatrix,
              eslaveresidual);
        }

        break;
      }
      case CORE::Conditions::S2IKineticsGrowth:
      {
        // equilibrium electric potential difference and its derivative w.r.t. concentration at
        // electrode surface
        const double epd = 0.0;
        const double epdderiv = matelectrode->compute_d_open_circuit_potential_d_concentration(
            eslavephiint, faraday, frt, detF);

        // compute exchange mass flux density of growth kinetics
        const double j0 = CalculateGrowthExchangeMassFluxDensity(
            kr, alphaa, emasterphiint, kineticmodel, s2iconditiontype);

        // compute mass flux density of growth kinetics via Newton-Raphson iteration
        const double j = CalculateGrowthMassFluxDensity(j0, frt, eslavepotint, emasterpotint, epd,
            eslaveresistanceint, eslavegrowthint, faraday, my::scatraparams_,
            my::scatraparamsboundary_);

        // continue with evaluation of linearizations and residual contributions only in case of
        // non-zero Butler-Volmer mass flux density to avoid unnecessary effort and to consistently
        // enforce the lithium plating condition
        if (std::abs(j) > 1.0e-16)
        {
          // calculate electrode-electrolyte overpotential at integration point and regularization
          // factor
          const double eta = eslavepotint - emasterpotint - epd - j * faraday * eslaveresistanceint;
          const double regfac =
              GetRegularizationFactor(eslavegrowthint, eta, my::scatraparamsboundary_);

          double dj_dc_slave(0.0), dj_dc_master(0.0), dj_dpot_slave(0.0), dj_dpot_master(0.0);
          CalculateS2IGrowthElchLinearizations(j0, frt, epdderiv, eta, eslaveresistanceint, regfac,
              emasterphiint, eslavephiint, cmax, my::scatraparamsboundary_, dj_dc_slave,
              dj_dc_master, dj_dpot_slave, dj_dpot_master);

          calculate_rhs_and_linearization(numelectrons, timefacfac, timefacrhsfac, j, dj_dc_slave,
              dj_dc_master, dj_dpot_slave, dj_dpot_master, eslavematrix, emastermatrix,
              eslaveresidual);
        }
        break;
      }
      default:
        FOUR_C_THROW("S2I condition type not recognized!");
    }
  }  // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::calculate_rhs_and_linearization(const int numelectrons, const double timefacfac,
    const double timefacrhsfac, const double j, const double dj_dc_slave, const double dj_dc_master,
    const double dj_dpot_slave, const double dj_dpot_master,
    CORE::LINALG::SerialDenseMatrix& eslavematrix, CORE::LINALG::SerialDenseMatrix& emastermatrix,
    CORE::LINALG::SerialDenseVector& eslaveresidual) const
{
  for (int irow = 0; irow < nen_; ++irow)
  {
    const int row_conc = irow * my::numdofpernode_;
    const int row_pot = row_conc + 1;
    const double funct_irow_timefacfac = my::funct_(irow) * timefacfac;

    for (int icol = 0; icol < nen_; ++icol)
    {
      const int col_conc = icol * my::numdofpernode_;
      const int col_pot = col_conc + 1;

      eslavematrix(row_conc, col_conc) += funct_irow_timefacfac * dj_dc_slave * my::funct_(icol);
      eslavematrix(row_conc, col_pot) += funct_irow_timefacfac * dj_dpot_slave * my::funct_(icol);
      eslavematrix(row_pot, col_conc) +=
          numelectrons * funct_irow_timefacfac * dj_dc_slave * my::funct_(icol);
      eslavematrix(row_pot, col_pot) +=
          numelectrons * funct_irow_timefacfac * dj_dpot_slave * my::funct_(icol);

      emastermatrix(row_conc, col_conc) += funct_irow_timefacfac * dj_dc_master * my::funct_(icol);
      emastermatrix(row_conc, col_pot) += funct_irow_timefacfac * dj_dpot_master * my::funct_(icol);
      emastermatrix(row_pot, col_conc) +=
          numelectrons * funct_irow_timefacfac * dj_dc_master * my::funct_(icol);
      emastermatrix(row_pot, col_pot) +=
          numelectrons * funct_irow_timefacfac * dj_dpot_master * my::funct_(icol);
    }

    eslaveresidual[row_conc] -= my::funct_(irow) * timefacrhsfac * j;
    eslaveresidual[row_pot] -= numelectrons * my::funct_(irow) * timefacrhsfac * j;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::evaluate_action(
    CORE::Elements::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, SCATRA::BoundaryAction action,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_s2icoupling_growthgrowth:
    {
      evaluate_s2_i_coupling_growth_growth(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_s2icoupling_growthscatra:
    {
      evaluate_s2_i_coupling_growth_scatra(
          ele, params, discretization, la, elemat1_epetra, elemat2_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_s2icoupling_scatragrowth:
    {
      evaluate_s2_i_coupling_scatra_growth(ele, params, discretization, la, elemat1_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_elch_minmax_overpotential:
    {
      evaluate_min_max_overpotential(ele, params, discretization, la);
      break;
    }

    default:
    {
      myelch::evaluate_action(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch action

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::evaluate_s2_i_coupling_scatra_growth(const CORE::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->parent_element()->Material());
  if (matelectrode == Teuchos::null)
    FOUR_C_THROW("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  extract_node_values(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::extract_node_values(emasterphinp, discretization, la, "imasterphinp");

  // extract condition type
  const CORE::Conditions::ConditionType& s2iconditiontype =
      my::scatraparamsboundary_->ConditionType();
  if (s2iconditiontype != CORE::Conditions::S2IKinetics and
      s2iconditiontype != CORE::Conditions::S2IKineticsGrowth)
    FOUR_C_THROW("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->charge_transfer_constant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_int_fac(intpoints, gpid);
    const double detF = my::calculate_det_f_of_parent_element(ele, intpoints.Point(gpid));

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) FOUR_C_THROW("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelch::elchparams_->FRT();

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint * resistivity;

    switch (s2iconditiontype)
    {
      case CORE::Conditions::S2IKinetics:
      {
        // equilibrium electric potential difference at electrode surface
        const double epd =
            matelectrode->compute_open_circuit_potential(eslavephiint, faraday, frt, detF);

        // compute exchange mass flux density
        const double j0 = CalculateButlerVolmerExchangeMassFluxDensity(
            kr, alphaa, alphac, cmax, eslavephiint, emasterphiint, kineticmodel, s2iconditiontype);

        // compute Butler-Volmer mass flux density via Newton-Raphson iteration
        const double j = CalculateModifiedButlerVolmerMassFluxDensity(j0, alphaa, alphac, frt,
            eslavepotint, emasterpotint, epd, eslaveresistanceint,
            my::scatraparams_->int_layer_growth_ite_max(),
            my::scatraparams_->int_layer_growth_conv_tol(), faraday);

        // continue with evaluation of linearizations only in case of non-zero Butler-Volmer mass
        // flux density to avoid unnecessary effort and to consistently enforce the lithium plating
        // condition
        if (std::abs(j) > 1.0e-16)
        {
          // calculate electrode-electrolyte overpotential at integration point, regularization
          // factor and derivative of regularization factor
          const double eta = eslavepotint - emasterpotint - epd - j * faraday * eslaveresistanceint;
          // no regfac required in this case
          const double regfac_dummy = 1.0;

          // exponential Butler-Volmer terms
          const double expterm1 = std::exp(alphaa * frt * eta);
          const double expterm2 = std::exp(-alphac * frt * eta);

          const double dj_dgrowth =
              CalculateS2IElchGrowthLinearizations(j0, j, frt, eslaveresistanceint, resistivity,
                  regfac_dummy, regfac_dummy, expterm1, expterm2, my::scatraparamsboundary_);

          // compute linearizations
          for (int irow = 0; irow < nen_; ++irow)
          {
            const int row_conc = irow * my::numdofpernode_;
            const int row_pot = row_conc + 1;
            const double funct_irow_timefacfac = my::funct_(irow) * timefacfac;

            for (int icol = 0; icol < nen_; ++icol)
            {
              eslavematrix(row_conc, icol) += funct_irow_timefacfac * dj_dgrowth * my::funct_(icol);
              eslavematrix(row_pot, icol) +=
                  numelectrons * funct_irow_timefacfac * dj_dgrowth * my::funct_(icol);
            }
          }
        }

        break;
      }
      case CORE::Conditions::S2IKineticsGrowth:
      {
        // equilibrium electric potential difference at electrode surface
        const double epd = 0.0;

        // compute exchange mass flux density
        const double j0 = CalculateGrowthExchangeMassFluxDensity(
            kr, alphaa, emasterphiint, kineticmodel, s2iconditiontype);

        // compute Butler-Volmer mass flux density via Newton-Raphson iteration
        const double j = CalculateGrowthMassFluxDensity(j0, frt, eslavepotint, emasterpotint, epd,
            eslaveresistanceint, eslavegrowthint, faraday, my::scatraparams_,
            my::scatraparamsboundary_);

        // continue with evaluation of linearizations only in case of non-zero Butler-Volmer mass
        // flux density to avoid unnecessary effort and to consistently enforce the lithium plating
        // condition
        if (std::abs(j) > 1.0e-16)
        {
          // calculate electrode-electrolyte overpotential at integration point, regularization
          // factor and derivative of regularization factor
          const double eta = eslavepotint - emasterpotint - epd - j * faraday * eslaveresistanceint;
          const double regfac =
              GetRegularizationFactor(eslavegrowthint, eta, my::scatraparamsboundary_);
          const double regfacderiv =
              GetRegularizationFactorDerivative(eslavegrowthint, eta, my::scatraparamsboundary_);

          // exponential Butler-Volmer terms
          const double expterm1 = std::exp(alphaa * frt * eta);
          const double expterm2 = std::exp(-alphac * frt * eta);

          const double dj_dgrowth =
              CalculateS2IElchGrowthLinearizations(j0, j, frt, eslaveresistanceint, resistivity,
                  regfac, regfacderiv, expterm1, expterm2, my::scatraparamsboundary_);

          // compute linearizations
          for (int irow = 0; irow < nen_; ++irow)
          {
            const int row_conc = irow * my::numdofpernode_;
            const int row_pot = row_conc + 1;
            const double funct_irow_timefacfac = my::funct_(irow) * timefacfac;

            for (int icol = 0; icol < nen_; ++icol)
            {
              eslavematrix(row_conc, icol) += funct_irow_timefacfac * dj_dgrowth * my::funct_(icol);
              eslavematrix(row_pot, icol) +=
                  numelectrons * funct_irow_timefacfac * dj_dgrowth * my::funct_(icol);
            }
          }
        }

        break;
      }
      default:
        FOUR_C_THROW("S2I condition type not recognized!");
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::evaluate_s2_i_coupling_growth_scatra(const CORE::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix,
    CORE::LINALG::SerialDenseMatrix& emastermatrix)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->parent_element()->Material());
  if (matelectrode == Teuchos::null)
    FOUR_C_THROW("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  extract_node_values(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::extract_node_values(emasterphinp, discretization, la, "imasterphinp");

  if (my::scatraparamsboundary_->ConditionType() != CORE::Conditions::S2IKineticsGrowth)
    FOUR_C_THROW("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    FOUR_C_THROW(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double kr = my::scatraparamsboundary_->charge_transfer_constant();
  if (kr < 0.0) FOUR_C_THROW("Charge transfer constant k_r is negative!");
  const double resistivity = my::scatraparamsboundary_->Resistivity();
  const double factor =
      my::scatraparamsboundary_->MolarMass() / (my::scatraparamsboundary_->Density());

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_int_fac(intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) FOUR_C_THROW("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelch::elchparams_->FRT();

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint * resistivity;

    // compute exchange mass flux density
    const double j0 = kr * std::pow(emasterphiint, alphaa);

    // compute mass flux density of growth kinetics via Newton-Raphson iteration
    const double j = CalculateGrowthMassFluxDensity(j0, frt, eslavepotint, emasterpotint, 0.0,
        eslaveresistanceint, eslavegrowthint, faraday, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations only in case of non-zero Butler-Volmer mass flux
    // density to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if (std::abs(j) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point and regularization
      // factor
      const double eta = eslavepotint - emasterpotint - j * faraday * eslaveresistanceint;
      const double regfac =
          GetRegularizationFactor(eslavegrowthint, eta, my::scatraparamsboundary_);

      double dummy(0.0), dj_dc_master(0.0), dj_dpot_slave(0.0), dj_dpot_master(0.0);
      CalculateS2IGrowthElchLinearizations(j0, frt, dummy, eta, eslaveresistanceint, regfac,
          emasterphiint, dummy, dummy, my::scatraparamsboundary_, dummy, dj_dc_master,
          dj_dpot_slave, dj_dpot_master);

      // compute linearizations associated with equation for scatra-scatra interface layer growth
      for (int irow = 0; irow < nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow) * factor * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
        {
          const int col_conc = icol * my::numdofpernode_;
          const int col_pot = col_conc + 1;

          eslavematrix(irow, col_pot) +=
              funct_irow_factor_timefacfac * dj_dpot_slave * my::funct_(icol);
          emastermatrix(irow, col_conc) +=
              funct_irow_factor_timefacfac * dj_dc_master * my::funct_(icol);
          emastermatrix(irow, col_pot) +=
              funct_irow_factor_timefacfac * dj_dpot_master * my::funct_(icol);
        }
      }
    }  // if(std::abs(j) > 1.e-16)
  }    // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::evaluate_s2_i_coupling_growth_growth(const CORE::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix,
    CORE::LINALG::SerialDenseVector& eslaveresidual)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->parent_element()->Material());
  if (matelectrode == Teuchos::null)
    FOUR_C_THROW("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  extract_node_values(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  CORE::LINALG::Matrix<nen_, 1> eslavegrowthhist(true);
  my::extract_node_values(emasterphinp, discretization, la, "imasterphinp");
  my::extract_node_values(
      eslavegrowthhist, discretization, la, "growthhist", my::scatraparams_->NdsGrowth());

  if (my::scatraparamsboundary_->ConditionType() != CORE::Conditions::S2IKineticsGrowth)
    FOUR_C_THROW("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    FOUR_C_THROW(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->charge_transfer_constant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();
  const double factor =
      my::scatraparamsboundary_->MolarMass() / (my::scatraparamsboundary_->Density());

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_int_fac(intpoints, gpid);

    // evaluate mass matrix
    for (int irow = 0; irow < nen_; ++irow)
      for (int icol = 0; icol < nen_; ++icol)
        eslavematrix(irow, icol) += my::funct_(irow) * my::funct_(icol) * fac;

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) FOUR_C_THROW("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelch::elchparams_->FRT();

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double eslavegrowthhistint = my::funct_.Dot(eslavegrowthhist);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint * resistivity;

    // compute exchange mass flux density
    const double j0 = kr * std::pow(emasterphiint, alphaa);

    // compute mass flux density of growth kinetics via Newton-Raphson iteration
    const double j = CalculateGrowthMassFluxDensity(j0, frt, eslavepotint, emasterpotint, 0.0,
        eslaveresistanceint, eslavegrowthint, faraday, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations and residual contributions only in case of
    // non-zero Butler-Volmer mass flux density to avoid unnecessary effort and to consistently
    // enforce the lithium plating condition. (If the plating condition is not fulfilled, we
    // manually set the Butler-Volmer mass flux density to zero, and thus we need to make sure that
    // all linearizations are also zero, i.e., that nothing is added to the element matrix.)
    if (std::abs(j) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point, regularization factor
      // and derivative of regularization factor
      const double eta = eslavepotint - emasterpotint - j * faraday * eslaveresistanceint;
      const double regfac =
          GetRegularizationFactor(eslavegrowthint, eta, my::scatraparamsboundary_);
      const double regfacderiv =
          GetRegularizationFactorDerivative(eslavegrowthint, eta, my::scatraparamsboundary_);

      // exponential Butler-Volmer terms
      const double expterm1 = std::exp(alphaa * frt * eta);
      const double expterm2 = std::exp(-alphac * frt * eta);

      const double dj_dgrowth =
          CalculateS2IElchGrowthLinearizations(j0, j, frt, eslaveresistanceint, resistivity, regfac,
              regfacderiv, expterm1, expterm2, my::scatraparamsboundary_);

      // compute linearizations and residual contributions associated with equation for
      // scatra-scatra interface layer growth
      for (int irow = 0; irow < nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow) * factor * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
          eslavematrix(irow, icol) += funct_irow_factor_timefacfac * dj_dgrowth * my::funct_(icol);

        eslaveresidual[irow] -= my::funct_(irow) * (eslavegrowthint - eslavegrowthhistint) * fac;
        eslaveresidual[irow] -= my::funct_(irow) * factor * j * timefacrhsfac;
      }
    }  // if(std::abs(i) > 1.e-16)
  }    // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::extract_node_values(
    const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la)
{
  // call base class routine
  my::extract_node_values(discretization, la);

  // extract nodal growth variables associated with boundary element
  my::extract_node_values(egrowth_, discretization, la, "growth", my::scatraparams_->NdsGrowth());
}

// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::quad4,
    3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::quad8,
    3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::quad9,
    3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::line2,
    2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::line2,
    3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::line3,
    2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::nurbs3,
    2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<CORE::FE::CellType::nurbs9,
    3>;

FOUR_C_NAMESPACE_CLOSE
