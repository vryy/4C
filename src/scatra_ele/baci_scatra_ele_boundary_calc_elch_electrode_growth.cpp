/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer
growth, e.g., lithium plating

\level 2

 */
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_boundary_calc_elch_electrode_growth.H"

#include "baci_scatra_ele_boundary_calc_elch_electrode_growth_utils.H"
#include "baci_scatra_ele_parameter_elch.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_scatra_ele_parameter_timint.H"
#include "baci_scatra_ele_parameter_boundary.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"

#include "baci_mat_electrode.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
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
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::ScaTraEleBoundaryCalcElchElectrodeGrowth(const int numdofpernode, const int numscal,
    const std::string& disname)
    : myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode, numscal, disname),
      egrowth_(true)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::EvaluateMinMaxOverpotential(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  if (my::scatraparamsboundary_->ConditionType() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    dserror(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    my::EvalShapeFuncAndIntFac(intpoints, gpid);

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

    // compute exchange current density
    double i0 = kr * faraday * pow(emasterphiint, alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = myelectrodegrowthutils::GetButlerVolmerCurrentDensity(i0, frt, eslavepotint,
        emasterpotint, 0.0, eslaveresistanceint, eslavegrowthint, my::scatraparams_,
        my::scatraparamsboundary_);

    // calculate electrode-electrolyte overpotential at integration point
    const double eta = eslavepotint - emasterpotint - i * eslaveresistanceint;

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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& eslavematrix, CORE::LINALG::SerialDenseMatrix& emastermatrix,
    CORE::LINALG::SerialDenseVector& eslaveresidual)
{
  // safety checks
  if (my::numscal_ != 1) dserror("Invalid number of transported scalars!");
  if (my::numdofpernode_ != 2) dserror("Invalid number of degrees of freedom per node!");
  if (myelch::elchparams_->EquPot() != INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // extract condition type
  const DRT::Condition::ConditionType& s2iconditiontype =
      my::scatraparamsboundary_->ConditionType();
  if (s2iconditiontype != DRT::Condition::S2IKinetics and
      s2iconditiontype != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if ((s2iconditiontype == DRT::Condition::S2IKinetics and
          kineticmodel != INPAR::S2I::kinetics_butlervolmer) or
      (s2iconditiontype == DRT::Condition::S2ICouplingGrowth and
          kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer))
  {
    dserror(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const double faraday = myelch::elchparams_->Faraday();
  const double invF = 1.0 / faraday;
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);
    const double detF = my::CalculateDetFOfParentElement(ele, intpoints.Point(gpid));

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");

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

    // equilibrium electric potential difference and its derivative w.r.t. concentration at
    // electrode surface
    const double epd =
        s2iconditiontype == DRT::Condition::S2IKinetics
            ? matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF)
            : 0.0;
    const double epdderiv =
        matelectrode->ComputeDOpenCircuitPotentialDConcentration(eslavephiint, faraday, frt, detF);

    // compute exchange current density
    double i0 = kr * faraday * pow(emasterphiint, alphaa);
    if (s2iconditiontype == DRT::Condition::S2IKinetics and not std::isinf(epd))
      i0 *= pow(cmax - eslavephiint, alphaa) * pow(eslavephiint, alphac);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = myelectrodegrowthutils::GetButlerVolmerCurrentDensity(i0, frt, eslavepotint,
        emasterpotint, epd, eslaveresistanceint, eslavegrowthint, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations and residual contributions only in case of
    // non-zero Butler-Volmer current density to avoid unnecessary effort and to consistently
    // enforce the lithium plating condition
    if (std::abs(i) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point and regularization
      // factor
      const double eta = eslavepotint - emasterpotint - epd - i * eslaveresistanceint;
      const double regfac = myelectrodegrowthutils::GetRegularizationFactor(
          eslavegrowthint, eta, my::scatraparamsboundary_);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);

      double di_dc_slave(0.0), di_dc_master(0.0), di_dpot_slave(0.0), di_dpot_master(0.0);

      myelectrodegrowthutils::CalculateS2IElchElchLinearizations(i0, frt, epdderiv,
          eslaveresistanceint, regfac, expterm1, expterm2, faraday, emasterphiint, eslavephiint,
          cmax, my::scatraparamsboundary_, di_dc_slave, di_dc_master, di_dpot_slave,
          di_dpot_master);

      // compute linearizations and residual contributions associated with equations for lithium
      // transport
      for (int irow = 0; irow < nen_; ++irow)
      {
        const int row_conc = irow * my::numdofpernode_;
        const double funct_irow_invF_timefacfac = my::funct_(irow) * invF * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
        {
          const int col_conc = icol * my::numdofpernode_;
          const int col_pot = col_conc + 1;

          eslavematrix(row_conc, col_conc) +=
              funct_irow_invF_timefacfac * di_dc_slave * my::funct_(icol);
          eslavematrix(row_conc, col_pot) +=
              funct_irow_invF_timefacfac * di_dpot_slave * my::funct_(icol);
          emastermatrix(row_conc, col_conc) +=
              funct_irow_invF_timefacfac * di_dc_master * my::funct_(icol);
          emastermatrix(row_conc, col_pot) +=
              funct_irow_invF_timefacfac * di_dpot_master * my::funct_(icol);
        }

        eslaveresidual[row_conc] -= my::funct_(irow) * invF * timefacrhsfac * i;
      }
    }  // if(std::abs(i) > 1.e-16)
  }    // loop over integration points

  // compute linearizations and residual contributions associated with closing equations for
  // electric potential
  for (int irow = 0; irow < nen_; ++irow)
  {
    const int row_conc = irow * my::numdofpernode_;
    const int row_pot = row_conc + 1;

    for (int icol = 0; icol < nen_; ++icol)
    {
      const int col_conc = icol * my::numdofpernode_;
      const int col_pot = col_conc + 1;

      eslavematrix(row_pot, col_conc) += numelectrons * eslavematrix(row_conc, col_conc);
      eslavematrix(row_pot, col_pot) += numelectrons * eslavematrix(row_conc, col_pot);
      emastermatrix(row_pot, col_conc) += numelectrons * emastermatrix(row_conc, col_conc);
      emastermatrix(row_pot, col_pot) += numelectrons * emastermatrix(row_conc, col_pot);
    }

    eslaveresidual[row_pot] += numelectrons * eslaveresidual[row_conc];
  }
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::EvaluateS2ICoupling

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
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
      EvaluateS2ICouplingGrowthGrowth(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_s2icoupling_growthscatra:
    {
      EvaluateS2ICouplingGrowthScatra(
          ele, params, discretization, la, elemat1_epetra, elemat2_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_s2icoupling_scatragrowth:
    {
      EvaluateS2ICouplingScatraGrowth(ele, params, discretization, la, elemat1_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_elch_minmax_overpotential:
    {
      EvaluateMinMaxOverpotential(ele, params, discretization, la);
      break;
    }

    default:
    {
      myelch::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch action

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::EvaluateS2ICouplingScatraGrowth(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // extract condition type
  const DRT::Condition::ConditionType& s2iconditiontype =
      my::scatraparamsboundary_->ConditionType();
  if (s2iconditiontype != DRT::Condition::S2IKinetics and
      s2iconditiontype != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if ((s2iconditiontype == DRT::Condition::S2IKinetics and
          kineticmodel != INPAR::S2I::kinetics_butlervolmer) or
      (s2iconditiontype == DRT::Condition::S2ICouplingGrowth and
          kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer))
  {
    dserror(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const double faraday = myelch::elchparams_->Faraday();
  const double invF = 1.0 / faraday;
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);
    const double detF = my::CalculateDetFOfParentElement(ele, intpoints.Point(gpid));

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");

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

    // equilibrium electric potential difference at electrode surface
    const double epd =
        s2iconditiontype == DRT::Condition::S2IKinetics
            ? matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF)
            : 0.0;

    // compute exchange current density
    double i0 = kr * faraday * pow(emasterphiint, alphaa);
    if (s2iconditiontype == DRT::Condition::S2IKinetics and not std::isinf(epd))
      i0 *= pow(cmax - eslavephiint, alphaa) * pow(eslavephiint, alphac);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = myelectrodegrowthutils::GetButlerVolmerCurrentDensity(i0, frt, eslavepotint,
        emasterpotint, epd, eslaveresistanceint, eslavegrowthint, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations only in case of non-zero Butler-Volmer current
    // density to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if (std::abs(i) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point, regularization factor
      // and derivative of regularization factor
      const double eta = eslavepotint - emasterpotint - epd - i * eslaveresistanceint;
      const double regfac = myelectrodegrowthutils::GetRegularizationFactor(
          eslavegrowthint, eta, my::scatraparamsboundary_);
      const double regfacderiv = myelectrodegrowthutils::GetRegularizationFactorDerivative(
          eslavegrowthint, eta, my::scatraparamsboundary_);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);

      const double di_dgrowth = myelectrodegrowthutils::CalculateS2IElchGrowthLinearizations(i0, i,
          frt, eslaveresistanceint, resistivity, regfac, regfacderiv, expterm1, expterm2,
          my::scatraparamsboundary_);

      // compute linearizations associated with equations for lithium transport
      for (int irow = 0; irow < nen_; ++irow)
      {
        const int row_conc = irow * my::numdofpernode_;
        const double funct_irow_invF_timefacfac = my::funct_(irow) * invF * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
          eslavematrix(row_conc, icol) +=
              funct_irow_invF_timefacfac * di_dgrowth * my::funct_(icol);
      }
    }  // if(std::abs(i) > 1.e-16)
  }    // loop over integration points

  // compute linearizations associated with closing equations for electric potential
  for (int irow = 0; irow < nen_; ++irow)
  {
    const int row_conc = irow * my::numdofpernode_;
    const int row_pot = row_conc + 1;

    for (int icol = 0; icol < nen_; ++icol)
      eslavematrix(row_pot, icol) += numelectrons * eslavematrix(row_conc, icol);
  }
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
   // probdim>::EvaluateS2ICouplingScatraGrowth

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::EvaluateS2ICouplingGrowthScatra(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix,
    CORE::LINALG::SerialDenseMatrix& emastermatrix)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  if (my::scatraparamsboundary_->ConditionType() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    dserror(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  if (kr < 0.0) dserror("Charge transfer constant k_r is negative!");
  const double resistivity = my::scatraparamsboundary_->Resistivity();
  const double factor =
      my::scatraparamsboundary_->MolarMass() / (my::scatraparamsboundary_->Density() * faraday);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");

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

    // compute exchange current density
    const double i0 = kr * faraday * pow(emasterphiint, alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = myelectrodegrowthutils::GetButlerVolmerCurrentDensity(i0, frt, eslavepotint,
        emasterpotint, 0.0, eslaveresistanceint, eslavegrowthint, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations only in case of non-zero Butler-Volmer current
    // density to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if (std::abs(i) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point and regularization
      // factor
      const double eta = eslavepotint - emasterpotint - i * eslaveresistanceint;
      const double regfac = myelectrodegrowthutils::GetRegularizationFactor(
          eslavegrowthint, eta, my::scatraparamsboundary_);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);

      double dummy(0.0), di_dc_master(0.0), di_dpot_slave(0.0), di_dpot_master(0.0);

      myelectrodegrowthutils::CalculateS2IElchElchLinearizations(i0, frt, dummy,
          eslaveresistanceint, regfac, expterm1, expterm2, faraday, emasterphiint, dummy, dummy,
          my::scatraparamsboundary_, dummy, di_dc_master, di_dpot_slave, di_dpot_master);

      // compute linearizations associated with equation for scatra-scatra interface layer growth
      for (int irow = 0; irow < nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow) * factor * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
        {
          const int col_conc = icol * my::numdofpernode_;
          const int col_pot = col_conc + 1;

          eslavematrix(irow, col_pot) +=
              funct_irow_factor_timefacfac * di_dpot_slave * my::funct_(icol);
          emastermatrix(irow, col_conc) +=
              funct_irow_factor_timefacfac * di_dc_master * my::funct_(icol);
          emastermatrix(irow, col_pot) +=
              funct_irow_factor_timefacfac * di_dpot_master * my::funct_(icol);
        }
      }
    }  // if(std::abs(i) > 1.e-16)
  }    // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
    probdim>::EvaluateS2ICouplingGrowthGrowth(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& eslavematrix,
    CORE::LINALG::SerialDenseVector& eslaveresidual)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  CORE::LINALG::Matrix<nen_, 1> eslavegrowthhist(true);
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");
  my::ExtractNodeValues(eslavegrowthhist, discretization, la, "growthhist", 2);

  if (my::scatraparamsboundary_->ConditionType() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  if (kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
  {
    dserror(
        "Received illegal kinetic model for scatra-scatra interface coupling involving interface "
        "layer growth!");
  }
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  const double resistivity = my::scatraparamsboundary_->Resistivity();
  const double factor =
      my::scatraparamsboundary_->MolarMass() / (my::scatraparamsboundary_->Density() * faraday);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate mass matrix
    for (int irow = 0; irow < nen_; ++irow)
      for (int icol = 0; icol < nen_; ++icol)
        eslavematrix(irow, icol) += my::funct_(irow) * my::funct_(icol) * fac;

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");

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

    // compute exchange current density
    const double i0 = kr * faraday * pow(emasterphiint, alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = myelectrodegrowthutils::GetButlerVolmerCurrentDensity(i0, frt, eslavepotint,
        emasterpotint, 0.0, eslaveresistanceint, eslavegrowthint, my::scatraparams_,
        my::scatraparamsboundary_);

    // continue with evaluation of linearizations and residual contributions only in case of
    // non-zero Butler-Volmer current density to avoid unnecessary effort and to consistently
    // enforce the lithium plating condition. (If the plating condition is not fulfilled, we
    // manually set the Butler-Volmer current density to zero, and thus we need to make sure that
    // all linearizations are also zero, i.e., that nothing is added to the element matrix.)
    if (std::abs(i) > 1.0e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point, regularization factor
      // and derivative of regularization factor
      const double eta = eslavepotint - emasterpotint - i * eslaveresistanceint;
      const double regfac = myelectrodegrowthutils::GetRegularizationFactor(
          eslavegrowthint, eta, my::scatraparamsboundary_);
      const double regfacderiv = myelectrodegrowthutils::GetRegularizationFactorDerivative(
          eslavegrowthint, eta, my::scatraparamsboundary_);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);

      const double di_dgrowth = myelectrodegrowthutils::CalculateS2IElchGrowthLinearizations(i0, i,
          frt, eslaveresistanceint, resistivity, regfac, regfacderiv, expterm1, expterm2,
          my::scatraparamsboundary_);

      // compute linearizations and residual contributions associated with equation for
      // scatra-scatra interface layer growth
      for (int irow = 0; irow < nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow) * factor * timefacfac;

        for (int icol = 0; icol < nen_; ++icol)
          eslavematrix(irow, icol) += funct_irow_factor_timefacfac * di_dgrowth * my::funct_(icol);

        eslaveresidual[irow] -= my::funct_(irow) * (eslavegrowthint - eslavegrowthhistint) * fac;
        eslaveresidual[irow] -= my::funct_(irow) * factor * i * timefacrhsfac;
      }
    }  // if(std::abs(i) > 1.e-16)
  }    // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::ExtractNodeValues(
    const DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // call base class routine
  my::ExtractNodeValues(discretization, la);

  // extract nodal growth variables associated with boundary element
  my::ExtractNodeValues(egrowth_, discretization, la, "growth", 2);
}

// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::nurbs9, 3>;
