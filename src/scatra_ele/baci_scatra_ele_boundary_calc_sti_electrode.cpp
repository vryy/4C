/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_boundary_calc_sti_electrode.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode_utils.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"

#include "baci_lib_discret.H"

#include "baci_mat_electrode.H"
#include "baci_mat_soret.H"

#include "baci_inpar_s2i.H"

#include "baci_scatra_ele_parameter_elch.H"
#include "baci_scatra_ele_parameter_timint.H"
#include "baci_scatra_ele_parameter_boundary.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_so3_utils.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>>(
            new ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype,
    probdim>::ScaTraEleBoundaryCalcSTIElectrode(const int numdofpernode, const int numscal,
    const std::string& disname)
    :  // constructor of base class
      my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname),

      // initialize member variable
      eelchnp_(2, CORE::LINALG::Matrix<nen_, 1>(true))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix,
    Epetra_SerialDenseVector& eslaveresidual)
{
  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret =
      Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::FourierIso> matfourier =
      Teuchos::rcp_dynamic_cast<const MAT::FourierIso>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if ((matsoret == Teuchos::null and matfourier == Teuchos::null) or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterscatra(2, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterscatra, discretization, la, "imasterscatra", 2);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();

  CORE::LINALG::Matrix<nen_, 1> emastertemp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
    my::ExtractNodeValues(emastertemp, discretization, la, "imastertemp", 3);

  // element slave mechanical stress tensor
  const bool is_pseudo_contact = my::scatraparamsboundary_->IsPseudoContact();
  std::vector<CORE::LINALG::Matrix<nen_, 1>> eslavestress_vector(
      6, CORE::LINALG::Matrix<nen_, 1>(true));
  if (is_pseudo_contact)
    my::ExtractNodeValues(eslavestress_vector, discretization, la, "mechanicalStressState",
        my::scatraparams_->NdsTwoTensorQuantity());

  CORE::LINALG::Matrix<nsd_, 1> normal;

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid, &normal);
    const double detF = my::CalculateDetFOfParentElement(ele, intpoints.Point(gpid));

    const double pseudo_contact_fac = my::CalculatePseudoContactFactor(
        is_pseudo_contact, eslavestress_vector, normal, my::funct_);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingAtIntegrationPoint<distype>(matelectrode, my::ephinp_[0], emastertemp,
        eelchnp_, emasterscatra, pseudo_contact_fac, my::funct_, my::funct_,
        my::scatraparamsboundary_, timefacfac, timefacrhsfac, detF, eslavematrix, emastermatrix,
        eslaveresidual);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype,
    probdim>::EvaluateS2ICouplingAtIntegrationPoint(const Teuchos::RCP<const MAT::Electrode>&
                                                        matelectrode,
    const CORE::LINALG::Matrix<nen_, 1>& eslavetempnp,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        emastertempnp,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& eslavephinp,
    const std::vector<CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const double pseudo_contact_fac, const CORE::LINALG::Matrix<nen_, 1>& funct_slave,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
    const double timefacfac, const double timefacrhsfac, const double detF,
    Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm, Epetra_SerialDenseVector& r_s)
{
  // get condition specific parameters
  const int kineticmodel = scatra_parameter_boundary->KineticModel();
  const double kr = scatra_parameter_boundary->ChargeTransferConstant();
  const double alphaa = scatra_parameter_boundary->AlphaA();
  const double alphac = scatra_parameter_boundary->AlphaC();
  const double peltier = scatra_parameter_boundary->Peltier();
  const double thermoperm = scatra_parameter_boundary->ThermoPerm();
  const double molar_heat_capacity = scatra_parameter_boundary->MolarHeatCapacity();

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  const double emastertempint = funct_master.Dot(emastertempnp);
  if (eslavetempint <= 0.) dserror("Temperature is non-positive!");
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  const int nen_master =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // access input parameters associated with current condition
  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // compute matrix and vector contributions according to kinetic model for current scatra-scatra
  // interface coupling condition
  switch (kineticmodel)
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      // extract saturation value of intercalated lithium concentration from electrode material
      const double cmax = matelectrode->CMax();

      // evaluate factor F/RT
      const double frt = faraday / (gasconstant * eslavetempint);

      // equilibrium electric potential difference at electrode surface
      const double epd =
          matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd;

      // Butler-Volmer exchange current density
      const double i0 = kr * faraday * pow(emasterphiint, alphaa) *
                        pow(cmax - eslavephiint, alphaa) * pow(eslavephiint, alphac);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
      const double expterm = expterm1 - expterm2;

      // safety check
      if (std::abs(expterm) > 1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // core residual term
      const double residual_timefacrhsfac =
          pseudo_contact_fac * timefacrhsfac * i0 * expterm * (eta + peltier);

      // core linearization w.r.t. temperature
      const double linearization_timefacfac =
          -pseudo_contact_fac * timefacfac * i0 * frt / eslavetempint * eta *
          (alphaa * expterm1 + alphac * expterm2) * (eta + peltier);

      // compute matrix and vector contributions
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int ui = 0; ui < nen_; ++ui)
          k_ss(vi, ui) -= funct_slave(vi) * linearization_timefacfac * funct_slave(ui);
        r_s[vi] += funct_slave(vi) * residual_timefacrhsfac;
      }

      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreducedthermoresistance:
    {
      // Flux of energy from mass flux and from difference in termperature
      double j_timefacrhsfac(0.0);
      double djdT_slave_timefacfac(0.0);
      double djdT_master_timefacfac(0.0);

      // Part 1: Energy flux from mass flux
      const double etempint = (eslavetempint + emastertempint) / 2.0;
      const double frt = faraday / (etempint * gasconstant);

      // equilibrium electric potential difference at electrode surface
      const double epd =
          matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);
      const double depddT = matelectrode->ComputeDOpenCircuitPotentialDTemperature(
          eslavephiint, faraday, gasconstant);

      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (std::isinf(epd)) break;

      // Butler-Volmer exchange mass flux density
      const double j0 = kr;

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd;

      // exponential Butler-Volmer terms
      const double expterm1 = std::exp(alphaa * frt * eta);
      const double expterm2 = std::exp(-alphac * frt * eta);
      const double expterm = expterm1 - expterm2;

      // safety check
      if (std::abs(expterm) > 1.0e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // core residual term associated with Butler-Volmer mass flux density
      const double j_mass = j0 * expterm;

      // temperature from source side of flux
      const double j_mass_energy = j_mass * molar_heat_capacity * etempint;

      j_timefacrhsfac += pseudo_contact_fac * timefacrhsfac * j_mass_energy;

      // forward declarations
      double dj_dT_slave(0.0);

      // calculate linearizations of Butler-Volmer kinetics w.r.t. temperature dofs
      myelectrodeutils::CalculateButlerVolmerTempLinearizations(
          alphaa, alphac, depddT, eta, etempint, faraday, frt, gasconstant, j0, dj_dT_slave);

      const double dj_mass_energydT_slave =
          dj_dT_slave * molar_heat_capacity * etempint * 0.5 + j0 * expterm * 0.5;
      const double dj_mass_energydT_master =
          dj_dT_slave * molar_heat_capacity * etempint * 0.5 - j0 * expterm * 0.5;

      djdT_slave_timefacfac += pseudo_contact_fac * dj_mass_energydT_slave * timefacfac;
      djdT_master_timefacfac += pseudo_contact_fac * dj_mass_energydT_master * timefacfac;

      // Part 2: Energy flux from temperature drop
      // core residual
      j_timefacrhsfac +=
          pseudo_contact_fac * timefacrhsfac * thermoperm * (eslavetempint - emastertempint);

      // core linearizations
      djdT_slave_timefacfac += pseudo_contact_fac * timefacfac * thermoperm;
      djdT_master_timefacfac -= pseudo_contact_fac * timefacfac * thermoperm;

      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int ui = 0; ui < nen_; ++ui)
          k_ss(vi, ui) += funct_slave(vi) * djdT_slave_timefacfac * funct_slave(ui);
        r_s[vi] -= funct_slave(vi) * j_timefacrhsfac;

        for (int ui = 0; ui < nen_master; ++ui)
          k_sm(vi, ui) += funct_slave(vi) * djdT_master_timefacfac * funct_master(ui);
      }

      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreduced:
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    case INPAR::S2I::kinetics_nointerfaceflux:
    {
      // do nothing
      break;
    }

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix)
{
  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret =
      Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::FourierIso> matfourier =
      Teuchos::rcp_dynamic_cast<const MAT::FourierIso>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if ((matsoret == Teuchos::null and matfourier == Teuchos::null) or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterscatra(2, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterscatra, discretization, la, "imasterscatra", 2);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();

  CORE::LINALG::Matrix<nen_, 1> emastertemp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
    my::ExtractNodeValues(emastertemp, discretization, la, "imastertemp", 3);

  // get primary variable to derive the linearization
  const auto differentiationtype =
      Teuchos::getIntegralValue<SCATRA::DifferentiationType>(params, "differentiationtype");

  // element slave mechanical stress tensor
  const bool is_pseudo_contact = my::scatraparamsboundary_->IsPseudoContact();
  std::vector<CORE::LINALG::Matrix<nen_, 1>> eslavestress_vector(
      6, CORE::LINALG::Matrix<nen_, 1>(true));
  if (is_pseudo_contact)
    my::ExtractNodeValues(eslavestress_vector, discretization, la, "mechanicalStressState",
        my::scatraparams_->NdsTwoTensorQuantity());

  CORE::LINALG::Matrix<nsd_, 1> normal;

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid, &normal);
    const double detF = my::CalculateDetFOfParentElement(ele, intpoints.Point(gpid));

    const double pseudo_contact_fac = my::CalculatePseudoContactFactor(
        is_pseudo_contact, eslavestress_vector, normal, my::funct_);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    if (timefacfac < 0.0) dserror("Integration factor is negative!");

    const double timefacwgt = my::scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];

    static CORE::LINALG::Matrix<nsd_, nen_> dsqrtdetg_dd;
    if (differentiationtype == SCATRA::DifferentiationType::disp)
    {
      static CORE::LINALG::Matrix<nen_, nsd_> xyze_transposed;
      xyze_transposed.UpdateT(my::xyze_);
      DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<distype, nsd_>(
          my::derxy_, my::deriv_, xyze_transposed, normal);
      my::EvaluateSpatialDerivativeOfAreaIntegrationFactor(intpoints, gpid, dsqrtdetg_dd);
    }

    EvaluateS2ICouplingODAtIntegrationPoint<distype>(matelectrode, my::ephinp_[0], emastertemp,
        eelchnp_, emasterscatra, pseudo_contact_fac, my::funct_, my::funct_,
        my::scatraparamsboundary_, timefacfac, timefacwgt, detF, differentiationtype, dsqrtdetg_dd,
        my::derxy_, eslavematrix, emastermatrix);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype,
    probdim>::EvaluateS2ICouplingODAtIntegrationPoint(const Teuchos::RCP<const MAT::Electrode>&
                                                          matelectrode,
    const CORE::LINALG::Matrix<nen_, 1>& eslavetempnp,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        emastertempnp,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& eslavephinp,
    const std::vector<CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const double pseudo_contact_fac, const CORE::LINALG::Matrix<nen_, 1>& funct_slave,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
    const double timefacfac, const double timefacwgt, const double detF,
    const SCATRA::DifferentiationType differentiationtype,
    const CORE::LINALG::Matrix<nsd_, nen_>& dsqrtdetg_dd,
    const CORE::LINALG::Matrix<nsd_, nen_>& shape_spatial_derivatives,
    Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm)
{
  // get condition specific parameters
  const int kineticmodel = scatra_parameter_boundary->KineticModel();
  const double kr = scatra_parameter_boundary->ChargeTransferConstant();
  const double alphaa = scatra_parameter_boundary->AlphaA();
  const double alphac = scatra_parameter_boundary->AlphaC();
  const double peltier = scatra_parameter_boundary->Peltier();
  const double thermoperm = scatra_parameter_boundary->ThermoPerm();
  const double molar_heat_capacity = scatra_parameter_boundary->MolarHeatCapacity();

  // number of nodes of master-side element
  const int nen_master =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  if (eslavetempint <= 0.0) dserror("Temperature is non-positive!");
  const double emastertempint = funct_master.Dot(emastertempnp);
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // compute derivatives of scatra-scatra interface coupling residuals w.r.t. concentration and
  // electric potential according to kinetic model for current thermo-thermo interface coupling
  // condition
  switch (kineticmodel)
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      switch (differentiationtype)
      {
        case SCATRA::DifferentiationType::elch:
        {
          // access input parameters associated with current condition
          const double faraday =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
          DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
          const double gasconstant =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

          // extract saturation value of intercalated lithium concentration from electrode material
          const double cmax = matelectrode->CMax();

          // evaluate factor F/RT
          const double frt = faraday / (gasconstant * eslavetempint);

          // equilibrium electric potential difference at electrode surface and its derivative
          // w.r.t. concentration at electrode surface
          const double epd =
              matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);
          const double epdderiv = matelectrode->ComputeDOpenCircuitPotentialDConcentration(
              eslavephiint, faraday, frt, detF);

          // electrode-electrolyte overpotential at integration point
          const double eta = eslavepotint - emasterpotint - epd;

          // Butler-Volmer exchange current density
          const double i0 = kr * faraday * std::pow(emasterphiint, alphaa) *
                            std::pow(cmax - eslavephiint, alphaa) * std::pow(eslavephiint, alphac);

          // exponential Butler-Volmer terms
          const double expterm1 = std::exp(alphaa * frt * eta);
          const double expterm2 = std::exp(-alphac * frt * eta);
          const double expterm = expterm1 - expterm2;

          // safety check
          if (std::abs(expterm) > 1.e5)
          {
            dserror(
                "Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
                expterm);
          }

          // core linearizations w.r.t. master-side and slave-side concentrations and electric
          // potentials
          const double dres_dc_slave =
              (kr * faraday * std::pow(emasterphiint, alphaa) *
                      std::pow(cmax - eslavephiint, alphaa - 1.) *
                      std::pow(eslavephiint, alphac - 1.) *
                      (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
                  i0 * (-alphaa * frt * epdderiv * expterm1 - alphac * frt * epdderiv * expterm2)) *
                  (eta + peltier) -
              i0 * expterm * epdderiv;
          const double dres_dc_slave_timefacfac = pseudo_contact_fac * timefacfac * dres_dc_slave;

          const double dres_dc_master = i0 * alphaa / emasterphiint * expterm * (eta + peltier);
          const double dres_dc_master_timefacfac = pseudo_contact_fac * timefacfac * dres_dc_master;

          const double dres_dpot_slave =
              i0 * frt * (alphaa * expterm1 + alphac * expterm2) * (eta + peltier) + i0 * expterm;
          const double dres_dpot_slave_timefacfac =
              pseudo_contact_fac * timefacfac * dres_dpot_slave;

          const double dres_dpot_master_timefacfac = -dres_dpot_slave_timefacfac;

          // compute matrix contributions associated with slave-side residuals
          for (int vi = 0; vi < nen_; ++vi)
          {
            for (int ui = 0; ui < nen_; ++ui)
            {
              // compute linearizations w.r.t. slave-side concentrations
              k_ss(vi, ui * 2) -= funct_slave(vi) * dres_dc_slave_timefacfac * funct_slave(ui);

              // compute linearizations w.r.t. slave-side electric potentials
              k_ss(vi, ui * 2 + 1) -=
                  funct_slave(vi) * dres_dpot_slave_timefacfac * funct_slave(ui);
            }

            for (int ui = 0; ui < nen_master; ++ui)
            {
              // compute linearizations w.r.t. master-side concentrations
              k_sm(vi, ui * 2) -= funct_slave(vi) * dres_dc_master_timefacfac * funct_master(ui);

              // compute linearizations w.r.t. master-side electric potentials
              k_sm(vi, ui * 2 + 1) -=
                  funct_slave(vi) * dres_dpot_master_timefacfac * funct_master(ui);
            }
          }
          break;
        }
        default:
        {
          dserror("Unknown primary quantity to calculate derivative");
        }
      }
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreducedthermoresistance:
    {
      switch (differentiationtype)
      {
        case SCATRA::DifferentiationType::disp:
        {
          // Part 1
          // compute factor F/(RT)
          const double gasconstant =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
          const double faraday =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
          const double etempint = (eslavetempint + emastertempint) / 2;
          const double frt = faraday / (etempint * gasconstant);

          // equilibrium electric potential difference at electrode surface
          const double epd =
              matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);

          // skip further computation in case equilibrium electric potential difference is outside
          // physically meaningful range
          if (not std::isinf(epd))
          {
            const double depd_ddetF =
                matelectrode->ComputeDOpenCircuitPotentialDDetF(eslavephiint, faraday, frt, detF);

            // electrode-electrolyte overpotential at integration point
            const double eta = eslavepotint - emasterpotint - epd;

            // Butler-Volmer exchange mass flux density
            const double j0 = kr;

            // exponential Butler-Volmer terms
            const double expterm1 = std::exp(alphaa * frt * eta);
            const double expterm2 = std::exp(-alphac * frt * eta);
            const double expterm = expterm1 - expterm2;

            // safety check
            if (std::abs(expterm) > 1.0e5)
            {
              dserror(
                  "Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
                  expterm);
            }

            // core linearization associated with Butler-Volmer mass flux density
            const double dj_dsqrtdetg_timefacwgt =
                pseudo_contact_fac * timefacwgt * j0 * expterm * molar_heat_capacity * etempint;
            const double dj_depd = -j0 * frt * (alphaa * expterm1 + alphac * expterm2) *
                                   molar_heat_capacity * etempint;
            const double dj_ddetF = dj_depd * depd_ddetF;
            const double dj_ddetF_timefacfac = pseudo_contact_fac * dj_ddetF * timefacfac;

            // loop over matrix columns
            for (int ui = 0; ui < nen_; ++ui)
            {
              const int fui = ui * 3;

              // loop over matrix rows
              for (int vi = 0; vi < nen_; ++vi)
              {
                const double vi_dj_dsqrtdetg = funct_slave(vi) * dj_dsqrtdetg_timefacwgt;
                const double vi_dj_ddetF = funct_slave(vi) * dj_ddetF_timefacfac;

                // loop over spatial dimensions
                for (int dim = 0; dim < 3; ++dim)
                {
                  // compute linearizations w.r.t. slave-side structural displacements
                  k_ss(vi, fui + dim) += vi_dj_dsqrtdetg * dsqrtdetg_dd(dim, ui);
                  k_ss(vi, fui + dim) += vi_dj_ddetF * detF * shape_spatial_derivatives(dim, ui);
                }
              }
            }
          }
          // Part 2
          const double dj_dsqrtdetg_timefacwgt =
              pseudo_contact_fac * timefacwgt * (eslavetempint - emastertempint) * thermoperm;

          // loop over matrix columns
          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * 3;

            // loop over matrix rows
            for (int vi = 0; vi < nen_; ++vi)
            {
              const double vi_dj_dsqrtdetg = funct_slave(vi) * dj_dsqrtdetg_timefacwgt;

              // loop over spatial dimensions
              for (int dim = 0; dim < 3; ++dim)
              {
                // finalize linearizations w.r.t. slave-side structural displacements
                k_ss(vi, fui + dim) += vi_dj_dsqrtdetg * dsqrtdetg_dd(dim, ui);
              }
            }
          }
          break;
        }
        case SCATRA::DifferentiationType::elch:
        {
          const double gasconstant =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
          const double faraday =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
          const double etempint = (eslavetempint + emastertempint) / 2;
          const double frt = faraday / (etempint * gasconstant);

          const double epd =
              matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);

          const double epdderiv = matelectrode->ComputeDOpenCircuitPotentialDConcentration(
              eslavephiint, faraday, frt, detF);

          const double cmax = matelectrode->CMax();

          if (not std::isinf(epd))
          {
            const double eta = eslavepotint - emasterpotint - epd;

            const double expterm1 = std::exp(alphaa * frt * eta);
            const double expterm2 = std::exp(-alphac * frt * eta);

            // Butler-Volmer exchange mass flux density
            const double j0 = kr;

            // forward declarations
            double dj_dc_slave(0.0);
            double dj_dc_master(0.0);
            double dj_dpot_slave(0.0);
            double dj_dpot_master(0.0);

            // calculate linearizations of Butler-Volmer kinetics w.r.t. elch dofs
            myelectrodeutils::CalculateButlerVolmerElchLinearizations(kineticmodel, j0, frt,
                epdderiv, alphaa, alphac, 0.0, expterm1, expterm2, kr, faraday, emasterphiint,
                eslavephiint, cmax, eta, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

            const double dj_energydc_slave =
                pseudo_contact_fac * dj_dc_slave * molar_heat_capacity * etempint;
            const double dj_energydpot_slave =
                pseudo_contact_fac * dj_dpot_slave * molar_heat_capacity * etempint;
            const double dj_energydc_master =
                pseudo_contact_fac * dj_dc_master * molar_heat_capacity * etempint;
            const double dj_energydpot_master =
                pseudo_contact_fac * dj_dpot_master * molar_heat_capacity * etempint;
            // compute matrix contributions associated with slave-side residuals
            for (int vi = 0; vi < nen_; ++vi)
            {
              for (int ui = 0; ui < nen_; ++ui)
              {
                k_ss(vi, ui * 2) -= funct_slave(vi) * dj_energydc_slave * funct_slave(ui);
                // compute linearizations w.r.t. slave-side electric potentials
                k_ss(vi, ui * 2 + 1) -= funct_slave(vi) * dj_energydpot_slave * funct_slave(ui);
              }

              for (int ui = 0; ui < nen_master; ++ui)
              {
                k_sm(vi, ui * 2) -= funct_slave(vi) * dj_energydc_master * funct_master(ui);
                // compute linearizations w.r.t. master-side electric potentials
                k_sm(vi, ui * 2 + 1) -= funct_slave(vi) * dj_energydpot_master * funct_master(ui);
              }
            }
          }
          break;
        }
        default:
        {
          dserror("Unknown type of primary variable");
        }
      }
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreduced:
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    case INPAR::S2I::kinetics_nointerfaceflux:
    {
      // do nothing
      break;
    }

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_s2icoupling:
    {
      EvaluateS2ICoupling(
          ele, params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra);
      break;
    }

    case SCATRA::BoundaryAction::calc_s2icoupling_od:
    {
      EvaluateS2ICouplingOD(ele, params, discretization, la, elemat1_epetra, elemat2_epetra);
      break;
    }

    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::ExtractNodeValues(
    const DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // call base class routine
  my::ExtractNodeValues(discretization, la);

  // extract nodal electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
  my::ExtractNodeValues(eelchnp_, discretization, la, "scatra", 2);
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs9, 3>;

// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, const SCATRA::DifferentiationType, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, const SCATRA::DifferentiationType, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, const SCATRA::DifferentiationType, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const double, const double,
        const double, const SCATRA::DifferentiationType, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
