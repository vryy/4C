/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_boundary.H"
#include "scatra_ele_boundary_calc_elch_electrode_utils.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_mat/electrode.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname,
    const ScaTraEleBoundaryCalcElchElectrode* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcElchElectrode<distype>*> instances;

  if (delete_me == nullptr)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleBoundaryCalcElchElectrode<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (auto& i : instances)
      if (i.second == delete_me)
      {
        delete i.second;
        instances.erase(i.first);
        return nullptr;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::ScaTraEleBoundaryCalcElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleBoundaryCalcElch(numdofpernode, numscal, disname)
{
}


/*-------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)   fang 04/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix,
    Epetra_SerialDenseVector& eslaveresidual)
{
  // safety check
  if (myelch::elchparams_->EquPot() != INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // get parameters for condition
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const double kr = my::scatraparamsboundary_->ChargeTransferConstant();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double resistance = my::scatraparamsboundary_->Resistance();
  const double itemaxmimplicitBV = my::scatraparamsboundary_->ItemaximplicitBV();
  const double convtolimplicitBV = my::scatraparamsboundary_->ConvtolimplicitBV();
  const bool isReducedBV =
      (kineticmodel == INPAR::S2I::kinetics_butlervolmerreduced or
          kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedwithresistance or
          kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance);

  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  this->ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterphinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  LINALG::Matrix<my::nen_, 1> eslavetempnp(true);
  LINALG::Matrix<my::nen_, 1> emastertempnp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
  {
    my::ExtractNodeValues(eslavetempnp, discretization, la, "islavetemp", 2);
    my::ExtractNodeValues(emastertempnp, discretization, la, "imastertemp", 2);
  }

  // dummy element matrix and vector
  Epetra_SerialDenseMatrix dummymatrix;
  Epetra_SerialDenseVector dummyvector;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
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

    EvaluateS2ICouplingAtIntegrationPoint<distype>(matelectrode, my::ephinp_, emasterphinp,
        eslavetempnp, emastertempnp, my::funct_, my::funct_, my::funct_, my::funct_, kineticmodel,
        numelectrons, kr, alphaa, alphac, resistance, itemaxmimplicitBV, convtolimplicitBV,
        timefacfac, timefacrhsfac, GetFRT(), isReducedBV, eslavematrix, emastermatrix, dummymatrix,
        dummymatrix, eslaveresidual, dummyvector);
  }  // loop over integration points
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling


/*---------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition at integration point   fang 05/16 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
    distype>::EvaluateS2ICouplingAtIntegrationPoint(const Teuchos::RCP<const MAT::Electrode>&
                                                        matelectrode,
    const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
    const std::vector<
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const LINALG::Matrix<my::nen_, 1>& eslavetempnp,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        emastertempnp,
    const LINALG::Matrix<my::nen_, 1>& funct_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const LINALG::Matrix<my::nen_, 1>& test_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const int kineticmodel, const int numelectrons, const double kr, const double alphaa,
    const double alphac, const double resistance, const double itemaxmimplicitBV,
    const double convtolimplicitBV, const double timefacfac, const double timefacrhsfac, double frt,
    const bool isReducedBV, Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm,
    Epetra_SerialDenseMatrix& k_ms, Epetra_SerialDenseMatrix& k_mm, Epetra_SerialDenseVector& r_s,
    Epetra_SerialDenseVector& r_m)
{
  // number of nodes of master-side mortar element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  const double emastertempint = funct_master.Dot(emastertempnp);

  const double etempint = 0.5 * (eslavetempint + emastertempint);

  // get faraday constant
  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();

  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
  {
    const double gasconstant =
        DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
    frt = faraday / (etempint * gasconstant);
  }

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // equilibrium electric potential difference at electrode surface
  const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);

  // derivative of equilibrium electric potential difference w.r.t. concentration at electrode
  // surface
  const double epdderiv =
      matelectrode->ComputeFirstDerivOpenCircuitPotentialConc(eslavephiint, faraday, frt);

  // Butler-Volmer exchange mass flux density
  const double j0 = isReducedBV ? kr
                                : kr * std::pow(emasterphiint, alphaa) *
                                      std::pow(cmax - eslavephiint, alphaa) *
                                      std::pow(eslavephiint, alphac);

  // compute matrix and vector contributions according to kinetic model for current scatra-scatra
  // interface coupling condition
  switch (kineticmodel)
  {
    // Butler-Volmer kinetics
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    case INPAR::S2I::kinetics_butlervolmerreducedthermoresistance:
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (std::isinf(epd)) break;

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
      const double j = j0 * expterm;

      // forward declarations
      double dj_dc_slave(0.0);
      double dj_dc_master(0.0);
      double dj_dpot_slave(0.0);
      double dj_dpot_master(0.0);

      // calculate linearizations of Butler-Volmer kinetics w.r.t. elch dofs
      myelectrodeutils::CalculateButlerVolmerElchLinearizations(kineticmodel, j0, frt, epdderiv,
          alphaa, alphac, resistance, expterm1, expterm2, kr, faraday, emasterphiint, eslavephiint,
          cmax, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

      // calculate RHS and linearizations of master and slave-side residuals
      CalculateRHSandGlobalSystem<distype_master>(funct_slave, funct_master, test_slave,
          test_master, numelectrons, nen_master, timefacfac, timefacrhsfac, dj_dc_slave,
          dj_dc_master, dj_dpot_slave, dj_dpot_master, j, k_ss, k_sm, k_ms, k_mm, r_s, r_m);

      break;
    }

    case INPAR::S2I::kinetics_butlervolmerresistance:
    case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
    {
      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (std::isinf(epd)) break;

      // compute Butler-Volmer mass flux density via Newton-Raphson method
      const double j = myelectrodeutils::CalculateModifiedButlerVolmerMassFluxDensity(j0, alphaa,
          alphac, frt, eslavepotint, emasterpotint, epd, resistance, itemaxmimplicitBV,
          convtolimplicitBV, faraday);

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd - j * faraday * resistance;

      // exponential Butler-Volmer terms
      const double expterm1 = std::exp(alphaa * frt * eta);
      const double expterm2 = std::exp(-alphac * frt * eta);
      const double expterm = expterm1 - expterm2;

      // safety check
      if (std::abs(expterm) > 1.0e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // forward declarations
      double dj_dc_slave(0.0);
      double dj_dc_master(0.0);
      double dj_dpot_slave(0.0);
      double dj_dpot_master(0.0);

      // calculate linearizations of Butler-Volmer kinetics w.r.t. elch dofs
      myelectrodeutils::CalculateButlerVolmerElchLinearizations(kineticmodel, j0, frt, epdderiv,
          alphaa, alphac, resistance, expterm1, expterm2, kr, faraday, emasterphiint, eslavephiint,
          cmax, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

      // calculate RHS and linearizations of master and slave-side residuals
      CalculateRHSandGlobalSystem<distype_master>(funct_slave, funct_master, test_slave,
          test_master, numelectrons, nen_master, timefacfac, timefacrhsfac, dj_dc_slave,
          dj_dc_master, dj_dpot_slave, dj_dpot_master, j, k_ss, k_sm, k_ms, k_mm, r_s, r_m);

      break;
    }  // case INPAR::S2I::kinetics_butlervolmerresistance:

    case INPAR::S2I::kinetics_constantinterfaceresistance:
    {
      // core residual
      const double inv_massfluxresistance = 1.0 / (resistance * faraday);
      const double jtimefacrhsfac =
          timefacrhsfac * (eslavepotint - emasterpotint) * inv_massfluxresistance;

      // calculate core linearizations
      const double dj_dpot_slave_timefacfac = timefacfac * inv_massfluxresistance;
      const double dj_dpot_master_timefacfac = -dj_dpot_slave_timefacfac;

      // calculate RHS and linearizations of master and slave-side residuals
      if (k_ss.M() and k_sm.M() and r_s.Length())
      {
        for (int vi = 0; vi < my::nen_; ++vi)
        {
          const int row_pot = vi * 2 + 1;

          for (int ui = 0; ui < my::nen_; ++ui)
          {
            const int col_pot = ui * 2 + 1;

            k_ss(row_pot, col_pot) += test_slave(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
          }

          for (int ui = 0; ui < nen_master; ++ui)
          {
            const int col_pot = ui * 2 + 1;

            k_sm(row_pot, col_pot) += test_slave(vi) * dj_dpot_master_timefacfac * funct_master(ui);
          }

          r_s[row_pot] -= test_slave(vi) * jtimefacrhsfac;
        }
      }
      else if (k_ss.M() or k_sm.M() or r_s.Length())
        dserror("Must provide both slave-side matrices and slave-side vector or none of them!");

      if (k_ms.M() and k_mm.M() and r_m.Length())
      {
        for (int vi = 0; vi < nen_master; ++vi)
        {
          const int row_pot = vi * 2 + 1;

          for (int ui = 0; ui < my::nen_; ++ui)
          {
            const int col_pot = ui * 2 + 1;

            k_ms(row_pot, col_pot) -= test_master(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
          }

          for (int ui = 0; ui < nen_master; ++ui)
          {
            const int col_pot = ui * 2 + 1;

            k_mm(row_pot, col_pot) -=
                test_master(vi) * dj_dpot_master_timefacfac * funct_master(ui);
          }

          r_m[row_pot] += test_master(vi) * jtimefacrhsfac;
        }
      }
      else if (k_ms.M() or k_mm.M() or r_m.Length())
        dserror("Must provide both master-side matrices and master-side vector or none of them!");

      break;
    }  // case INPAR::S2I::kinetics_constantinterfaceresistance

    case INPAR::S2I::kinetics_nointerfaceflux:
    {
      // do nothing
      break;
    }  // case INPAR::S2I::kinetics_nointerfaceflux

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
      break;
    }
  }  // switch(kineticmodel)
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 11/17 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& eslavematrix)
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // get condition specific parameters
  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int differentiationtype =
      params.get<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::none));
  const bool isReducedBV =
      (kineticmodel == INPAR::S2I::kinetics_butlervolmerreduced or
          kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedwithresistance or
          kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance);

  // extract local nodal values on present and opposite side of scatra-scatra interface
  this->ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterphinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  LINALG::Matrix<my::nen_, 1> eslavetempnp(true);
  LINALG::Matrix<my::nen_, 1> emastertempnp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
  {
    my::ExtractNodeValues(eslavetempnp, discretization, la, "islavetemp", 2);
    my::ExtractNodeValues(emastertempnp, discretization, la, "imastertemp", 2);
  }

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate shape derivatives
    static LINALG::Matrix<my::nsd_ + 1, my::nen_> shapederivatives;
    my::EvalShapeDerivatives(shapederivatives);

    // evaluate overall integration factors
    const double timefacwgt = my::scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];
    if (timefacwgt < 0.0) dserror("Integration factor is negative!");
    const double timefacfac =
        my::scatraparamstimint_->TimeFac() * my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);
    const double etempslaveint = my::funct_.Dot(eslavetempnp);
    const double etempmasterint = my::funct_.Dot(emastertempnp);

    const double etempint = 0.5 * (etempslaveint + etempmasterint);

    // compute matrix and vector contributions according to kinetic
    // model for current scatra-scatra interface coupling condition
    switch (kineticmodel)
    {
      // Butler-Volmer kinetics
      case INPAR::S2I::kinetics_butlervolmer:
      case INPAR::S2I::kinetics_butlervolmerreduced:
      case INPAR::S2I::kinetics_butlervolmerreducedthermoresistance:
      {
        // access input parameters associated with current condition
        const int numelectrons = my::scatraparamsboundary_->NumElectrons();
        const double faraday = myelch::elchparams_->Faraday();
        const double gasconstant = myelch::elchparams_->GasConstant();
        const double alphaa = my::scatraparamsboundary_->AlphaA();
        const double alphac = my::scatraparamsboundary_->AlphaC();
        const double kr = my::scatraparamsboundary_->ChargeTransferConstant();

        // extract saturation value of intercalated lithium concentration from electrode
        // material
        const double cmax = matelectrode->CMax();

        // compute factor F/(RT)
        const double frt =
            (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
                ? faraday / (etempint * gasconstant)
                : myelch::elchparams_->FRT();

        // equilibrium electric potential difference at electrode surface
        const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);

        // skip further computation in case equilibrium electric potential difference is
        // outside physically meaningful range
        if (std::isinf(epd)) break;

        // Butler-Volmer exchange mass flux density
        const double j0 = isReducedBV ? kr
                                      : kr * std::pow(emasterphiint, alphaa) *
                                            std::pow(cmax - eslavephiint, alphaa) *
                                            std::pow(eslavephiint, alphac);

        // electrode-electrolyte overpotential at integration point
        const double eta = eslavepotint - emasterpotint - epd;

        // dervivative of interface flux w.r.t. displacement
        switch (differentiationtype)
        {
          case static_cast<int>(SCATRA::DifferentiationType::disp):
          {
            // exponential Butler-Volmer terms
            const double expterm1 = std::exp(alphaa * frt * eta);
            const double expterm2 = std::exp(-alphac * frt * eta);
            const double expterm = expterm1 - expterm2;

            // safety check
            if (std::abs(expterm) > 1.0e5)
            {
              dserror(
                  "Overflow of exponential term in Butler-Volmer formulation detected! Value: "
                  "%lf",
                  expterm);
            }

            // core linearization associated with Butler-Volmer mass flux density
            const double dj_dd_slave_timefacwgt = timefacwgt * j0 * expterm;

            // loop over matrix columns
            for (int ui = 0; ui < my::nen_; ++ui)
            {
              const int fui = ui * 3;

              // loop over matrix rows
              for (int vi = 0; vi < my::nen_; ++vi)
              {
                const int row_conc = vi * 2;
                const int row_pot = row_conc + 1;
                const double vi_dj_dd_slave = my::funct_(vi) * dj_dd_slave_timefacwgt;

                // loop over spatial dimensions
                for (int dim = 0; dim < 3; ++dim)
                {
                  // compute linearizations w.r.t. slave-side structural displacements
                  eslavematrix(row_conc, fui + dim) += vi_dj_dd_slave * shapederivatives(dim, ui);
                  eslavematrix(row_pot, fui + dim) +=
                      numelectrons * vi_dj_dd_slave * shapederivatives(dim, ui);
                }
              }
            }
            break;
          }
          case static_cast<int>(SCATRA::DifferentiationType::temp):
          {
            if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
            {
              // derivative of epd w.r.t temperature
              const double depddT = matelectrode->ComputeFirstDerivOpenCircuitPotentialTemp(
                  eslavephiint, faraday, gasconstant);

              // forward declarations
              double dj_dT_slave(0.0);

              // calculate linearizations of Butler-Volmer kinetics w.r.t. tmperature dofs
              myelectrodeutils::CalculateButlerVolmerTempLinearizations(alphaa, alphac, depddT, eta,
                  etempint, faraday, frt, gasconstant, j0, dj_dT_slave);

              const double djdT_slave_timefacfac = dj_dT_slave * timefacfac;

              // loop over matrix columns
              for (int ui = 0; ui < my::nen_; ++ui)
              {
                // loop over matrix rows
                for (int vi = 0; vi < my::nen_; ++vi)
                {
                  const int row_conc = vi * 2;
                  const int row_pot = row_conc + 1;
                  const double vi_dj_dT_slave = my::funct_(vi) * djdT_slave_timefacfac;

                  // compute linearizations w.r.t. temperature
                  eslavematrix(row_conc, ui) += vi_dj_dT_slave * my::funct_(ui);
                  eslavematrix(row_pot, ui) += numelectrons * vi_dj_dT_slave * my::funct_(ui);
                }
              }
            }
            else
              dserror("Unknown kinetics type");
            break;
          }
          default:
          {
            dserror("Unknown differentiation type");
            break;
          }
        }
        break;
      }
      case INPAR::S2I::kinetics_constantinterfaceresistance:
      {
        switch (differentiationtype)
        {
          case static_cast<int>(SCATRA::DifferentiationType::disp):
          {
            // calculate linearizations
            const double inv_massfluxresistance =
                1.0 / (my::scatraparamsboundary_->Resistance() * myelch::elchparams_->Faraday());
            const double dj_dd_slave_timefacwgt =
                timefacwgt * (eslavepotint - emasterpotint) * inv_massfluxresistance;

            // loop over matrix columns
            for (int ui = 0; ui < my::nen_; ++ui)
            {
              const int fui = ui * 3;

              // loop over matrix rows
              for (int vi = 0; vi < my::nen_; ++vi)
              {
                const int row_pot = vi * 2 + 1;
                const double vi_dj_dd_slave = my::funct_(vi) * dj_dd_slave_timefacwgt;

                // loop over spatial dimensions
                for (int dim = 0; dim < 3; ++dim)
                {
                  // finalize linearizations w.r.t. slave-side structural displacements
                  eslavematrix(row_pot, fui + dim) += vi_dj_dd_slave * shapederivatives(dim, ui);
                }
              }
            }
            break;
          }
          case static_cast<int>(SCATRA::DifferentiationType::temp):
            break;
          default:
          {
            dserror("Unknown primary quantity to calculate derivative");
            break;
          }
        }

        break;
      }
      case INPAR::S2I::kinetics_nointerfaceflux:
      {
        // nothing to do
        break;
      }
      default:
      {
        dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
        break;
      }
    }  // switch(kineticmodel)
  }    // loop over integration points
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>& material, const int k) const
{
  // valence cannot be computed for electrode material
  dserror("Valence cannot be computed for electrode material!");

  return 0.0;
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetFRT() const
{
  // fetch factor F/RT from electrochemistry parameter list in isothermal case
  return myelch::elchparams_->FRT();
}

/*------------------------------------------------------------------------------------*
 | calculate RHS and global system                                      civaner 09/19 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::CalculateRHSandGlobalSystem(
    const LINALG::Matrix<my::nen_, 1>& funct_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const LINALG::Matrix<my::nen_, 1>& test_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const double numelectrons, const int nen_master, const double timefacfac,
    const double timefacrhsfac, const double dj_dc_slave, const double dj_dc_master,
    const double dj_dpot_slave, const double dj_dpot_master, const double j,
    Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm, Epetra_SerialDenseMatrix& k_ms,
    Epetra_SerialDenseMatrix& k_mm, Epetra_SerialDenseVector& r_s, Epetra_SerialDenseVector& r_m)
{
  // pre calculate integrand values
  const double jtimefacrhsfac = j * timefacrhsfac;
  const double dj_dc_slave_timefacfac = dj_dc_slave * timefacfac;
  const double dj_dc_master_timefacfac = dj_dc_master * timefacfac;
  const double dj_dpot_slave_timefacfac = dj_dpot_slave * timefacfac;
  const double dj_dpot_master_timefacfac = dj_dpot_master * timefacfac;

  // assemble slave side element rhs and linearizations
  if (k_ss.M() and k_sm.M() and r_s.Length())
  {
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      const int row_conc = vi * 2;
      const int row_pot = row_conc + 1;

      for (int ui = 0; ui < my::nen_; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_ss(row_conc, col_conc) += test_slave(vi) * dj_dc_slave_timefacfac * funct_slave(ui);
        k_ss(row_conc, col_pot) += test_slave(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
        k_ss(row_pot, col_conc) +=
            numelectrons * test_slave(vi) * dj_dc_slave_timefacfac * funct_slave(ui);
        k_ss(row_pot, col_pot) +=
            numelectrons * test_slave(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
      }

      for (int ui = 0; ui < nen_master; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_sm(row_conc, col_conc) += test_slave(vi) * dj_dc_master_timefacfac * funct_master(ui);
        k_sm(row_conc, col_pot) += test_slave(vi) * dj_dpot_master_timefacfac * funct_master(ui);
        k_sm(row_pot, col_conc) +=
            numelectrons * test_slave(vi) * dj_dc_master_timefacfac * funct_master(ui);
        k_sm(row_pot, col_pot) +=
            numelectrons * test_slave(vi) * dj_dpot_master_timefacfac * funct_master(ui);
      }

      r_s[row_conc] -= test_slave(vi) * jtimefacrhsfac;
      r_s[row_pot] -= numelectrons * test_slave(vi) * jtimefacrhsfac;
    }
  }
  else if (k_ss.M() or k_sm.M() or r_s.Length())
    dserror("Must provide both slave-side matrices and slave-side vector or none of them!");

  // assemble master side element rhs and linearizations
  if (k_ms.M() and k_mm.M() and r_m.Length())
  {
    for (int vi = 0; vi < nen_master; ++vi)
    {
      const int row_conc = vi * 2;
      const int row_pot = row_conc + 1;

      for (int ui = 0; ui < my::nen_; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_ms(row_conc, col_conc) -= test_master(vi) * dj_dc_slave_timefacfac * funct_slave(ui);
        k_ms(row_conc, col_pot) -= test_master(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
        k_ms(row_pot, col_conc) -=
            numelectrons * test_master(vi) * dj_dc_slave_timefacfac * funct_slave(ui);
        k_ms(row_pot, col_pot) -=
            numelectrons * test_master(vi) * dj_dpot_slave_timefacfac * funct_slave(ui);
      }

      for (int ui = 0; ui < nen_master; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_mm(row_conc, col_conc) -= test_master(vi) * dj_dc_master_timefacfac * funct_master(ui);
        k_mm(row_conc, col_pot) -= test_master(vi) * dj_dpot_master_timefacfac * funct_master(ui);
        k_mm(row_pot, col_conc) -=
            numelectrons * test_master(vi) * dj_dc_master_timefacfac * funct_master(ui);
        k_mm(row_pot, col_pot) -=
            numelectrons * test_master(vi) * dj_dpot_master_timefacfac * funct_master(ui);
      }

      r_m[row_conc] += test_master(vi) * jtimefacrhsfac;
      r_m[row_pot] += numelectrons * test_master(vi) * jtimefacrhsfac;
    }
  }
  else if (k_ms.M() or k_mm.M() or r_m.Length())
    dserror("Must provide both master-side matrices and master-side vector or none of them!");
}

// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>;
// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const int, const int, const double, const double, const double, const double, const double,
        const double, const double, const double, const double, const bool,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const int, const int, const double, const double, const double, const double, const double,
        const double, const double, const double, const double, const bool,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const int, const int, const double, const double, const double, const double, const double,
        const double, const double, const double, const double, const bool,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const int, const int, const double, const double, const double, const double, const double,
        const double, const double, const double, const double, const bool,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&);
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs9>;
