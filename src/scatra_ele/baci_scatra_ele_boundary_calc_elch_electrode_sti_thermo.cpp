/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"

#include "baci_scatra_ele_boundary_calc_elch_electrode_utils.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"

#include "baci_mat_electrode.H"

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
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
    probdim>::ScaTraEleBoundaryCalcElchElectrodeSTIThermo(const int numdofpernode,
    const int numscal, const std::string& disname)
    :  // constructor of base class
      myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
    probdim>::EvaluateS2ICouplingOD(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& eslavematrix)
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

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const auto differentiationtype =
      Teuchos::getIntegralValue<SCATRA::DifferentiationType>(params, "differentiationtype");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> emasterphinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  CORE::LINALG::Matrix<nen_, 1> emastertempnp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
    my::ExtractNodeValues(emastertempnp, discretization, la, "imastertemp", 2);

  // element slave mechanical stress tensor
  const bool is_pseudo_contact = my::scatraparamsboundary_->IsPseudoContact();
  std::vector<CORE::LINALG::Matrix<nen_, 1>> eslavestress_vector(
      6, CORE::LINALG::Matrix<nen_, 1>(true));
  if (is_pseudo_contact)
    my::ExtractNodeValues(eslavestress_vector, discretization, la, "mechanicalStressState",
        my::scatraparams_->NdsTwoTensorQuantity());

  CORE::LINALG::Matrix<nsd_, 1> normal;

  // dummy element matrix
  Epetra_SerialDenseMatrix dummymatrix;

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  static CORE::LINALG::Matrix<nsd_, nen_> dsqrtdetg_dd;

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid, &normal);
    const double detF = my::CalculateDetFOfParentElement(ele, intpoints.Point(gpid));

    const double pseudo_contact_fac = my::CalculatePseudoContactFactor(
        is_pseudo_contact, eslavestress_vector, normal, my::funct_);

    if (differentiationtype == SCATRA::DifferentiationType::disp)
    {
      static CORE::LINALG::Matrix<nen_, nsd_> xyze_transposed;
      xyze_transposed.UpdateT(my::xyze_);
      DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<distype, nsd_>(
          my::derxy_, my::deriv_, xyze_transposed, normal);
      my::EvaluateSpatialDerivativeOfAreaIntegrationFactor(intpoints, gpid, dsqrtdetg_dd);
    }

    // evaluate overall integration factor
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    if (timefacfac < 0.0) dserror("Integration factor is negative!");

    const double timefacwgt = my::scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];
    if (timefacwgt < 0.0) dserror("Integration factor is negative!");

    EvaluateS2ICouplingODAtIntegrationPoint<distype>(matelectrode, my::ephinp_, etempnp_,
        emastertempnp, emasterphinp, pseudo_contact_fac, my::funct_, my::funct_, my::funct_,
        my::funct_, dsqrtdetg_dd, my::derxy_, my::scatraparamsboundary_, differentiationtype,
        timefacfac, timefacwgt, detF, my::numdofpernode_, eslavematrix, dummymatrix);
  }  // loop over integration points
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
   // probdim>::EvaluateS2ICouplingOD

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
    probdim>::EvaluateS2ICouplingODAtIntegrationPoint(const Teuchos::RCP<const MAT::Electrode>&
                                                          matelectrode,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& eslavephinp,
    const CORE::LINALG::Matrix<nen_, 1>& eslavetempnp,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        emastertempnp,
    const std::vector<CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const double pseudo_contact_fac, const CORE::LINALG::Matrix<nen_, 1>& funct_slave,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const CORE::LINALG::Matrix<nen_, 1>& test_slave,
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const CORE::LINALG::Matrix<nsd_, nen_>& dsqrtdetg_dd,
    const CORE::LINALG::Matrix<nsd_, nen_>& shape_spatial_derivatives,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
    const SCATRA::DifferentiationType differentiationtype, const double timefacfac,
    const double timefacwgt, const double detF, const int num_dof_per_node,
    Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_ms)
{
  // get condition specific parameters
  const int kineticmodel = scatra_parameter_boundary->KineticModel();
  const int numelectrons = scatra_parameter_boundary->NumElectrons();
  const double kr = scatra_parameter_boundary->ChargeTransferConstant();
  const double alphaa = scatra_parameter_boundary->AlphaA();
  const double alphac = scatra_parameter_boundary->AlphaC();

  // number of nodes of master-side element
  const int nen_master =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  const double emastertempint = funct_master.Dot(emastertempnp);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();

  // compute derivatives of scatra-scatra interface coupling residuals w.r.t. thermo dofs according
  // to kinetic model for current scatra-scatra interface coupling condition
  switch (kineticmodel)
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      switch (differentiationtype)
      {
        case SCATRA::DifferentiationType::temp:
        {
          // evaluate factor F/RT
          const double frt = faraday / (gasconstant * eslavetempint);

          // equilibrium electric potential difference at electrode surface
          const double epd =
              matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);

          // electrode-electrolyte overpotential at integration point
          const double eta = eslavepotint - emasterpotint - epd;

          // Butler-Volmer exchange mass flux density
          const double j0 = kr * std::pow(emasterphiint, alphaa) *
                            std::pow(cmax - eslavephiint, alphaa) * std::pow(eslavephiint, alphac);

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

          // linearization of Butler-Volmer mass flux density w.r.t. temperature
          const double dj_dT_timefacfac = -pseudo_contact_fac * timefacfac * j0 * frt /
                                          eslavetempint * eta *
                                          (alphaa * expterm1 + alphac * expterm2);

          // compute matrix contributions associated with slave-side residuals
          for (int vi = 0; vi < nen_; ++vi)
          {
            const int row_conc = vi * num_dof_per_node;

            for (int ui = 0; ui < nen_; ++ui)
            {
              // compute linearizations associated with slave-side equations for lithium transport
              k_ss(row_conc, ui) += test_slave(vi) * dj_dT_timefacfac * funct_slave(ui);

              // compute linearizations associated with slave-side closing equations for electric
              // potential
              k_ss(row_conc + 1, ui) +=
                  numelectrons * test_slave(vi) * dj_dT_timefacfac * funct_slave(ui);
            }
          }

          // compute matrix contributions associated with master-side residuals if necessary
          if (k_ms.M())
          {
            for (int vi = 0; vi < nen_master; ++vi)
            {
              const int row_conc = vi * num_dof_per_node;

              for (int ui = 0; ui < nen_; ++ui)
              {
                // compute linearizations associated with master-side equations for lithium
                // transport
                k_ms(row_conc, ui) -= test_master(vi) * dj_dT_timefacfac * funct_slave(ui);

                // compute linearizations associated with master-side closing equations for electric
                // potential
                k_ms(row_conc + 1, ui) -=
                    numelectrons * test_master(vi) * dj_dT_timefacfac * funct_slave(ui);
              }
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
      // average temperature at interface
      const double etempint = 0.5 * (eslavetempint + emastertempint);

      // factor FRT
      const double frt = faraday / (etempint * gasconstant);

      // equilibrium electric potential difference at electrode surface
      const double epd =
          matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt, detF);

      // skip further computation in case equilibrium electric potential difference is
      // outside physically meaningful range
      if (std::isinf(epd)) break;

      const double depd_ddetF =
          matelectrode->ComputeDOpenCircuitPotentialDDetF(eslavephiint, faraday, frt, detF);

      // Butler-Volmer exchange mass flux density
      const double j0 = kr;

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd;

      // dervivative of interface flux w.r.t. displacement
      switch (differentiationtype)
      {
        case SCATRA::DifferentiationType::disp:
        {
          double dj_dsqrtdetg(0.0), dj_ddetF(0.0);
          myelectrodeutils::CalculateButlerVolmerDispLinearizations(
              kineticmodel, alphaa, alphac, frt, j0, eta, depd_ddetF, dj_dsqrtdetg, dj_ddetF);

          const double dj_dsqrtdetg_timefacwgt = pseudo_contact_fac * dj_dsqrtdetg * timefacwgt;
          const double dj_ddetF_timefacfac = pseudo_contact_fac * dj_ddetF * timefacfac;

          // loop over matrix columns
          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * 3;

            // loop over matrix rows
            for (int vi = 0; vi < nen_; ++vi)
            {
              const int row_conc = vi * num_dof_per_node;
              const int row_pot = row_conc + 1;
              const double vi_dj_dsqrtdetg = test_slave(vi) * dj_dsqrtdetg_timefacwgt;
              const double vi_dj_ddetF = test_slave(vi) * dj_ddetF_timefacfac;

              // loop over spatial dimensions
              for (int dim = 0; dim < 3; ++dim)
              {
                // compute linearizations w.r.t. slave-side structural displacements
                k_ss(row_conc, fui + dim) += vi_dj_dsqrtdetg * dsqrtdetg_dd(dim, ui);
                k_ss(row_conc, fui + dim) +=
                    vi_dj_ddetF * detF * shape_spatial_derivatives(dim, ui);
                k_ss(row_pot, fui + dim) += numelectrons * vi_dj_dsqrtdetg * dsqrtdetg_dd(dim, ui);
                k_ss(row_pot, fui + dim) +=
                    numelectrons * vi_dj_ddetF * detF * shape_spatial_derivatives(dim, ui);
              }
            }
          }
          break;
        }
        case SCATRA::DifferentiationType::temp:
        {
          // derivative of epd w.r.t temperature
          const double depddT = matelectrode->ComputeDOpenCircuitPotentialDTemperature(
              eslavephiint, faraday, gasconstant);

          // forward declarations
          double dj_dT_slave(0.0);

          // calculate linearizations of Butler-Volmer kinetics w.r.t. tmperature dofs
          myelectrodeutils::CalculateButlerVolmerTempLinearizations(
              alphaa, alphac, depddT, eta, etempint, faraday, frt, gasconstant, j0, dj_dT_slave);

          const double djdT_slave_timefacfac = dj_dT_slave * timefacfac;

          // loop over matrix columns
          for (int ui = 0; ui < nen_; ++ui)
          {
            // loop over matrix rows
            for (int vi = 0; vi < nen_; ++vi)
            {
              const int row_conc = vi * num_dof_per_node;
              const int row_pot = row_conc + 1;
              const double vi_dj_dT_slave =
                  test_slave(vi) * pseudo_contact_fac * djdT_slave_timefacfac;

              // compute linearizations w.r.t. temperature
              k_ss(row_conc, ui) += vi_dj_dT_slave * funct_slave(ui);
              k_ss(row_pot, ui) += numelectrons * vi_dj_dT_slave * funct_slave(ui);
            }
          }
          break;
        }
        default:
        {
          dserror("Unknown differentiation type");
        }
      }
      break;
    }
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    case INPAR::S2I::kinetics_nointerfaceflux:
      break;

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
    }
  }  // select kinetic model
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
   // probdim>::EvaluateS2ICouplingODAtIntegrationPoint

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_s2icoupling_od:
    {
      EvaluateS2ICouplingOD(ele, params, discretization, la, elemat1_epetra);
      break;
    }

    default:
    {
      myelectrode::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch action

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
    probdim>::ExtractNodeValues(const DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // call base class routine
  my::ExtractNodeValues(discretization, la);

  // extract nodal temperature variables associated with time t_{n+1} or t_{n+alpha_f}
  my::ExtractNodeValues(etempnp_, discretization, la, "thermo", 2);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>::GetFRT() const
{
  // evaluate local temperature value
  const double temperature = my::funct_.Dot(etempnp_);

  // safety check
  if (temperature <= 0.) dserror("Temperature is non-positive!");

  const double faraday = myelch::elchparams_->Faraday();
  const double gasconstant = myelch::elchparams_->GasConstant();

  // evaluate factor F/RT
  return faraday / (gasconstant * temperature);
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs9, 3>;

// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, const int, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, const int, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, const int, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&,
        const std::vector<CORE::LINALG::Matrix<nen_, 1>>&, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nen_, 1>&,
        const CORE::LINALG::Matrix<
            CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const CORE::LINALG::Matrix<nsd_, nen_>&, const CORE::LINALG::Matrix<nsd_, nen_>&,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, const int, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
