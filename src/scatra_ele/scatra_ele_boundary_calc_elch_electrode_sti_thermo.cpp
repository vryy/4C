/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"

#include "scatra_ele_boundary_calc_elch_electrode_utils.H"

#include "fem_general_utils_boundary_integration.H"

#include "mat_electrode.H"

#include "inpar_s2i.H"

#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_boundary.H"
#include "scatra_ele_parameter_std.H"
#include "headers_singleton_owner.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>>(
            new ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      ::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<
    distype>::ScaTraEleBoundaryCalcElchElectrodeSTIThermo(const int numdofpernode,
    const int numscal, const std::string& disname)
    :  // constructor of base class
      myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode, numscal, disname)
{
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 08/15 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele,            ///< current boundary element
    Teuchos::ParameterList& params,         ///< parameter list
    DRT::Discretization& discretization,    ///< discretization
    DRT::Element::LocationArray& la,        ///< location array
    Epetra_SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
)
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
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterphinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  LINALG::Matrix<my::nen_, 1> emastertempnp(true);
  if (kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)
    my::ExtractNodeValues(emastertempnp, discretization, la, "imastertemp", 2);

  // element slave mechanical stress tensor
  const bool is_pseudo_contact = my::scatraparamsboundary_->IsPseudoContact();
  std::vector<LINALG::Matrix<my::nen_, 1>> eslavestress_vector(
      6, LINALG::Matrix<my::nen_, 1>(true));
  if (is_pseudo_contact)
    my::ExtractNodeValues(eslavestress_vector, discretization, la, "mechanicalStressState",
        my::scatraparams_->NdsTwoTensorQuantity());

  LINALG::Matrix<my::nsd_ + 1, 1> normal;

  // dummy element matrix
  Epetra_SerialDenseMatrix dummymatrix;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // get shape derivatives when linearization w.r.t. displacement is calculated
  static LINALG::Matrix<my::nsd_ + 1, my::nen_> shapederivatives;
  if (differentiationtype == SCATRA::DifferentiationType::disp)
    my::EvalShapeDerivatives(shapederivatives);

  // get the node coordinates in material configuration (we have a nsd_+1 dimensional domain!)
  LINALG::Matrix<my::nsd_ + 1, my::nen_> XYZe;
  GEO::fillInitialPositionArray<distype, my::nsd_ + 1, LINALG::Matrix<my::nsd_ + 1, my::nen_>>(
      ele, XYZe);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid, &normal);
    const double detg =
        my::EvaluateSquareRootOfDeterminantOfMetricTensorAtIntPoint(intpoints, gpid, XYZe);
    const double detF = fac / detg;

    const double pseudo_contact_fac = my::CalculatePseudoContactFactor(
        is_pseudo_contact, eslavestress_vector, normal, my::funct_);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    const double timefacwgt = my::scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];
    if (timefacwgt < 0.0) dserror("Integration factor is negative!");

    EvaluateS2ICouplingODAtIntegrationPoint<distype>(matelectrode, my::ephinp_, etempnp_,
        emastertempnp, emasterphinp, pseudo_contact_fac, my::funct_, my::funct_, my::funct_,
        my::funct_, shapederivatives, my::scatraparamsboundary_, differentiationtype, timefacfac,
        timefacwgt, detF, eslavematrix, dummymatrix);
  }  // loop over integration points
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateS2ICouplingOD


/*------------------------------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition at integration point   fang 01/17 |
 *------------------------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<
    distype>::EvaluateS2ICouplingODAtIntegrationPoint(const Teuchos::RCP<const MAT::Electrode>&
                                                          matelectrode,
    const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
    const LINALG::Matrix<my::nen_, 1>& eslavetempnp,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        emastertempnp,
    const std::vector<
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const double pseudo_contact_fac, const LINALG::Matrix<my::nen_, 1>& funct_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const LINALG::Matrix<my::nen_, 1>& test_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const LINALG::Matrix<my::nsd_ + 1, my::nen_> shapederivatives,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
    const SCATRA::DifferentiationType differentiationtype, const double timefacfac,
    const double timefacwgt, const double detF, Epetra_SerialDenseMatrix& k_ss,
    Epetra_SerialDenseMatrix& k_ms)
{
  // get condition specific parameters
  const int kineticmodel = scatra_parameter_boundary->KineticModel();
  const int numelectrons = scatra_parameter_boundary->NumElectrons();
  const double kr = scatra_parameter_boundary->ChargeTransferConstant();
  const double alphaa = scatra_parameter_boundary->AlphaA();
  const double alphac = scatra_parameter_boundary->AlphaC();

  // number of nodes of master-side element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

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
          for (int vi = 0; vi < my::nen_; ++vi)
          {
            const int row_conc = vi * 2;

            for (int ui = 0; ui < my::nen_; ++ui)
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
              const int row_conc = vi * 2;

              for (int ui = 0; ui < my::nen_; ++ui)
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
          break;
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

      // Butler-Volmer exchange mass flux density
      const double j0 = kr;

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd;

      // dervivative of interface flux w.r.t. displacement
      switch (differentiationtype)
      {
        case SCATRA::DifferentiationType::disp:
        {
          double dj_dd_slave_timefacwgt(0.0);
          myelectrodeutils::CalculateButlerVolmerDispLinearizations(
              kineticmodel, alphaa, alphac, frt, j0, eta, timefacwgt, dj_dd_slave_timefacwgt);

          // loop over matrix columns
          for (int ui = 0; ui < my::nen_; ++ui)
          {
            const int fui = ui * 3;

            // loop over matrix rows
            for (int vi = 0; vi < my::nen_; ++vi)
            {
              const int row_conc = vi * 2;
              const int row_pot = row_conc + 1;
              const double vi_dj_dd_slave =
                  test_slave(vi) * pseudo_contact_fac * dj_dd_slave_timefacwgt;

              // loop over spatial dimensions
              for (int dim = 0; dim < 3; ++dim)
              {
                // compute linearizations w.r.t. slave-side structural displacements
                k_ss(row_conc, fui + dim) += vi_dj_dd_slave * shapederivatives(dim, ui);
                k_ss(row_pot, fui + dim) +=
                    numelectrons * vi_dj_dd_slave * shapederivatives(dim, ui);
              }
            }
          }
          break;
        }
        case SCATRA::DifferentiationType::temp:
        {
          // derivative of epd w.r.t temperature
          const double depddT = matelectrode->ComputeFirstDerivOpenCircuitPotentialTemp(
              eslavephiint, faraday, gasconstant);

          // forward declarations
          double dj_dT_slave(0.0);

          // calculate linearizations of Butler-Volmer kinetics w.r.t. tmperature dofs
          myelectrodeutils::CalculateButlerVolmerTempLinearizations(
              alphaa, alphac, depddT, eta, etempint, faraday, frt, gasconstant, j0, dj_dT_slave);

          const double djdT_slave_timefacfac = dj_dT_slave * timefacfac;

          // loop over matrix columns
          for (int ui = 0; ui < my::nen_; ++ui)
          {
            // loop over matrix rows
            for (int vi = 0; vi < my::nen_; ++vi)
            {
              const int row_conc = vi * 2;
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
          break;
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
      break;
    }
  }  // select kinetic model
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateS2ICouplingODAtIntegrationPoint


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateAction(
    DRT::FaceElement* ele,                     //!< boundary element
    Teuchos::ParameterList& params,            //!< parameter list
    DRT::Discretization& discretization,       //!< discretization
    SCATRA::BoundaryAction action,             //!< action
    DRT::Element::LocationArray& la,           //!< location array
    Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
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


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::ExtractNodeValues(
    const DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la             //!< location array
)
{
  // call base class routine
  my::ExtractNodeValues(discretization, la);

  // extract nodal temperature variables associated with time t_{n+1} or t_{n+alpha_f}
  my::ExtractNodeValues(etempnp_, discretization, la, "thermo", 2);
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::GetFRT() const
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
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs9>;

// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nsd_ + 1, my::nen_>,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nsd_ + 1, my::nen_>,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const double, const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nsd_ + 1, my::nen_>,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const double, const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nsd_ + 1, my::nen_>,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const, const SCATRA::DifferentiationType,
        const double, const double, const double, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&);
