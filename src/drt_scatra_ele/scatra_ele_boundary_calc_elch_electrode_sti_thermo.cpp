/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/soret.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../drt_scatra_ele/scatra_ele_parameter_boundary.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname,
    const ScaTraEleBoundaryCalcElchElectrodeSTIThermo* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string,
             ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
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
  return;
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

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterphinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // dummy element matrix
  Epetra_SerialDenseMatrix dummymatrix;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingODAtIntegrationPoint<distype>(matelectrode, my::ephinp_, etempnp_,
        emasterphinp, my::funct_, my::funct_, my::funct_, my::funct_, kineticmodel, numelectrons,
        stoichiometries, kr, alphaa, alphac, timefacfac, eslavematrix, dummymatrix);
  }  // loop over integration points

  return;
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
    const std::vector<
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,
    const LINALG::Matrix<my::nen_, 1>& funct_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const LINALG::Matrix<my::nen_, 1>& test_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const int kineticmodel, const int numelectrons, const std::vector<int>* stoichiometries,
    const double kr, const double alphaa, const double alphac, const double timefacfac,
    Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_ms)
{
  // number of nodes of master-side element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  if (eslavetempint <= 0.) dserror("Temperature is non-positive!");
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // compute derivatives of scatra-scatra interface coupling residuals w.r.t. thermo dofs according
  // to kinetic model for current scatra-scatra interface coupling condition
  switch (kineticmodel)
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      // get or check input parameters associated with current condition
      if (numelectrons != 1)
        dserror(
            "Invalid number of electrons involved in charge transfer at electrode-electrolyte "
            "interface!");
      if (stoichiometries == NULL)
        dserror(
            "Cannot access vector of stoichiometric coefficients for scatra-scatra interface "
            "coupling!");
      if (stoichiometries->size() != 1)
        dserror("Number of stoichiometric coefficients does not match number of scalars!");
      if ((*stoichiometries)[0] != -1) dserror("Invalid stoichiometric coefficient!");
      const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
      const double gasconstant =
          DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
      if (kr < 0.) dserror("Charge transfer constant k_r is negative!");

      // extract saturation value of intercalated lithium concentration from electrode material
      const double cmax = matelectrode->CMax();
      if (cmax < 1.e-12)
        dserror("Saturation value c_max of intercalated lithium concentration is too small!");

      // evaluate factor F/RT
      const double frt = faraday / (gasconstant * eslavetempint);

      // equilibrium electric potential difference at electrode surface
      const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint - emasterpotint - epd;

      // Butler-Volmer exchange mass flux density
      const double j0 = kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa) *
                        pow(eslavephiint, alphac);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
      const double expterm = expterm1 - expterm2;

      // safety check
      if (abs(expterm) > 1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // linearization of Butler-Volmer mass flux density w.r.t. temperature
      const double dj_dT =
          -timefacfac * j0 * frt / eslavetempint * eta * (alphaa * expterm1 + alphac * expterm2);

      // compute matrix contributions associated with slave-side residuals
      for (int vi = 0; vi < my::nen_; ++vi)
      {
        const int row_conc = vi * 2;

        for (int ui = 0; ui < my::nen_; ++ui)
        {
          // compute linearizations associated with slave-side equations for lithium transport
          k_ss(row_conc, ui) += test_slave(vi) * dj_dT * funct_slave(ui);

          // compute linearizations associated with slave-side closing equations for electric
          // potential
          k_ss(row_conc + 1, ui) += numelectrons * test_slave(vi) * dj_dT * funct_slave(ui);
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
            // compute linearizations associated with master-side equations for lithium transport
            k_ms(row_conc, ui) -= test_master(vi) * dj_dT * funct_slave(ui);

            // compute linearizations associated with master-side closing equations for electric
            // potential
            k_ms(row_conc + 1, ui) -= numelectrons * test_master(vi) * dj_dT * funct_slave(ui);
          }
        }
      }

      break;
    }

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
      break;
    }
  }  // select kinetic model

  return;
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
    case SCATRA::bd_calc_s2icoupling_od:
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

  return;
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
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(
        const Teuchos::RCP<const MAT::Electrode>&, const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
