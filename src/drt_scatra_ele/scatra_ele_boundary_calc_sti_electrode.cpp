/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_sti_electrode.cpp

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_sti_electrode.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/soret.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname,
    const ScaTraEleBoundaryCalcSTIElectrode* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcSTIElectrode<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleBoundaryCalcSTIElectrode<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleBoundaryCalcSTIElectrode<distype>*>::iterator i =
             instances.begin();
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
 | singleton destruction                                     fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::ScaTraEleBoundaryCalcSTIElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname),

      // initialize member variable
      eelchnp_(2, LINALG::Matrix<my::nen_, 1>(true))
{
  return;
}


/*----------------------------------------------------------------------------------------------------------------------------*
 | evaluate main-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 08/15 |
 *----------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele,              ///< current boundary element
    Teuchos::ParameterList& params,           ///< parameter list
    DRT::Discretization& discretization,      ///< discretization
    DRT::Element::LocationArray& la,          ///< location array
    Epetra_SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
    Epetra_SerialDenseVector& eslaveresidual  ///< element residual for slave side
)
{
  // safety check
  if (my::numscal_ != 1 or my::numdofpernode_ != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret =
      Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if (matsoret == Teuchos::null or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterscatra(2, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterscatra, discretization, la, "imasterscatra", 2);

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

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
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingAtIntegrationPoint<distype>(*s2icondition, matelectrode, my::ephinp_[0],
        eelchnp_, emasterscatra, my::funct_, my::funct_, timefacfac, timefacrhsfac, eslavematrix,
        eslaveresidual);
  }  // loop over integration points

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICoupling


/*---------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition at integration point   fang 01/17 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::
    EvaluateS2ICouplingAtIntegrationPoint(
        DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
        const Teuchos::RCP<const MAT::Electrode>& matelectrode,  //!< electrode material
        const LINALG::Matrix<my::nen_, 1>&
            eslavetempnp,  //!< thermo state variables at slave-side nodes
        const std::vector<LINALG::Matrix<my::nen_, 1>>&
            eslavephinp,  //!< scatra state variables at slave-side nodes
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
            emasterphinp,  //!< scatra state variables at master-side nodes
        const LINALG::Matrix<my::nen_, 1>& funct_slave,  //!< slave-side shape function values
        const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,
            1>& funct_master,        //!< master-side shape function values
        const double timefacfac,     //!< time-integration factor times domain-integration factor
        const double timefacrhsfac,  //!< time-integration factor for right-hand side times
                                     //!< domain-integration factor
        Epetra_SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseVector& r_s  //!< slave-side residual vector
    )
{
  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  if (eslavetempint <= 0.) dserror("Temperature is non-positive!");
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // compute matrix and vector contributions according to kinetic model for current scatra-scatra
  // interface coupling condition
  switch (s2icondition.GetInt("kinetic model"))
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      // access input parameters associated with current condition
      const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
      const double gasconstant =
          DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
      const double alphaa = s2icondition.GetDouble("alpha_a");
      const double alphac = s2icondition.GetDouble("alpha_c");
      const double kr = s2icondition.GetDouble("k_r");
      if (kr < 0.) dserror("Charge transfer constant k_r is negative!");
      const double peltier = s2icondition.GetDouble("peltier");

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

      // Butler-Volmer exchange current density
      const double i0 = kr * faraday * pow(emasterphiint, alphaa) *
                        pow(cmax - eslavephiint, alphaa) * pow(eslavephiint, alphac);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
      const double expterm = expterm1 - expterm2;

      // safety check
      if (abs(expterm) > 1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // core residual term
      const double residual = timefacrhsfac * i0 * expterm * (eta + peltier);

      // core linearization w.r.t. temperature
      const double linearization = -timefacfac * i0 * frt / eslavetempint * eta *
                                   (alphaa * expterm1 + alphac * expterm2) * (eta + peltier);

      // compute matrix and vector contributions
      for (int vi = 0; vi < my::nen_; ++vi)
      {
        for (int ui = 0; ui < my::nen_; ++ui)
          k_ss(vi, ui) -= funct_slave(vi) * linearization * funct_slave(ui);
        r_s[vi] += funct_slave(vi) * residual;
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
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICouplingAtIntegrationPoint


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 08/15 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele,             ///< current boundary element
    Teuchos::ParameterList& params,          ///< parameter list
    DRT::Discretization& discretization,     ///< discretization
    DRT::Element::LocationArray& la,         ///< location array
    Epetra_SerialDenseMatrix& eslavematrix,  ///< element matrix for slave side
    Epetra_SerialDenseMatrix& emastermatrix  ///< element matrix for master side
)
{
  // safety check
  if (my::numscal_ != 1 or my::numdofpernode_ != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret =
      Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if (matsoret == Teuchos::null or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterscatra(2, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterscatra, discretization, la, "imasterscatra", 2);

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingODAtIntegrationPoint<distype>(*s2icondition, matelectrode, my::ephinp_[0],
        eelchnp_, emasterscatra, my::funct_, my::funct_, timefacfac, eslavematrix, emastermatrix);
  }  // loop over integration points

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition at integration point   fang 01/17 |
 *------------------------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::
    EvaluateS2ICouplingODAtIntegrationPoint(
        DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
        const Teuchos::RCP<const MAT::Electrode>& matelectrode,  //!< electrode material
        const LINALG::Matrix<my::nen_, 1>&
            eslavetempnp,  //!< thermo state variables at slave-side nodes
        const std::vector<LINALG::Matrix<my::nen_, 1>>&
            eslavephinp,  //!< scatra state variables at slave-side nodes
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
            emasterphinp,  //!< scatra state variables at master-side nodes
        const LINALG::Matrix<my::nen_, 1>& funct_slave,  //!< slave-side shape function values
        const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,
            1>& funct_master,     //!< master-side shape function values
        const double timefacfac,  //!< time-integration factor times domain-integration factor
        Epetra_SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&
            k_sm  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    )
{
  // number of nodes of master-side element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavetempint = funct_slave.Dot(eslavetempnp);
  if (eslavetempint <= 0.) dserror("Temperature is non-positive!");
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // compute derivatives of scatra-scatra interface coupling residuals w.r.t. concentration and
  // electric potential according to kinetic model for current thermo-thermo interface coupling
  // condition
  switch (s2icondition.GetInt("kinetic model"))
  {
    // Butler-Volmer-Peltier kinetics
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      // access input parameters associated with current condition
      const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
      const double gasconstant =
          DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
      const double alphaa = s2icondition.GetDouble("alpha_a");
      const double alphac = s2icondition.GetDouble("alpha_c");
      const double kr = s2icondition.GetDouble("k_r");
      if (kr < 0.) dserror("Charge transfer constant k_r is negative!");
      const double peltier = s2icondition.GetDouble("peltier");

      // extract saturation value of intercalated lithium concentration from electrode material
      const double cmax = matelectrode->CMax();
      if (cmax < 1.e-12)
        dserror("Saturation value c_max of intercalated lithium concentration is too small!");

      // evaluate factor F/RT
      const double frt = faraday / (gasconstant * eslavetempint);

      // equilibrium electric potential difference at electrode surface and its derivative w.r.t.
      // concentration at electrode surface
      const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);
      const double epdderiv =
          matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint, faraday, frt);

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
      if (abs(expterm) > 1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // core linearizations w.r.t. master-side and slave-side concentrations and electric
      // potentials
      const double dres_dc_slave =
          timefacfac *
          ((kr * faraday * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.) *
                   pow(eslavephiint, alphac - 1.) *
                   (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
               i0 * (-alphaa * frt * epdderiv * expterm1 - alphac * frt * epdderiv * expterm2)) *
                  (eta + peltier) -
              i0 * expterm * epdderiv);
      const double dres_dc_master =
          timefacfac * i0 * alphaa / emasterphiint * expterm * (eta + peltier);
      const double dres_dpot_slave =
          timefacfac *
          (i0 * frt * (alphaa * expterm1 + alphac * expterm2) * (eta + peltier) + i0 * expterm);
      const double dres_dpot_master = -dres_dpot_slave;

      // compute matrix contributions associated with slave-side residuals
      for (int vi = 0; vi < my::nen_; ++vi)
      {
        for (int ui = 0; ui < my::nen_; ++ui)
        {
          // compute linearizations w.r.t. slave-side concentrations
          k_ss(vi, ui * 2) -= funct_slave(vi) * dres_dc_slave * funct_slave(ui);

          // compute linearizations w.r.t. slave-side electric potentials
          k_ss(vi, ui * 2 + 1) -= funct_slave(vi) * dres_dpot_slave * funct_slave(ui);
        }

        for (int ui = 0; ui < nen_master; ++ui)
        {
          // compute linearizations w.r.t. master-side concentrations
          k_sm(vi, ui * 2) -= funct_slave(vi) * dres_dc_master * funct_master(ui);

          // compute linearizations w.r.t. master-side electric potentials
          k_sm(vi, ui * 2 + 1) -= funct_slave(vi) * dres_dpot_master * funct_master(ui);
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
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICouplingODAtIntegrationPoint


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateAction(
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
    case SCATRA::bd_calc_s2icoupling:
    {
      EvaluateS2ICoupling(ele, params, discretization, la, elemat1_epetra, elevec1_epetra);
      break;
    }

    case SCATRA::bd_calc_s2icoupling_od:
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
  }  // switch action

  return 0;
}


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::ExtractNodeValues(
    const DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la             //!< location array
)
{
  // call base class routine
  my::ExtractNodeValues(discretization, la);

  // extract nodal electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
  my::ExtractNodeValues(eelchnp_, discretization, la, "scatra", 2);

  return;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs9>;

// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,
        const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
    DRT::Element::quad4>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
    const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
    const std::vector<LINALG::Matrix<my::nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
    const LINALG::Matrix<my::nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
    DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,
    const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
    const std::vector<LINALG::Matrix<my::nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
    const LINALG::Matrix<my::nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,
        1>&,
    const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
    DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
    const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
    const std::vector<LINALG::Matrix<my::nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
    const LINALG::Matrix<my::nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,
        const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
        const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,
        const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>&,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>::
    EvaluateS2ICouplingODAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
        const Teuchos::RCP<const MAT::Electrode>&, const LINALG::Matrix<my::nen_, 1>&,
        const std::vector<LINALG::Matrix<my::nen_, 1>>&,
        const std::vector<LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
        const LINALG::Matrix<my::nen_, 1>&,
        const LINALG::Matrix<
            DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>&,
        const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&);
