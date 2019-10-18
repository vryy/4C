/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_boundary.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

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

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleBoundaryCalcElchElectrode<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleBoundaryCalcElchElectrode<distype>*>::iterator i =
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
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::ScaTraEleBoundaryCalcElchElectrode(
    const int numdofpernode, const int numscal,
    const std::string& disname)
    :  // constructor of base class
      myelch::ScaTraEleBoundaryCalcElch(numdofpernode, numscal, disname)
{
  return;
}


/*-------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)   fang 04/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele,              ///< current boundary element
    Teuchos::ParameterList& params,           ///< parameter list
    DRT::Discretization& discretization,      ///< discretization
    DRT::Element::LocationArray& la,          ///< location array
    Epetra_SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
    Epetra_SerialDenseMatrix& emastermatrix,  ///< element matrix for master side
    Epetra_SerialDenseVector& eslaveresidual  ///< element residual for slave side
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
  this->ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<my::nen_, 1>> emasterphinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  my::ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // dummy element matrix and vector
  Epetra_SerialDenseMatrix dummymatrix;
  Epetra_SerialDenseVector dummyvector;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double resistance = my::scatraparamsboundary_->Resistance();
  const double itemaxmimplicitBV = my::scatraparamsboundary_->ItemaximplicitBV();
  const double convtolimplicitBV = my::scatraparamsboundary_->ConvtolimplicitBV();

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingAtIntegrationPoint<distype>(matelectrode, my::ephinp_, emasterphinp,
        my::funct_, my::funct_, my::funct_, my::funct_, kineticmodel, numelectrons, stoichiometries,
        kr, alphaa, alphac, resistance, itemaxmimplicitBV, convtolimplicitBV, timefacfac,
        timefacrhsfac, GetFRT(), eslavematrix, emastermatrix, dummymatrix, dummymatrix,
        eslaveresidual, dummyvector);
  }  // loop over integration points

  return;
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
    const LINALG::Matrix<my::nen_, 1>& funct_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,
    const LINALG::Matrix<my::nen_, 1>& test_slave,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,
    const int kineticmodel, const int numelectrons, const std::vector<int>* stoichiometries,
    const double kr, const double alphaa, const double alphac, const double resistance,
    const double itemaxmimplicitBV, const double convtolimplicitBV, const double timefacfac,
    const double timefacrhsfac, const double frt, Epetra_SerialDenseMatrix& k_ss,
    Epetra_SerialDenseMatrix& k_sm, Epetra_SerialDenseMatrix& k_ms, Epetra_SerialDenseMatrix& k_mm,
    Epetra_SerialDenseVector& r_s, Epetra_SerialDenseVector& r_m)
{
  // safety checks
  // access input parameters associated with current condition
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
  if (kr < 0.) dserror("Charge transfer constant k_r is negative!");

  // number of nodes of master-side mortar element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra
  // interface
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // get faraday constant
  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();
  if (cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated lithium concentration is too small!");

  // equilibrium electric potential difference at electrode surface
  const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);

  // derivative of equilibrium electric potential difference w.r.t. concentration at electrode
  // surface
  const double epdderiv =
      matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint, faraday, frt);

  // Butler-Volmer exchange mass flux density
  const double j0(kineticmodel == INPAR::S2I::kinetics_butlervolmerreduced or
                          kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedwithresistance
                      ? kr
                      : kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa) *
                            pow(eslavephiint, alphac));

  // compute matrix and vector contributions according to kinetic model for current scatra-scatra
  // interface coupling condition
  switch (kineticmodel)
  {
    // Butler-Volmer kinetics
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (not std::isinf(epd))
      {
        // electrode-electrolyte overpotential at integration point
        const double eta = eslavepotint - emasterpotint - epd;

        // exponential Butler-Volmer terms
        const double expterm1 = exp(alphaa * frt * eta);
        const double expterm2 = exp(-alphac * frt * eta);
        const double expterm = expterm1 - expterm2;

        // safety check
        if (abs(expterm) > 1.e5)
          dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
              expterm);

        // core residual term associated with Butler-Volmer mass flux density
        const double jfacrhsfac = j0 * expterm * timefacrhsfac;

        // forward declarations
        double dj_dc_slave(.0);
        double dj_dc_master(.0);
        double dj_dpot_slave(.0);
        double dj_dpot_master(.0);

        // calculate core linearizations
        CalculateCoreLinearizations(kineticmodel, timefacfac, timefacrhsfac, j0, frt, epdderiv,
            alphaa, alphac, resistance, expterm1, expterm2, kr, faraday, emasterphiint,
            eslavephiint, cmax, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

        // calculate RHS and linearizations of master and slave-side residuals
        CalculateRHSandGlobalSystem<distype_master>(funct_slave, funct_master, test_slave,
            test_master, numelectrons, nen_master, dj_dc_slave, dj_dc_master, dj_dpot_slave,
            dj_dpot_master, jfacrhsfac, k_ss, k_sm, k_ms, k_mm, r_s, r_m);
      }

      break;
    }

    case INPAR::S2I::kinetics_butlervolmerresistance:
    case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
    {
      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (not std::isinf(epd))
      {
        // compute Butler-Volmer mass flux density via Newton-Raphson method
        const double j =
            CalculateModifiedButlerVolmerMassFluxDensity(j0, alphaa, alphac, frt, eslavepotint,
                emasterpotint, epd, resistance, itemaxmimplicitBV, convtolimplicitBV, faraday);

        // electrode-electrolyte overpotential at integration point
        const double eta = eslavepotint - emasterpotint - epd - j * faraday * resistance;

        // exponential Butler-Volmer terms
        const double expterm1 = exp(alphaa * frt * eta);
        const double expterm2 = exp(-alphac * frt * eta);
        const double expterm = expterm1 - expterm2;

        // safety check
        if (abs(expterm) > 1.e5)
          dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
              expterm);

        // forward declarations
        double dj_dc_slave(.0);
        double dj_dc_master(.0);
        double dj_dpot_slave(.0);
        double dj_dpot_master(.0);

        // calculate core residual term associated with Butler-Volmer mass flux density
        const double jfacrhsfac = j * timefacrhsfac;

        // calculate core linearizations
        CalculateCoreLinearizations(kineticmodel, timefacfac, timefacrhsfac, j0, frt, epdderiv,
            alphaa, alphac, resistance, expterm1, expterm2, kr, faraday, emasterphiint,
            eslavephiint, cmax, dj_dc_slave, dj_dc_master, dj_dpot_slave, dj_dpot_master);

        // calculate RHS and linearizations of master and slave-side residuals
        CalculateRHSandGlobalSystem<distype_master>(funct_slave, funct_master, test_slave,
            test_master, numelectrons, nen_master, dj_dc_slave, dj_dc_master, dj_dpot_slave,
            dj_dpot_master, jfacrhsfac, k_ss, k_sm, k_ms, k_mm, r_s, r_m);
      }

      break;
    }  // case INPAR::S2I::kinetics_butlervolmerresistance:

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
      break;
    }
  }  // switch(kineticmodel)

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 11/17 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele,            ///< current boundary element
    Teuchos::ParameterList& params,         ///< parameter list
    DRT::Discretization& discretization,    ///< discretization
    DRT::Element::LocationArray& la,        ///< location array
    Epetra_SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
)
{
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

    // evaluate overall integration factor
    const double timefacwgt = my::scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];
    if (timefacwgt < 0.) dserror("Integration factor is negative!");

    // evaluate dof values at current integration point on present and opposite side of
    // scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // compute matrix and vector contributions according to kinetic model for current scatra-scatra
    // interface coupling condition
    const int kineticmodel = my::scatraparamsboundary_->KineticModel();
    switch (kineticmodel)
    {
      // Butler-Volmer kinetics
      case INPAR::S2I::kinetics_butlervolmer:
      case INPAR::S2I::kinetics_butlervolmerreduced:
      {
        // access input parameters associated with current condition
        const int numelectrons = my::scatraparamsboundary_->NumElectrons();
        if (numelectrons != 1)
          dserror(
              "Invalid number of electrons involved in charge transfer at electrode-electrolyte "
              "interface: %i",
              numelectrons);
        const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
        if (stoichiometries == NULL)
          dserror(
              "Cannot access vector of stoichiometric coefficients for scatra-scatra interface "
              "coupling!");
        if (stoichiometries->size() != 1)
          dserror("Number of stoichiometric coefficients does not match number of scalars!");
        if ((*stoichiometries)[0] != -1) dserror("Invalid stoichiometric coefficient!");
        const double faraday = myelch::elchparams_->Faraday();
        const double alphaa = my::scatraparamsboundary_->AlphaA();
        const double alphac = my::scatraparamsboundary_->AlphaC();
        const double kr = my::scatraparamsboundary_->Kr();
        if (kr < 0.) dserror("Charge transfer constant k_r is negative!");

        // extract saturation value of intercalated lithium concentration from electrode material
        const double cmax = matelectrode->CMax();
        if (cmax < 1.e-12)
          dserror("Saturation value c_max of intercalated lithium concentration is too small!");

        // compute factor F/(RT)
        const double frt = GetFRT();

        // equilibrium electric potential difference at electrode surface
        const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint, faraday, frt);

        // skip further computation in case equilibrium electric potential difference is outside
        // physically meaningful range
        if (not std::isinf(epd))
        {
          // electrode-electrolyte overpotential at integration point
          const double eta = eslavepotint - emasterpotint - epd;

          // Butler-Volmer exchange mass flux density
          const double j0(kineticmodel == INPAR::S2I::kinetics_butlervolmerreduced
                              ? kr
                              : kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa) *
                                    pow(eslavephiint, alphac));

          // exponential Butler-Volmer terms
          const double expterm1 = exp(alphaa * frt * eta);
          const double expterm2 = exp(-alphac * frt * eta);
          const double expterm = expterm1 - expterm2;

          // safety check
          if (abs(expterm) > 1.e5)
            dserror(
                "Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
                expterm);

          // core linearization associated with Butler-Volmer mass flux density
          const double dj_dd_slave = timefacwgt * j0 * expterm;

          // loop over matrix columns
          for (int ui = 0; ui < my::nen_; ++ui)
          {
            const int fui = ui * 3;

            // loop over matrix rows
            for (int vi = 0; vi < my::nen_; ++vi)
            {
              const int row_conc = vi * 2;
              const int row_pot = row_conc + 1;
              const double vi_dj_dd_slave = my::funct_(vi) * dj_dd_slave;

              // loop over spatial dimensions
              for (unsigned dim = 0; dim < 3; ++dim)
              {
                // compute linearizations w.r.t. slave-side structural displacements
                eslavematrix(row_conc, fui + dim) += vi_dj_dd_slave * shapederivatives(dim, ui);
                eslavematrix(row_pot, fui + dim) +=
                    numelectrons * vi_dj_dd_slave * shapederivatives(dim, ui);
              }
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
    }  // switch(kineticmodel)
  }    // loop over integration points

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>& material,  // element material
    const int k                                         // species number
    ) const
{
  // valence cannot be computed for electrode material
  dserror("Valence cannot be computed for electrode material!");

  return 0.;
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetFRT() const
{
  // fetch factor F/RT from electrochemistry parameter list in isothermal case
  return myelch::elchparams_->FRT();
};

/*------------------------------------------------------------------------------------*
 | compute Butler-Volmer current density via Newton-Raphson iteration   civaner 08/19 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
    distype>::CalculateModifiedButlerVolmerMassFluxDensity(const double j0, const double alphaa,
    const double alphac, const double frt, const double pot_ed, const double pot_el,
    const double epd, const double resistance, const double itemax, const double convtol,
    const double faraday)
{
  // Iterations are conducted over current density i which is scaled down to mass flux
  // density j by j = i / faraday at the end of the function in order to reduce the effect
  // of numerical error introduced into global problem since i is roughly 10^5 times bigger than j

  // initialize Butler-Volmer current density
  double i(0.);

  const double i0 = j0 * faraday;
  // compute Butler-Volmer current density in case of physically reasonable half-cell open-circuit
  // potential
  if (not std::isinf(epd))
  {
    // initialize Newton-Raphson iteration counter
    unsigned iternum(0);

    // apply Newton-Raphson method to compute Butler-Volmer current density, involving overpotential
    // due to scatra-scatra interface layer resistance
    while (true)
    {
      // increment counter
      ++iternum;

      // compute current Newton-Raphson residual
      const double eta = pot_ed - pot_el - epd - resistance * i;
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
      const double residual = i0 * (expterm1 - expterm2) - i;

      // convergence check
      if (std::abs(residual) < convtol)
        break;

      else if (iternum == itemax)
        dserror(
            "Local Newton-Raphson iteration for Butler-Volmer current density did not converge!");

      // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current
      // density
      const double linearization =
          -i0 * resistance * frt * (alphaa * expterm1 + alphac * expterm2) - 1.;

      // update Butler-Volmer current density
      i -= residual / linearization;
    }
  }
  // final scaling
  return i / faraday;
}

/*------------------------------------------------------------------------------------*
 | calculate core linearizations associated with Butler-Volmer                        |
 | mass flux density                                                    civaner 09/19 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::CalculateCoreLinearizations(
    const int kineticmodel, const double timefacfac, const double timefacrhsfac, const double j0,
    const double frt, const double epdderiv, const double alphaa, const double alphac,
    const double resistance, const double expterm1, const double expterm2, const double kr,
    const double faraday, const double emasterphiint, const double eslavephiint, const double cmax,
    double& dj_dc_slave, double& dj_dc_master, double& dj_dpot_slave, double& dj_dpot_master)
{
  const double expterm = expterm1 - expterm2;
  // core linearizations associated with Butler-Volmer mass flux density
  switch (kineticmodel)
  {
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      dj_dc_slave = timefacfac * j0 * frt * epdderiv * (-alphaa * expterm1 - alphac * expterm2);
      dj_dc_master = 0.0;
      dj_dpot_slave = timefacfac * j0 * (alphaa * frt * expterm1 + alphac * frt * expterm2);
      dj_dpot_master = -dj_dpot_slave;
      break;
    }
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      dj_dc_slave =
          timefacfac * (kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.) *
                               pow(eslavephiint, alphac - 1.) *
                               (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
                           j0 * frt * epdderiv * (-alphaa * expterm1 - alphac * expterm2));
      dj_dc_master = timefacfac * j0 * alphaa / emasterphiint * expterm;
      dj_dpot_slave = timefacfac * j0 * (alphaa * frt * expterm1 + alphac * frt * expterm2);
      dj_dpot_master = -dj_dpot_slave;
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerresistance:
    {
      // core linearizations associated with Butler-Volmer current density according to MA Schmidt
      // 2016 via implicit differentiation where F(x,i) = i - i0 * expterm
      const double dF_di_inverse =
          1. / (1.0 + j0 * faraday * resistance * frt * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_slave =
          timefacfac * (-kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.) *
                               pow(eslavephiint, alphac - 1.) *
                               (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
                           j0 * frt * epdderiv * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_master = -timefacfac * j0 * alphaa / emasterphiint * expterm;
      const double dF_dpot_slave = -timefacfac * j0 * frt * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dpot_master = -dF_dpot_slave;

      // rule of implicit differentiation
      dj_dc_slave = -dF_dc_slave * dF_di_inverse;
      dj_dc_master = -dF_dc_master * dF_di_inverse;
      dj_dpot_slave = -dF_dpot_slave * dF_di_inverse;
      dj_dpot_master = -dF_dpot_master * dF_di_inverse;
      break;
    }  // case INPAR::S2I::kinetics_butlervolmerresistance
    case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
    {
      // core linearizations associated with Butler-Volmer current density according to MA Schmidt
      // 2016 via implicit differentiation where F(x,i) = i - i0 * expterm
      const double dF_di_inverse =
          1. / (1.0 + j0 * faraday * resistance * frt * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_slave =
          timefacfac * j0 * frt * epdderiv * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dc_master = 0.0;
      const double dF_dpot_slave = -timefacfac * j0 * frt * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dpot_master = -dF_dpot_slave;

      // rule of implicit differentiation
      dj_dc_slave = -dF_dc_slave * dF_di_inverse;
      dj_dc_master = -dF_dc_master * dF_di_inverse;
      dj_dpot_slave = -dF_dpot_slave * dF_di_inverse;
      dj_dpot_master = -dF_dpot_master * dF_di_inverse;
      break;
    }  // case INPAR::S2I::kinetics_butlervolmerreducedwithresistance
  }    // switch(kineticmodel)
  return;
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
    const double numelectrons, const int nen_master, const double dj_dc_slave,
    const double dj_dc_master, const double dj_dpot_slave, const double dj_dpot_master,
    const double j, Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm,
    Epetra_SerialDenseMatrix& k_ms, Epetra_SerialDenseMatrix& k_mm, Epetra_SerialDenseVector& r_s,
    Epetra_SerialDenseVector& r_m)
{
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

        k_ss(row_conc, col_conc) += test_slave(vi) * dj_dc_slave * funct_slave(ui);
        k_ss(row_conc, col_pot) += test_slave(vi) * dj_dpot_slave * funct_slave(ui);
        k_ss(row_pot, col_conc) += numelectrons * test_slave(vi) * dj_dc_slave * funct_slave(ui);
        k_ss(row_pot, col_pot) += numelectrons * test_slave(vi) * dj_dpot_slave * funct_slave(ui);
      }

      for (int ui = 0; ui < nen_master; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_sm(row_conc, col_conc) += test_slave(vi) * dj_dc_master * funct_master(ui);
        k_sm(row_conc, col_pot) += test_slave(vi) * dj_dpot_master * funct_master(ui);
        k_sm(row_pot, col_conc) += numelectrons * test_slave(vi) * dj_dc_master * funct_master(ui);
        k_sm(row_pot, col_pot) += numelectrons * test_slave(vi) * dj_dpot_master * funct_master(ui);
      }

      r_s[row_conc] -= test_slave(vi) * j;
      r_s[row_pot] -= numelectrons * test_slave(vi) * j;
    }
  }
  else if (k_ss.M() or k_sm.M() or r_s.Length())
    dserror("Must provide both slave-side matrices and slave-side vector or none of them!");


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

        k_ms(row_conc, col_conc) -= test_master(vi) * dj_dc_slave * funct_slave(ui);
        k_ms(row_conc, col_pot) -= test_master(vi) * dj_dpot_slave * funct_slave(ui);
        k_ms(row_pot, col_conc) -= numelectrons * test_master(vi) * dj_dc_slave * funct_slave(ui);
        k_ms(row_pot, col_pot) -= numelectrons * test_master(vi) * dj_dpot_slave * funct_slave(ui);
      }

      for (int ui = 0; ui < nen_master; ++ui)
      {
        const int col_conc = ui * 2;
        const int col_pot = col_conc + 1;

        k_mm(row_conc, col_conc) -= test_master(vi) * dj_dc_master * funct_master(ui);
        k_mm(row_conc, col_pot) -= test_master(vi) * dj_dpot_master * funct_master(ui);
        k_mm(row_pot, col_conc) -= numelectrons * test_master(vi) * dj_dc_master * funct_master(ui);
        k_mm(row_pot, col_pot) -=
            numelectrons * test_master(vi) * dj_dpot_master * funct_master(ui);
      }

      r_m[row_conc] += test_master(vi) * j;
      r_m[row_pot] += numelectrons * test_master(vi) * j;
    }
  }
  else if (k_ms.M() or k_mm.M() or r_m.Length())
    dserror("Must provide both master-side matrices and master-side vector or none of them!");
  return;
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
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, const double, const double, const double, const double, const double,
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
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, const double, const double, const double, const double, const double,
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
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, const double, const double, const double, const double, const double,
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
        const int, const int, const std::vector<int>*, const double, const double, const double,
        const double, const double, const double, const double, const double, const double,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
        Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&);
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs9>;
