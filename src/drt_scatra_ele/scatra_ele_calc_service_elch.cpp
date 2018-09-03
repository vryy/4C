/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch.cpp

\brief evaluation of scatra elements for elch

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_utils_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"  // consistency check of formulation and material

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

#include "scatra_ele_calc_elch.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::check_scatra_element_parameter:
    {
      CheckElchElementParameter(ele);
      break;
    }

    case SCATRA::calc_flux_domain:
    {
      //--------------------------------------------------------------------------------
      // extract element based or nodal values
      //--------------------------------------------------------------------------------

      // get number of dofset associated with velocity related dofs
      const int ndsvel = params.get<int>("ndsvel");

      // get velocity values at nodes
      const Teuchos::RCP<const Epetra_Vector> convel =
          discretization.GetState(ndsvel, "convective velocity field");
      const Teuchos::RCP<const Epetra_Vector> vel =
          discretization.GetState(ndsvel, "velocity field");

      // safety check
      if (convel == Teuchos::null or vel == Teuchos::null) dserror("Cannot get state vector");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / my::nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel(my::nsd_ * my::nen_, -1);
      for (unsigned inode = 0; inode < my::nen_; ++inode)
        for (unsigned idim = 0; idim < my::nsd_; ++idim)
          lmvel[inode * my::nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // extract local values of (convective) velocity field from global state vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nsd_, my::nen_>>(
          *convel, my::econvelnp_, lmvel);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nsd_, my::nen_>>(*vel, my::evelnp_, lmvel);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      my::rotsymmpbc_->RotateMyValuesIfNecessary(my::econvelnp_);
      my::rotsymmpbc_->RotateMyValuesIfNecessary(my::evelnp_);

      // need current values of transported scalar
      // -> extract local values from global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, la[0].lm_);

      //----------------------------------------------------------------------
      // calculation of element volume both for tau at ele. cent. and int. pt.
      //----------------------------------------------------------------------
      my::EvalShapeFuncAndDerivsAtEleCenter();

      //----------------------------------------------------------------------
      // get material and stabilization parameters (evaluation at element center)
      //----------------------------------------------------------------------
      // density at t_(n)
      std::vector<double> densn(my::numscal_, 1.0);
      // density at t_(n+1) or t_(n+alpha_F)
      std::vector<double> densnp(my::numscal_, 1.0);
      // density at t_(n+alpha_M)
      std::vector<double> densam(my::numscal_, 1.0);

      // fluid viscosity
      double visc(0.0);

      // material parameter at the element center
      if (not my::scatrapara_->MatGP())
      {
        SetInternalVariablesForMatAndRHS();

        GetMaterialParams(ele, densn, densnp, densam, visc);
      }

      //----------------------------------------------------------------------
      // integration loop for one element
      //----------------------------------------------------------------------
      // integration points and weights
      const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // set internal variables
        SetInternalVariablesForMatAndRHS();

        //----------------------------------------------------------------------
        // get material parameters (evaluation at integration point)
        //----------------------------------------------------------------------
        if (my::scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

        // access control parameter for flux calculation
        INPAR::SCATRA::FluxType fluxtype = my::scatrapara_->CalcFluxDomain();
        Teuchos::RCP<std::vector<int>> writefluxids = my::scatrapara_->WriteFluxIds();

        // do a loop for systems of transported scalars
        for (std::vector<int>::iterator it = writefluxids->begin(); it != writefluxids->end(); ++it)
        {
          int k = 0;
          // Actually, we compute here a weighted (and integrated) form of the fluxes!
          // On time integration level, these contributions are then used to calculate
          // an L2-projected representation of fluxes.
          // Thus, this method here DOES NOT YET provide flux values that are ready to use!!

          // allocate and initialize!
          LINALG::Matrix<my::nsd_, 1> q(true);

          if ((*it) != my::numdofpernode_)
          {
            k = (*it) - 1;
            CalculateFlux(q, fluxtype, k);
          }
          else if ((*it) == my::numdofpernode_)
          {
            k = my::numdofpernode_ - 1;
            CalculateCurrent(q, fluxtype, fac);
          }
          else
            dserror("Flux id, which should be calculated, does not exit in the dof set.");

          // integrate and assemble everything into the "flux" vector
          for (unsigned vi = 0; vi < my::nen_; vi++)
          {
            const int fvi = vi * my::numdofpernode_ + k;
            elevec1_epetra[fvi] += fac * my::funct_(vi) * q(0);
            elevec2_epetra[fvi] += fac * my::funct_(vi) * q(1);
            if (my::nsd_ < 3)
              elevec3_epetra[fvi] = 0.0;
            else
              elevec3_epetra[fvi] += fac * my::funct_(vi) * q(2);
          }  // vi
        }
      }

      break;
    }
    case SCATRA::calc_error:
    {
      // check if length suffices
      if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

      // need current solution
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, la[0].lm_);

      CalErrorComparedToAnalytSolution(ele, params, elevec1_epetra);

      break;
    }

    case SCATRA::calc_elch_conductivity:
    {
      // get flag if effective conductivity should be calculated
      bool effCond = params.get<bool>("effCond");
      // get flag if the inverse of the conductivity should be calculated -> specific resistance
      bool specresist = params.get<bool>("specresist");

      // extract quantities for element evaluation
      this->ExtractElementAndNodeValues(ele, params, discretization, la);

      // elevec1_epetra(0):          conductivity of ionic species 0
      // elevec1_epetra(numscal_-1): conductivity of ionic species (numscal_-1)
      // elevec1_epetra(numscal_):   conductivity of the electrolyte solution (sum_k sigma(k))
      // elevec1_epetra(numscal_+1): domain integral
      CalculateConductivity(ele, elchparams_->EquPot(), elevec1_epetra, effCond, specresist);
      break;
    }

    case SCATRA::calc_elch_boundary_kinetics_point:
    {
      // process electrode boundary kinetics point condition
      CalcElchBoundaryKineticsPoint(
          ele, params, discretization, la[0].lm_, elemat1_epetra, elevec1_epetra, 1.);

      break;
    }

    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------------------------*
 | calculate error of numerical solution with respect to analytical solution   fang 10/16 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele,          //!< element
    Teuchos::ParameterList& params,   //!< parameter list
    Epetra_SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
)
{
  // call base class routine
  my::CalErrorComparedToAnalytSolution(ele, params, errors);

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalErrorComparedToAnalytSolution


/*----------------------------------------------------------------------*
  |  Calculate conductivity (ELCH) (private)                   gjb 07/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateConductivity(const DRT::Element* ele,
    const enum INPAR::ELCH::EquPot equpot, Epetra_SerialDenseVector& sigma_domint, bool effCond,
    bool specresist)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    //----------------------------------------------------------------------
    // get material and stabilization parameters (evaluation at element center)
    //----------------------------------------------------------------------
    // density at t_(n)
    std::vector<double> densn(my::numscal_, 1.0);
    // density at t_(n+1) or t_(n+alpha_F)
    std::vector<double> densnp(my::numscal_, 1.0);
    // density at t_(n+alpha_M)
    std::vector<double> densam(my::numscal_, 1.0);

    // fluid viscosity
    double visc(0.0);

    // set internal variables at integration point
    SetInternalVariablesForMatAndRHS();

    // material parameter at integration point
    GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // calculate integrals of (inverted) scalar(s) and domain
    for (unsigned i = 0; i < my::nen_; i++)
    {
      double sigma_all(0.0);
      std::vector<double> sigma(my::numscal_, 0.0);
      // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
      GetConductivity(equpot, sigma_all, sigma, effCond);

      const double fac_funct_i = fac * my::funct_(i);

      // sigma_domint(0):          conductivity of ionic species 0
      // sigma_domint(numscal_-1): conductivity of ionic species (numscal_-1)
      // sigma_domint(numscal_):   conductivity of the electrolyte solution (sum_k sigma(k))
      // sigma_domint(numscal_+1): domain integral
      for (int k = 0; k < my::numscal_; k++)
      {
        sigma_domint[k] += sigma[k] * fac_funct_i;
      }

      // calculation of conductivity or specific resistance of electrolyte solution
      if (specresist == false)
        sigma_domint[my::numscal_] += sigma_all * fac_funct_i;
      else
        sigma_domint[my::numscal_] += 1 / sigma_all * fac_funct_i;

      // domain integral
      sigma_domint[my::numscal_ + 1] += fac_funct_i;
    }
  }  // loop over integration points

  return;

}  // ScaTraEleCalcElch()


/*----------------------------------------------------------------------*
 | process an electrode boundary kinetics point condition    fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcElchBoundaryKineticsPoint(
    DRT::Element* ele,                         ///< current element
    Teuchos::ParameterList& params,            ///< parameter list
    DRT::Discretization& discretization,       ///< discretization
    std::vector<int>& lm,                      ///< location vector
    Epetra_SerialDenseMatrix& elemat1_epetra,  ///< element matrix
    Epetra_SerialDenseVector& elevec1_epetra,  ///< element right-hand side vector
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // get actual values of transported scalars
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from the global vector
  std::vector<LINALG::Matrix<my::nen_, 1>> ephinp(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, ephinp, lm);

  // get history variable (needed for double layer modeling)
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'hist'");

  // extract local values from the global vector
  std::vector<LINALG::Matrix<my::nen_, 1>> ehist(
      my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*hist, ehist, lm);

  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) dserror("Cannot access condition 'ElchBoundaryKineticsPoint'!");

  // access parameters of the condition
  const int kinetics = cond->GetInt("kinetic model");
  double pot0 = cond->GetDouble("pot");
  const int functnum = cond->GetInt("funct");
  const int nume = cond->GetInt("e-");
  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman
  // condition) but the electrode status is evaluated
  const int zerocur = cond->GetInt("zero_cur");
  if (nume < 0)
    dserror(
        "The convention for electrochemical reactions at the electrodes does not allow \n"
        "a negative number of transferred electrons");

  // convention for stoichiometric coefficients s_i:
  // Sum_i (s_i  M_i^(z_i)) -> n e- (n needs to be positive)
  const std::vector<int>* stoich = cond->GetMutable<std::vector<int>>("stoich");
  if ((unsigned int)my::numscal_ != (*stoich).size())
    dserror(
        "Electrode kinetics: number of stoichiometry coefficients %u does not match"
        " the number of ionic species %d",
        (*stoich).size(), my::numscal_);

  // the classical implementations of kinetic electrode models does not support
  // more than one reagent or product!! There are alternative formulations
  // as e.g. Newman (2004), pp. 205, eq. 8.6 with 8.10
  {
    int reactspecies = 0;
    for (int kk = 0; kk < my::numscal_; ++kk) reactspecies += abs((*stoich)[kk]);

    if (reactspecies > 1 and (kinetics == INPAR::ELCH::butler_volmer or
                                 kinetics == INPAR::ELCH::butler_volmer_yang1997 or
                                 kinetics == INPAR::ELCH::tafel or kinetics == INPAR::ELCH::linear))
      dserror(
          "Kinetic model Butler-Volmer / Butler-Volmer-Yang / Tafel and Linear: \n"
          "Only one educt and no product is allowed in the implemented version");
  }

  // access input parameter
  const double frt = elchparams_->FRT();
  if (frt <= 0.0) dserror("A negative factor frt is not possible by definition");

  // get control parameter from parameter list
  const bool is_stationary = my::scatraparatimint_->IsStationary();
  const double time = my::scatraparatimint_->Time();
  double timefac = 1.0;
  double rhsfac = 1.0;
  // find out whether we shell use a time curve and get the factor
  // this feature can be also used for stationary "pseudo time loops"
  if (functnum >= 0)
  {
    const double functfac = DRT::Problem::Instance()->Funct(functnum).EvaluateTime(time);

    // adjust potential at metal side accordingly
    pot0 *= functfac;
  }

  // TODO: Ist es besser über parameter?
  if (!(params.get<bool>("calc_status", false)))
  {
    if (not is_stationary)
    {
      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
      timefac = my::scatraparatimint_->TimeFac();
      if (timefac < 0.0) dserror("time factor is negative.");
      // for correct scaling of rhs contribution (see below)
      rhsfac = 1 / my::scatraparatimint_->AlphaF();
    }

    if (zerocur == 0)
    {
      EvaluateElchBoundaryKineticsPoint(ele, elemat1_epetra, elevec1_epetra, ephinp, ehist, timefac,
          cond, nume, *stoich, kinetics, pot0, frt, scalar);
    }

    // realize correct scaling of rhs contribution for gen.alpha case
    // with dt*(gamma/alpha_M) = timefac/alpha_F
    // matrix contributions are already scaled correctly with
    // timefac=dt*(gamma*alpha_F/alpha_M)
    elevec1_epetra.Scale(rhsfac);
  }
  else
  {
    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
    if (phidtnp == Teuchos::null) dserror("Cannot get state vector 'ephidtnp'");
    // extract local values from the global vector
    std::vector<LINALG::Matrix<my::nen_, 1>> ephidtnp(
        my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true));
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phidtnp, ephidtnp, lm);

    if (not is_stationary)
    {
      // One-step-Theta:    timefacrhs = theta*dt
      // BDF2:              timefacrhs = 2/3 * dt
      // generalized-alpha: timefacrhs = (gamma/alpha_M) * dt
      timefac = my::scatraparatimint_->TimeFacRhs();
      if (timefac < 0.) dserror("time factor is negative.");
    }

    EvaluateElectrodeStatusPoint(ele, elevec1_epetra, params, cond, ephinp, ephidtnp, kinetics,
        *stoich, nume, pot0, frt, timefac, scalar);
  }

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcElchBoundaryKineticsPoint


/*----------------------------------------------------------------------*
 | evaluate an electrode boundary kinetics point condition   fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateElchBoundaryKineticsPoint(
    const DRT::Element* ele,                                 ///< current element
    Epetra_SerialDenseMatrix& emat,                          ///< element matrix
    Epetra_SerialDenseVector& erhs,                          ///< element right-hand side vector
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ephinp,  ///< state variables at element nodes
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ehist,   ///< history variables at element nodes
    double timefac,                                          ///< time factor
    Teuchos::RCP<DRT::Condition> cond,  ///< electrode kinetics boundary condition
    const int nume,                     ///< number of transferred electrons
    const std::vector<int> stoich,      ///< stoichiometry of the reaction
    const int kinetics,                 ///< desired electrode kinetics model
    const double pot0,                  ///< electrode potential on metal side
    const double frt,                   ///< factor F/RT
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // get boundary porosity from condition if available, or set equal to volume porosity otherwise
  double epsilon = cond->GetDouble("epsilon");
  if (epsilon == -1)
    epsilon = scalar;
  else if (epsilon <= 0 or epsilon > 1)
    dserror("Boundary porosity has to be between 0 and 1, or -1 by default!");

  // extract nodal cloud of current condition
  const std::vector<int>* nodeids = cond->Nodes();

  // safety checks
  if (!nodeids) dserror("Electrode kinetics point boundary condition doesn't have nodal cloud!");
  if (nodeids->size() != 1)
    dserror(
        "Electrode kinetics point boundary condition must be associated with exactly one node!");
  if (my::nsd_ele_ != 1)
    dserror(
        "Electrode kinetics point boundary conditions are applicable to one-dimensional problems "
        "only!");

  // extract global ID of conditioned node
  const int nodeid = (*nodeids)[0];

  // find out whether conditioned node is the leftmost (position 0) or rightmost (position 1) node
  // of the current line element
  int position(-1);
  if (nodeid == ele->Nodes()[0]->Id())
    position = 0;
  else if (nodeid == ele->Nodes()[1]->Id())
    position = 1;
  else
    dserror(
        "Electrode kinetics point boundary condition must be imposed either on the leftmost or on "
        "the rightmost node of a line element!");

  // manipulate shape function array according to node position
  my::funct_.Clear();
  my::funct_(position) = 1.;

  // loop over all scalars
  for (int k = 0; k < my::numscal_; ++k)
  {
    if (stoich[k] == 0) continue;

    //(- N^(d+m)*n) = j = s_k / (nume * faraday * z_e-) * i
    //                  = s_k / (nume * faraday * (-1)) * i
    //                    |_______fns_________________|
    // see, e.g. in Ehrl et al., "A computational approach for the simulation of natural convection
    // in electrochemical cells", JCP, 2012
    double fns = -1.0 / elchparams_->Faraday() / nume;
    // stoichiometry as a consequence of the reaction convention
    fns *= stoich[k];

    // get valence of the single reactant
    const double valence_k = DiffManager()->GetValence(k);

    // call utility class for evaluation of electrode boundary kinetics point condition
    utils_->EvaluateElchKineticsAtIntegrationPoint(ele, emat, erhs, ephinp, ehist, timefac, 1.,
        my::funct_, cond, nume, stoich, valence_k, kinetics, pot0, frt, fns, epsilon, k);
  }  // loop over all scalars

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateElchBoundaryKineticsPoint


/*----------------------------------------------------------------------*
 | evaluate status information on point electrode            fang 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateElectrodeStatusPoint(
    const DRT::Element* ele,                                   ///< current element
    Epetra_SerialDenseVector& scalars,                         ///< scalars to be integrated
    Teuchos::ParameterList& params,                            ///< parameter list
    Teuchos::RCP<DRT::Condition> cond,                         ///< condition
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ephinp,    ///< state variables at element nodes
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ephidtnp,  ///< nodal time derivative vector
    const int kinetics,                                        ///< desired electrode kinetics model
    const std::vector<int> stoich,                             ///< stoichiometry of the reaction
    const int nume,                                            ///< number of transferred electrons
    const double pot0,     ///< electrode potential on metal side
    const double frt,      ///< factor F/RT
    const double timefac,  ///< time factor
    const double scalar    ///< scaling factor for current related quantities
)
{
  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at
  // the electrode. A different approach is not possible (without major hacks) since the
  // time-integration scheme is necessary to perform galvanostatic simulations, for instance. Think
  // about: double layer effects for genalpha time-integration scheme

  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neumann
  // condition) but the electrode status is evaluated
  const int zerocur = cond->GetInt("zero_cur");

  // get boundary porosity from condition if available, or set equal to volume porosity otherwise
  double epsilon = cond->GetDouble("epsilon");
  if (epsilon == -1)
    epsilon = scalar;
  else if (epsilon <= 0 or epsilon > 1)
    dserror("Boundary porosity has to be between 0 and 1, or -1 by default!");

  bool statistics = false;

  // extract nodal cloud of current condition
  const std::vector<int>* nodeids = cond->Nodes();

  // safety checks
  if (!nodeids) dserror("Electrode kinetics point boundary condition doesn't have nodal cloud!");
  if (nodeids->size() != 1)
    dserror(
        "Electrode kinetics point boundary condition must be associated with exactly one node!");
  if (my::nsd_ele_ != 1)
    dserror(
        "Electrode kinetics point boundary conditions are applicable to one-dimensional problems "
        "only!");

  // extract global ID of conditioned node
  const int nodeid = (*nodeids)[0];

  // find out whether conditioned node is the leftmost (position 0) or rightmost (position 1) node
  // of the current line element
  int position(-1);
  if (nodeid == ele->Nodes()[0]->Id())
    position = 0;
  else if (nodeid == ele->Nodes()[1]->Id())
    position = 1;
  else
    dserror(
        "Electrode kinetics point boundary condition must be imposed either on the leftmost or on "
        "the rightmost node of a line element!");

  // manipulate shape function array according to node position
  my::funct_.Clear();
  my::funct_(position) = 1.;

  // index of reactive species (starting from zero)
  for (int k = 0; k < my::numscal_; ++k)
  {
    // only the first oxidized species O is considered for statistics
    // statistics of other species result directly from the oxidized species (current density, ...)
    // or need to be implemented (surface concentration, OCV, ...)
    if (stoich[k] >= 0) continue;

    statistics = true;

    // call utility class for element evaluation
    utils_->EvaluateElectrodeStatusAtIntegrationPoint(ele, scalars, params, cond, ephinp, ephidtnp,
        my::funct_, zerocur, kinetics, stoich, nume, pot0, frt, timefac, 1., epsilon, k);

    // stop loop over ionic species after one evaluation (see also comment above)
    break;
  }  // loop over scalars

  // safety check
  if (statistics == false)
    dserror(
        "There is no oxidized species O (stoich<0) defined in your input file!! \n"
        " Statistics could not be evaluated");

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateElectrodeStatusPoint


/*----------------------------------------------------------------------------------------*
 | finite difference check on element level (for debugging only) (protected)   fang 10/14 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::FDCheck(DRT::Element* ele,
    Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs,
    Epetra_SerialDenseVector& subgrdiff)
{
  // screen output
  std::cout << "FINITE DIFFERENCE CHECK FOR ELEMENT " << ele->Id();

  // make a copy of state variables to undo perturbations later
  std::vector<LINALG::Matrix<my::nen_, 1>> ephinp_original(my::numdofpernode_);
  for (int k = 0; k < my::numdofpernode_; ++k)
    for (unsigned i = 0; i < my::nen_; ++i) ephinp_original[k](i, 0) = my::ephinp_[k](i, 0);

  // generalized-alpha time integration requires a copy of history variables as well
  std::vector<LINALG::Matrix<my::nen_, 1>> ehist_original(my::numscal_);
  if (my::scatraparatimint_->IsGenAlpha())
  {
    for (int k = 0; k < my::numscal_; ++k)
      for (unsigned i = 0; i < my::nen_; ++i) ehist_original[k](i, 0) = my::ehist_[k](i, 0);
  }

  // initialize element matrix and vectors for perturbed state
  Epetra_SerialDenseMatrix emat_dummy(emat);
  Epetra_SerialDenseVector erhs_perturbed(erhs);
  Epetra_SerialDenseVector subgrdiff_dummy(subgrdiff);

  // initialize counter for failed finite difference checks
  unsigned counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over columns of element matrix by first looping over nodes and then over dofs at each node
  for (unsigned inode = 0; inode < my::nen_; ++inode)
  {
    for (int idof = 0; idof < my::numdofpernode_; ++idof)
    {
      // number of current column of element matrix
      unsigned col = inode * my::numdofpernode_ + idof;

      // clear element matrix and vectors for perturbed state
      emat_dummy.Scale(0.0);
      erhs_perturbed.Scale(0.0);
      subgrdiff_dummy.Scale(0.0);

      // fill state vectors with original state variables
      for (int k = 0; k < my::numdofpernode_; ++k)
        for (unsigned i = 0; i < my::nen_; ++i) my::ephinp_[k](i, 0) = ephinp_original[k](i, 0);
      if (my::scatraparatimint_->IsGenAlpha())
        for (int k = 0; k < my::numscal_; ++k)
          for (unsigned i = 0; i < my::nen_; ++i) my::ehist_[k](i, 0) = ehist_original[k](i, 0);

      // impose perturbation
      if (my::scatraparatimint_->IsGenAlpha())
      {
        // perturbation of phi(n+alphaF), not of phi(n+1) => scale epsilon by factor alphaF
        my::ephinp_[idof](inode, 0) +=
            my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();

        // perturbation of phi(n+alphaF) by alphaF*epsilon corresponds to perturbation of phidtam
        // (stored in ehist_) by alphaM*epsilon/(gamma*dt); note: alphaF/timefac = alphaM/(gamma*dt)
        if (idof < my::numscal_)
          my::ehist_[idof](inode, 0) += my::scatraparatimint_->AlphaF() /
                                        my::scatraparatimint_->TimeFac() *
                                        my::scatrapara_->FDCheckEps();
      }
      else
        my::ephinp_[idof](inode, 0) += my::scatrapara_->FDCheckEps();

      // calculate element right-hand side vector for perturbed state
      Sysmat(ele, emat_dummy, erhs_perturbed, subgrdiff_dummy);

      // Now we compare the difference between the current entries in the element matrix
      // and their finite difference approximations according to
      // entries ?= (-erhs_perturbed + erhs_original) / epsilon

      // Note that the element right-hand side equals the negative element residual.
      // To account for errors due to numerical cancellation, we additionally consider
      // entries - erhs_original / epsilon ?= -erhs_perturbed / epsilon

      // Note that we still need to evaluate the first comparison as well. For small entries in the
      // element matrix, the second comparison might yield good agreement in spite of the entries
      // being wrong!
      for (unsigned row = 0; row < static_cast<unsigned>(my::numdofpernode_ * my::nen_); ++row)
      {
        // get current entry in original element matrix
        const double entry = emat(row, col);

        // finite difference suggestion (first divide by epsilon and then subtract for better
        // conditioning)
        const double fdval = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps() +
                             erhs(row) / my::scatrapara_->FDCheckEps();

        // confirm accuracy of first comparison
        if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in first comparison
        const double abserr1 = entry - fdval;
        if (abs(abserr1) > abs(maxabserr)) maxabserr = abserr1;
        double relerr1(0.);
        if (abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if (abs(fdval) > 1.e-17)
          relerr1 = abserr1 / abs(fdval);
        if (abs(relerr1) > abs(maxrelerr)) maxrelerr = relerr1;

        // evaluate first comparison
        if (abs(relerr1) > my::scatrapara_->FDCheckTol())
        {
          if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
          std::cout << "emat[" << row << "," << col << "]:  " << entry << "   ";
          std::cout << "finite difference suggestion:  " << fdval << "   ";
          std::cout << "absolute error:  " << abserr1 << "   ";
          std::cout << "relative error:  " << relerr1 << std::endl;

          counter++;
        }

        // first comparison OK
        else
        {
          // left-hand side in second comparison
          const double left = entry - erhs(row) / my::scatrapara_->FDCheckEps();

          // right-hand side in second comparison
          const double right = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps();

          // confirm accuracy of second comparison
          if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
            dserror("Finite difference check involves values too close to numerical zero!");

          // absolute and relative errors in second comparison
          const double abserr2 = left - right;
          if (abs(abserr2) > abs(maxabserr)) maxabserr = abserr2;
          double relerr2(0.);
          if (abs(left) > 1.e-17)
            relerr2 = abserr2 / abs(left);
          else if (abs(right) > 1.e-17)
            relerr2 = abserr2 / abs(right);
          if (abs(relerr2) > abs(maxrelerr)) maxrelerr = relerr2;

          // evaluate second comparison
          if (abs(relerr2) > my::scatrapara_->FDCheckTol())
          {
            if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
            std::cout << "emat[" << row << "," << col << "]-erhs[" << row << "]/eps:  " << left
                      << "   ";
            std::cout << "-erhs_perturbed[" << row << "]/eps:  " << right << "   ";
            std::cout << "absolute error:  " << abserr2 << "   ";
            std::cout << "relative error:  " << relerr2 << std::endl;

            counter++;
          }
        }
      }
    }
  }

  // screen output in case finite difference check is passed
  if (!counter)
    std::cout << " --> PASSED WITH MAXIMUM ABSOLUTE ERROR " << maxabserr
              << " AND MAXIMUM RELATIVE ERROR " << maxrelerr << std::endl;

  // undo perturbations of state variables
  for (int k = 0; k < my::numdofpernode_; ++k)
    for (unsigned i = 0; i < my::nen_; ++i) my::ephinp_[k](i, 0) = ephinp_original[k](i, 0);
  if (my::scatraparatimint_->IsGenAlpha())
    for (int k = 0; k < my::numscal_; ++k)
      for (unsigned i = 0; i < my::nen_; ++i) my::ehist_[k](i, 0) = ehist_original[k](i, 0);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
