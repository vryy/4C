/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for elch

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_scatra_ele_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch_diffcond.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CheckElchElementParameter(
    DRT::Element* ele)
{
  // 1) Check material specific options
  // 2) Check if numdofpernode, numscal is set correctly
  if (ele->Material()->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());

    int numphase = actmat->NumPhase();

    // access mat_elchmat: container material for porous structures in elch
    if (numphase != 1) FOUR_C_THROW("In the moment a single phase is only allowed.");

    // 1) loop over single phases
    for (int iphase = 0; iphase < actmat->NumPhase(); ++iphase)
    {
      // access phase material
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlephase = actmat->PhaseById(phaseid);

      // dynmic cast: get access to mat_phase
      const Teuchos::RCP<const MAT::ElchPhase>& actphase =
          Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(singlephase);

      // Check if numdofpernode, numscal is set correctly
      int nummat = actphase->NumMat();
      // enough materials defined
      if (nummat != my::numscal_)
      {
        FOUR_C_THROW(
            "The number of scalars defined in the material ElchMat does not correspond with "
            "the number of materials defined in the material MatPhase.");
      }

      int numdofpernode = 0;
      if (diffcondparams_->CurSolVar())
        numdofpernode = nummat + GLOBAL::Problem::Instance()->NDim() + numphase;
      else if (actphase->MatById(actphase->MatID(0))->MaterialType() ==
               INPAR::MAT::m_newman_multiscale)
        numdofpernode = 3;
      else
        numdofpernode = nummat + numphase;

      if (numdofpernode != my::numdofpernode_)
      {
        FOUR_C_THROW(
            "The chosen element formulation (e.g. current as solution variable) "
            "does not correspond with the number of dof's defined in your material");
      }

      // 2) loop over materials of the single phase
      for (int imat = 0; imat < actphase->NumMat(); ++imat)
      {
        const int matid = actphase->MatID(imat);
        Teuchos::RCP<const MAT::Material> singlemat = actphase->MatById(matid);

        if (singlemat->MaterialType() == INPAR::MAT::m_newman)
        {
          // Newman material must be combined with divi closing equation for electric potential
          if (myelch::elchparams_->EquPot() != INPAR::ELCH::equpot_divi)
          {
            FOUR_C_THROW(
                "Newman material must be combined with divi closing equation for electric "
                "potential!");
          }

          // Material Newman is derived for a binary electrolyte utilizing the ENC to condense the
          // non-reacting species
          if (my::numscal_ > 1)
          {
            FOUR_C_THROW(
                "Material Newman is only valid for one scalar (binary electrolyte utilizing the "
                "ENC)");
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CalcInitialTimeDerivative(
    DRT::Element* ele, CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // call base class routine
  myelch::CalcInitialTimeDerivative(ele, emat, erhs, params, discretization, la);

  // In the moment the diffusion manager contains the porosity at the last Gauss point (previous
  // call my::CalcInitialTimeDerivative()) Since the whole approach is valid only for constant
  // porosities, we do not fill the diffusion manager again at the element center The solution
  // variable is the initial time derivative. Therefore, we have to correct emat by the initial
  // porosity Attention: this procedure is only valid for a constant porosity in the beginning
  emat.scale(DiffManager()->GetPhasePoro(0));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CorrectRHSFromCalcRHSLinMass(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double fac, const double densnp,
    const double phinp)
{
  // fac->-fac to change sign of rhs
  if (my::scatraparatimint_->IsIncremental())
    my::CalcRHSLinMass(erhs, k, 0.0, -fac, 0.0, DiffManager()->GetPhasePoro(0));
  else
    FOUR_C_THROW("Must be incremental!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::Action::calc_elch_boundary_kinetics_point:
    {
      // access material of parent element
      Teuchos::RCP<MAT::Material> material = ele->Material();

      // extract porosity from material and store in diffusion manager
      if (material->MaterialType() == INPAR::MAT::m_elchmat)
      {
        const auto* elchmat = static_cast<const MAT::ElchMat*>(material.get());

        for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
        {
          Teuchos::RCP<const MAT::Material> phase = elchmat->PhaseById(elchmat->PhaseID(iphase));

          if (phase->MaterialType() == INPAR::MAT::m_elchphase)
          {
            DiffManager()->SetPhasePoro(
                (static_cast<const MAT::ElchPhase*>(phase.get()))->Epsilon(), iphase);
          }
          else
            FOUR_C_THROW("Invalid material!");
        }
      }

      else
        FOUR_C_THROW("Invalid material!");

      // process electrode boundary kinetics point condition
      myelch::CalcElchBoundaryKineticsPoint(ele, params, discretization, la[0].lm_, elemat1_epetra,
          elevec1_epetra, DiffManager()->GetPhasePoro(0));

      break;
    }

    case SCATRA::Action::calc_elch_domain_kinetics:
    {
      CalcElchDomainKinetics(
          ele, params, discretization, la[0].lm_, elemat1_epetra, elevec1_epetra);

      break;
    }

    default:
    {
      myelectrode::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CalcElchDomainKinetics(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra)
{
  // from scatra_ele_boundary_calc_elch_diffcond.cpp
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const auto* elchmat = static_cast<const MAT::ElchMat*>(material.get());

    for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
    {
      Teuchos::RCP<const MAT::Material> phase = elchmat->PhaseById(elchmat->PhaseID(iphase));

      if (phase->MaterialType() == INPAR::MAT::m_elchphase)
      {
        DiffManager()->SetPhasePoro(
            (static_cast<const MAT::ElchPhase*>(phase.get()))->Epsilon(), iphase);
      }
      else
        FOUR_C_THROW("Invalid material!");
    }
  }

  else
    FOUR_C_THROW("Invalid material!");

  // get actual values of transported scalars
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phinp'");

  // get history variable (needed for double layer modeling)
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  if (phinp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'hist'");

  // state and history variables at element nodes
  std::vector<CORE::LINALG::Matrix<nen_, 1>> ephinp(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  std::vector<CORE::LINALG::Matrix<nen_, 1>> ehist(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);
  CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*hist, ehist, lm);

  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) FOUR_C_THROW("Cannot access condition 'ElchDomainKinetics'");

  // access parameters of the condition
  const int kinetics = *cond->Get<int>("kinetic model");
  double pot0 = *cond->Get<double>("pot");
  const int curvenum = *cond->Get<int>("funct");
  const int nume = *cond->Get<int>("e-");
  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman
  // condition) but the electrode status is evaluated
  const int zerocur = *cond->Get<int>("zero_cur");
  if (nume < 0)
  {
    FOUR_C_THROW(
        "The convention for electrochemical reactions at the electrodes does not allow \n"
        "a negative number of transferred electrons");
  }

  // convention for stoichiometric coefficients s_i:
  // Sum_i (s_i  M_i^(z_i)) -> n e- (n needs to be positive)
  const auto* stoich = cond->Get<std::vector<int>>("stoich");
  if ((unsigned int)my::numscal_ != (*stoich).size())
  {
    FOUR_C_THROW(
        "Electrode kinetics: number of stoichiometry coefficients %u does not match"
        " the number of ionic species %d",
        (*stoich).size(), my::numscal_);
  }

  // the classical implementations of kinetic electrode models does not support
  // more than one reagent or product!! There are alternative formulations
  // as e.g. Newman (2004), pp. 205, eq. 8.6 with 8.10
  {
    int reactspecies = 0;
    for (int kk = 0; kk < my::numscal_; ++kk) reactspecies += abs((*stoich)[kk]);

    if (reactspecies > 1 and (kinetics == INPAR::ELCH::butler_volmer or
                                 kinetics == INPAR::ELCH::butler_volmer_yang1997 or
                                 kinetics == INPAR::ELCH::tafel or kinetics == INPAR::ELCH::linear))
    {
      FOUR_C_THROW(
          "Kinetic model Butler-Volmer / Butler-Volmer-Yang / Tafel and Linear: \n"
          "Only one educt and no product is allowed in the implemented version");
    }
  }

  // get control parameter from parameter list
  const bool is_stationary = my::scatraparatimint_->IsStationary();
  const double time = my::scatraparatimint_->Time();
  double timefac = 1.0;
  double rhsfac = 1.0;
  // find out whether we shell use a time curve and get the factor
  // this feature can be also used for stationary "pseudo time loops"
  if (curvenum >= 0)
  {
    const double curvefac =
        GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfTime>(curvenum).Evaluate(
            time);
    // adjust potential at metal side accordingly
    pot0 *= curvefac;
  }

  if (!(params.get<bool>("calc_status", false)))
  {
    if (not is_stationary)
    {
      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
      timefac = my::scatraparatimint_->TimeFac();
      if (timefac < 0.0) FOUR_C_THROW("time factor is negative.");
      // for correct scaling of rhs contribution (see below)
      rhsfac = 1 / my::scatraparatimint_->AlphaF();
    }

    if (zerocur == 0)
    {
      EvaluateElchDomainKinetics(ele, elemat1_epetra, elevec1_epetra, ephinp, ehist, timefac, cond,
          nume, *stoich, kinetics, pot0);
    }

    // realize correct scaling of rhs contribution for gen.alpha case
    // with dt*(gamma/alpha_M) = timefac/alpha_F
    // matrix contributions are already scaled correctly with
    // timefac=dt*(gamma*alpha_F/alpha_M)
    elevec1_epetra.scale(rhsfac);
  }

  else
  {
    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
    if (phidtnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'ephidtnp'");
    std::vector<CORE::LINALG::Matrix<nen_, 1>> ephidtnp(
        my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
    CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phidtnp, ephidtnp, lm);

    if (not is_stationary)
    {
      // One-step-Theta:    timefacrhs = theta*dt
      // BDF2:              timefacrhs = 2/3 * dt
      // generalized-alpha: timefacrhs = (gamma/alpha_M) * dt
      timefac = my::scatraparatimint_->TimeFacRhs();
      if (timefac < 0.) FOUR_C_THROW("time factor is negative.");
    }

    EvaluateElectrodeStatus(ele, elevec1_epetra, params, cond, ephinp, ephidtnp, kinetics, *stoich,
        nume, pot0, timefac);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::EvaluateElchBoundaryKineticsPoint(
    const DRT::Element* ele, CORE::LINALG::SerialDenseMatrix& emat,
    CORE::LINALG::SerialDenseVector& erhs, const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist, double timefac,
    Teuchos::RCP<DRT::Condition> cond, const int nume, const std::vector<int> stoich,
    const int kinetics, const double pot0, const double frt, const double scalar)
{
  // call base class routine
  myelch::EvaluateElchBoundaryKineticsPoint(
      ele, emat, erhs, ephinp, ehist, timefac, cond, nume, stoich, kinetics, pot0, frt, scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case INPAR::ELCH::equpot_divi:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] += nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown closing equation for electric potential!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::EvaluateElchDomainKinetics(
    const DRT::Element* ele, CORE::LINALG::SerialDenseMatrix& emat,
    CORE::LINALG::SerialDenseVector& erhs, const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist, double timefac,
    Teuchos::RCP<DRT::Condition> cond, const int nume, const std::vector<int> stoich,
    const int kinetics, const double pot0)
{
  // for pre-multiplication of i0 with 1/(F z_k)
  const double faraday = myelch::elchparams_->Faraday();

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < my::numscal_; ++k)
  {
    if (stoich[k] == 0) continue;

    //(- N^(d+m)*n) = j = s_k / (nume * faraday * z_e-) * i
    //                  = s_k / (nume * faraday * (-1)) * i
    //                    |_______fns_________________|
    // see, e.g. in Ehrl et al., "A computational approach for the simulation of natural convection
    // in electrochemical cells", JCP, 2012
    double fns = -1.0 / faraday / nume;
    // stoichiometry as a consequence of the reaction convention
    fns *= stoich[k];

    // get valence of the single reactant
    const double valence_k = myelch::DiffManager()->GetValence(k);

    /*----------------------------------------------------------------------*
     |               start loop over integration points                     |
     *----------------------------------------------------------------------*/
    for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
    {
      SetInternalVariablesForMatAndRHS();

      // access input parameter
      const double frt = VarManager()->FRT();
      if (frt <= 0.0) FOUR_C_THROW("A negative factor frt is not possible by definition");

      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

      // extract specific electrode surface area A_s from condition
      double A_s = *cond->Get<double>("A_s");

      // call utility class for element evaluation
      Utils()->EvaluateElchKineticsAtIntegrationPoint(ele, emat, erhs, ephinp, ehist, timefac, fac,
          my::funct_, cond, nume, stoich, valence_k, kinetics, pot0, frt, fns, A_s, k);
    }  // end of loop over integration points gpid
  }    // end loop over scalars

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case INPAR::ELCH::equpot_divi:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] += nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown closing equation for electric potential!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::EvaluateElectrodeStatus(
    const DRT::Element* ele, CORE::LINALG::SerialDenseVector& scalars,
    Teuchos::ParameterList& params, Teuchos::RCP<DRT::Condition> cond,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephidtnp, const int kinetics,
    const std::vector<int> stoich, const int nume, const double pot0, const double timefac)
{
  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at
  // the electrode. A different approach is not possible (without major hacks) since the
  // time-integration scheme is necessary to perform galvanostatic simulations, for instance. Think
  // about: double layer effects for genalpha time-integration scheme

  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman
  // condition) but the electrode status is evaluated
  const int zerocur = *cond->Get<int>("zero_cur");

  // extract volumetric electrode surface area A_s from condition
  double A_s = *cond->Get<double>("A_s");

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  bool statistics = false;

  // index of reactive species (starting from zero)
  for (int k = 0; k < my::numscal_; ++k)
  {
    // only the first oxidized species O is considered for statistics
    // statistics of other species result directly from the oxidized species (current density, ...)
    // or need to be implemented (surface concentration, OCV, ...)
    if (stoich[k] >= 0) continue;

    statistics = true;

    // loop over integration points
    for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

      SetInternalVariablesForMatAndRHS();

      // access input parameter
      const double frt = VarManager()->FRT();
      if (frt <= 0.0) FOUR_C_THROW("A negative factor frt is not possible by definition");

      // call utility class for element evaluation
      Utils()->EvaluateElectrodeStatusAtIntegrationPoint(ele, scalars, params, cond, ephinp,
          ephidtnp, my::funct_, zerocur, kinetics, stoich, nume, pot0, frt, timefac, fac, A_s, k);
    }  // loop over integration points

    // stop loop over ionic species after one evaluation (see also comment above)
    break;
  }  // loop over scalars

  // safety check
  if (!statistics)
  {
    FOUR_C_THROW(
        "There is no oxidized species O (stoich<0) defined in your input file!! \n"
        " Statistics could not be evaluated");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CalculateFlux(
    CORE::LINALG::Matrix<nsd_, 1>& q, const INPAR::SCATRA::FluxType fluxtype, const int k)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
    case INPAR::SCATRA::flux_total:
      // convective flux contribution
      q.Update(VarManager()->Phinp(k), VarManager()->ConVel(k));

      [[fallthrough]];
    case INPAR::SCATRA::flux_diffusive:
      // diffusive flux contribution
      q.Update(-DiffManager()->GetIsotropicDiff(k) * DiffManager()->GetPhasePoroTort(0),
          VarManager()->GradPhi(k), 1.0);
      // flux due to ohmic overpotential
      q.Update(-DiffManager()->GetTransNum(k) * DiffManager()->InvFVal(k) *
                   DiffManager()->GetCond() * DiffManager()->GetPhasePoroTort(0),
          VarManager()->GradPot(), 1.0);
      // flux due to concentration overpotential
      q.Update(-DiffManager()->GetTransNum(k) * VarManager()->RTFFC() /
                   DiffManager()->GetValence(k) * DiffManager()->GetCond() *
                   DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetThermFac() *
                   (diffcondparams_->NewmanConstA() +
                       (diffcondparams_->NewmanConstB() * DiffManager()->GetTransNum(k))) *
                   VarManager()->ConIntInv(k),
          VarManager()->GradPhi(k), 1.0);
      break;
    default:
      FOUR_C_THROW("received illegal flag inside flux evaluation for whole domain");
      break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CalculateCurrent(
    CORE::LINALG::Matrix<nsd_, 1>& q, const INPAR::SCATRA::FluxType fluxtype, const double fac)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
    case INPAR::SCATRA::flux_total:
    case INPAR::SCATRA::flux_diffusive:
      // ohmic flux contribution
      q.Update(-DiffManager()->GetCond(), VarManager()->GradPot());
      // diffusion overpotential flux contribution
      for (int k = 0; k < my::numscal_; ++k)
      {
        q.Update(-VarManager()->RTF() / diffcondparams_->NewmanConstC() * DiffManager()->GetCond() *
                     DiffManager()->GetThermFac() *
                     (diffcondparams_->NewmanConstA() +
                         (diffcondparams_->NewmanConstB() * DiffManager()->GetTransNum(k))) *
                     VarManager()->ConIntInv(k),
            VarManager()->GradPhi(k), 1.0);
      }

      break;
    default:
      FOUR_C_THROW("received illegal flag inside flux evaluation for whole domain");
      break;
  }
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& errors)
{
  switch (CORE::UTILS::GetAsEnum<INPAR::SCATRA::CalcError>(params, "calcerrorflag"))
  {
    case INPAR::SCATRA::calcerror_Kwok_Wu:
    {
      //   References:
      //   Kwok, Yue-Kuen and Wu, Charles C. K.
      //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
      //   Numerical Methods for Partial Differential Equations
      //   1995, Vol 11, 389-397

      //   G. Bauer, V. Gravemeier, W.A. Wall,
      //   A 3D finite element approach for the coupled numerical simulation of
      //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339-1359.

      // safety checks
      if (Teuchos::getIntegralValue<SCATRA::Action>(params, "action") != SCATRA::Action::calc_error)
        FOUR_C_THROW("How did you get here?");
      if (my::scatrapara_->IsAle()) FOUR_C_THROW("No ALE for Kwok & Wu error calculation allowed.");
      if (my::numscal_ != 1) FOUR_C_THROW("Numscal_ != 1 for desired error calculation.");

      // set constants for analytical solution
      const double t =
          my::scatraparatimint_->Time() +
          (1 - my::scatraparatimint_->AlphaF()) * my::scatraparatimint_->Dt();  //-(1-alphaF_)*dta_
      const double frt = VarManager()->FRT();

      // integration points and weights
      // more GP than usual due to (possible) cos/exp fcts in analytical solutions
      const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

      // working arrays
      double potint(0.0);
      CORE::LINALG::Matrix<1, 1> conint(true);
      CORE::LINALG::Matrix<nsd_, 1> xint(true);
      CORE::LINALG::Matrix<1, 1> c(true);
      double deltapot(0.0);
      CORE::LINALG::Matrix<1, 1> deltacon(true);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // density at t_(n)
        std::vector<double> densn(my::numscal_, 1.0);
        // density at t_(n+1) or t_(n+alpha_F)
        std::vector<double> densnp(my::numscal_, 1.0);
        // density at t_(n+alpha_M)
        std::vector<double> densam(my::numscal_, 1.0);

        // fluid viscosity
        double visc(0.0);

        // get material parameter (constants values)
        SetInternalVariablesForMatAndRHS();
        GetMaterialParams(ele, densn, densnp, densam, visc);

        // get values of all transported scalars at integration point
        conint(0) = my::funct_.Dot(my::ephinp_[0]);

        // get el. potential solution at integration point
        potint = my::funct_.Dot(my::ephinp_[my::numscal_]);

        // get global coordinate of integration point
        xint.Multiply(my::xyze_, my::funct_);

        const double D = DiffManager()->GetIsotropicDiff(0);

        // compute analytical solution for cation and anion concentrations
        const double A0 = 2.0;
        const double m = 1.0;
        const double n = 2.0;
        const double k = 3.0;
        const double A_mnk = 1.0;
        double expterm;
        double c_0_0_0_t;

        if (nsd_ == 3)
        {
          expterm = exp((-D) * (m * m + n * n + k * k) * t * M_PI * M_PI);
          c(0) = A0 +
                 (A_mnk *
                     (cos(m * M_PI * xint(0)) * cos(n * M_PI * xint(1)) * cos(k * M_PI * xint(2))) *
                     expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m + n * n + k * k) * t * M_PI * M_PI));
        }
        else if (nsd_ == 2)
        {
          expterm = exp((-D) * (m * m + n * n) * t * M_PI * M_PI);
          c(0) = A0 + (A_mnk * (cos(m * M_PI * xint(0)) * cos(n * M_PI * xint(1))) * expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m + n * n) * t * M_PI * M_PI));
        }
        else if (nsd_ == 1)
        {
          expterm = exp((-D) * (m * m) * t * M_PI * M_PI);
          c(0) = A0 + (A_mnk * (cos(m * M_PI * xint(0))) * expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m) * t * M_PI * M_PI));
        }
        else
          FOUR_C_THROW("Illegal number of space dimensions for analyt. solution: %d", nsd_);

        // compute analytical solution for el. potential
        // const double pot =
        // ((DiffManager()->GetIsotropicDiff(1)-DiffManager()->GetIsotropicDiff(0))/d) *
        // log(c(0)/c_0_0_0_t);
        const double pot = -1 / frt *
                           (diffcondparams_->NewmanConstA() +
                               (diffcondparams_->NewmanConstB() * DiffManager()->GetTransNum(0))) /
                           diffcondparams_->NewmanConstC() * log(c(0) / c_0_0_0_t);

        // compute differences between analytical solution and numerical solution
        deltapot = potint - pot;
        deltacon.Update(1.0, conint, -1.0, c);

        // add square to L2 error
        errors[0] += deltacon(0) * deltacon(0) * fac;  // cation concentration
        errors[1] += 0;                                // anion concentration
        errors[2] += deltapot * deltapot * fac;        // electric potential in electrolyte solution

      }  // end of loop over integration points
    }    // Kwok and Wu
    break;
    default:
    {
      myelectrode::CalErrorComparedToAnalytSolution(ele, params, errors);
      break;
    }
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchDiffCond(
      my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);
}

template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype,
    probdim>::CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
    const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& conc)
{
  // for complete 1D simulation of battery:
  // Micro state must exist for electrolyte -> set value to 0.0
  for (int node = 0; node < static_cast<int>(nen_); ++node) conc(node) = 0.0;
}
// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::hex8, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tet10, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::pyramid5, 3>;

FOUR_C_NAMESPACE_CLOSE
