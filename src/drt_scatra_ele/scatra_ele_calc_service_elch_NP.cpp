/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_NP.cpp

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
#include "scatra_ele_calc_elch_NP.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------------------------------------------*
 | validity check with respect to input parameters, degrees of freedom, number of scalars etc. fang
 02/15 |
 *----------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CheckElchElementParameter(
    DRT::Element* ele  //!< current element
)
{
  // check material
  if (ele->Material()->MaterialType() != INPAR::MAT::m_matlist) dserror("Invalid material type!");

  // check type of closing equation
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc:
    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    case INPAR::ELCH::equpot_poisson:
    case INPAR::ELCH::equpot_laplace:
    {
      // valid closing equations for electric potential
      break;
    }
    default:
    {
      dserror("Invalid closing equation for electric potential!");
      break;
    }
  }

  // check stabilization
  if (my::scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization and
      my::scatrapara_->StabType() != INPAR::SCATRA::stabtype_SUPG)
    dserror(
        "Only SUPG-type stabilization available for electrochemistry problems governed by "
        "Nernst-Planck formulation!");

  return;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode boundary kinetics point condition   fang 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::EvaluateElchBoundaryKineticsPoint(
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

    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (unsigned vi = 0; vi < my::nen_; ++vi)
        {
          for (unsigned ui = 0; ui < my::nen_; ++ui)
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

    // need special treatment for Laplace equation due to missing scaling with inverse of Faraday
    // constant
    case INPAR::ELCH::equpot_laplace:
    {
      const double faraday = myelch::elchparams_->Faraday();
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (unsigned vi = 0; vi < my::nen_; ++vi)
        {
          for (unsigned ui = 0; ui < my::nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                faraday * nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                faraday * nume *
                emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] +=
              faraday * nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    case INPAR::ELCH::equpot_poisson:
    {
      dserror("Poisson equation combined with electrode boundary conditions not implemented!");
      break;
    }

    default:
    {
      dserror("Unknown closing equation for electric potential!");
      break;
    }
  }  // switch(myelch::elchparams_->EquPot())

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::EvaluateElchBoundaryKineticsPoint


/*----------------------------------------------------------------------*
 * Get Conductivity
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetConductivity(
    const enum INPAR::ELCH::EquPot equpot, double& sigma_all, std::vector<double>& sigma,
    bool effCond  // the bool effCond is not used for the NP formulation since the volume averaging
                  // is not implemented
)
{
  // calculate conductivity of electrolyte solution
  const double frt = VarManager()->FRT();
  const double factor = frt * myelch::elchparams_->Faraday();  // = F^2/RT

  // Dilute solution theory:
  // Conductivity is computed by
  // sigma = F^2/RT*Sum(z_k^2 D_k c_k)
  for (int k = 0; k < my::numscal_; k++)
  {
    double sigma_k = factor * myelch::DiffManager()->GetValence(k) *
                     myelch::DiffManager()->GetIsotropicDiff(k) *
                     myelch::DiffManager()->GetValence(k) * VarManager()->Phinp(k);
    sigma[k] += sigma_k;  // insert value for this ionic species
    sigma_all += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if (equpot == INPAR::ELCH::equpot_enc_pde_elim)
    {
      sigma_all += factor * myelch::DiffManager()->GetIsotropicDiff(my::numscal_) *
                   myelch::DiffManager()->GetValence(my::numscal_) *
                   myelch::DiffManager()->GetValence(k) * (-VarManager()->Phinp(k));
    }
  }

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetConductivity


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_, 1>& q,          //!< flux of species k
    const INPAR::SCATRA::FluxType fluxtype,  //!< type fo flux
    const int k                              //!< index of current scalar
)
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

      // no break statement here!
    case INPAR::SCATRA::flux_diffusive:
      // diffusive flux contribution
      q.Update(-myelch::DiffManager()->GetIsotropicDiff(k), VarManager()->GradPhi(k), 1.0);

      q.Update(-VarManager()->FRT() * myelch::DiffManager()->GetIsotropicDiff(k) *
                   myelch::DiffManager()->GetValence(k) * VarManager()->Phinp(k),
          VarManager()->GradPot(), 1.0);

      break;

    default:
    {
      dserror("received illegal flag inside flux evaluation for whole domain");
      break;
    }
  };

  return;
}  // ScaTraCalc::CalculateFlux

/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& errors)
{
  // at the moment, there is only one analytical test problem available!
  if (DRT::INPUT::get<SCATRA::Action>(params, "action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // in the ALE case add nodal displacements
  if (my::scatrapara_->IsAle()) dserror("No ALE for Kwok & Wu error calculation allowed.");

  // set constants for analytical solution
  const double t =
      my::scatraparatimint_->Time() +
      (1 - my::scatraparatimint_->AlphaF()) * my::scatraparatimint_->Dt();  //-(1-alphaF_)*dta_
  const double frt = VarManager()->FRT();

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

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype =
      DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
  switch (errortype)
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
      //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

      if (my::numscal_ != 2) dserror("Numscal_ != 2 for desired error calculation.");

      // working arrays
      double potint(0.0);
      LINALG::Matrix<2, 1> conint(true);
      LINALG::Matrix<my::nsd_, 1> xint(true);
      LINALG::Matrix<2, 1> c(true);
      double deltapot(0.0);
      LINALG::Matrix<2, 1> deltacon(true);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get values of all transported scalars at integration point
        for (int k = 0; k < my::numscal_; ++k) conint(k) = my::funct_.Dot(my::ephinp_[k]);

        // get el. potential solution at integration point
        potint = my::funct_.Dot(my::ephinp_[my::numscal_]);

        // get global coordinate of integration point
        xint.Multiply(my::xyze_, my::funct_);

        // compute various constants
        const double d =
            frt *
            ((myelch::DiffManager()->GetIsotropicDiff(0) * myelch::DiffManager()->GetValence(0)) -
                (myelch::DiffManager()->GetIsotropicDiff(1) *
                    myelch::DiffManager()->GetValence(1)));
        if (abs(d) == 0.0) dserror("division by zero");
        const double D =
            frt *
            ((myelch::DiffManager()->GetValence(0) * myelch::DiffManager()->GetIsotropicDiff(0) *
                 myelch::DiffManager()->GetIsotropicDiff(1)) -
                (myelch::DiffManager()->GetValence(1) * myelch::DiffManager()->GetIsotropicDiff(1) *
                    myelch::DiffManager()->GetIsotropicDiff(0))) /
            d;

        // compute analytical solution for cation and anion concentrations
        const double A0 = 2.0;
        const double m = 1.0;
        const double n = 2.0;
        const double k = 3.0;
        const double A_mnk = 1.0;
        double expterm;
        double c_0_0_0_t;

        if (my::nsd_ == 3)
        {
          expterm = exp((-D) * (m * m + n * n + k * k) * t * PI * PI);
          c(0) = A0 +
                 (A_mnk * (cos(m * PI * xint(0)) * cos(n * PI * xint(1)) * cos(k * PI * xint(2))) *
                     expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m + n * n + k * k) * t * PI * PI));
        }
        else if (my::nsd_ == 2)
        {
          expterm = exp((-D) * (m * m + n * n) * t * PI * PI);
          c(0) = A0 + (A_mnk * (cos(m * PI * xint(0)) * cos(n * PI * xint(1))) * expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m + n * n) * t * PI * PI));
        }
        else if (my::nsd_ == 1)
        {
          expterm = exp((-D) * (m * m) * t * PI * PI);
          c(0) = A0 + (A_mnk * (cos(m * PI * xint(0))) * expterm);
          c_0_0_0_t = A0 + (A_mnk * exp((-D) * (m * m) * t * PI * PI));
        }
        else
          dserror("Illegal number of space dimensions for analyt. solution: %d", my::nsd_);

        // compute analytical solution for anion concentration
        c(1) =
            (-myelch::DiffManager()->GetValence(0) / myelch::DiffManager()->GetValence(1)) * c(0);
        // compute analytical solution for el. potential
        const double pot = ((myelch::DiffManager()->GetIsotropicDiff(1) -
                                myelch::DiffManager()->GetIsotropicDiff(0)) /
                               d) *
                           log(c(0) / c_0_0_0_t);

        // compute differences between analytical solution and numerical solution
        deltapot = potint - pot;
        deltacon.Update(1.0, conint, -1.0, c);

        // add square to L2 error
        errors[0] += deltacon(0) * deltacon(0) * fac;  // cation concentration
        errors[1] += deltacon(1) * deltacon(1) * fac;  // anion concentration
        errors[2] += deltapot * deltapot * fac;        // electric potential in electrolyte solution

      }  // end of loop over integration points
    }    // Kwok and Wu
    break;
    case INPAR::SCATRA::calcerror_cylinder:
    {
      // two-ion system with Butler-Volmer kinetics between two concentric cylinders
      //   G. Bauer, V. Gravemeier, W.A. Wall,
      //   A 3D finite element approach for the coupled numerical simulation of
      //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

      if (my::numscal_ != 2) dserror("Numscal_ != 2 for desired error calculation.");

      // working arrays
      LINALG::Matrix<2, 1> conint(true);
      LINALG::Matrix<my::nsd_, 1> xint(true);
      LINALG::Matrix<2, 1> c(true);
      LINALG::Matrix<2, 1> deltacon(true);

      // some constants that are needed
      const double c0_inner = 0.6147737641011396;
      const double c0_outer = 1.244249192148809;
      const double r_inner = 1.0;
      const double r_outer = 2.0;
      const double pot_inner = 2.758240847314454;
      const double b = log(r_outer / r_inner);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get values of all transported scalars at integration point
        for (int k = 0; k < my::numscal_; ++k) conint(k) = my::funct_.Dot(my::ephinp_[k]);

        // get el. potential solution at integration point
        const double potint = my::funct_.Dot(my::ephinp_[my::numscal_]);

        // get global coordinate of integration point
        xint.Multiply(my::xyze_, my::funct_);

        // evaluate analytical solution for cation concentration at radial position r
        if (my::nsd_ == 3)
        {
          const double r = sqrt(xint(0) * xint(0) + xint(1) * xint(1));
          c(0) = c0_inner + ((c0_outer - c0_inner) * (log(r) - log(r_inner)) / b);
        }
        else
          dserror("Illegal number of space dimensions for analyt. solution: %d", my::nsd_);

        // compute analytical solution for anion concentration
        c(1) =
            (-myelch::DiffManager()->GetValence(0) / myelch::DiffManager()->GetValence(1)) * c(0);
        // compute analytical solution for el. potential
        const double d =
            frt *
            ((myelch::DiffManager()->GetIsotropicDiff(0) * myelch::DiffManager()->GetValence(0)) -
                (myelch::DiffManager()->GetIsotropicDiff(1) *
                    myelch::DiffManager()->GetValence(1)));
        if (abs(d) == 0.0) dserror("division by zero");
        // reference value + ohmic resistance + concentration potential
        const double pot =
            pot_inner +
            log(c(0) / c0_inner);  // + (((diffus_[1]-diffus_[0])/d) * log(c(0)/c0_inner));

        // compute differences between analytical solution and numerical solution
        double deltapot = potint - pot;
        deltacon.Update(1.0, conint, -1.0, c);

        // add square to L2 error
        errors[0] += deltacon(0) * deltacon(0) * fac;  // cation concentration
        errors[1] += deltacon(1) * deltacon(1) * fac;  // anion concentration
        errors[2] += deltapot * deltapot * fac;        // electric potential in electrolyte solution

      }  // end of loop over integration points
    }    // concentric cylinders
    break;
    case INPAR::SCATRA::calcerror_electroneutrality:
    {
      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get values of transported scalars at integration point
        // and compute electroneutrality
        double deviation(0.0);
        for (int k = 0; k < my::numscal_; ++k)
        {
          const double conint_k = my::funct_.Dot(my::ephinp_[k]);
          deviation += myelch::DiffManager()->GetValence(k) * conint_k;
        }

        // add square to L2 error
        errors[0] += deviation * deviation * fac;
      }  // loop over integration points
    }
    break;
    default:
      dserror("Unknown analytical solution!");
      break;
  }  // switch(errortype)

  return;
}  // CalErrorComparedToAnalytSolution


/*------------------------------------------------------------------------------*
 | set internal variables for Nernst-Planck formulation              fang 02/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchNP(
      my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs27>;
