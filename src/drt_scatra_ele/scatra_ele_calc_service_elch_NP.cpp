/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_NP.cpp

\brief evaluation of scatra elements for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "../drt_mat/material.H"

#include "scatra_ele.H"
#include "scatra_ele_calc_elch_NP.H"


/*----------------------------------------------------------------------------------------------------------*
 | validity check with respect to input parameters, degrees of freedom, number of scalars etc.   fang 02/15 |
 *----------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CheckElchElementParameter(
    DRT::Element*              ele   //!< current element
)
{
  // safety checks
  if(ele->Material()->MaterialType() != INPAR::MAT::m_matlist)
    dserror("Invalid material type!");

  switch(myelch::ElchPara()->EquPot())
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

  return;
}


/*----------------------------------------------------------------------*
 * Get Conductivity
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetConductivity(
  const enum INPAR::ELCH::EquPot    equpot,
  double&                           sigma_all,
  Epetra_SerialDenseVector&         sigma
)
{
  // calculate conductivity of electrolyte solution
  const double frt = myelch::ElchPara()->FRT();
  const double factor = frt*INPAR::ELCH::faraday_const; // = F^2/RT

  // get concentration of transported scalar k at integration point
  std::vector<double> conint(my::numscal_);
  for (int k = 0;k<my::numscal_;++k)
    conint[k] = my::funct_.Dot(my::ephinp_[k]);

  // Dilute solution theory:
  // Conductivity is computed by
  // sigma = F^2/RT*Sum(z_k^2 D_k c_k)
  for(int k=0; k < my::numscal_; k++)
  {
    double sigma_k = factor*myelch::DiffManager()->GetValence(k)*myelch::DiffManager()->GetIsotropicDiff(k)*myelch::DiffManager()->GetValence(k)*conint[k];
    sigma[k] += sigma_k; // insert value for this ionic species
    sigma_all += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if(equpot==INPAR::ELCH::equpot_enc_pde_elim)
    {
      sigma_all += factor*myelch::DiffManager()->GetIsotropicDiff(my::numscal_)*myelch::DiffManager()->GetValence(my::numscal_)*myelch::DiffManager()->GetValence(k)*(-conint[k]);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac         //!< integration factor
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
  case INPAR::SCATRA::flux_total_domain:
    // convective flux contribution
    q.Update(VarManager()->ConInt(k),VarManager()->ConVelInt());

    // no break statement here!
  case INPAR::SCATRA::flux_diffusive_domain:
    // diffusive flux contribution
    q.Update(-myelch::DiffManager()->GetIsotropicDiff(k),VarManager()->GradPhi(k),1.0);

    q.Update(-myelch::ElchPara()->FRT()*myelch::DiffManager()->GetIsotropicDiff(k)*myelch::DiffManager()->GetValence(k)*VarManager()->ConInt(k),VarManager()->GradPot(),1.0);

    break;

  default:
  {
    dserror("received illegal flag inside flux evaluation for whole domain");
    break;
  }
  };

  return;
} // ScaTraCalc::CalculateFlux


/*------------------------------------------------------------------------------*
 | set internal variables for Nernst-Planck formulation              fang 02/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchNP(my::funct_,my::derxy_,my::ephinp_,my::ephin_,myelch::epotnp_,my::econvelnp_,my::ehist_);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs27>;
