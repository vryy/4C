/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch.cpp

\brief evaluation of scatra elements for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch_NP.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
//TODO: SCATRA_ELE_CLEANING: Wie bekommen wir das sonst?
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_globalproblem.H"  // consistency check of formulation and material

#include "../drt_inpar/inpar_elch.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

//#include "scatra_ele_parameter_timint.H"
//
//#include "../drt_lib/drt_utils.H"

//#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 * Add dummy mass matrix to sysmat
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::PrepMatAndRhsInitialTimeDerivative(
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra
)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  /*----------------------------------------------------------------------*/
  // element integration loop
  /*----------------------------------------------------------------------*/
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // loop starts at k=numscal_ !!
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac*my::funct_(vi); // no density required here
      const int fvi = vi*my::numdofpernode_+my::numscal_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+my::numscal_;

        elemat1_epetra(fvi,fui) += v*my::funct_(ui);
      }
    }
  }

  // set zero for the rhs of the potential
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+my::numscal_;

    elevec1_epetra[fvi] = 0.0; // zero out!
  }

  return;
}

/*----------------------------------------------------------------------*
 * Get Conductivity
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetConductivity(
  const enum INPAR::ELCH::ElchType  elchtype,
  double&                           sigma_all,
  Epetra_SerialDenseVector&         sigma
)
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();
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
    double sigma_k = factor*dme->GetValence(k)*dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k];
    sigma[k] += sigma_k; // insert value for this ionic species
    sigma_all += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
    {
      sigma_all += factor*dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*dme->GetValence(k)*(-conint[k]);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * Calculate Mat and Rhs for electric potential field        ehrl 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalMatAndRhsElectricPotentialField(
  Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
  const enum INPAR::ELCH::ElchType  elchtype,
  Epetra_SerialDenseMatrix&         emat,
  Epetra_SerialDenseVector&         erhs,
  const double                      fac,
  Teuchos::RCP<ScaTraEleDiffManagerElch>& dme
)
{
  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();

  double sigmaint(0.0);
  for (int k=0; k<my::numscal_; ++k)
  {
    double sigma_k = frt*dme->GetValence(k)*dme->GetIsotropicDiff(k)*dme->GetValence(k)*vm->ConInt(k);
    sigmaint += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
      sigmaint += frt*dme->GetValence(k)*dme->GetIsotropicDiff(k)*dme->GetValence(my::numscal_)*(-vm->ConInt(k));

    // diffusive terms on rhs
    const double vrhs = fac*dme->GetIsotropicDiff(k)*dme->GetValence(k);
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+my::numscal_;
      double laplawf(0.0);
      my::GetLaplacianWeakFormRHS(laplawf,vm->GradPhi(k),vi);
      erhs[fvi] -= vrhs*laplawf;
      // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
      if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
        erhs[fvi] -= -fac*dme->GetValence(k)*dme->GetIsotropicDiff(my::numscal_)*laplawf;
    }

    // provide something for conc. dofs: a standard mass matrix
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int    fvi = vi*my::numdofpernode_+k;
      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;
        emat(fvi,fui) += fac*my::funct_(vi)*my::funct_(ui);
      }
    }
  } // for k

  // ----------------------------------------matrix entries
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int    fvi = vi*my::numdofpernode_+my::numscal_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+my::numscal_;
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi);
      emat(fvi,fui) += fac*sigmaint*laplawf;
    }

    double laplawf(0.0);
    my::GetLaplacianWeakFormRHS(laplawf,vm->GradPot(),vi);
    erhs[fvi] -= fac*sigmaint*laplawf;
  }

  return;
}


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs27>;
