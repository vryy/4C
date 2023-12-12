/*----------------------------------------------------------------------*/
/*! \file

\brief quadrature rules for wall-enrichment element

\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele.H"
#include "baci_fluid_ele_calc_xwall.H"
#include "baci_fluid_ele_parameter_std.H"
#include "baci_fluid_ele_parameter_timint.H"
#include "baci_lib_condition_utils.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_newtonianfluid.H"

BACI_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------*
 | Prepare custom (direction-dependent) Gauss rule                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::PrepareGaussRule()
{
  // which is the wall-normal element direction?
  // calculate jacobian at element center
  my::is_higher_order_ele_ = false;
  my::EvalShapeFuncAndDerivsAtEleCenter();
  my::is_higher_order_ele_ = true;

  // the derivative of the wall distance with respect to the local coordinates
  // shows how the local axes are oriented with respect to the wall-normal vector
  CORE::LINALG::Matrix<nsd_, 1> normwallrst(true);
  normwallrst.Multiply(deriv_, ewdist_);
  double normwallrstnorm2 = normwallrst.Norm2();
  const double dot1 = abs(normwallrst(0) / normwallrstnorm2);
  const double dot2 = abs(normwallrst(1) / normwallrstnorm2);
  const double dot3 = abs(normwallrst(2) / normwallrstnorm2);

  double minyp = 1e9;
  for (int inode = 0; inode < enren_; inode++)
  {
    double utaunode = sqrt(etauw_(inode) * densinv_);
    double yp = ewdist_(inode) * viscinv_ * utaunode;
    if (yp < minyp) minyp = yp;
  }

  if (minyp > 15.0) numgpnorm_ = numgpnormow_;

  if (distype == CORE::FE::CellType::tet4)
  {
    CORE::DRT::UTILS::GaussIntegration intpointsplane(CORE::FE::CellType::tet4, 2 * numgpnorm_ - 1);
    my::intpoints_ = intpointsplane;
  }
  else  // hex8
  {
    cgp_ = Teuchos::rcp(
        new CORE::DRT::UTILS::CollectedGaussPoints(numgpnorm_ * numgpplane_ * numgpplane_));
    // get the quad9 gaussrule for the in plane integration
    CORE::DRT::UTILS::GaussIntegration intpointsplane(
        CORE::FE::CellType::quad8, 2 * numgpplane_ - 1);
    // get the quad9 gaussrule for the in normal integration
    CORE::DRT::UTILS::GaussIntegration intpointsnormal(
        CORE::FE::CellType::line3, 2 * numgpnorm_ - 1);

    // 0.9 corresponds to an angle of 25.8 deg
    if (dot1 < 0.90 && dot2 < 0.90 && dot3 < 0.90)
    {  // element, where the wall normal direction does not point in one specific element direction,
       // e.g. in corners
      cgp_->IncreaseReserved(
          (numgpnorm_ * numgpnorm_ * numgpnorm_) - (numgpnorm_ * numgpplane_ * numgpplane_));
      CORE::DRT::UTILS::GaussIntegration intpointsplane(
          CORE::FE::CellType::quad8, 2 * numgpnorm_ - 1);
      // start loop over integration points in layer
      for (CORE::DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadnorm.Point()[0], iquadplane.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else if (dot1 > dot2 && dot1 > dot3)
    {
      // start loop over integration points in layer
      for (CORE::DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadnorm.Point()[0], iquadplane.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else if (dot2 > dot3)
    {
      // start loop over integration points in layer
      for (CORE::DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadplane.Point()[0], iquadnorm.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else
    {
      // start loop over integration points in layer
      for (CORE::DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadplane.Point()[0], iquadplane.Point()[1], iquadnorm.Point()[0],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    CORE::DRT::UTILS::GaussIntegration grule(cgp_);
    my::intpoints_ = grule;
  }

  return;
}

template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Sysmat(
    const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& ebofon,
    const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgn,
    const CORE::LINALG::Matrix<nsd_, nen_>& evelaf, const CORE::LINALG::Matrix<nsd_, nen_>& evelam,
    const CORE::LINALG::Matrix<nsd_, nen_>& eveln, const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
    const CORE::LINALG::Matrix<nsd_, nen_>& fsevelaf, const CORE::LINALG::Matrix<nen_, 1>& fsescaaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& evel_hat,
    const CORE::LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
    const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& epream,
    const CORE::LINALG::Matrix<nen_, 1>& epren, const CORE::LINALG::Matrix<nen_, 1>& eprenp,
    const CORE::LINALG::Matrix<nsd_, nen_>& eaccam, const CORE::LINALG::Matrix<nen_, 1>& escaaf,
    const CORE::LINALG::Matrix<nen_, 1>& escaam, const CORE::LINALG::Matrix<nen_, 1>& escadtam,
    const CORE::LINALG::Matrix<nsd_, nen_>& eveldtam, const CORE::LINALG::Matrix<nen_, 1>& epredtam,
    const CORE::LINALG::Matrix<nen_, 1>& escabofoaf, const CORE::LINALG::Matrix<nen_, 1>& escabofon,
    const CORE::LINALG::Matrix<nsd_, nen_>& emhist, const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
    const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const CORE::LINALG::Matrix<nen_, 1>& eporo,
    const CORE::LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
    const CORE::LINALG::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
    double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
    const CORE::DRT::UTILS::GaussIntegration& intpoints)
{
  my::Sysmat(ebofoaf, eprescpgaf, ebofon, eprescpgn, evelaf, evelam, eveln, evelnp, fsevelaf,
      fsescaaf, evel_hat, ereynoldsstress_hat, epreaf, epream, epren, eprenp, eaccam, escaaf,
      escaam, escadtam, eveldtam, epredtam, escabofoaf, escabofon, emhist, edispnp, egridv, estif,
      emesh,  // -> emesh
      eforce, eporo, egradphi, ecurvature, thermpressaf, thermpressam, thermpressdtaf,
      thermpressdtam, material, Cs_delta_sq, Ci_delta_sq, Cv, isale, saccn, sveln, svelnp,
      intpoints);
  return;
}

template class DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::tet4,
    DRT::ELEMENTS::Fluid::xwall>;

BACI_NAMESPACE_CLOSE
