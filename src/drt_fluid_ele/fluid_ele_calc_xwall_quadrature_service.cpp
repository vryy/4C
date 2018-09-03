/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xwall_quadrature_service.cpp

\brief quadrature rules for xwall

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_xwall.H"

#include "fluid_ele.H"
#include "fluid_ele_xwall.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_action.H"

#include "../drt_geometry/position_array.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_mat/newtonianfluid.H"

#include "../drt_lib/drt_condition_utils.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"


/*-----------------------------------------------------------------------------*
 | Prepare custom (direction-dependent) Gauss rule                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::PrepareGaussRule()
{
  // which is the wall-normal element direction?
  // calculate jacobian at element center
  my::is_higher_order_ele_ = false;
  my::EvalShapeFuncAndDerivsAtEleCenter();
  my::is_higher_order_ele_ = true;

  // the derivative of the wall distance with respect to the local coordinates
  // shows how the local axes are oriented with respect to the wall-normal vector
  LINALG::Matrix<my::nsd_, 1> normwallrst(true);
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

  if (distype == DRT::Element::tet4)
  {
    DRT::UTILS::GaussIntegration intpointsplane(DRT::Element::tet4, 2 * numgpnorm_ - 1);
    my::intpoints_ = intpointsplane;
  }
  else  // hex8
  {
    cgp_ =
        Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(numgpnorm_ * numgpplane_ * numgpplane_));
    // get the quad9 gaussrule for the in plane integration
    DRT::UTILS::GaussIntegration intpointsplane(DRT::Element::quad8, 2 * numgpplane_ - 1);
    // get the quad9 gaussrule for the in normal integration
    DRT::UTILS::GaussIntegration intpointsnormal(DRT::Element::line3, 2 * numgpnorm_ - 1);

    // 0.9 corresponds to an angle of 25.8 deg
    if (dot1 < 0.90 && dot2 < 0.90 && dot3 < 0.90)
    {  // element, where the wall normal direction does not point in one specific element direction,
       // e.g. in corners
      cgp_->IncreaseReserved(
          (numgpnorm_ * numgpnorm_ * numgpnorm_) - (numgpnorm_ * numgpplane_ * numgpplane_));
      DRT::UTILS::GaussIntegration intpointsplane(DRT::Element::quad8, 2 * numgpnorm_ - 1);
      // start loop over integration points in layer
      for (DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
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
      for (DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
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
      for (DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
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
      for (DRT::UTILS::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (DRT::UTILS::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadplane.Point()[0], iquadplane.Point()[1], iquadnorm.Point()[0],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    DRT::UTILS::GaussIntegration grule(cgp_);
    my::intpoints_ = grule;
  }

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Sysmat(
    const LINALG::Matrix<my::nsd_, my::nen_>& ebofoaf,
    const LINALG::Matrix<my::nsd_, my::nen_>& eprescpgaf,
    const LINALG::Matrix<my::nsd_, my::nen_>& ebofon,
    const LINALG::Matrix<my::nsd_, my::nen_>& eprescpgn,
    const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>& evelam,
    const LINALG::Matrix<my::nsd_, my::nen_>& eveln,
    const LINALG::Matrix<my::nsd_, my::nen_>& evelnp,
    const LINALG::Matrix<my::nsd_, my::nen_>& fsevelaf, const LINALG::Matrix<my::nen_, 1>& fsescaaf,
    const LINALG::Matrix<my::nsd_, my::nen_>& evel_hat,
    const LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>& ereynoldsstress_hat,
    const LINALG::Matrix<my::nen_, 1>& epreaf, const LINALG::Matrix<my::nen_, 1>& epream,
    const LINALG::Matrix<my::nen_, 1>& epren, const LINALG::Matrix<my::nen_, 1>& eprenp,
    const LINALG::Matrix<my::nsd_, my::nen_>& eaccam, const LINALG::Matrix<my::nen_, 1>& escaaf,
    const LINALG::Matrix<my::nen_, 1>& escaam, const LINALG::Matrix<my::nen_, 1>& escadtam,
    const LINALG::Matrix<my::nsd_, my::nen_>& eveldtam, const LINALG::Matrix<my::nen_, 1>& epredtam,
    const LINALG::Matrix<my::nen_, 1>& escabofoaf, const LINALG::Matrix<my::nen_, 1>& escabofon,
    const LINALG::Matrix<my::nsd_, my::nen_>& emhist,
    const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_>& egridv,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_>& estif,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_>& emesh,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>& eforce, const LINALG::Matrix<my::nen_, 1>& eporo,
    const LINALG::Matrix<my::nsd_, 2 * my::nen_>& egradphi,
    const LINALG::Matrix<my::nen_, 2 * 1>& ecurvature, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
    double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
    const DRT::UTILS::GaussIntegration& intpoints)
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

template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::hex8, DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::tet4, DRT::ELEMENTS::Fluid::xwall>;
