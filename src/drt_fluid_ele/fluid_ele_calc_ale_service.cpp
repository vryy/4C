/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_ale_service.cpp

\brief ALE service routines for calculation of fluid element

\level 2
<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/



#include "fluid_ele_calc.H"
#include "fluid_ele_parameter_timint.H"


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::LinMeshMotion_2D(
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const double& press, const double& timefac,
    const double& timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_

  const double fac0 =
      densam_ * velint_(0) - rhsmom_(0) * fldparatimint_->Dt() * fldparatimint_->Theta();
  const double fac1 =
      densam_ * velint_(1) - rhsmom_(1) * fldparatimint_->Dt() * fldparatimint_->Theta();

  // mass + rhs
  for (int vi = 0; vi < nen_; ++vi)
  {
    const int tvi = 3 * vi;
    const int tvip = tvi + 1;

    const double v = fac_ * funct_(vi);
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int tui = 3 * ui;
      const int tuip = tui + 1;

      emesh(tvi, tui) += v * fac0 * derxy_(0, ui);
      emesh(tvi, tuip) += v * fac0 * derxy_(1, ui);

      emesh(tvip, tui) += v * fac1 * derxy_(0, ui);
      emesh(tvip, tuip) += v * fac1 * derxy_(1, ui);
    }
  }

  vderiv_.MultiplyNT(evelaf, deriv_);

  const double vderiv_0_0 = vderiv_(0, 0);
  const double vderiv_0_1 = vderiv_(0, 1);
  const double vderiv_1_0 = vderiv_(1, 0);
  const double vderiv_1_1 = vderiv_(1, 1);

  {
    const double convvelint_0 = convvelint_(0);
    const double convvelint_1 = convvelint_(1);
    const double densaftimefacfac_det = densaf_ * timefacfac / det_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int tvi = 3 * vi;
      const int tvip = tvi + 1;
      const double v = densaftimefacfac_det * funct_(vi);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int tui = 3 * ui;
        const int tuip = tui + 1;

        emesh(tvi, tui) +=
            v * (+convvelint_1 * (-vderiv_0_0 * deriv_(1, ui) + vderiv_0_1 * deriv_(0, ui)));

        emesh(tvi, tuip) +=
            v * (+convvelint_0 * (-vderiv_0_0 * deriv_(1, ui) + vderiv_0_1 * deriv_(0, ui)));

        emesh(tvip, tui) +=
            v * (+convvelint_1 * (-vderiv_1_0 * deriv_(1, ui) + vderiv_1_1 * deriv_(0, ui)));

        emesh(tvip, tuip) +=
            v * (+convvelint_0 * (-vderiv_1_0 * deriv_(1, ui) + vderiv_1_1 * deriv_(0, ui)));
      }
    }
  }

  // pressure
  const double v = press * timefacfac / det_;
  for (int vi = 0; vi < nen_; ++vi)
  {
    const int tvi = 3 * vi;
    const int tvip = tvi + 1;
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int tui = 3 * ui;
      emesh(tvi, tui + 1) += v * (deriv_(0, vi) * deriv_(1, ui) - deriv_(0, ui) * deriv_(1, vi));
      emesh(tvip, tui) += v * (deriv_(0, vi) * deriv_(1, ui) - deriv_(0, ui) * deriv_(1, vi));
    }
  }

  // div u
  const double timefacfac_det = timefacfac / det_;
  for (int vi = 0; vi < nen_; ++vi)
  {
    const int tvipp = 3 * vi + 2;
    const double v = timefacfac_det * funct_(vi);
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int tui = 3 * ui;
      emesh(tvipp, tui) += v * (deriv_(0, ui) * vderiv_1_1 - deriv_(1, ui) * vderiv_1_0);

      emesh(tvipp, tui + 1) += v * (deriv_(0, ui) * vderiv_0_1 - deriv_(1, ui) * vderiv_0_0);
    }
  }


  return;
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::LinMeshMotion_3D(
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const double& press, const double& timefac,
    const double& timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_

  // mass + rhs
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = fac_ * funct_(vi, 0);
    const double fac0 =
        v * (densam_ * velint_(0) - rhsmom_(0) * fldparatimint_->Dt() * fldparatimint_->Theta());
    for (int ui = 0; ui < nen_; ++ui)
    {
      emesh(vi * 4, ui * 4) += fac0 * derxy_(0, ui);
      emesh(vi * 4, ui * 4 + 1) += fac0 * derxy_(1, ui);
      emesh(vi * 4, ui * 4 + 2) += fac0 * derxy_(2, ui);
    }

    const double fac1 =
        v * (densam_ * velint_(1) - rhsmom_(1) * fldparatimint_->Dt() * fldparatimint_->Theta());
    for (int ui = 0; ui < nen_; ++ui)
    {
      emesh(vi * 4 + 1, ui * 4) += fac1 * derxy_(0, ui);
      emesh(vi * 4 + 1, ui * 4 + 1) += fac1 * derxy_(1, ui);
      emesh(vi * 4 + 1, ui * 4 + 2) += fac1 * derxy_(2, ui);
    }

    const double fac2 =
        v * (densam_ * velint_(2) - rhsmom_(2) * fldparatimint_->Dt() * fldparatimint_->Theta());
    for (int ui = 0; ui < nen_; ++ui)
    {
      emesh(vi * 4 + 2, ui * 4) += fac2 * derxy_(0, ui);
      emesh(vi * 4 + 2, ui * 4 + 1) += fac2 * derxy_(1, ui);
      emesh(vi * 4 + 2, ui * 4 + 2) += fac2 * derxy_(2, ui);
    }
  }

  // vderiv_  = sum(evelaf(i,k) * deriv_(j,k), k);
  vderiv_.MultiplyNT(evelaf, deriv_);

#define derxjm_(r, c, d, i) derxjm_##r##c##d(i)

#define derxjm_001(ui) (deriv_(2, ui) * xjm_1_2 - deriv_(1, ui) * xjm_2_2)
#define derxjm_002(ui) (deriv_(1, ui) * xjm_2_1 - deriv_(2, ui) * xjm_1_1)

#define derxjm_100(ui) (deriv_(1, ui) * xjm_2_2 - deriv_(2, ui) * xjm_1_2)
#define derxjm_102(ui) (deriv_(2, ui) * xjm_1_0 - deriv_(1, ui) * xjm_2_0)

#define derxjm_200(ui) (deriv_(2, ui) * xjm_1_1 - deriv_(1, ui) * xjm_2_1)
#define derxjm_201(ui) (deriv_(1, ui) * xjm_2_0 - deriv_(2, ui) * xjm_1_0)

#define derxjm_011(ui) (deriv_(0, ui) * xjm_2_2 - deriv_(2, ui) * xjm_0_2)
#define derxjm_012(ui) (deriv_(2, ui) * xjm_0_1 - deriv_(0, ui) * xjm_2_1)

#define derxjm_110(ui) (deriv_(2, ui) * xjm_0_2 - deriv_(0, ui) * xjm_2_2)
#define derxjm_112(ui) (deriv_(0, ui) * xjm_2_0 - deriv_(2, ui) * xjm_0_0)

#define derxjm_210(ui) (deriv_(0, ui) * xjm_2_1 - deriv_(2, ui) * xjm_0_1)
#define derxjm_211(ui) (deriv_(2, ui) * xjm_0_0 - deriv_(0, ui) * xjm_2_0)

#define derxjm_021(ui) (deriv_(1, ui) * xjm_0_2 - deriv_(0, ui) * xjm_1_2)
#define derxjm_022(ui) (deriv_(0, ui) * xjm_1_1 - deriv_(1, ui) * xjm_0_1)

#define derxjm_120(ui) (deriv_(0, ui) * xjm_1_2 - deriv_(1, ui) * xjm_0_2)
#define derxjm_122(ui) (deriv_(1, ui) * xjm_0_0 - deriv_(0, ui) * xjm_1_0)

#define derxjm_220(ui) (deriv_(1, ui) * xjm_0_1 - deriv_(0, ui) * xjm_1_1)
#define derxjm_221(ui) (deriv_(0, ui) * xjm_1_0 - deriv_(1, ui) * xjm_0_0)

  const double vderiv_0_0 = vderiv_(0, 0);
  const double vderiv_0_1 = vderiv_(0, 1);
  const double vderiv_0_2 = vderiv_(0, 2);
  const double vderiv_1_0 = vderiv_(1, 0);
  const double vderiv_1_1 = vderiv_(1, 1);
  const double vderiv_1_2 = vderiv_(1, 2);
  const double vderiv_2_0 = vderiv_(2, 0);
  const double vderiv_2_1 = vderiv_(2, 1);
  const double vderiv_2_2 = vderiv_(2, 2);

  const double xjm_0_0 = xjm_(0, 0);
  const double xjm_0_1 = xjm_(0, 1);
  const double xjm_0_2 = xjm_(0, 2);
  const double xjm_1_0 = xjm_(1, 0);
  const double xjm_1_1 = xjm_(1, 1);
  const double xjm_1_2 = xjm_(1, 2);
  const double xjm_2_0 = xjm_(2, 0);
  const double xjm_2_1 = xjm_(2, 1);
  const double xjm_2_2 = xjm_(2, 2);

  {
    const double convvelint_0 = convvelint_(0);
    const double convvelint_1 = convvelint_(1);
    const double convvelint_2 = convvelint_(2);
    const double denstimefacfac_det = densaf_ * timefacfac / det_;

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v00 =
          +convvelint_1 * (vderiv_0_0 * derxjm_(0, 0, 1, ui) + vderiv_0_1 * derxjm_(0, 1, 1, ui) +
                              vderiv_0_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (vderiv_0_0 * derxjm_(0, 0, 2, ui) + vderiv_0_1 * derxjm_(0, 1, 2, ui) +
                             vderiv_0_2 * derxjm_(0, 2, 2, ui));
      const double v01 =
          +convvelint_0 * (vderiv_0_0 * derxjm_(1, 0, 0, ui) + vderiv_0_1 * derxjm_(1, 1, 0, ui) +
                              vderiv_0_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (vderiv_0_0 * derxjm_(1, 0, 2, ui) + vderiv_0_1 * derxjm_(1, 1, 2, ui) +
                             vderiv_0_2 * derxjm_(1, 2, 2, ui));
      const double v02 =
          +convvelint_0 * (vderiv_0_0 * derxjm_(2, 0, 0, ui) + vderiv_0_1 * derxjm_(2, 1, 0, ui) +
                              vderiv_0_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (vderiv_0_0 * derxjm_(2, 0, 1, ui) + vderiv_0_1 * derxjm_(2, 1, 1, ui) +
                             vderiv_0_2 * derxjm_(2, 2, 1, ui));
      const double v10 =
          +convvelint_1 * (vderiv_1_0 * derxjm_(0, 0, 1, ui) + vderiv_1_1 * derxjm_(0, 1, 1, ui) +
                              vderiv_1_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (vderiv_1_0 * derxjm_(0, 0, 2, ui) + vderiv_1_1 * derxjm_(0, 1, 2, ui) +
                             vderiv_1_2 * derxjm_(0, 2, 2, ui));
      const double v11 =
          +convvelint_0 * (vderiv_1_0 * derxjm_(1, 0, 0, ui) + vderiv_1_1 * derxjm_(1, 1, 0, ui) +
                              vderiv_1_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (vderiv_1_0 * derxjm_(1, 0, 2, ui) + vderiv_1_1 * derxjm_(1, 1, 2, ui) +
                             vderiv_1_2 * derxjm_(1, 2, 2, ui));
      const double v12 =
          +convvelint_0 * (vderiv_1_0 * derxjm_(2, 0, 0, ui) + vderiv_1_1 * derxjm_(2, 1, 0, ui) +
                              vderiv_1_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (vderiv_1_0 * derxjm_(2, 0, 1, ui) + vderiv_1_1 * derxjm_(2, 1, 1, ui) +
                             vderiv_1_2 * derxjm_(2, 2, 1, ui));
      const double v20 =
          +convvelint_1 * (vderiv_2_0 * derxjm_(0, 0, 1, ui) + vderiv_2_1 * derxjm_(0, 1, 1, ui) +
                              vderiv_2_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (vderiv_2_0 * derxjm_(0, 0, 2, ui) + vderiv_2_1 * derxjm_(0, 1, 2, ui) +
                             vderiv_2_2 * derxjm_(0, 2, 2, ui));
      const double v21 =
          +convvelint_0 * (vderiv_2_0 * derxjm_(1, 0, 0, ui) + vderiv_2_1 * derxjm_(1, 1, 0, ui) +
                              vderiv_2_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (vderiv_2_0 * derxjm_(1, 0, 2, ui) + vderiv_2_1 * derxjm_(1, 1, 2, ui) +
                             vderiv_2_2 * derxjm_(1, 2, 2, ui));
      const double v22 =
          +convvelint_0 * (vderiv_2_0 * derxjm_(2, 0, 0, ui) + vderiv_2_1 * derxjm_(2, 1, 0, ui) +
                              vderiv_2_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (vderiv_2_0 * derxjm_(2, 0, 1, ui) + vderiv_2_1 * derxjm_(2, 1, 1, ui) +
                             vderiv_2_2 * derxjm_(2, 2, 1, ui));

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = denstimefacfac_det * funct_(vi);

        emesh(vi * 4 + 0, ui * 4 + 0) += v * v00;
        emesh(vi * 4 + 0, ui * 4 + 1) += v * v01;
        emesh(vi * 4 + 0, ui * 4 + 2) += v * v02;

        emesh(vi * 4 + 1, ui * 4 + 0) += v * v10;
        emesh(vi * 4 + 1, ui * 4 + 1) += v * v11;
        emesh(vi * 4 + 1, ui * 4 + 2) += v * v12;

        emesh(vi * 4 + 2, ui * 4 + 0) += v * v20;
        emesh(vi * 4 + 2, ui * 4 + 1) += v * v21;
        emesh(vi * 4 + 2, ui * 4 + 2) += v * v22;
      }
    }
  }

  // viscosity

  const double xji_00 = xji_(0, 0);
  const double xji_01 = xji_(0, 1);
  const double xji_02 = xji_(0, 2);
  const double xji_10 = xji_(1, 0);
  const double xji_11 = xji_(1, 1);
  const double xji_12 = xji_(1, 2);
  const double xji_20 = xji_(2, 0);
  const double xji_21 = xji_(2, 1);
  const double xji_22 = xji_(2, 2);

  // part 1: derivative of 1/det
  {
    const double vderxy_0_0 = 2.0 * vderxy_(0, 0);
    const double vderxy_1_1 = 2.0 * vderxy_(1, 1);
    const double vderxy_2_2 = 2.0 * vderxy_(2, 2);
    const double vderxy_0_1 = vderxy_(0, 1) + vderxy_(1, 0);
    const double vderxy_0_2 = vderxy_(0, 2) + vderxy_(2, 0);
    const double vderxy_1_2 = vderxy_(1, 2) + vderxy_(2, 1);

    const double v = visceff_ * timefac * fac_;
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double derinvJ0 =
          -v * (deriv_(0, ui) * xji_00 + deriv_(1, ui) * xji_01 + deriv_(2, ui) * xji_02);
      const double derinvJ1 =
          -v * (deriv_(0, ui) * xji_10 + deriv_(1, ui) * xji_11 + deriv_(2, ui) * xji_12);
      const double derinvJ2 =
          -v * (deriv_(0, ui) * xji_20 + deriv_(1, ui) * xji_21 + deriv_(2, ui) * xji_22);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double visres0 =
            derxy_(0, vi) * vderxy_0_0 + derxy_(1, vi) * vderxy_0_1 + derxy_(2, vi) * vderxy_0_2;
        const double visres1 =
            derxy_(0, vi) * vderxy_0_1 + derxy_(1, vi) * vderxy_1_1 + derxy_(2, vi) * vderxy_1_2;
        const double visres2 =
            derxy_(0, vi) * vderxy_0_2 + derxy_(1, vi) * vderxy_1_2 + derxy_(2, vi) * vderxy_2_2;
        emesh(vi * 4 + 0, ui * 4 + 0) += derinvJ0 * visres0;
        emesh(vi * 4 + 1, ui * 4 + 0) += derinvJ0 * visres1;
        emesh(vi * 4 + 2, ui * 4 + 0) += derinvJ0 * visres2;

        emesh(vi * 4 + 0, ui * 4 + 1) += derinvJ1 * visres0;
        emesh(vi * 4 + 1, ui * 4 + 1) += derinvJ1 * visres1;
        emesh(vi * 4 + 2, ui * 4 + 1) += derinvJ1 * visres2;

        emesh(vi * 4 + 0, ui * 4 + 2) += derinvJ2 * visres0;
        emesh(vi * 4 + 1, ui * 4 + 2) += derinvJ2 * visres1;
        emesh(vi * 4 + 2, ui * 4 + 2) += derinvJ2 * visres2;
      }
    }
  }

  // part 2: derivative of viscosity residual

  {
    const double v = timefacfac * visceff_ / det_;
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double ui_derxjm_001 = (deriv_(2, ui) * xjm_1_2 - deriv_(1, ui) * xjm_2_2);
      const double ui_derxjm_002 = (deriv_(1, ui) * xjm_2_1 - deriv_(2, ui) * xjm_1_1);
      const double ui_derxjm_100 = (deriv_(1, ui) * xjm_2_2 - deriv_(2, ui) * xjm_1_2);
      const double ui_derxjm_102 = (deriv_(2, ui) * xjm_1_0 - deriv_(1, ui) * xjm_2_0);
      const double ui_derxjm_200 = (deriv_(2, ui) * xjm_1_1 - deriv_(1, ui) * xjm_2_1);
      const double ui_derxjm_201 = (deriv_(1, ui) * xjm_2_0 - deriv_(2, ui) * xjm_1_0);
      const double ui_derxjm_011 = (deriv_(0, ui) * xjm_2_2 - deriv_(2, ui) * xjm_0_2);
      const double ui_derxjm_012 = (deriv_(2, ui) * xjm_0_1 - deriv_(0, ui) * xjm_2_1);
      const double ui_derxjm_110 = (deriv_(2, ui) * xjm_0_2 - deriv_(0, ui) * xjm_2_2);
      const double ui_derxjm_112 = (deriv_(0, ui) * xjm_2_0 - deriv_(2, ui) * xjm_0_0);
      const double ui_derxjm_210 = (deriv_(0, ui) * xjm_2_1 - deriv_(2, ui) * xjm_0_1);
      const double ui_derxjm_211 = (deriv_(2, ui) * xjm_0_0 - deriv_(0, ui) * xjm_2_0);
      const double ui_derxjm_021 = (deriv_(1, ui) * xjm_0_2 - deriv_(0, ui) * xjm_1_2);
      const double ui_derxjm_022 = (deriv_(0, ui) * xjm_1_1 - deriv_(1, ui) * xjm_0_1);
      const double ui_derxjm_120 = (deriv_(0, ui) * xjm_1_2 - deriv_(1, ui) * xjm_0_2);
      const double ui_derxjm_122 = (deriv_(1, ui) * xjm_0_0 - deriv_(0, ui) * xjm_1_0);
      const double ui_derxjm_220 = (deriv_(1, ui) * xjm_0_1 - deriv_(0, ui) * xjm_1_1);
      const double ui_derxjm_221 = (deriv_(0, ui) * xjm_1_0 - deriv_(1, ui) * xjm_0_0);

      {
        const double v0 =
            -vderiv_0_0 * (xji_10 * ui_derxjm_100 + xji_10 * ui_derxjm_100 +
                              xji_20 * ui_derxjm_200 + xji_20 * ui_derxjm_200) -
            vderiv_0_1 * (xji_11 * ui_derxjm_100 + xji_10 * ui_derxjm_110 + xji_21 * ui_derxjm_200 +
                             xji_20 * ui_derxjm_210) -
            vderiv_0_2 * (xji_12 * ui_derxjm_100 + xji_10 * ui_derxjm_120 + xji_22 * ui_derxjm_200 +
                             xji_20 * ui_derxjm_220) -
            vderiv_1_0 * (ui_derxjm_100 * xji_00) - vderiv_1_1 * (ui_derxjm_100 * xji_01) -
            vderiv_1_2 * (ui_derxjm_100 * xji_02) - vderiv_2_0 * (ui_derxjm_200 * xji_00) -
            vderiv_2_1 * (ui_derxjm_200 * xji_01) - vderiv_2_2 * (ui_derxjm_200 * xji_02);
        const double v1 =
            -vderiv_0_0 * (xji_10 * ui_derxjm_110 + xji_11 * ui_derxjm_100 +
                              xji_20 * ui_derxjm_210 + xji_21 * ui_derxjm_200) -
            vderiv_0_1 * (xji_11 * ui_derxjm_110 + xji_11 * ui_derxjm_110 + xji_21 * ui_derxjm_210 +
                             xji_21 * ui_derxjm_210) -
            vderiv_0_2 * (xji_12 * ui_derxjm_110 + xji_11 * ui_derxjm_120 + xji_22 * ui_derxjm_210 +
                             xji_21 * ui_derxjm_220) -
            vderiv_1_0 * (ui_derxjm_110 * xji_00) - vderiv_1_1 * (ui_derxjm_110 * xji_01) -
            vderiv_1_2 * (ui_derxjm_110 * xji_02) - vderiv_2_0 * (ui_derxjm_210 * xji_00) -
            vderiv_2_1 * (ui_derxjm_210 * xji_01) - vderiv_2_2 * (ui_derxjm_210 * xji_02);
        const double v2 =
            -vderiv_0_0 * (xji_10 * ui_derxjm_120 + xji_12 * ui_derxjm_100 +
                              xji_20 * ui_derxjm_220 + xji_22 * ui_derxjm_200) -
            vderiv_0_1 * (xji_11 * ui_derxjm_120 + xji_12 * ui_derxjm_110 + xji_21 * ui_derxjm_220 +
                             xji_22 * ui_derxjm_210) -
            vderiv_0_2 * (xji_12 * ui_derxjm_120 + xji_12 * ui_derxjm_120 + xji_22 * ui_derxjm_220 +
                             xji_22 * ui_derxjm_220) -
            vderiv_1_0 * (ui_derxjm_120 * xji_00) - vderiv_1_1 * (ui_derxjm_120 * xji_01) -
            vderiv_1_2 * (ui_derxjm_120 * xji_02) - vderiv_2_0 * (ui_derxjm_220 * xji_00) -
            vderiv_2_1 * (ui_derxjm_220 * xji_01) - vderiv_2_2 * (ui_derxjm_220 * xji_02);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 0, ui * 4 + 0) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 = -vderiv_0_0 * (2 * ui_derxjm_001 * xji_00 + 2 * ui_derxjm_001 * xji_00 +
                                            xji_20 * ui_derxjm_201 + xji_20 * ui_derxjm_201) -
                          vderiv_0_1 * (2 * ui_derxjm_011 * xji_00 + 2 * ui_derxjm_001 * xji_01 +
                                           xji_21 * ui_derxjm_201 + xji_20 * ui_derxjm_211) -
                          vderiv_0_2 * (2 * ui_derxjm_021 * xji_00 + 2 * ui_derxjm_001 * xji_02 +
                                           xji_22 * ui_derxjm_201 + xji_20 * ui_derxjm_221) -
                          vderiv_1_0 * (ui_derxjm_001 * xji_10) -
                          vderiv_1_1 * (ui_derxjm_011 * xji_10) -
                          vderiv_1_2 * (ui_derxjm_021 * xji_10) -
                          vderiv_2_0 * (ui_derxjm_201 * xji_00 + ui_derxjm_001 * xji_20) -
                          vderiv_2_1 * (ui_derxjm_201 * xji_01 + ui_derxjm_011 * xji_20) -
                          vderiv_2_2 * (ui_derxjm_201 * xji_02 + ui_derxjm_021 * xji_20);
        const double v1 = -vderiv_0_0 * (2 * ui_derxjm_011 * xji_00 + 2 * ui_derxjm_001 * xji_01 +
                                            xji_21 * ui_derxjm_201 + xji_20 * ui_derxjm_211) -
                          vderiv_0_1 * (2 * ui_derxjm_011 * xji_01 + 2 * ui_derxjm_011 * xji_01 +
                                           xji_21 * ui_derxjm_211 + xji_21 * ui_derxjm_211) -
                          vderiv_0_2 * (2 * ui_derxjm_011 * xji_02 + 2 * ui_derxjm_021 * xji_01 +
                                           xji_21 * ui_derxjm_221 + xji_22 * ui_derxjm_211) -
                          vderiv_1_0 * (ui_derxjm_001 * xji_11) -
                          vderiv_1_1 * (ui_derxjm_011 * xji_11) -
                          vderiv_1_2 * (ui_derxjm_021 * xji_11) -
                          vderiv_2_0 * (ui_derxjm_211 * xji_00 + ui_derxjm_001 * xji_21) -
                          vderiv_2_1 * (ui_derxjm_211 * xji_01 + ui_derxjm_011 * xji_21) -
                          vderiv_2_2 * (ui_derxjm_211 * xji_02 + ui_derxjm_021 * xji_21);
        const double v2 = -vderiv_0_0 * (2 * ui_derxjm_021 * xji_00 + 2 * ui_derxjm_001 * xji_02 +
                                            xji_22 * ui_derxjm_201 + xji_20 * ui_derxjm_221) -
                          vderiv_0_1 * (2 * ui_derxjm_011 * xji_02 + 2 * ui_derxjm_021 * xji_01 +
                                           xji_21 * ui_derxjm_221 + xji_22 * ui_derxjm_211) -
                          vderiv_0_2 * (2 * ui_derxjm_021 * xji_02 + 2 * ui_derxjm_021 * xji_02 +
                                           xji_22 * ui_derxjm_221 + xji_22 * ui_derxjm_221) -
                          vderiv_1_0 * (ui_derxjm_001 * xji_12) -
                          vderiv_1_1 * (ui_derxjm_011 * xji_12) -
                          vderiv_1_2 * (ui_derxjm_021 * xji_12) -
                          vderiv_2_0 * (ui_derxjm_221 * xji_00 + ui_derxjm_001 * xji_22) -
                          vderiv_2_1 * (ui_derxjm_221 * xji_01 + ui_derxjm_011 * xji_22) -
                          vderiv_2_2 * (ui_derxjm_221 * xji_02 + ui_derxjm_021 * xji_22);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 0, ui * 4 + 1) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 = -vderiv_0_0 * (2 * ui_derxjm_002 * xji_00 + 2 * ui_derxjm_002 * xji_00 +
                                            xji_10 * ui_derxjm_102 + xji_10 * ui_derxjm_102) -
                          vderiv_0_1 * (2 * ui_derxjm_012 * xji_00 + 2 * ui_derxjm_002 * xji_01 +
                                           xji_11 * ui_derxjm_102 + xji_10 * ui_derxjm_112) -
                          vderiv_0_2 * (2 * ui_derxjm_022 * xji_00 + 2 * ui_derxjm_002 * xji_02 +
                                           xji_12 * ui_derxjm_102 + xji_10 * ui_derxjm_122) -
                          vderiv_1_0 * (ui_derxjm_002 * xji_10 + ui_derxjm_102 * xji_00) -
                          vderiv_1_1 * (ui_derxjm_012 * xji_10 + ui_derxjm_102 * xji_01) -
                          vderiv_1_2 * (ui_derxjm_022 * xji_10 + ui_derxjm_102 * xji_02) -
                          vderiv_2_0 * (ui_derxjm_002 * xji_20) -
                          vderiv_2_1 * (ui_derxjm_012 * xji_20) -
                          vderiv_2_2 * (ui_derxjm_022 * xji_20);
        const double v1 = -vderiv_0_0 * (2 * ui_derxjm_012 * xji_00 + 2 * ui_derxjm_002 * xji_01 +
                                            xji_11 * ui_derxjm_102 + xji_10 * ui_derxjm_112) -
                          vderiv_0_1 * (2 * ui_derxjm_012 * xji_01 + 2 * ui_derxjm_012 * xji_01 +
                                           xji_11 * ui_derxjm_112 + xji_11 * ui_derxjm_112) -
                          vderiv_0_2 * (2 * ui_derxjm_012 * xji_02 + 2 * ui_derxjm_022 * xji_01 +
                                           xji_11 * ui_derxjm_122 + xji_12 * ui_derxjm_112) -
                          vderiv_1_0 * (ui_derxjm_002 * xji_11 + ui_derxjm_112 * xji_00) -
                          vderiv_1_1 * (ui_derxjm_012 * xji_11 + ui_derxjm_112 * xji_01) -
                          vderiv_1_2 * (ui_derxjm_022 * xji_11 + ui_derxjm_112 * xji_02) -
                          vderiv_2_0 * (ui_derxjm_002 * xji_21) -
                          vderiv_2_1 * (ui_derxjm_012 * xji_21) -
                          vderiv_2_2 * (ui_derxjm_022 * xji_21);
        const double v2 = -vderiv_0_0 * (2 * ui_derxjm_022 * xji_00 + 2 * ui_derxjm_002 * xji_02 +
                                            xji_12 * ui_derxjm_102 + xji_10 * ui_derxjm_122) -
                          vderiv_0_1 * (2 * ui_derxjm_012 * xji_02 + 2 * ui_derxjm_022 * xji_01 +
                                           xji_11 * ui_derxjm_122 + xji_12 * ui_derxjm_112) -
                          vderiv_0_2 * (2 * ui_derxjm_022 * xji_02 + 2 * ui_derxjm_022 * xji_02 +
                                           xji_12 * ui_derxjm_122 + xji_12 * ui_derxjm_122) -
                          vderiv_1_0 * (ui_derxjm_002 * xji_12 + ui_derxjm_122 * xji_00) -
                          vderiv_1_1 * (ui_derxjm_012 * xji_12 + ui_derxjm_122 * xji_01) -
                          vderiv_1_2 * (ui_derxjm_022 * xji_12 + ui_derxjm_122 * xji_02) -
                          vderiv_2_0 * (ui_derxjm_002 * xji_22) -
                          vderiv_2_1 * (ui_derxjm_012 * xji_22) -
                          vderiv_2_2 * (ui_derxjm_022 * xji_22);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 0, ui * 4 + 2) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 = -vderiv_0_0 * (ui_derxjm_100 * xji_00) -
                          vderiv_0_1 * (ui_derxjm_110 * xji_00) -
                          vderiv_0_2 * (ui_derxjm_120 * xji_00) -
                          vderiv_1_0 * (2 * xji_10 * ui_derxjm_100 + 2 * xji_10 * ui_derxjm_100 +
                                           xji_20 * ui_derxjm_200 + xji_20 * ui_derxjm_200) -
                          vderiv_1_1 * (2 * xji_11 * ui_derxjm_100 + 2 * xji_10 * ui_derxjm_110 +
                                           xji_21 * ui_derxjm_200 + xji_20 * ui_derxjm_210) -
                          vderiv_1_2 * (2 * xji_12 * ui_derxjm_100 + 2 * xji_10 * ui_derxjm_120 +
                                           xji_22 * ui_derxjm_200 + xji_20 * ui_derxjm_220) -
                          vderiv_2_0 * (ui_derxjm_200 * xji_10 + ui_derxjm_100 * xji_20) -
                          vderiv_2_1 * (ui_derxjm_200 * xji_11 + ui_derxjm_110 * xji_20) -
                          vderiv_2_2 * (ui_derxjm_200 * xji_12 + ui_derxjm_120 * xji_20);
        const double v1 = -vderiv_0_0 * (ui_derxjm_100 * xji_01) -
                          vderiv_0_1 * (ui_derxjm_110 * xji_01) -
                          vderiv_0_2 * (ui_derxjm_120 * xji_01) -
                          vderiv_1_0 * (2 * xji_10 * ui_derxjm_110 + 2 * xji_11 * ui_derxjm_100 +
                                           xji_20 * ui_derxjm_210 + xji_21 * ui_derxjm_200) -
                          vderiv_1_1 * (2 * xji_11 * ui_derxjm_110 + 2 * xji_11 * ui_derxjm_110 +
                                           xji_21 * ui_derxjm_210 + xji_21 * ui_derxjm_210) -
                          vderiv_1_2 * (2 * xji_12 * ui_derxjm_110 + 2 * xji_11 * ui_derxjm_120 +
                                           xji_22 * ui_derxjm_210 + xji_21 * ui_derxjm_220) -
                          vderiv_2_0 * (ui_derxjm_210 * xji_10 + ui_derxjm_100 * xji_21) -
                          vderiv_2_1 * (ui_derxjm_210 * xji_11 + ui_derxjm_110 * xji_21) -
                          vderiv_2_2 * (ui_derxjm_210 * xji_12 + ui_derxjm_120 * xji_21);
        const double v2 = -vderiv_0_0 * (ui_derxjm_100 * xji_02) -
                          vderiv_0_1 * (ui_derxjm_110 * xji_02) -
                          vderiv_0_2 * (ui_derxjm_120 * xji_02) -
                          vderiv_1_0 * (2 * xji_10 * ui_derxjm_120 + 2 * xji_12 * ui_derxjm_100 +
                                           xji_20 * ui_derxjm_220 + xji_22 * ui_derxjm_200) -
                          vderiv_1_1 * (2 * xji_11 * ui_derxjm_120 + 2 * xji_12 * ui_derxjm_110 +
                                           xji_21 * ui_derxjm_220 + xji_22 * ui_derxjm_210) -
                          vderiv_1_2 * (2 * xji_12 * ui_derxjm_120 + 2 * xji_12 * ui_derxjm_120 +
                                           xji_22 * ui_derxjm_220 + xji_22 * ui_derxjm_220) -
                          vderiv_2_0 * (ui_derxjm_220 * xji_10 + ui_derxjm_100 * xji_22) -
                          vderiv_2_1 * (ui_derxjm_220 * xji_11 + ui_derxjm_110 * xji_22) -
                          vderiv_2_2 * (ui_derxjm_220 * xji_12 + ui_derxjm_120 * xji_22);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 1, ui * 4 + 0) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 =
            -vderiv_0_0 * (ui_derxjm_001 * xji_10) - vderiv_0_1 * (ui_derxjm_001 * xji_11) -
            vderiv_0_2 * (ui_derxjm_001 * xji_12) -
            vderiv_1_0 * (xji_00 * ui_derxjm_001 + xji_00 * ui_derxjm_001 + xji_20 * ui_derxjm_201 +
                             xji_20 * ui_derxjm_201) -
            vderiv_1_1 * (xji_01 * ui_derxjm_001 + xji_00 * ui_derxjm_011 + xji_21 * ui_derxjm_201 +
                             xji_20 * ui_derxjm_211) -
            vderiv_1_2 * (xji_02 * ui_derxjm_001 + xji_00 * ui_derxjm_021 + xji_22 * ui_derxjm_201 +
                             xji_20 * ui_derxjm_221) -
            vderiv_2_0 * (ui_derxjm_201 * xji_10) - vderiv_2_1 * (ui_derxjm_201 * xji_11) -
            vderiv_2_2 * (ui_derxjm_201 * xji_12);
        const double v1 =
            -vderiv_0_0 * (ui_derxjm_011 * xji_10) - vderiv_0_1 * (ui_derxjm_011 * xji_11) -
            vderiv_0_2 * (ui_derxjm_011 * xji_12) -
            vderiv_1_0 * (xji_00 * ui_derxjm_011 + xji_01 * ui_derxjm_001 + xji_20 * ui_derxjm_211 +
                             xji_21 * ui_derxjm_201) -
            vderiv_1_1 * (xji_01 * ui_derxjm_011 + xji_01 * ui_derxjm_011 + xji_21 * ui_derxjm_211 +
                             xji_21 * ui_derxjm_211) -
            vderiv_1_2 * (xji_02 * ui_derxjm_011 + xji_01 * ui_derxjm_021 + xji_22 * ui_derxjm_211 +
                             xji_21 * ui_derxjm_221) -
            vderiv_2_0 * (ui_derxjm_211 * xji_10) - vderiv_2_1 * (ui_derxjm_211 * xji_11) -
            vderiv_2_2 * (ui_derxjm_211 * xji_12);
        const double v2 =
            -vderiv_0_0 * (ui_derxjm_021 * xji_10) - vderiv_0_1 * (ui_derxjm_021 * xji_11) -
            vderiv_0_2 * (ui_derxjm_021 * xji_12) -
            vderiv_1_0 * (xji_00 * ui_derxjm_021 + xji_02 * ui_derxjm_001 + xji_20 * ui_derxjm_221 +
                             xji_22 * ui_derxjm_201) -
            vderiv_1_1 * (xji_01 * ui_derxjm_021 + xji_02 * ui_derxjm_011 + xji_21 * ui_derxjm_221 +
                             xji_22 * ui_derxjm_211) -
            vderiv_1_2 * (xji_02 * ui_derxjm_021 + xji_02 * ui_derxjm_021 + xji_22 * ui_derxjm_221 +
                             xji_22 * ui_derxjm_221) -
            vderiv_2_0 * (ui_derxjm_221 * xji_10) - vderiv_2_1 * (ui_derxjm_221 * xji_11) -
            vderiv_2_2 * (ui_derxjm_221 * xji_12);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 1, ui * 4 + 1) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 =
            -vderiv_0_0 * (ui_derxjm_002 * xji_10 + ui_derxjm_102 * xji_00) -
            vderiv_0_1 * (ui_derxjm_002 * xji_11 + ui_derxjm_112 * xji_00) -
            vderiv_0_2 * (ui_derxjm_002 * xji_12 + ui_derxjm_122 * xji_00) -
            vderiv_1_0 * (xji_00 * ui_derxjm_002 + xji_00 * ui_derxjm_002 +
                             2 * xji_10 * ui_derxjm_102 + 2 * xji_10 * ui_derxjm_102) -
            vderiv_1_1 * (xji_01 * ui_derxjm_002 + xji_00 * ui_derxjm_012 +
                             2 * xji_11 * ui_derxjm_102 + 2 * xji_10 * ui_derxjm_112) -
            vderiv_1_2 * (xji_02 * ui_derxjm_002 + xji_00 * ui_derxjm_022 +
                             2 * xji_12 * ui_derxjm_102 + 2 * xji_10 * ui_derxjm_122) -
            vderiv_2_0 * (ui_derxjm_102 * xji_20) - vderiv_2_1 * (ui_derxjm_112 * xji_20) -
            vderiv_2_2 * (ui_derxjm_122 * xji_20);
        const double v1 =
            -vderiv_0_0 * (ui_derxjm_012 * xji_10 + ui_derxjm_102 * xji_01) -
            vderiv_0_1 * (ui_derxjm_012 * xji_11 + ui_derxjm_112 * xji_01) -
            vderiv_0_2 * (ui_derxjm_012 * xji_12 + ui_derxjm_122 * xji_01) -
            vderiv_1_0 * (xji_00 * ui_derxjm_012 + xji_01 * ui_derxjm_002 +
                             2 * xji_10 * ui_derxjm_112 + 2 * xji_11 * ui_derxjm_102) -
            vderiv_1_1 * (xji_01 * ui_derxjm_012 + xji_01 * ui_derxjm_012 +
                             2 * xji_11 * ui_derxjm_112 + 2 * xji_11 * ui_derxjm_112) -
            vderiv_1_2 * (xji_02 * ui_derxjm_012 + xji_01 * ui_derxjm_022 +
                             2 * xji_12 * ui_derxjm_112 + 2 * xji_11 * ui_derxjm_122) -
            vderiv_2_0 * (ui_derxjm_102 * xji_21) - vderiv_2_1 * (ui_derxjm_112 * xji_21) -
            vderiv_2_2 * (ui_derxjm_122 * xji_21);
        const double v2 =
            -vderiv_0_0 * (ui_derxjm_022 * xji_10 + ui_derxjm_102 * xji_02) -
            vderiv_0_1 * (ui_derxjm_022 * xji_11 + ui_derxjm_112 * xji_02) -
            vderiv_0_2 * (ui_derxjm_022 * xji_12 + ui_derxjm_122 * xji_02) -
            vderiv_1_0 * (xji_00 * ui_derxjm_022 + xji_02 * ui_derxjm_002 +
                             2 * xji_10 * ui_derxjm_122 + 2 * xji_12 * ui_derxjm_102) -
            vderiv_1_1 * (xji_01 * ui_derxjm_022 + xji_02 * ui_derxjm_012 +
                             2 * xji_11 * ui_derxjm_122 + 2 * xji_12 * ui_derxjm_112) -
            vderiv_1_2 * (xji_02 * ui_derxjm_022 + xji_02 * ui_derxjm_022 +
                             2 * xji_12 * ui_derxjm_122 + 2 * xji_12 * ui_derxjm_122) -
            vderiv_2_0 * (ui_derxjm_102 * xji_22) - vderiv_2_1 * (ui_derxjm_112 * xji_22) -
            vderiv_2_2 * (ui_derxjm_122 * xji_22);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 1, ui * 4 + 2) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 =
            -vderiv_0_0 * (ui_derxjm_200 * xji_00) - vderiv_0_1 * (ui_derxjm_210 * xji_00) -
            vderiv_0_2 * (ui_derxjm_220 * xji_00) -
            vderiv_1_0 * (ui_derxjm_200 * xji_10 + ui_derxjm_100 * xji_20) -
            vderiv_1_1 * (ui_derxjm_210 * xji_10 + ui_derxjm_100 * xji_21) -
            vderiv_1_2 * (ui_derxjm_220 * xji_10 + ui_derxjm_100 * xji_22) -
            vderiv_2_0 * (xji_10 * ui_derxjm_100 + xji_10 * ui_derxjm_100 +
                             2 * xji_20 * ui_derxjm_200 + 2 * xji_20 * ui_derxjm_200) -
            vderiv_2_1 * (xji_11 * ui_derxjm_100 + xji_10 * ui_derxjm_110 +
                             2 * xji_21 * ui_derxjm_200 + 2 * xji_20 * ui_derxjm_210) -
            vderiv_2_2 * (xji_12 * ui_derxjm_100 + xji_10 * ui_derxjm_120 +
                             2 * xji_22 * ui_derxjm_200 + 2 * xji_20 * ui_derxjm_220);
        const double v1 =
            -vderiv_0_0 * (ui_derxjm_200 * xji_01) - vderiv_0_1 * (ui_derxjm_210 * xji_01) -
            vderiv_0_2 * (ui_derxjm_220 * xji_01) -
            vderiv_1_0 * (ui_derxjm_200 * xji_11 + ui_derxjm_110 * xji_20) -
            vderiv_1_1 * (ui_derxjm_210 * xji_11 + ui_derxjm_110 * xji_21) -
            vderiv_1_2 * (ui_derxjm_220 * xji_11 + ui_derxjm_110 * xji_22) -
            vderiv_2_0 * (xji_10 * ui_derxjm_110 + xji_11 * ui_derxjm_100 +
                             2 * xji_20 * ui_derxjm_210 + 2 * xji_21 * ui_derxjm_200) -
            vderiv_2_1 * (xji_11 * ui_derxjm_110 + xji_11 * ui_derxjm_110 +
                             2 * xji_21 * ui_derxjm_210 + 2 * xji_21 * ui_derxjm_210) -
            vderiv_2_2 * (xji_12 * ui_derxjm_110 + xji_11 * ui_derxjm_120 +
                             2 * xji_22 * ui_derxjm_210 + 2 * xji_21 * ui_derxjm_220);
        const double v2 =
            -vderiv_0_0 * (ui_derxjm_200 * xji_02) - vderiv_0_1 * (ui_derxjm_210 * xji_02) -
            vderiv_0_2 * (ui_derxjm_220 * xji_02) -
            vderiv_1_0 * (ui_derxjm_200 * xji_12 + ui_derxjm_120 * xji_20) -
            vderiv_1_1 * (ui_derxjm_210 * xji_12 + ui_derxjm_120 * xji_21) -
            vderiv_1_2 * (ui_derxjm_220 * xji_12 + ui_derxjm_120 * xji_22) -
            vderiv_2_0 * (xji_10 * ui_derxjm_120 + xji_12 * ui_derxjm_100 +
                             2 * xji_20 * ui_derxjm_220 + 2 * xji_22 * ui_derxjm_200) -
            vderiv_2_1 * (xji_11 * ui_derxjm_120 + xji_12 * ui_derxjm_110 +
                             2 * xji_21 * ui_derxjm_220 + 2 * xji_22 * ui_derxjm_210) -
            vderiv_2_2 * (xji_12 * ui_derxjm_120 + xji_12 * ui_derxjm_120 +
                             2 * xji_22 * ui_derxjm_220 + 2 * xji_22 * ui_derxjm_220);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 2, ui * 4 + 0) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 =
            -vderiv_0_0 * (ui_derxjm_201 * xji_00 + ui_derxjm_001 * xji_20) -
            vderiv_0_1 * (ui_derxjm_211 * xji_00 + ui_derxjm_001 * xji_21) -
            vderiv_0_2 * (ui_derxjm_221 * xji_00 + ui_derxjm_001 * xji_22) -
            vderiv_1_0 * (ui_derxjm_201 * xji_10) - vderiv_1_1 * (ui_derxjm_211 * xji_10) -
            vderiv_1_2 * (ui_derxjm_221 * xji_10) -
            vderiv_2_0 * (xji_00 * ui_derxjm_001 + xji_00 * ui_derxjm_001 +
                             2 * xji_20 * ui_derxjm_201 + 2 * xji_20 * ui_derxjm_201) -
            vderiv_2_1 * (xji_01 * ui_derxjm_001 + xji_00 * ui_derxjm_011 +
                             2 * xji_21 * ui_derxjm_201 + 2 * xji_20 * ui_derxjm_211) -
            vderiv_2_2 * (xji_02 * ui_derxjm_001 + xji_00 * ui_derxjm_021 +
                             2 * xji_22 * ui_derxjm_201 + 2 * xji_20 * ui_derxjm_221);
        const double v1 =
            -vderiv_0_0 * (ui_derxjm_201 * xji_01 + ui_derxjm_011 * xji_20) -
            vderiv_0_1 * (ui_derxjm_211 * xji_01 + ui_derxjm_011 * xji_21) -
            vderiv_0_2 * (ui_derxjm_221 * xji_01 + ui_derxjm_011 * xji_22) -
            vderiv_1_0 * (ui_derxjm_201 * xji_11) - vderiv_1_1 * (ui_derxjm_211 * xji_11) -
            vderiv_1_2 * (ui_derxjm_221 * xji_11) -
            vderiv_2_0 * (xji_00 * ui_derxjm_011 + xji_01 * ui_derxjm_001 +
                             2 * xji_20 * ui_derxjm_211 + 2 * xji_21 * ui_derxjm_201) -
            vderiv_2_1 * (xji_01 * ui_derxjm_011 + xji_01 * ui_derxjm_011 +
                             2 * xji_21 * ui_derxjm_211 + 2 * xji_21 * ui_derxjm_211) -
            vderiv_2_2 * (xji_02 * ui_derxjm_011 + xji_01 * ui_derxjm_021 +
                             2 * xji_22 * ui_derxjm_211 + 2 * xji_21 * ui_derxjm_221);
        const double v2 =
            -vderiv_0_0 * (ui_derxjm_201 * xji_02 + ui_derxjm_021 * xji_20) -
            vderiv_0_1 * (ui_derxjm_211 * xji_02 + ui_derxjm_021 * xji_21) -
            vderiv_0_2 * (ui_derxjm_221 * xji_02 + ui_derxjm_021 * xji_22) -
            vderiv_1_0 * (ui_derxjm_201 * xji_12) - vderiv_1_1 * (ui_derxjm_211 * xji_12) -
            vderiv_1_2 * (ui_derxjm_221 * xji_12) -
            vderiv_2_0 * (xji_00 * ui_derxjm_021 + xji_02 * ui_derxjm_001 +
                             2 * xji_20 * ui_derxjm_221 + 2 * xji_22 * ui_derxjm_201) -
            vderiv_2_1 * (xji_01 * ui_derxjm_021 + xji_02 * ui_derxjm_011 +
                             2 * xji_21 * ui_derxjm_221 + 2 * xji_22 * ui_derxjm_211) -
            vderiv_2_2 * (xji_02 * ui_derxjm_021 + xji_02 * ui_derxjm_021 +
                             2 * xji_22 * ui_derxjm_221 + 2 * xji_22 * ui_derxjm_221);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 2, ui * 4 + 1) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }

      ////////////////////////////////////////////////////////////////

      {
        const double v0 =
            -vderiv_0_0 * (ui_derxjm_002 * xji_20) - vderiv_0_1 * (ui_derxjm_002 * xji_21) -
            vderiv_0_2 * (ui_derxjm_002 * xji_22) - vderiv_1_0 * (ui_derxjm_102 * xji_20) -
            vderiv_1_1 * (ui_derxjm_102 * xji_21) - vderiv_1_2 * (ui_derxjm_102 * xji_22) -
            vderiv_2_0 * (xji_00 * ui_derxjm_002 + xji_00 * ui_derxjm_002 + xji_10 * ui_derxjm_102 +
                             xji_10 * ui_derxjm_102) -
            vderiv_2_1 * (xji_01 * ui_derxjm_002 + xji_00 * ui_derxjm_012 + xji_11 * ui_derxjm_102 +
                             xji_10 * ui_derxjm_112) -
            vderiv_2_2 * (xji_02 * ui_derxjm_002 + xji_00 * ui_derxjm_022 + xji_12 * ui_derxjm_102 +
                             xji_10 * ui_derxjm_122);
        const double v1 =
            -vderiv_0_0 * (ui_derxjm_012 * xji_20) - vderiv_0_1 * (ui_derxjm_012 * xji_21) -
            vderiv_0_2 * (ui_derxjm_012 * xji_22) - vderiv_1_0 * (ui_derxjm_112 * xji_20) -
            vderiv_1_1 * (ui_derxjm_112 * xji_21) - vderiv_1_2 * (ui_derxjm_112 * xji_22) -
            vderiv_2_0 * (xji_00 * ui_derxjm_012 + xji_01 * ui_derxjm_002 + xji_10 * ui_derxjm_112 +
                             xji_11 * ui_derxjm_102) -
            vderiv_2_1 * (xji_01 * ui_derxjm_012 + xji_01 * ui_derxjm_012 + xji_11 * ui_derxjm_112 +
                             xji_11 * ui_derxjm_112) -
            vderiv_2_2 * (xji_02 * ui_derxjm_012 + xji_01 * ui_derxjm_022 + xji_12 * ui_derxjm_112 +
                             xji_11 * ui_derxjm_122);
        const double v2 =
            -vderiv_0_0 * (ui_derxjm_022 * xji_20) - vderiv_0_1 * (ui_derxjm_022 * xji_21) -
            vderiv_0_2 * (ui_derxjm_022 * xji_22) - vderiv_1_0 * (ui_derxjm_122 * xji_20) -
            vderiv_1_1 * (ui_derxjm_122 * xji_21) - vderiv_1_2 * (ui_derxjm_122 * xji_22) -
            vderiv_2_0 * (xji_00 * ui_derxjm_022 + xji_02 * ui_derxjm_002 + xji_10 * ui_derxjm_122 +
                             xji_12 * ui_derxjm_102) -
            vderiv_2_1 * (xji_01 * ui_derxjm_022 + xji_02 * ui_derxjm_012 + xji_11 * ui_derxjm_122 +
                             xji_12 * ui_derxjm_112) -
            vderiv_2_2 * (xji_02 * ui_derxjm_022 + xji_02 * ui_derxjm_022 + xji_12 * ui_derxjm_122 +
                             xji_12 * ui_derxjm_122);

        for (int vi = 0; vi < nen_; ++vi)
        {
          emesh(vi * 4 + 2, ui * 4 + 2) +=
              v * (deriv_(0, vi) * v0 + deriv_(1, vi) * v1 + deriv_(2, vi) * v2);
        }
      }
    }
  }


  // pressure
  const double v = press * timefacfac / det_;
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      emesh(vi * 4, ui * 4 + 1) +=
          v * (deriv_(0, vi) * derxjm_(0, 0, 1, ui) + deriv_(1, vi) * derxjm_(0, 1, 1, ui) +
                  deriv_(2, vi) * derxjm_(0, 2, 1, ui));
      emesh(vi * 4, ui * 4 + 2) +=
          v * (deriv_(0, vi) * derxjm_(0, 0, 2, ui) + deriv_(1, vi) * derxjm_(0, 1, 2, ui) +
                  deriv_(2, vi) * derxjm_(0, 2, 2, ui));

      emesh(vi * 4 + 1, ui * 4 + 0) +=
          v * (deriv_(0, vi) * derxjm_(1, 0, 0, ui) + deriv_(1, vi) * derxjm_(1, 1, 0, ui) +
                  deriv_(2, vi) * derxjm_(1, 2, 0, ui));
      emesh(vi * 4 + 1, ui * 4 + 2) +=
          v * (deriv_(0, vi) * derxjm_(1, 0, 2, ui) + deriv_(1, vi) * derxjm_(1, 1, 2, ui) +
                  deriv_(2, vi) * derxjm_(1, 2, 2, ui));

      emesh(vi * 4 + 2, ui * 4 + 0) +=
          v * (deriv_(0, vi) * derxjm_(2, 0, 0, ui) + deriv_(1, vi) * derxjm_(2, 1, 0, ui) +
                  deriv_(2, vi) * derxjm_(2, 2, 0, ui));
      emesh(vi * 4 + 2, ui * 4 + 1) +=
          v * (deriv_(0, vi) * derxjm_(2, 0, 1, ui) + deriv_(1, vi) * derxjm_(2, 1, 1, ui) +
                  deriv_(2, vi) * derxjm_(2, 2, 1, ui));
    }
  }

  // div u
  const double timefacfac_det = timefacfac / det_;
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = timefacfac_det * funct_(vi, 0);
    for (int ui = 0; ui < nen_; ++ui)
    {
      emesh(vi * 4 + 3, ui * 4 + 0) +=
          v * (+vderiv_1_0 * derxjm_(0, 0, 1, ui) + vderiv_1_1 * derxjm_(0, 1, 1, ui) +
                  vderiv_1_2 * derxjm_(0, 2, 1, ui) + vderiv_2_0 * derxjm_(0, 0, 2, ui) +
                  vderiv_2_1 * derxjm_(0, 1, 2, ui) + vderiv_2_2 * derxjm_(0, 2, 2, ui));

      emesh(vi * 4 + 3, ui * 4 + 1) +=
          v * (+vderiv_0_0 * derxjm_(1, 0, 0, ui) + vderiv_0_1 * derxjm_(1, 1, 0, ui) +
                  vderiv_0_2 * derxjm_(1, 2, 0, ui) + vderiv_2_0 * derxjm_(1, 0, 2, ui) +
                  vderiv_2_1 * derxjm_(1, 1, 2, ui) + vderiv_2_2 * derxjm_(1, 2, 2, ui));

      emesh(vi * 4 + 3, ui * 4 + 2) +=
          v * (+vderiv_0_0 * derxjm_(2, 0, 0, ui) + vderiv_0_1 * derxjm_(2, 1, 0, ui) +
                  vderiv_0_2 * derxjm_(2, 2, 0, ui) + vderiv_1_0 * derxjm_(2, 0, 1, ui) +
                  vderiv_1_1 * derxjm_(2, 1, 1, ui) + vderiv_1_2 * derxjm_(2, 2, 1, ui));
    }
  }

  return;
}


// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8, DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4, DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex20, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex27, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet10, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge6, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge15, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::pyramid5, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad4, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad8, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad9, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri3, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri6, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs9, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs27, DRT::ELEMENTS::Fluid::none>;
