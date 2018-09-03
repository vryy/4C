/*!------------------------------------------------------------------------------------------------*
\file fluid_ele_calc_topopt_service.cpp

\brief

\level 2

<pre>
\maintainer Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"

#include "../drt_inpar/inpar_topopt.H"


/*----------------------------------------------------------------------*
 |  compute porosity = reacoeff at gausspoint          winklmaier 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::GetPorosityAtGP(
    const LINALG::Matrix<nen_, 1>& eporo)
{
  // the eporo actually has density values here
  double densint = funct_.Dot(eporo);
  const double* params = fldpara_->TopoptParams();

  if (params[2] > -1.0e-15)  // >=0 as it should be -> standard case
  {
    reacoeff_ =
        params[1] + (params[0] - params[1]) * densint * (1 + params[2]) / (densint + params[2]);
  }
  else  // special cases
  {
    INPAR::TOPOPT::OptiCase testcase = (INPAR::TOPOPT::OptiCase)round(-params[2]);

    switch (testcase)
    {
      case INPAR::TOPOPT::optitest_channel:
      {
        LINALG::Matrix<nsd_, 1> gp(true);
        gp.Multiply(xyze_, funct_);

        if (gp(1) < -0.2 || gp(1) > 0.2)  // wall area
          reacoeff_ = params[1];
        else
          reacoeff_ = params[0];
        break;
      }
      case INPAR::TOPOPT::optitest_channel_with_step:
      {
        LINALG::Matrix<nsd_, 1> gp(true);
        gp.Multiply(xyze_, funct_);

        if ((gp(0) > 1.5 and gp(0) < 1.9) and (gp(1) < 0.4))  // step -> wall
          reacoeff_ = params[1];
        else
          reacoeff_ = params[0];
        break;
      }
      case INPAR::TOPOPT::optitest_cornerflow:
      {
        LINALG::Matrix<nsd_, 1> gp(true);
        gp.Multiply(xyze_, funct_);

        if ((gp(0) > 0.875) or (gp(1) > 0.875) or (gp(0) < 0.625 and gp(1) < 0.625))  // wall
          reacoeff_ = params[1];
        else
          reacoeff_ = params[0];
        break;
      }
      case INPAR::TOPOPT::optitest_lin_poro:
      {
        double diff = params[1] - params[0];
        double pmax = params[1];

        reacoeff_ = -diff * densint + pmax;
        break;
      }
      case INPAR::TOPOPT::optitest_quad_poro:
      {
        double diff = params[1] - params[0];
        double pmax = params[1];

        double k = 0.1;

        reacoeff_ = (diff - k * pmax) * densint * densint + (-2 * diff + k * pmax) * densint + pmax;
        break;
      }
      case INPAR::TOPOPT::optitest_cub_poro:
      {
        double diff = params[1] - params[0];
        double pmax = params[1];

        double k1 = -50000.0;
        double k2 = 0.1;

        reacoeff_ = (2 * diff + k1 - k2 * pmax) * densint * densint * densint +
                    (-3 * diff - 2 * k1 + k2 * pmax) * densint * densint + k1 * densint + pmax;
        break;
      }
      default:
      {
        dserror("you should not be here with a testcase not handled above");
        break;
      }
    }
  }
}


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
