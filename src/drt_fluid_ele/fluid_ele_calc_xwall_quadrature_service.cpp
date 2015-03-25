/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xwall_quadrature_service.cpp

\brief quadrature rules for xwall

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
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
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::PrepareGaussRule()
{
  cgp_ = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints( numgpnorm_*numgpplane_*numgpplane_) );
  //which is the wall-normal element direction?
  //calculate jacobian at element center
  my::is_higher_order_ele_=false;
  my::EvalShapeFuncAndDerivsAtEleCenter();
  my::is_higher_order_ele_=true;
  //test the three element directions:
  LINALG::Matrix<my::nsd_,1> lvec1(true);
  LINALG::Matrix<my::nsd_,1> lvec2(true);
  LINALG::Matrix<my::nsd_,1> lvec3(true);
  lvec1(0)=1.0;
  lvec2(1)=1.0;
  lvec3(2)=1.0;
  LINALG::Matrix<my::nsd_,1> normv1(true);
  LINALG::Matrix<my::nsd_,1> normv2(true);
  LINALG::Matrix<my::nsd_,1> normv3(true);
  normv1.Multiply(my::xji_,lvec1);
  normv2.Multiply(my::xji_,lvec2);
  normv3.Multiply(my::xji_,lvec3);

  LINALG::Matrix<my::nsd_,1> normwall(true);
  normwall.Multiply(derxy_,ewdist_);
  const double dot1=abs(normwall.Dot(normv1)/normwall.Norm2()/normv1.Norm2());
  const double dot2=abs(normwall.Dot(normv2)/normwall.Norm2()/normv2.Norm2());
  const double dot3=abs(normwall.Dot(normv3)/normwall.Norm2()/normv3.Norm2());

  // get the quad9 gaussrule for the in plane integration
  DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::quad8 ,2*numgpplane_-1);
  // get the quad9 gaussrule for the in normal integration
  DRT::UTILS::GaussIntegration intpointsnormal( DRT::Element::line3 ,2*numgpnorm_-1);

  //0.9 corresponds to an angle of 25.8 deg
  if(dot1<0.90&&dot2<0.90&&dot3<0.90)
  { //element, where the wall normal direction does not point in one specific element direction, e.g. in corners
    cgp_->IncreaseReserved((numgpnorm_*numgpnorm_*numgpnorm_)-(numgpnorm_*numgpplane_*numgpplane_) );
    DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::quad8 ,2*numgpnorm_-1);
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadnorm.Point()[0],iquadplane.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else if(dot1>dot2&&dot1>dot3)
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadnorm.Point()[0],iquadplane.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else if(dot2>dot3)
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadplane.Point()[0],iquadnorm.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadplane.Point()[0],iquadplane.Point()[1],iquadnorm.Point()[0],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  DRT::UTILS::GaussIntegration grule(cgp_);
  my::intpoints_=grule;

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Sysmat(
    const LINALG::Matrix<my::nsd_,my::nen_>&              ebofoaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eprescpgaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              ebofon,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eprescpgn,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evelaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eveln,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evelnp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              fsevelaf,
    const LINALG::Matrix<my::nen_,1>&                 fsescaaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evel_hat,
    const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&         ereynoldsstress_hat,
    const LINALG::Matrix<my::nen_,1>&                 epreaf,
    const LINALG::Matrix<my::nen_,1>&                 epren,
    const LINALG::Matrix<my::nen_,1>&                 eprenp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eaccam,
    const LINALG::Matrix<my::nen_,1>&                 escaaf,
    const LINALG::Matrix<my::nen_,1>&                 escaam,
    const LINALG::Matrix<my::nen_,1>&                 escadtam,
    const LINALG::Matrix<my::nen_,1>&                 escabofoaf,
    const LINALG::Matrix<my::nen_,1>&                 escabofon,
    const LINALG::Matrix<my::nsd_,my::nen_>&              emhist,
    const LINALG::Matrix<my::nsd_,my::nen_>&              edispnp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              egridv,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&              eforce,
    const LINALG::Matrix<my::nen_,1> &                eporo,
    const LINALG::Matrix<my::nsd_,2*my::nen_> &             egradphi,
    const LINALG::Matrix<my::nen_,2*1> &                ecurvature,
    const double                                  thermpressaf,
    const double                                  thermpressam,
    const double                                  thermpressdtaf,
    const double                                  thermpressdtam,
    Teuchos::RCP<const MAT::Material>             material,
    double&                                       Cs_delta_sq,
    double&                                       Ci_delta_sq,
    double&                                       Cv,
    bool                                          isale,
    double * saccn,
    double * sveln,
    double * svelnp,
    const DRT::UTILS::GaussIntegration & intpoints
    )
{

  if(quadraturetol_>0.0 && quadraturetol_<1.0)
  {
    //this saves a lot of calculations
    my::is_higher_order_ele_=false;

    //dummy
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>  estifdummy(true);
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>  exactestifdummy(true);
    my::EvalShapeFuncAndDerivsAtEleCenter();

    //test the three element directions:
    LINALG::Matrix<my::nsd_,1> lvec1(true);
    LINALG::Matrix<my::nsd_,1> lvec2(true);
    LINALG::Matrix<my::nsd_,1> lvec3(true);
    lvec1(0)=1.0;
    lvec2(1)=1.0;
    lvec3(2)=1.0;
    LINALG::Matrix<my::nsd_,1> normv1(true);
    LINALG::Matrix<my::nsd_,1> normv2(true);
    LINALG::Matrix<my::nsd_,1> normv3(true);
    normv1.Multiply(my::xji_,lvec1);
    normv2.Multiply(my::xji_,lvec2);
    normv3.Multiply(my::xji_,lvec3);

    LINALG::Matrix<my::nsd_,1> normwall(true);
    normwall.Multiply(derxy_,ewdist_);
    double dot1=abs(normwall.Dot(normv1)/normwall.Norm2()/normv1.Norm2());
    double dot2=abs(normwall.Dot(normv2)/normwall.Norm2()/normv2.Norm2());
    double dot3=abs(normwall.Dot(normv3)/normwall.Norm2()/normv3.Norm2());

    //test in wall-normal direction
    LINALG::Matrix<(my::nsd_+1)*my::nen_,1> newtestforce(true);
    LINALG::Matrix<(my::nsd_+1)*my::nen_,1> oldtestforce(true);

    {
      int gp=31;

      SysmatForErrorEstimation(ebofoaf,
             eprescpgaf,
             ebofon,
             eprescpgn,
             evelaf,
             eveln,
             evelnp,
             fsevelaf,
             fsescaaf,
             evel_hat,
             ereynoldsstress_hat,
             epreaf,
             epren,
             eprenp,
             eaccam,
             escaaf,
             escaam,
             escadtam,
             escabofoaf,
             escabofon,
             emhist,
             edispnp,
             egridv,
             exactestifdummy,
             emesh,  // -> emesh
             oldtestforce,
             eporo,
             egradphi,
             ecurvature,
             thermpressaf,
             thermpressam,
             thermpressdtaf,
             thermpressdtam,
             material,
             Cs_delta_sq,
             Ci_delta_sq,
             Cv,
             isale,
             saccn,
             sveln,
             svelnp,
             intpoints,
             dot1,
             dot2,
             dot3,
             gp);
    }

    double err=10.0;
    int gp=7;

    while((err>quadraturetol_)&&gp<31)
    {

      gp++;

      SysmatForErrorEstimation(ebofoaf,
             eprescpgaf,
             ebofon,
             eprescpgn,
             evelaf,
             eveln,
             evelnp,
             fsevelaf,
             fsescaaf,
             evel_hat,
             ereynoldsstress_hat,
             epreaf,
             epren,
             eprenp,
             eaccam,
             escaaf,
             escaam,
             escadtam,
             escabofoaf,
             escabofon,
             emhist,
             edispnp,
             egridv,
             estifdummy,
             emesh,  // -> emesh
             oldtestforce,
             eporo,
             egradphi,
             ecurvature,
             thermpressaf,
             thermpressam,
             thermpressdtaf,
             thermpressdtam,
             material,
             Cs_delta_sq,
             Ci_delta_sq,
             Cv,
             isale,
             saccn,
             sveln,
             svelnp,
             intpoints,
             dot1,
             dot2,
             dot3,
             gp);
      err=0.0;
      double count=0.0;
      for(int idof=0 ; (my::nsd_+1)*my::nen_!=idof; idof++)
        for(int jdof=0 ; (my::nsd_+1)*my::nen_!=jdof; jdof++)
        {
          if(abs(exactestifdummy(idof,jdof))>1.0e-12)
          {
            err+=abs(estifdummy(idof,jdof)-exactestifdummy(idof,jdof))/abs(exactestifdummy(idof,jdof));
            count++;
//               if(err<newerr)
//                 err=newerr;
          }
          estifdummy(idof,jdof)=0.0;
        }
      err/=count;


    }
    numgpnorm_=gp;

    //for the in-plane direction, use the direction of the flow
    my::EvalShapeFuncAndDerivsAtEleCenter();
    //normwall.Multiply(evelaf,my::funct_);
    //use the direction of the highest gradient in tauw
    //this gradient is wall-parallel because tauw is constant in wall-normal direction
    normwall.Multiply(derxy_,etauw_);
    //if there is no gradient, e.g. if aggregated, use streamwise direction
    if(normwall.Norm2()<1.0e-8)
      normwall.Multiply(evelaf,my::funct_);
    numgpplane_=15;
    if(normwall.Norm2()>1.0e-8)
    {
      dot1=abs(normwall.Dot(normv1)/normwall.Norm2()/normv1.Norm2());
      dot2=abs(normwall.Dot(normv2)/normwall.Norm2()/normv2.Norm2());
      dot3=abs(normwall.Dot(normv3)/normwall.Norm2()/normv3.Norm2());
      for(int idof=0 ; (my::nsd_+1)*my::nen_!=idof; idof++)
        for(int jdof=0 ; (my::nsd_+1)*my::nen_!=jdof; jdof++)
        {
          exactestifdummy(idof,jdof)=0.0;
        }

      {
        int gp=31;

        SysmatForErrorEstimation(ebofoaf,
               eprescpgaf,
               ebofon,
               eprescpgn,
               evelaf,
               eveln,
               evelnp,
               fsevelaf,
               fsescaaf,
               evel_hat,
               ereynoldsstress_hat,
               epreaf,
               epren,
               eprenp,
               eaccam,
               escaaf,
               escaam,
               escadtam,
               escabofoaf,
               escabofon,
               emhist,
               edispnp,
               egridv,
               exactestifdummy,
               emesh,  // -> emesh
               oldtestforce,
               eporo,
               egradphi,
               ecurvature,
               thermpressaf,
               thermpressam,
               thermpressdtaf,
               thermpressdtam,
               material,
               Cs_delta_sq,
               Ci_delta_sq,
               Cv,
               isale,
               saccn,
               sveln,
               svelnp,
               intpoints,
               dot1,
               dot2,
               dot3,
               gp);
      }

      err=10.0;
      gp=3;
      while((err>quadraturetol_)&&gp<31)
      {
        gp++;

        SysmatForErrorEstimation(ebofoaf,
               eprescpgaf,
               ebofon,
               eprescpgn,
               evelaf,
               eveln,
               evelnp,
               fsevelaf,
               fsescaaf,
               evel_hat,
               ereynoldsstress_hat,
               epreaf,
               epren,
               eprenp,
               eaccam,
               escaaf,
               escaam,
               escadtam,
               escabofoaf,
               escabofon,
               emhist,
               edispnp,
               egridv,
               estifdummy,
               emesh,  // -> emesh
               oldtestforce,
               eporo,
               egradphi,
               ecurvature,
               thermpressaf,
               thermpressam,
               thermpressdtaf,
               thermpressdtam,
               material,
               Cs_delta_sq,
               Ci_delta_sq,
               Cv,
               isale,
               saccn,
               sveln,
               svelnp,
               intpoints,
               dot1,
               dot2,
               dot3,
               gp);
        err=0.0;
        double count=0.0;
        for(int idof=0 ; (my::nsd_+1)*my::nen_!=idof; idof++)
          for(int jdof=0 ; (my::nsd_+1)*my::nen_!=jdof; jdof++)
          {
            if(abs(exactestifdummy(idof,jdof))>1.0e-12)
            {
              err+=abs(estifdummy(idof,jdof)-exactestifdummy(idof,jdof))/abs(exactestifdummy(idof,jdof));
              count++;
//               if(err<newerr)
//                 err=newerr;
            }
            estifdummy(idof,jdof)=0.0;
          }
        err/=count;
        //std::cout << "order  " << gp << "  err  " << err << std::endl;
      }
      numgpplane_=gp;
    }
//int cost=0;
//for (int i=3;i<=numgpnorm_; i++)
//  cost+=i;
//for (int i=2;i<=numgpplane_; i++)
//  cost+=i;
//cost+=62;
//std::cout << "evaluation cost:  " << cost << "/" << numgpnorm_*numgpplane_*numgpplane_ << std::endl;
    cgp_ = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints(numgpnorm_*numgpplane_*numgpplane_) );
    // get the quad9 gaussrule for the in plane integration
    DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::quad8 ,2*numgpplane_-1);
    // get the quad9 gaussrule for the in normal integration
    DRT::UTILS::GaussIntegration intpointsnormal( DRT::Element::line3 ,2*numgpnorm_-1);
    my::EvalShapeFuncAndDerivsAtEleCenter();
    normwall.Multiply(derxy_,ewdist_);
    dot1=abs(normwall.Dot(normv1)/normwall.Norm2()/normv1.Norm2());
    dot2=abs(normwall.Dot(normv2)/normwall.Norm2()/normv2.Norm2());
    dot3=abs(normwall.Dot(normv3)/normwall.Norm2()/normv3.Norm2());

    if(dot1>dot2&&dot1>dot3)
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
      {
        // start loop over integration points in layer
        for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
        {
          cgp_->Append(iquadnorm.Point()[0],iquadplane.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
        }
      }
    }
    else if(dot2>dot3)
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
      {
        // start loop over integration points in layer
        for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
        {
          cgp_->Append(iquadplane.Point()[0],iquadnorm.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
        }
      }
    }
    else
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
      {
        // start loop over integration points in layer
        for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
        {
          cgp_->Append(iquadplane.Point()[0],iquadplane.Point()[1],iquadnorm.Point()[0],iquadplane.Weight()*iquadnorm.Weight());
        }
      }
    }

    DRT::UTILS::GaussIntegration intpointsnew(cgp_);

    my::is_higher_order_ele_=true;

    my::Sysmat(ebofoaf,
           eprescpgaf,
           ebofon,
           eprescpgn,
           evelaf,
           eveln,
           evelnp,
           fsevelaf,
           fsescaaf,
           evel_hat,
           ereynoldsstress_hat,
           epreaf,
           epren,
           eprenp,
           eaccam,
           escaaf,
           escaam,
           escadtam,
           escabofoaf,
           escabofon,
           emhist,
           edispnp,
           egridv,
           estif,
           emesh,  // -> emesh
           eforce,
           eporo,
           egradphi,
           ecurvature,
           thermpressaf,
           thermpressam,
           thermpressdtaf,
           thermpressdtam,
           material,
           Cs_delta_sq,
           Ci_delta_sq,
           Cv,
           isale,
           saccn,
           sveln,
           svelnp,
           intpointsnew);
  }
  else
  {

    my::Sysmat(ebofoaf,
           eprescpgaf,
           ebofon,
           eprescpgn,
           evelaf,
           eveln,
           evelnp,
           fsevelaf,
           fsescaaf,
           evel_hat,
           ereynoldsstress_hat,
           epreaf,
           epren,
           eprenp,
           eaccam,
           escaaf,
           escaam,
           escadtam,
           escabofoaf,
           escabofon,
           emhist,
           edispnp,
           egridv,
           estif,
           emesh,  // -> emesh
           eforce,
           eporo,
           egradphi,
           ecurvature,
           thermpressaf,
           thermpressam,
           thermpressdtaf,
           thermpressdtam,
           material,
           Cs_delta_sq,
           Ci_delta_sq,
           Cv,
           isale,
           saccn,
           sveln,
           svelnp,
           intpoints);
    return;
  }
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::SysmatForErrorEstimation(
    const LINALG::Matrix<my::nsd_,my::nen_>&              ebofoaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eprescpgaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              ebofon,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eprescpgn,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evelaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eveln,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evelnp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              fsevelaf,
    const LINALG::Matrix<my::nen_,1>&                 fsescaaf,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evel_hat,
    const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&         ereynoldsstress_hat,
    const LINALG::Matrix<my::nen_,1>&                 epreaf,
    const LINALG::Matrix<my::nen_,1>&                 epren,
    const LINALG::Matrix<my::nen_,1>&                 eprenp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              eaccam,
    const LINALG::Matrix<my::nen_,1>&                 escaaf,
    const LINALG::Matrix<my::nen_,1>&                 escaam,
    const LINALG::Matrix<my::nen_,1>&                 escadtam,
    const LINALG::Matrix<my::nen_,1>&                 escabofoaf,
    const LINALG::Matrix<my::nen_,1>&                 escabofon,
    const LINALG::Matrix<my::nsd_,my::nen_>&              emhist,
    const LINALG::Matrix<my::nsd_,my::nen_>&              edispnp,
    const LINALG::Matrix<my::nsd_,my::nen_>&              egridv,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
    LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&              eforce,
    const LINALG::Matrix<my::nen_,1> &                eporo,
    const LINALG::Matrix<my::nsd_,2*my::nen_> &             egradphi,
    const LINALG::Matrix<my::nen_,2*1> &                ecurvature,
    const double                                  thermpressaf,
    const double                                  thermpressam,
    const double                                  thermpressdtaf,
    const double                                  thermpressdtam,
    Teuchos::RCP<const MAT::Material>             material,
    double&                                       Cs_delta_sq,
    double&                                       Ci_delta_sq,
    double&                                       Cv,
    bool                                          isale,
    double * saccn,
    double * sveln,
    double * svelnp,
    const DRT::UTILS::GaussIntegration & intpoints,
    double dot1,
    double dot2,
    double dot3,
    int gp
    )
{
cgp_ = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints(gp) );
// get the quad9 gaussrule for the in normal integration
DRT::UTILS::GaussIntegration intpointsnormal( DRT::Element::line3 ,2*gp-1);
if(dot1>dot2&&dot1>dot3)
{
  // start loop over integration points in layer
  for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
  {
    cgp_->Append(iquadnorm.Point()[0],0.0,0.0,iquadnorm.Weight());
  }
}
else if(dot2>dot3)
{
  // start loop over integration points in layer
  for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
  {
    cgp_->Append(0.0,iquadnorm.Point()[0],0.0,iquadnorm.Weight());
  }
}
else
{
  // start loop over integration points in layer
  for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
  {
    cgp_->Append(0.0,0.0,iquadnorm.Point()[0],iquadnorm.Weight());
  }
}
DRT::UTILS::GaussIntegration grule(cgp_);
my::Sysmat(ebofoaf,
       eprescpgaf,
       ebofon,
       eprescpgn,
       evelaf,
       eveln,
       evelnp,
       fsevelaf,
       fsescaaf,
       evel_hat,
       ereynoldsstress_hat,
       epreaf,
       epren,
       eprenp,
       eaccam,
       escaaf,
       escaam,
       escadtam,
       escabofoaf,
       escabofon,
       emhist,
       edispnp,
       egridv,
       estif,
       emesh,  // -> emesh
       eforce,
       eporo,
       egradphi,
       ecurvature,
       thermpressaf,
       thermpressam,
       thermpressdtaf,
       thermpressdtam,
       material,
       Cs_delta_sq,
       Ci_delta_sq,
       Cv,
       isale,
       saccn,
       sveln,
       svelnp,
       grule);
  return;
}

template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>;
