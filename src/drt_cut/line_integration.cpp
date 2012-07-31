#include "line_integration.H"
#include "base.H"
#include "base_boundarycell.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include<cmath>
#include<iostream>

/*---------------------------------------------------------------------------------------------------------------------*
*      Compute normal vector for the line. if normal(0)==0, the line need not be integrated (Divergence theorem)       *
*----------------------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<2,1> LineIntegration::compute_normal()
{
  LINALG::Matrix<2,1> normal;
  double dy = end_pts_(1,1)-end_pts_(1,0);
  double dx = -end_pts_(0,1)+end_pts_(0,0);
  double modd = sqrt(dx*dx+dy*dy);

  normal(0,0) = dy/modd;
  normal(1,0) = dx/modd;

  return normal;
}

/*-----------------------------------------------------------------------------*
*                performs integration over the given line                      *
*------------------------------------------------------------------------------*/
double LineIntegration::integrate_line()
{
  LINALG::Matrix<2,1> normal;
  normal = compute_normal();

  if (fabs(normal(0,0))<0.000000001)
    return 0.0;

      double inte = 0.0;

    //8 is the order of Gauss integration used in the line integration
    //since we integrate 6th order polynomial in volume, 8th order must be used for line
    DRT::UTILS::GaussIntegration gi( DRT::Element::line2, 8 );

    for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
    {
      const LINALG::Matrix<1,1> eta( iquad.Point() );
      double weight = iquad.Weight();
      LINALG::Matrix<2,1>normaltemp,actCoord;
      double drs=0.0;
      Transform<DRT::Element::line2>(end_pts_,eta(0,0),actCoord,normaltemp,drs);

      if(bcellInt_==false)  //integration over volumecell
      {
        double linein = base_func_line_int(actCoord, inte_num_,alpha_);
        inte = inte+weight*linein*drs;
      }
      else                 // integration over boundarycell
      {
        double linein=0.0;
        if(intType_=="x")
          linein = base_func_surfX(actCoord,inte_num_,alpha_);
        else if(intType_=="y")
          linein = base_func_surfY(actCoord,inte_num_,alpha_);
        else if(intType_=="z")
          linein = base_func_surfZ(actCoord,inte_num_,alpha_);
        else
          dserror("Integration type unspecified");
        inte = inte+weight*linein*drs;
      }
    }
    inte = inte*normal(0,0);

  return inte;
}

/*---------------------------------------------------------------------------------------------------------------------*
*     Transform the Gaussian point coordinates and weight from (-1,1) interval to actual coordinate of the lines       *
*----------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType celldistype>
void LineIntegration::Transform( const LINALG::Matrix<2,2> &xyze,
                const double              & eta,
                LINALG::Matrix<2,1>       & x_gp_lin,
                LINALG::Matrix<2,1>       & normal,
                double                    & drs)
{
    const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
    LINALG::Matrix<numnodes,1> funct;
    LINALG::Matrix<1,numnodes> deriv;
    LINALG::Matrix<1,1> metrictensor;

    DRT::UTILS::shape_function_1D( funct, eta, celldistype );
    DRT::UTILS::shape_function_1D_deriv1( deriv, eta, celldistype );
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<celldistype>( xyze, deriv, metrictensor, drs, &normal );

    x_gp_lin.Multiply(xyze, funct);

    return;
}
