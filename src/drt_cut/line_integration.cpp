#include "line_integration.H"
#include "base.H"
#include "base_boundarycell.H"

#include<cmath>
#include<iostream>


//Compute normal vector for the line
LINALG::Matrix<2,1> LineIntegration::compute_normal()
{
  LINALG::Matrix<2,1> normal;
  double dy = point_end_[1]-point_begin_[1];
  double dx = -point_end_[0]+point_begin_[0];
  double modd = sqrt(dx*dx+dy*dy);

  normal(0,0) = dy/modd;
  normal(1,0) = dx/modd;

  return normal;
}

//returns the weights of 8-point gaussian quadrture in (-1,1) interval
std::vector<double> LineIntegration::get_Gauss_weights()
{
  std::vector<double> line_wei(8);
  line_wei[0] = 0.1012285362903762591525314;
  line_wei[1] = 0.2223810344533744705443560;
  line_wei[2] = 0.3137066458778872873379622;
  line_wei[3] = 0.3626837833783619829651504;
  line_wei[4] = 0.3626837833783619829651504;
  line_wei[5] = 0.3137066458778872873379622;
  line_wei[6] = 0.2223810344533744705443560;
  line_wei[7] = 0.1012285362903762591525314;

  return line_wei;
}

//returns the positions of 8-point gaussian quadrture in (-1,1) interval
std::vector<double> LineIntegration::get_Gauss_line_pts()
{
  std::vector<double> line_tau(8);
  line_tau[0] = -0.9602898564975362316835609;
  line_tau[1] = -0.7966664774136267395915539;
  line_tau[2] = -0.5255324099163289858177390;
  line_tau[3] = -0.1834346424956498049394761;
  line_tau[4] = 0.1834346424956498049394761;
  line_tau[5] = 0.5255324099163289858177390;
  line_tau[6] = 0.7966664774136267395915539;
  line_tau[7] = 0.9602898564975362316835609;

  return line_tau;
}

//Obtain the actual integration points from the points available in (-1,1) interval
LINALG::Matrix<2,8> LineIntegration::find_line_integration_pts()
{
  std::vector<double> line_tau = get_Gauss_line_pts();

  LINALG::Matrix<2,8> line_int_pts;
  double xmid[2];

  //middle point in all 2 coordinates
  for(int i=0;i<2;i++)
  {
    xmid[i] = 0.5*(point_begin_[i]+point_end_[i]);
  }

  //Finding the 8 integration points used in Gaussian quadrature
  for(int i=0;i<8;i++)
  {
    double tau_fac = line_tau[i];
    for(int j=0;j<2;j++)
    {
      line_int_pts(j,i) = (xmid[j]-point_begin_[j])*tau_fac+xmid[j];
    }
  }

  return line_int_pts;
}

//The half length of the line appears in the Gaussian quadrature when the integration
//is scaled from (a,b) to (-1,1) for which the weights are available 
double LineIntegration::half_length()
{
  double x1 = point_begin_[0];
  double y1 = point_begin_[1];

  double x2 = point_end_[0];
  double y2 = point_end_[1];

  double hal_len = 0.5*sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  return hal_len;
}

//performs integration over the given line
double LineIntegration::integrate_line()
{
  LINALG::Matrix<2,1> normal;
  normal = compute_normal();

  if (fabs(normal(0,0))<0.000000001)
    return 0.0;

  LINALG::Matrix<2,8> line_int_pts;
  line_int_pts = find_line_integration_pts();

  std::vector<double> line_wei = get_Gauss_weights();

  double half_len = half_length();

  double inte = 0.0;

  for (int i=0;i<8;i++)
  {
    std::vector<double> pt(2);
    pt[0] = line_int_pts(0,i);
    pt[1] = line_int_pts(1,i);

    if(bcellInt_==false)
    {
      double linein = base_func_line_int(pt, inte_num_,alpha_);
      inte = inte+line_wei[i]*linein;
    }
    else
    {
      double linein=0.0;
      if(intType_=="x")
        linein = base_func_surfX(pt,inte_num_,alpha_);
      if(intType_=="y")
        linein = base_func_surfY(pt,inte_num_,alpha_);
      if(intType_=="z")
        linein = base_func_surfZ(pt,inte_num_,alpha_);
      inte = inte+line_wei[i]*linein;
    }

  }
  inte = inte*normal(0,0)*half_len;

  return inte;
}
