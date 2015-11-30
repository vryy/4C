/*!-----------------------------------------------------------------------------------------------*
 \file xfluid_functions.cpp

 \brief Managing and evaluating of spatial functions for Xfluid problems

  detailed description in header file combust_interface.H

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "xfluid_functions.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/standardtypes_cpp.H"


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 04/10 |
 *----------------------------------------------------------------------*/
DRT::UTILS::GerstenbergerForwardfacingStep::GerstenbergerForwardfacingStep() :
Function()
{
}

/*----------------------------------------------------------------------*
 | evaluation of xfluid level set test case             winter    09/14 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::GerstenbergerForwardfacingStep::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  //  //cube_Gerstenberger:
  //  //1.6x1.6
  //  //=========================================================

  double xp_corner[2];
  double distance = 0.0;

  xp_corner[0]=0.75;
  xp_corner[1]=0.75;

  double xtmp=-(xp_corner[0]-xp[0]);
  double ytmp=-(xp[1]-xp_corner[1]);

  if(xtmp<ytmp)
    distance = -xtmp;
  else
    distance = -ytmp;

  //  //1.6x1.6 (collapsing water column from the right)
  //  //=========================================================
  //  radius = 0.375; //0.375
  //
  //  xp_corner[0]=0.75;
  //  xp_corner[1]=0.75;
  //
  ////  double totallength_x=1.6;
  //
  ////  xp_center[0]= (xp_corner[0]+radius) + totallength_x;
  //  xp_center[0]= xp_corner[0]+radius;
  //  xp_center[1]= xp_corner[1]-radius;
  //
  ////  xp_center[0]=xp_corner[0]-radius;
  ////  xp_center[1]=xp_corner[1]-radius;
  //
  //  if (xp[0] >=xp_center[0] and xp[1] >= xp_center[1])
  //  {
  //     distance = xp[1]-xp_corner[1];
  //  }
  //  else if (xp[0] <=xp_center[0] and xp[1] <= xp_center[1] and !(xp[0]==xp_center[0] and xp[1]==xp_center[1]))
  //  {
  //      distance= -(xp[0]-xp_corner[0]);
  //  }
  //  else if (xp[0] > xp_center[0] and xp[1] < xp_center[1])
  //  {
  //      if(xp[1]>(xp_corner[1]+(xp[0]-xp_corner[0])))
  //      {
  //          distance = - fabs(xp_corner[1] - xp[1]);
  //      }
  //      else
  //      {
  //          distance = - fabs(xp_corner[0] - xp[0]);
  //      }
  //  }
  //  else
  //  {
  //      distance = sqrt(DSQR(xp[0]-xp_center[0])+DSQR(xp[1]-xp_center[1]))-radius;
  //  }
    //=========================================================

  return distance;
}

/*----------------------------------------------------------------------*
 | constructor                                          winter    11/15 |
 *----------------------------------------------------------------------*/
DRT::UTILS::SlipLengthLevelSetManipulator::SlipLengthLevelSetManipulator() :
Function()
{
}

/*----------------------------------------------------------------------*
 | evaluation of xfluid level set test case             winter    09/14 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::SlipLengthLevelSetManipulator::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{

  double x_cut = 0.4;
  double scalar;
//  if(xp[0]<x_cut)
//    scalar = 0.0;
//  else
//  {
//    scalar = 1000*(std::pow(1.6,xp[0]-x_cut));
//  }

    if(xp[0]>x_cut)
      scalar = 0.0;
    else
      scalar = 1000*(std::pow(1.6,x_cut-xp[0]));

  return scalar;
}

/*----------------------------------------------------------------------*
 | constructor                                            winter 10/15  |
 *----------------------------------------------------------------------*/
DRT::UTILS::MovingLevelSetCylinder::MovingLevelSetCylinder(std::vector<double>* origin, double radius,
    std::vector<double>* direction, double distance, double maxspeed) :
Function()
{
  // Origin of the geometry
  origin_ = *origin;

  // Radius
  radius_=radius;

  // Orientation of the geometry (symmetry axis)
  direction_= *direction; //needs to be normalized!

  double direction_norm = sqrt(direction_[0]*direction_[0]+direction_[1]*direction_[1]+direction_[2]*direction_[2]);

  direction_[0] /= direction_norm;
  direction_[1] /= direction_norm;
  direction_[2] /= direction_norm;

  // Distance traveled
  distance_= distance;

  // Distance traveled
  maxspeed_= maxspeed;

  //Initialize vector
  midpoint_trajectory_ = origin_;

  //Midpoint of trajectory
  midpoint_trajectory_[0] = 0.5*distance_*direction_[0]+origin_[0];
  midpoint_trajectory_[1] = 0.5*distance_*direction_[1]+origin_[1];
  midpoint_trajectory_[2] = 0.5*distance_*direction_[2]+origin_[2];

  //std::cout << "T_halfcycle: " << (L/2)/maxspeed * PI; //i.e. -d/2 -> d/2

}

/*----------------------------------------------------------------------*
 | evaluation of xfluid level set test case             winter    10/15 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::MovingLevelSetCylinder::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{

  // d = L/2 * sin(f*t-PI/2)
  // v = L*f/2 * cos(f*(t-PI/2)) = maxspeed * cos( (maxspeed*2/L)*t-PI/2 )

  //maxspeed =  L*f/2 -> f = 2*maxspeed/L
  //T_half   = (L/2)/maxspeed * PI/2

  //coefficient for sinus.
  double sin_coeff = 2.0 * maxspeed_*(1.0/distance_);
  //distance of cylinder viewed from the midpoint
  double dist = 0.5*distance_*sin(sin_coeff*t-PI*0.5);

  double x0_t = midpoint_trajectory_[0] + direction_[0]*dist;
  double x1_t = midpoint_trajectory_[1] + direction_[1]*dist;
//  double x2_t = midpoint_trajectory_[2] + direction_[2]*dist;

  double lsvalue = sqrt( (xp[0]-x0_t)*(xp[0]-x0_t) + (xp[1]-x1_t)*(xp[1]-x1_t) )-radius_;

  return lsvalue;
}
