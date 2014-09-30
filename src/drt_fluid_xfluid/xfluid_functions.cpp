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
