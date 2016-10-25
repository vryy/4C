/*!-------------------------------------------------------------------------*
\file xfluid_functions.H

\brief Managing and evaluating of spatial functions for xfluid cases

\level 3

\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef XFLUID_FUNCTION_H
#define XFLUID_FUNCTION_H

namespace DRT
{
  class Discretization;

  namespace UTILS
  {

    /// special implementation level set cut utilizing xfluid
    class GerstenbergerForwardfacingStep : public Function
    {
    public:

      /// ctor
      GerstenbergerForwardfacingStep();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };


    /* Moves a 2D cylinder back and forth along a predetermined direction in a sinus like fashion,

       x_pos = midpoint_trajectory_ + direction_*(L/2)*sin(sin_coeff*t-PI*0.5)

       The cylinder reaches its turning point distance at the time t_dist = n*PI/2   n=0,1,2,3,...
                                           and maximum velocity at t_max  = n*PI/4   n=1,2,3,4,...
    */
    class MovingLevelSetCylinder : public Function
    {
    public:

      /// ctor
      MovingLevelSetCylinder(std::vector<double>* origin, double radius,
          std::vector<double>* direction, double distance, double maxspeed);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    private:

      // Origin of the geometry
      std::vector<double> origin_;

      // Radius
      double radius_;

      // Orientation of the geometry (symmetry axis)
      std::vector<double> direction_;

      // Midpoint of trajectory
      std::vector<double> midpoint_trajectory_;

      // Distance traveled
      double distance_;

      // Maximum speed
      double maxspeed_;
    };


    /* Stationary Taylor-Couette flow with Navier-Slip type boundary condition at inner cylinder.

    */
    class TaylorCouetteFlow : public Function
    {
    public:

      /// ctor (for NavSlip at both boundaries)
      TaylorCouetteFlow(
          double radius_inner,
          double radius_outer,
          double vel_theta_inner,
          double vel_theta_outer,
          double epsilon_inner,
          double epsilon_outer,
          double traction_theta_inner,
          double traction_theta_outer,
          double viscosity);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

      std::vector<double> FctDer(int index, const double* x, const double t, DRT::Discretization* dis);

    private:

      // Origin of the geometry
//      std::vector<double> origin_; // Assume (0,0,0)

      // Radius
      double radius_inner_;
      double radius_outer_;

      // std::vector<double> prescribed_zero_pressure_point;

      // Traction prescribed (might be added at some point)

      double c1_;
      double c2_;
      double c3_;

    };


    /* SlipLengthManipulator. Useful for testing different slip-length functions.
    */
    class SlipLengthLevelSetManipulator : public Function
    {
    public:

      /// ctor
      SlipLengthLevelSetManipulator();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };


  }
}

#endif
