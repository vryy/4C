/*!-------------------------------------------------------------------------*
\file xfluid_functions.H

\brief Managing and evaluating of spatial functions for xfluid cases

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
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

       x_pos = maxspeed_ * sin(sin_coeff* t) * direction_ + origin_

       The cylinder reaches its destined distance at the time t_dist = 2.0*maxspeed_/(PI*distance_)
                                      and maximum velocity at t_max  = PI/2.
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

      // Distance traveled
      double maxspeed_;
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
