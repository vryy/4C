/*----------------------------------------------------------------------*/
/*!
\file xfluid_functions.H

\brief Managing and evaluating of spatial functions for Xfluid problems

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"
//#include "../drt_lib/drt_linedefinition.H"

#ifndef XFLUID_FUNCTION_H
#define XFLUID_FUNCTION_H

namespace DRT
{
  class DiscretizationInterface;

  namespace INPUT
  {
    class LineDefinition;
  }

  namespace UTILS
  {
    void XfluidValidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);
    bool XfluidFunctionHaveNamed(Teuchos::RCP<DRT::INPUT::LineDefinition> functionnamed,
        std::vector<Teuchos::RCP<Function>>* function_);

    /// special implementation level set cut utilizing xfluid
    class GerstenbergerForwardfacingStep : public Function
    {
     public:
      /// ctor
      GerstenbergerForwardfacingStep();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
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
      double Evaluate(int index, const double* x, double t);

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

    /* Moves a 3D torus back and forth along a predetermined direction in a sine like fashion with
      constant rotation,

       +++ LEVEL SET VALUE +++

       TORUS
       The torus is oriented in the orientationvec_torus direction and its midpoint is at origin
      (for t=0). The radius of the torus is radius, where as the radius for the tube is radius_tube.

       TRANSLATION
       x_pos = midpoint_trajectory_ + direction_*(L/2)*sin(sin_coeff*t-PI*0.5)

       The cylinder reaches its turning point distance at the time t_dist = n*PI/2   n=0,1,2,3,...
                                           and maximum velocity at t_max  = n*PI/4   n=1,2,3,4,...

       ROTATION
       The torus rotates around a given rotation vector rotvector, at the constant rotation speed
      rotspeed.

          R(theta) = cos(theta) I + sin(theta) * [w]_x +(1-cos(theta))*w(x)w
          theta = rotspeed*t

          x_rot_pos = R(j,i)*( x(j) - x_pos(j) )   (MAYBE TRANSPOSE R(j,i)!!!!)

       LEVEL SET

          x_orto = R(i,j) * ( x(j) - x_pos(j) ) * orientationvec_torus(i)

          x_base = x_rot_pos(i) - x_orto * orientationvec_torus(i)

          phi(theta,t) = sqrt( ( radius - sqrt( x_base*x_base ) )^2 + x_orto*x_orto ) - radius_tube

       #########################################################

       +++ VELOCITY +++

       TRANSLATION

          v_trans(i) = direction_(i)*sin_coeff(L/2)*cos(sin_coeff*t-PI*0.5)

       ROTATION

          v_rot(i) = rotspeed * [w]_x (i,j) * (x(j) - x_pos(j))   %(or should it be x_rot_pos?

       TOTAL

          v_tot = v_trans + v_rot

       #########################################################

      +++ SLIPLENGTH +++

      0)
      1) Increasing sliplength in shape of sphere.

    */
    class MovingLSTorus : public Function
    {
     public:
      /// ctor
      MovingLSTorus(std::vector<double>* origin, std::vector<double>* orientationvec_torus,
          double radius, double radius_tube, std::vector<double>* direction, double distance,
          double maxspeed, std::vector<double>* rotvector, double rotspeed, double rotramptime);

     protected:
      // Origin of the geometry
      std::vector<double> origin_;

      // Orientation vector
      std::vector<double> orientationvec_torus_;

      // Radius
      double radius_;

      // Radius of torus tube
      double radius_tube_;

      // Orientation of the geometry (symmetry axis)
      std::vector<double> direction_;

      // Midpoint of trajectory
      std::vector<double> midpoint_trajectory_;

      // Distance traveled
      double distance_;

      // Maximum speed
      double maxspeed_;

      // Rotation vector
      std::vector<double> rotvector_;

      // Rotation speed (revolutions per second)
      double rotspeed_;

      // Ramp time for the rotation ( ramping done as: 0.5*cos(PI*t/T_ramp) )
      double ramptime_;

      // MATRICES!!
      std::vector<std::vector<double>> eye_;
      std::vector<std::vector<double>> rot_joint_;
      std::vector<std::vector<double>> rot_cross_;
    };

    class MovingLevelSetTorus : public MovingLSTorus
    {
     public:
      /// ctor
      MovingLevelSetTorus(std::vector<double>* origin, std::vector<double>* orientationvec_torus,
          double radius, double radius_tube, std::vector<double>* direction, double distance,
          double maxspeed, std::vector<double>* rotvector, double rotspeed, double rotramptime);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    class MovingLevelSetTorusVelocity : public MovingLSTorus
    {
     public:
      /// ctor
      MovingLevelSetTorusVelocity(std::vector<double>* origin,
          std::vector<double>* orientationvec_torus, double radius, double radius_tube,
          std::vector<double>* direction, double distance, double maxspeed,
          std::vector<double>* rotvector, double rotspeed, double rotramptime);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

      /*
      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w)

      */
      virtual int NumberComponents() { return (3); };
    };

    class MovingLevelSetTorusSliplength : public MovingLSTorus
    {
     public:
      /// ctor
      MovingLevelSetTorusSliplength(std::vector<double>* origin,
          std::vector<double>* orientationvec_torus, double radius, double radius_tube,
          std::vector<double>* direction, double distance, double maxspeed,
          std::vector<double>* rotvector, double rotspeed, double rotramptime, int slipfunct);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

     private:
      int slipfunct_;
    };

    /* Stationary Taylor-Couette flow with Navier-Slip type boundary condition at inner cylinder.

    */
    class TaylorCouetteFlow : public Function
    {
     public:
      /// ctor (for NavSlip at both boundaries)
      TaylorCouetteFlow(double radius_inner, double radius_outer, double vel_theta_inner,
          double vel_theta_outer, double epsilon_inner, double epsilon_outer,
          double traction_theta_inner, double traction_theta_outer, double viscosity);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

      std::vector<double> EvaluateSpatialDerivative(int index, const double* x, const double t);

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

    /* Urquiza-box flow
     * Oseen test case (can be solved with Nav-Stokes as well as the advective velocity
     * is chosen as the analytic solution to the Nav-Stokes equations).
     */
    class UrquizaBoxFlow : public Function
    {
     public:
      /// ctor
      UrquizaBoxFlow(double lengthx, double lengthy, double rotation, double viscosity,
          double density, int functno, std::vector<double> lincomb);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

      // Derivative
      std::vector<double> EvaluateSpatialDerivative(int index, const double* x, const double t);

     protected:
      double lengthx_;   // Length lx -> [-lx,lx]
      double lengthy_;   // Length ly ->[-ly,ly]
      double rotation_;  // rot    [0,2*pi]
      std::vector<std::vector<double>> rotvector_;
      double kinvisc_;  // kinvisc = visc/rho
      int functno_;     // function number

      double c1_;
      double c2_;
    };


    /* Urquiza-box test case
     */
    class UrquizaBoxFlowForce : public UrquizaBoxFlow  // public Function
    {
     public:
      /// ctor
      UrquizaBoxFlowForce(double lengthx, double lengthy, double rotation, double viscosity,
          double density, int functno, std::vector<double> lincomb);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

      std::vector<double> EvaluateSpatialDerivative(int index, const double* x, const double t)
      {
        dserror("Derivative not implemented for UrquizaBoxFlowForce");
        return std::vector<double>();
      }
    };

    /* SlipLengthManipulator. Useful for testing different slip-length functions.
     */
    class UrquizaBoxFlowTraction : public UrquizaBoxFlow  // public Function
    {
     public:
      /// ctor
      UrquizaBoxFlowTraction(double lengthx, double lengthy, double rotation, double viscosity,
          double density, int functno, std::vector<double> lincomb);

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);

      std::vector<double> EvaluateSpatialDerivative(int index, const double* x, const double t)
      {
        dserror("Derivative not implemented for UrquizaBoxFlowTraction");
        return std::vector<double>();
      }
    };

    /* SlipLengthManipulator. Useful for testing different slip-length functions.
     */
    class SlipLengthLevelSetManipulator : public Function
    {
     public:
      /// ctor
      SlipLengthLevelSetManipulator();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };


  }  // namespace UTILS
}  // namespace DRT

#endif
