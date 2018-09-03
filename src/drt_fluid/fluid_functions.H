/*!-------------------------------------------------------------------------*
\file fluid_functions.H

\brief Managing and evaluating of spatial functions for fluid problems

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef FLUID_FUNCTION_H
#define FLUID_FUNCTION_H

namespace DRT
{
  class Discretization;
  namespace INPUT
  {
    class LineDefinition;
  }  // namespace INPUT
}  // namespace DRT


namespace FLD
{
  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public DRT::UTILS::Function
  {
   public:
    BeltramiUP(int mat_id);

    BeltramiUP(Teuchos::RCP<MAT::Material>& mat);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    virtual double Evaluate(const int index, const double* x, double t);

    /*!

      \brief evaluate first or second time derivative

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the time derivative will
                       be evaluated
      \param t     (i) The point in time in which the time derivative will
                       be evaluated
      \param deg   (i) The degree of the time derivative

      */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents() { return (4); };

   private:
    double density_;
    double kinviscosity_;
  };



  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public DRT::UTILS::Function
  {
   public:
    BeltramiGradU(int mat_id);

    BeltramiGradU(Teuchos::RCP<MAT::Material>& mat);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t);

    /*!

      \brief evaluate first or second time derivative

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the time derivative will
                       be evaluated
      \param t     (i) The point in time in which the time derivative will
                       be evaluated
      \param deg   (i) The degree of the time derivative

      */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents() { return (9); };

   private:
    double kinviscosity_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity,
  /// pressure)
  class KimMoinUP : public DRT::UTILS::Function
  {
   public:
    KimMoinUP(int mat_id, bool is_stationary);

    KimMoinUP(Teuchos::RCP<MAT::Material>& mat, bool is_stationary);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t);

    /*!

      \brief evaluate first or second time derivative

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the time derivative will
                       be evaluated
      \param t     (i) The point in time in which the time derivative will
                       be evaluated
      \param deg   (i) The degree of the time derivative

     */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents() { return (4); };

   private:
    double density_;
    double kinviscosity_;
    bool is_stationary_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of
  /// velocity)
  class KimMoinGradU : public DRT::UTILS::Function
  {
   public:
    KimMoinGradU(int mat_id, bool is_stationary);

    KimMoinGradU(Teuchos::RCP<MAT::Material>& mat, bool is_stationary);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated
      \param t     (i) The point in time in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t);

    /*!

      \brief evaluate first or second time derivative

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the time derivative will
                       be evaluated
      \param t     (i) The point in time in which the time derivative will
                       be evaluated
      \param deg   (i) The degree of the time derivative

     */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents() { return (9); };

   private:
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public DRT::UTILS::Function
  {
   public:
    BeltramiFunction(double c1);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (4); };

   private:
    double c1_;
  };

  /// special implementation for weakly compressible flow in a channel
  class ChannelWeaklyCompressibleFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (3); };
  };

  /// correction term for weakly compressible flow in a channel
  class CorrectionTermChannelWeaklyCompressibleFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (1); };
  };

  /// special implementation for weakly compressible flow in a channel with Fourier expansion up to
  /// 3rd order
  class ChannelWeaklyCompressibleFourier3Function : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (3); };
  };

  /// body force for weakly compressible flow in a channel with Fourier expansion up to 3rd order
  class BodyForceChannelWeaklyCompressibleFourier3Function : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (1); };
  };

  /// correction term for weakly compressible flow in a channel with Fourier expansion up to 3rd
  /// order
  class CorrectionTermChannelWeaklyCompressibleFourier3Function : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (1); };
  };


  /// special implementation for 2d Kim-Moin flow
  class KimMoinFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

      \brief evaluate first or second time derivative

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the time derivative will
                       be evaluated
      \param t     (i) The point in time in which the time derivative will
                       be evaluated
      \param deg   (i) The degree of the time derivative

     */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (this is a vector-valued functions)

    \return number of components (u,v,p)

    */
    virtual int NumberComponents() { return (3); };
  };


  /// special implementation for 2d Bochev test case (velocity and pressure)
  class BochevUPFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,p)

    */
    virtual int NumberComponents() { return (3); };
  };


  /// special implementation for 2d Bochev test case (rhs function)
  class BochevRHSFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (f1,f2)

    */
    virtual int NumberComponents() { return (2); };
  };



  /// special implementation for beltrami flow (rhs)
  class BeltramiRHS : public DRT::UTILS::Function
  {
   public:
    BeltramiRHS(int mat_id, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                       be evaluated
    \param x     (i) The point in space in which the time derivative will
                       be evaluated
    \param t     (i) The point in time in which the time derivative will
                       be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    virtual int NumberComponents() { return (3); };

   private:
    double kinviscosity_;
    bool is_stokes_;
  };



  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public DRT::UTILS::Function
  {
   public:
    KimMoinRHS(int mat_id, bool is_stationary, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                       be evaluated
    \param x     (i) The point in space in which the time derivative will
                       be evaluated
    \param t     (i) The point in time in which the time derivative will
                       be evaluated
    \param deg   (i) The degree of the time derivative

     */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    virtual int NumberComponents() { return (3); };

   private:
    double kinviscosity_;
    bool is_stationary_;
    bool is_stokes_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (analytical
  /// stress)
  class KimMoinStress : public DRT::UTILS::Function
  {
   public:
    KimMoinStress(int mat_id, bool is_stationary, double amplitude = 1.0);

    KimMoinStress(Teuchos::RCP<MAT::Material>& mat, bool is_stationary, double amplitude = 1.0);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents() { return (6); };

   private:
    double kinviscosity_;
    double density_;
    bool is_stationary_;
    double amplitude_;
  };


  /// special implementation for (randomly) disturbed 3d turbulent
  /// boundary-layer profile
  /// (currently fixed for low-Mach-number flow through a backward-facing step,
  ///  but may be easily manipulated to fit other profiles in other geometries)
  class TurbBouLayerFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (4); };
  };


  /// special implementation for (randomly) disturbed 3d turbulent
  /// boundary-layer profile
  /// (incompressible flow over backward-facing step,
  ///  corresponding to geometry of DNS by Le, Moin and Kim)
  class TurbBouLayerFunctionBFS : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (4); };
  };


  /// special implementation for (randomly) disturbed 3d turbulent boundary-layer profile
  /// (incompressible flow in the ORACLES test rig)
  class TurbBouLayerFunctionORACLES : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);
  };


  /// special implementation for stationary 2d Jeffery-Hamel flow
  class JefferyHamelFlowFunction : public DRT::UTILS::Function
  {
   public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents() { return (4); };

    /*!

    compute radial flow u as a function of the angle alpha

    \return radial velocity u

    \note this function is static such that we can get the radial velocity
          without creating the Function object

    */
    static double RadialVelocity(
        const double& theta  ///< angle between 0 and PI/4 (range is checked in debug mode)
    );
  };

}  // namespace FLD

#endif
