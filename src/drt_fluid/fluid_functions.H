/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "../drt_lib/drt_function.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/fluid_weakly_compressible.H"
#include "../drt_mat/stvenantkirchhoff.H"

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
  // SEBASTIAN
  // using FunctioMap = std::unordered_map<std::string,
  // std::function<Teuchos::RCP<DRT::Function>(DRT::LineDefinition)>>;

  // FunctionMap GetFunctions()
  //{
  //  return {
  //    {"Waklaynfkjgkjagj", BeltramiUP::from}
  //  };
  //}

  /// Valid fluid function lines
  void FluidValidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);


  /// Return fluid function for given line definition
  Teuchos::RCP<DRT::UTILS::Function> FluidTryCreateFunction(
      Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def);


  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public DRT::UTILS::Function
  {
   public:
    BeltramiUP(const MAT::PAR::NewtonianFluid& fparams);

    // SEBASTIAN
    // static Teuchos::RCP<DRT::Function> from(LineDefinition lineDef);

    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    std::size_t NumberComponents() override { return (4); };

   private:
    double density_;
    double kinviscosity_;
  };


  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public DRT::UTILS::Function
  {
   public:
    BeltramiGradU(const MAT::PAR::NewtonianFluid& fparams);

    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    std::size_t NumberComponents() override { return (9); };

   private:
    double kinviscosity_;
  };


  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity,
  /// pressure)
  class KimMoinUP : public DRT::UTILS::Function
  {
   public:
    KimMoinUP(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    std::size_t NumberComponents() override { return (4); };

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
    KimMoinGradU(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated
      \param t     (i) The point in time in which the function will be
                       evaluated

     */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    std::size_t NumberComponents() override { return (9); };

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
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    std::size_t NumberComponents() override { return (4); };

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
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    std::size_t NumberComponents() override { return (3); };
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
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    std::size_t NumberComponents() override { return (1); };
  };


  /// special implementation for weakly compressible Poiseuille flow
  class WeaklyCompressiblePoiseuilleFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressiblePoiseuilleFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (6); };

   private:
    double length_;
    double halfheight_;
    double meanvelocityexit_;
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible Poiseuille flow (force)
  class WeaklyCompressiblePoiseuilleForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressiblePoiseuilleForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (3); };

   private:
    double length_;
    double halfheight_;
    double meanvelocityexit_;
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow with manufactured solution
  class WeaklyCompressibleManufacturedFlowFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleManufacturedFlowFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (6); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow with manufactured solution (force)
  class WeaklyCompressibleManufacturedFlowForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleManufacturedFlowForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (3); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow - Etienne CFD problem
  class WeaklyCompressibleEtienneCFDFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneCFDFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow - Etienne CFD problem (force)
  class WeaklyCompressibleEtienneCFDForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneCFDForceFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (3); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow - Etienne CFD problem (viscosity)
  class WeaklyCompressibleEtienneCFDViscosityFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneCFDViscosityFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (1); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };


  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIFluidFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIFluidFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };


  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIFluidForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIFluidForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (3); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };


  /// special implementation for weakly compressible flow - Etienne FSI problem (viscosity)
  class WeaklyCompressibleEtienneFSIFluidViscosityFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIFluidViscosityFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    std::size_t NumberComponents() override { return (1); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };


  /// special implementation for beltrami flow (rhs)
  class BeltramiRHS : public DRT::UTILS::Function
  {
   public:
    BeltramiRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    std::size_t NumberComponents() override { return (3); };

   private:
    double kinviscosity_;
    bool is_stokes_;
  };


  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public DRT::UTILS::Function
  {
   public:
    KimMoinRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated
    \param t     (i) The point in time in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t) override;

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
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    std::size_t NumberComponents() override { return (3); };

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
    KimMoinStress(
        const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, double amplitude = 1.0);


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t) override;

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    std::size_t NumberComponents() override { return (6); };

   private:
    double kinviscosity_;
    double density_;
    bool is_stationary_;
    double amplitude_;
  };

}  // namespace FLD

#endif
