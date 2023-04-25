/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "lib_function.H"
#include "mat_newtonianfluid.H"
#include "mat_fluid_weakly_compressible.H"
#include "mat_stvenantkirchhoff.H"

#ifndef FLUID_FUNCTIONS_H
#define FLUID_FUNCTIONS_H

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
  /// add valid fluid-specific function lines
  void AddValidFluidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

  /// try to create fluid-specific functions from a given line definition
  Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> TryCreateFluidFunction(
      const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);

  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiUP(const MAT::PAR::NewtonianFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() override { return (4); };

   private:
    double density_;
    double kinviscosity_;
  };

  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiGradU(const MAT::PAR::NewtonianFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t NumberComponents() override { return (9); };

   private:
    double kinviscosity_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity,
  /// pressure)
  class KimMoinUP : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinUP(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() override { return (4); };

   private:
    double density_;
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of
  /// velocity)
  class KimMoinGradU : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinGradU(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t NumberComponents() override { return (9); };

   private:
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiFunction(double c1);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() override { return (4); };

   private:
    double c1_;
  };

  /// special implementation for weakly compressible flow in a channel
  class ChannelWeaklyCompressibleFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() override { return (3); };
  };

  /// correction term for weakly compressible flow in a channel
  class CorrectionTermChannelWeaklyCompressibleFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (1); };
  };

  /// special implementation for weakly compressible Poiseuille flow
  class WeaklyCompressiblePoiseuilleFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

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
  class WeaklyCompressiblePoiseuilleForceFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

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
  class WeaklyCompressibleManufacturedFlowFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (6); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow with manufactured solution (force)
  class WeaklyCompressibleManufacturedFlowForceFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (3); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem
  class WeaklyCompressibleEtienneCFDFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (force)
  class WeaklyCompressibleEtienneCFDForceFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDForceFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (3); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (viscosity)
  class WeaklyCompressibleEtienneCFDViscosityFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDViscosityFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (1); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIFluidFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

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
  class WeaklyCompressibleEtienneFSIFluidForceFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

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
  class WeaklyCompressibleEtienneFSIFluidViscosityFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidViscosityFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

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
  class BeltramiRHS : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stokes);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t NumberComponents() override { return (3); };

   private:
    double kinviscosity_;
    bool is_stokes_;
  };

  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, bool is_stokes);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t NumberComponents() override { return (3); };

   private:
    double kinviscosity_;
    bool is_stationary_;
    bool is_stokes_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (analytical
  /// stress)
  class KimMoinStress : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinStress(
        const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, double amplitude = 1.0);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::size_t NumberComponents() override { return (6); };

   private:
    double kinviscosity_;
    double density_;
    bool is_stationary_;
    double amplitude_;
  };
}  // namespace FLD

#endif
