/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_config.hpp"

#include "baci_mat_fluid_weakly_compressible.hpp"
#include "baci_mat_newtonianfluid.hpp"
#include "baci_mat_stvenantkirchhoff.hpp"
#include "baci_utils_function.hpp"

#ifndef FOUR_C_FLUID_FUNCTIONS_HPP
#define FOUR_C_FLUID_FUNCTIONS_HPP

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::UTILS
{
  class FunctionManager;
}


namespace FLD
{
  /// add valid fluid-specific function lines
  void AddValidFluidFunctions(CORE::UTILS::FunctionManager& function_manager);

  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiUP(const MAT::PAR::NewtonianFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() const override { return (4); };

   private:
    double density_;
  };

  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiGradU(const MAT::PAR::NewtonianFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t NumberComponents() const override { return (9); };
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity,
  /// pressure)
  class KimMoinUP : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinUP(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() const override { return (4); };

   private:
    double density_;
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of
  /// velocity)
  class KimMoinGradU : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinGradU(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t NumberComponents() const override { return (9); };

   private:
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiFunction(double c1);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() const override { return (4); };

   private:
    double c1_;
  };

  /// special implementation for weakly compressible flow in a channel
  class ChannelWeaklyCompressibleFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t NumberComponents() const override { return (3); };
  };

  /// correction term for weakly compressible flow in a channel
  class CorrectionTermChannelWeaklyCompressibleFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (1); };
  };

  /// special implementation for weakly compressible Poiseuille flow
  class WeaklyCompressiblePoiseuilleFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (6); };

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
  class WeaklyCompressiblePoiseuilleForceFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (3); };

   private:
    double length_;
    double halfheight_;
    double meanvelocityexit_;
    double viscosity_;
    double refdensity_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow with manufactured solution
  class WeaklyCompressibleManufacturedFlowFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (6); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow with manufactured solution (force)
  class WeaklyCompressibleManufacturedFlowForceFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (3); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem
  class WeaklyCompressibleEtienneCFDFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (force)
  class WeaklyCompressibleEtienneCFDForceFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDForceFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (3); };

   private:
    double refdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (viscosity)
  class WeaklyCompressibleEtienneCFDViscosityFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDViscosityFunction(const MAT::PAR::WeaklyCompressibleFluid& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (1); };
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIFluidFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIFluidForceFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (3); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (viscosity)
  class WeaklyCompressibleEtienneFSIFluidViscosityFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidViscosityFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (1); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for beltrami flow (rhs)
  class BeltramiRHS : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stokes);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t NumberComponents() const override { return (3); };

   private:
    double kinviscosity_;
    bool is_stokes_;
  };

  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, bool is_stokes);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t NumberComponents() const override { return (3); };

   private:
    double kinviscosity_;
    bool is_stationary_;
    bool is_stokes_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (analytical
  /// stress)
  class KimMoinStress : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinStress(
        const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, double amplitude = 1.0);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (6); };

   private:
    double kinviscosity_;
    double density_;
    bool is_stationary_;
    double amplitude_;
  };
}  // namespace FLD

BACI_NAMESPACE_CLOSE

#endif
