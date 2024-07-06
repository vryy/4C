/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_mat_fluid_weakly_compressible.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_utils_function.hpp"

#ifndef FOUR_C_FLUID_FUNCTIONS_HPP
#define FOUR_C_FLUID_FUNCTIONS_HPP

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::UTILS
{
  class FunctionManager;
}


namespace FLD
{
  /// add valid fluid-specific function lines
  void AddValidFluidFunctions(Core::UTILS::FunctionManager& function_manager);

  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiUP(const Mat::PAR::NewtonianFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t number_components() const override { return (4); };

   private:
    double density_;
  };

  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiGradU(const Mat::PAR::NewtonianFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t number_components() const override { return (9); };
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity,
  /// pressure)
  class KimMoinUP : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinUP(const Mat::PAR::NewtonianFluid& fparams, bool is_stationary);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,v,w,p)
     */
    std::size_t number_components() const override { return (4); };

   private:
    double density_;
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of
  /// velocity)
  class KimMoinGradU : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinGradU(const Mat::PAR::NewtonianFluid& fparams, bool is_stationary);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )
     */
    std::size_t number_components() const override { return (9); };

   private:
    double kinviscosity_;
    bool is_stationary_;
  };

  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiFunction(double c1);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t number_components() const override { return (4); };

   private:
    double c1_;
  };

  /// special implementation for weakly compressible flow in a channel
  class ChannelWeaklyCompressibleFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w,p)
     */
    std::size_t number_components() const override { return (3); };
  };

  /// correction term for weakly compressible flow in a channel
  class CorrectionTermChannelWeaklyCompressibleFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (1); };
  };

  /// special implementation for weakly compressible Poiseuille flow
  class WeaklyCompressiblePoiseuilleFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (6); };

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
  class WeaklyCompressiblePoiseuilleForceFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressiblePoiseuilleForceFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (3); };

   private:
    double length_;
    double halfheight_;
    double meanvelocityexit_;
    double viscosity_;
    double refdensity_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow with manufactured solution
  class WeaklyCompressibleManufacturedFlowFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowFunction(const Mat::PAR::WeaklyCompressibleFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (6); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow with manufactured solution (force)
  class WeaklyCompressibleManufacturedFlowForceFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleManufacturedFlowForceFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (3); };

   private:
    double viscosity_;
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem
  class WeaklyCompressibleEtienneCFDFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDFunction(const Mat::PAR::WeaklyCompressibleFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (force)
  class WeaklyCompressibleEtienneCFDForceFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDForceFunction(const Mat::PAR::WeaklyCompressibleFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (3); };

   private:
    double refdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne CFD problem (viscosity)
  class WeaklyCompressibleEtienneCFDViscosityFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneCFDViscosityFunction(const Mat::PAR::WeaklyCompressibleFluid& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (1); };
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIFluidFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const Mat::PAR::StVenantKirchhoff& fparams_struc);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (6); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIFluidForceFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidForceFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const Mat::PAR::StVenantKirchhoff& fparams_struc);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (3); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (viscosity)
  class WeaklyCompressibleEtienneFSIFluidViscosityFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIFluidViscosityFunction(
        const Mat::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const Mat::PAR::StVenantKirchhoff& fparams_struc);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t number_components() const override { return (1); };

   private:
    double refdensity_;
    double refpressure_;
    double comprcoeff_;
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

  /// special implementation for beltrami flow (rhs)
  class BeltramiRHS : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    BeltramiRHS(const Mat::PAR::NewtonianFluid& fparams, bool is_stokes);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t number_components() const override { return (3); };

   private:
    double kinviscosity_;
    bool is_stokes_;
  };

  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinRHS(const Mat::PAR::NewtonianFluid& fparams, bool is_stationary, bool is_stokes);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    /*!
     * \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     *  \return number of components (u,v,w)
     */
    std::size_t number_components() const override { return (3); };

   private:
    double kinviscosity_;
    bool is_stationary_;
    bool is_stokes_;
  };

  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (analytical
  /// stress)
  class KimMoinStress : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    KimMoinStress(
        const Mat::PAR::NewtonianFluid& fparams, bool is_stationary, double amplitude = 1.0);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::size_t number_components() const override { return (6); };

   private:
    double kinviscosity_;
    double density_;
    bool is_stationary_;
    double amplitude_;
  };
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
