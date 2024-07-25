/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_utils_function.hpp"

#ifndef FOUR_C_STRUCTURE_NEW_FUNCTIONS_HPP
#define FOUR_C_STRUCTURE_NEW_FUNCTIONS_HPP

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::UTILS
{
  class FunctionManager;
}

namespace Solid
{
  /// add valid structure-specific function lines
  void AddValidStructureFunctions(Core::UTILS::FunctionManager& function_manager);

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIStructureFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureFunction(const Mat::PAR::StVenantKirchhoff& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    [[nodiscard]] std::size_t number_components() const override { return (2); };
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIStructureForceFunction : public Core::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureForceFunction(const Mat::PAR::StVenantKirchhoff& fparams);

    double evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> evaluate_time_derivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    [[nodiscard]] std::size_t number_components() const override { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
