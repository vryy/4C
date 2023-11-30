/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_mat_stvenantkirchhoff.H"
#include "baci_utils_function.H"

#ifndef BACI_STRUCTURE_NEW_FUNCTIONS_H
#define BACI_STRUCTURE_NEW_FUNCTIONS_H

namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace CORE::UTILS
{
  class FunctionManager;
}


namespace STR
{
  /// add valid structure-specific function lines
  void AddValidStructureFunctions(CORE::UTILS::FunctionManager& function_manager);

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIStructureFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (2); };

   private:
    double poissonratio_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIStructureForceFunction : public CORE::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureForceFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

}  // namespace STR

#endif
