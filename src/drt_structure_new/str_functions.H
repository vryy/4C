/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "../drt_lib/drt_function.H"
#include "../drt_mat/stvenantkirchhoff.H"

#ifndef STR_FUNCTION_H
#define STR_FUNCTION_H

namespace DRT
{
  class Discretization;
  namespace INPUT
  {
    class LineDefinition;
  }  // namespace INPUT
}  // namespace DRT


namespace STR
{
  /// add valid structure-specific function lines
  void AddValidStructureFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

  /// try to create structure-specific functions from a given line definition
  Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> TryCreateStructureFunction(
      Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def,
      DRT::UTILS::FunctionManager& manager, const int index_current_funct_in_manager);

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIStructureFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (2); };

   private:
    double poissonratio_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIStructureForceFunction : public DRT::UTILS::FunctionOfSpaceTime
  {
   public:
    WeaklyCompressibleEtienneFSIStructureForceFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    double Evaluate(const double* x, double t, std::size_t component) override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) override;

    std::size_t NumberComponents() override { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

}  // namespace STR

#endif
