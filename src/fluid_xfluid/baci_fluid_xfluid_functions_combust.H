/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

*----------------------------------------------------------------------*/

#include "baci_config.H"

#include "baci_io_linedefinition.H"
#include "baci_utils_function.H"

#ifndef BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H
#define BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H

BACI_NAMESPACE_OPEN

namespace CORE::UTILS
{
  class FunctionManager;
}

namespace DRT
{
  class Discretization;

  // namespace UTILS

  namespace UTILS
  {
    /// add valid combustion-specific function lines
    void AddValidCombustFunctions(CORE::UTILS::FunctionManager& function_manager);

    /// special implementation for a level set test function "Zalesak's disk"
    class ZalesaksDiskFunction : public CORE::UTILS::FunctionOfSpaceTime
    {
     public:
      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) const override;
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public CORE::UTILS::FunctionOfSpaceTime
    {
     public:
      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) const override;
    };
  }  // namespace UTILS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
