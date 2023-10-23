/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_lib_function.H"

#ifndef BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H
#define BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H


namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    class FunctionManager;
  }  // namespace UTILS

  namespace UTILS
  {
    /// add valid combustion-specific function lines
    void AddValidCombustFunctions(DRT::UTILS::FunctionManager& function_manager);

    /// special implementation for a level set test function "Zalesak's disk"
    class ZalesaksDiskFunction : public DRT::UTILS::FunctionOfSpaceTime
    {
     public:
      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) const override;
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public DRT::UTILS::FunctionOfSpaceTime
    {
     public:
      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) const override;
    };
  }  // namespace UTILS
}  // namespace DRT

#endif
