/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

*----------------------------------------------------------------------*/

#include "baci_lib_function.H"
#include "baci_lib_linedefinition.H"

#ifndef BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H
#define BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H


namespace DRT
{
  class Discretization;

  namespace INPUT
  {
    class LineDefinition;
    class Lines;
  }  // namespace INPUT

  namespace UTILS
  {
    /// add valid combustion-specific function lines
    void AddValidCombustFunctionLines(DRT::INPUT::Lines& lines);

    /// try to create combustion-specific functions from a given line definition
    Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> TryCreateCombustFunction(
        const std::vector<DRT::INPUT::LineDefinition>& function_line_defs);

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
