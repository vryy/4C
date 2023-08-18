/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

*----------------------------------------------------------------------*/

#include "baci_lib_function.H"

#ifndef BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H
#define BACI_FLUID_XFLUID_FUNCTIONS_COMBUST_H


namespace DRT
{
  class DiscretizationInterface;

  namespace INPUT
  {
    class LineDefinition;
  }

  namespace UTILS
  {
    /// add valid combustion-specific function lines
    void AddValidCombustFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

    /// try to create combustion-specific functions from a given line definition
    Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> TryCreateCombustFunction(
        const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);

    /// special implementation for a level set test function "Zalesak's disk"
    class ZalesaksDiskFunction : public DRT::UTILS::FunctionOfSpaceTime
    {
     public:
      /// ctor
      ZalesaksDiskFunction();

      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) override;
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public DRT::UTILS::FunctionOfSpaceTime
    {
     public:
      /// ctor
      CollapsingWaterColumnFunction();

      /// evaluate function at given position in space
      double Evaluate(const double* x, double t, std::size_t component) override;
    };
  }  // namespace UTILS
}  // namespace DRT

#endif
