/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef XFLUID_FUNCTIONS_COMBUST_H
#define XFLUID_FUNCTIONS_COMBUST_H


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
    void CombustValidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

    /// try to create combustion-specific functions from a given line definition
    Teuchos::RCP<DRT::UTILS::Function> CombustTryCreateFunction(
        Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def);

    /// special implementation for a level set test function "Zalesak's disk"
    class ZalesaksDiskFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      ZalesaksDiskFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t) override;
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CollapsingWaterColumnFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t) override;
    };
  }  // namespace UTILS
}  // namespace DRT

#endif
