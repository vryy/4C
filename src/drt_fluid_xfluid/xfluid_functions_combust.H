/*!-------------------------------------------------------------------------*
\file xfluid_functions_combust.H

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

\level 2

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
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
    void CombustValidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);
    bool CombustFunctionHaveNamed(Teuchos::RCP<DRT::INPUT::LineDefinition> functionnamed,
        std::vector<Teuchos::RCP<DRT::UTILS::Function>>* function_);

    /// initial level-set field for two merging bubbles
    class BubbleFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      BubbleFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a level set test function
    class ZalesaksDiskFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      ZalesaksDiskFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a combustion test function
    class CircularFlame2Function : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CircularFlame2Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a combustion test function
    class CircularFlame3Function : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CircularFlame3Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a combustion test function
    class CircularFlame4Function : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CircularFlame4Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// initial level-set field for dam break with obstacle
    class DamBreakObstacle : public DRT::UTILS::Function
    {
     public:
      /// ctor
      DamBreakObstacle();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CollapsingWaterColumnFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunctionCoarse : public DRT::UTILS::Function
    {
     public:
      /// ctor
      CollapsingWaterColumnFunctionCoarse();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for the G-function in the ORACLES problem
    class ORACLESGFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      ORACLESGFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a level set test function
    class RotatingConeFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      RotatingConeFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// special implementation for a xfem test function
    class LevelSetCutTestFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      LevelSetCutTestFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };

    /// initial level-set field for impact of drop
    class ImpactFunction : public DRT::UTILS::Function
    {
     public:
      /// ctor
      ImpactFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t);
    };


  }  // namespace UTILS
}  // namespace DRT

#endif
