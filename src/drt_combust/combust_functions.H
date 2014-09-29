/*!-------------------------------------------------------------------------*
\file combust_functions.H

\brief Managing and evaluating of spatial functions for combustion and two-phase flow problems

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef COMBUST_FUNCTION_H
#define COMBUST_FUNCTION_H

namespace DRT
{
  class Discretization;

  namespace UTILS
  {

    /// initial level-set field for two merging bubbles
    class BubbleFunction : public Function
    {
      public:

      /// ctor
      BubbleFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a level set test function
    class ZalesaksDiskFunction : public Function
    {
    public:

      /// ctor
      ZalesaksDiskFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a combustion test function
    class CircularFlame2Function : public Function
    {
    public:

      /// ctor
      CircularFlame2Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a combustion test function
    class CircularFlame3Function : public Function
    {
    public:

      /// ctor
      CircularFlame3Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a combustion test function
    class CircularFlame4Function : public Function
    {
    public:

      /// ctor
      CircularFlame4Function();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// initial level-set field for dam break with obstacle
    class DamBreakObstacle : public Function
    {
      public:

      /// ctor
      DamBreakObstacle();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public Function
    {
    public:

      /// ctor
      CollapsingWaterColumnFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunctionCoarse : public Function
    {
    public:

      /// ctor
      CollapsingWaterColumnFunctionCoarse();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation level set cut utilizing xfluid
    class GerstenbergerForwardfacingStep : public Function
    {
    public:

      /// ctor
      GerstenbergerForwardfacingStep();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for the G-function in the ORACLES problem
    class ORACLESGFunction : public Function
    {
    public:

      /// ctor
      ORACLESGFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a level set test function
    class RotatingConeFunction : public Function
    {
    public:

      /// ctor
      RotatingConeFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

    /// special implementation for a xfem test function
    class LevelSetCutTestFunction : public Function
    {
    public:

      /// ctor
      LevelSetCutTestFunction();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

  }
}

#endif
