/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci nonlinear truss tutorial


\level 2

*/
/*---------------------------------------------------------------------*/

// All header files should have #define guards to prevent multiple inclusion.
// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef TUTORIAL_STRUCT_ALE_COUPLING_H
#define TUTORIAL_STRUCT_ALE_COUPLING_H

#include "tutorial_nln_truss.H"

namespace TUTORIAL
{
  class StructAleCoupling : public NonlinearTruss
  {
   public:
    /// constructor (ctor)
    StructAleCoupling();

    /// constructor (ctor)
    virtual ~StructAleCoupling(){};

   protected:
    /// print the type of the tutorial to screen
    virtual void PrintTutorialType();

    /// post evaluate operator ( e.g. for manipulation of rhs and stiffness matrix after usual
    /// evaluate)
    virtual void PostEvaluate(LINALG::Matrix<numele + 1, 1> rhs,
        LINALG::Matrix<numele + 1, numele + 1> stiff, bool eval_rhs = true, bool eval_stiff = true);

  };  // class StructAleCoupling
}  // namespace TUTORIAL

#endif /* TUTORIAL_STRUCT_ALE_COUPLING_H */
