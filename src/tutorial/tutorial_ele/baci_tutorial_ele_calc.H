/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial element evaluation


\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_TUTORIAL_ELE_CALC_H
#define BACI_TUTORIAL_ELE_CALC_H

#include "baci_linalg_fixedsizematrix.H"

#include "baci_utils_singleton_owner.H"

namespace TUTORIAL
{
  namespace ELEMENTS
  {
    /// FORWARD DECLARATION
    /// Here we declare the class 'TutorialElement'.
    /// This is needed because an pointer of this type is
    /// provided to Evaluate(...). We could also include
    /// the tutorial_ele.H-file but it should be avoided
    /// to include headers in headers. This increases
    /// dependencies in the code and slows down compilation
    /// and linking.
    class TutorialElement;

    /// class definition of class 'TutorialEleCalc'
    class TutorialEleCalc
    {
     public:
      /// Singleton access method
      static TutorialEleCalc* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Empty destructor
      virtual ~TutorialEleCalc(){};

      /// evaluates the desired matrix and/or vector
      void Evaluate(TutorialElement* element, CORE::LINALG::Matrix<7, 1>* state,
          CORE::LINALG::Matrix<2, 2>* elematrix, CORE::LINALG::Matrix<2, 1>* elevector)
      {
        Evaluate(element, state, elematrix, elevector, true, true);
        return;
      };

      /// evaluates the desired matrix and/or vector
      void Evaluate(TutorialElement* element, CORE::LINALG::Matrix<7, 1>* state,
          CORE::LINALG::Matrix<2, 2>* elematrix, CORE::LINALG::Matrix<2, 1>* elevector,
          bool eval_mat)
      {
        Evaluate(element, state, elematrix, elevector, eval_mat, true);
        return;
      };

      /// evaluates the desired matrix and/or vector
      void Evaluate(TutorialElement* element, CORE::LINALG::Matrix<7, 1>* state,
          CORE::LINALG::Matrix<2, 2>* elematrix, CORE::LINALG::Matrix<2, 1>* elevector,
          bool eval_mat, bool eval_vec);


     protected:
      /// protected constructor since we are singleton
      TutorialEleCalc();

    };  // class TutorialEleCalc

  }  // namespace ELEMENTS
}  // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef BACI_TUTORIAL_ELE_CALC_H*/
