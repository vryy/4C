/*!----------------------------------------------------------------------
\file tutorial_ele_calc.H

\brief student's c++/baci tutorial element evaluation


\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef TUTORIAL_ELE_CALC_H
#define TUTORIAL_ELE_CALC_H

#include "../../linalg/linalg_fixedsizematrix.H"


namespace TUTORIAL {

namespace ELEMENTS {

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
  static TutorialEleCalc* Instance( bool create=true );

  /// Empty destructor
  virtual ~TutorialEleCalc() {};

  /// called upon destruction
  virtual void Done();

  /// evaluates the desired matrix and/or vector
  void Evaluate(TutorialElement* element, LINALG::Matrix<7,1>* state, LINALG::Matrix<2,2>* elematrix, LINALG::Matrix<2,1>* elevector)
  {Evaluate( element, state, elematrix, elevector, true, true); return;};

  /// evaluates the desired matrix and/or vector
  void Evaluate(TutorialElement* element, LINALG::Matrix<7,1>* state, LINALG::Matrix<2,2>* elematrix, LINALG::Matrix<2,1>* elevector, bool eval_mat)
  {Evaluate( element, state, elematrix, elevector, eval_mat, true); return;};

  /// evaluates the desired matrix and/or vector
  void Evaluate(TutorialElement* element, LINALG::Matrix<7,1>* state, LINALG::Matrix<2,2>* elematrix, LINALG::Matrix<2,1>* elevector, bool eval_mat, bool eval_vec);


protected:

  /// protected constructor since we are singleton
  TutorialEleCalc();

}; // class TutorialEleCalc

} // namespace ELEMENTS
} // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef TUTORIAL_ELE_CALC_H*/
