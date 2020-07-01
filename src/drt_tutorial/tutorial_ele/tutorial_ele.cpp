/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial element implementation


\level 2

*/
/*---------------------------------------------------------------------*/

#include "../tutorial_ele/tutorial_ele.H"
#include "../tutorial_ele/tutorial_ele_calc.H"
#include "../tutorial_material/tutorial_mat.H"


/// ctor
TUTORIAL::ELEMENTS::TutorialElement::TutorialElement(int id, double E, double X0, double X1)
    : id_(id), numnode_(2), numdof_(2)
{
  // fill vector with initial node positions
  X_.push_back(X0);
  X_.push_back(X1);

  // assign the dof gids to the nodes of the element
  // so element 0 has dof gids 0 and 1
  // element 1 has dof gids 1 and 2, and so on ...
  // you see, that the discretization is built from left to right.
  //
  dofs_.push_back(id);
  dofs_.push_back(id + 1);

  // build the material for this element
  // we only have the St. Venant-Kirchhoff for the tutorial so far
  material_ = Teuchos::rcp(new TUTORIAL::MATERIAL::TutorialMatStVenantKirchhoff(E));
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::ELEMENTS::TutorialElement::Evaluate(LINALG::Matrix<7, 1>* state,
    LINALG::Matrix<2, 2>* elematrix, LINALG::Matrix<2, 1>* elevector, bool eval_mat, bool eval_vec)
{
  // get calc instance (singleton object)
  TUTORIAL::ELEMENTS::TutorialEleCalc* calculator_object = TutorialEleCalc::Instance();

  // evaluate 'this' element
  calculator_object->Evaluate(this, state, elematrix, elevector, eval_mat, eval_vec);

  return;
}
