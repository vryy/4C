/*!----------------------------------------------------------------------
\file tutorial_utils.H

\brief student's c++/baci tutorial utility methods

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

 *----------------------------------------------------------------------*/

/// All header files should have #define guards to prevent multiple inclusion.
/// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef TUTORIAL_UTILS_H
#define TUTORIAL_UTILS_H


#include <map>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_tutorial/tutorial_ele/tutorial_ele.H"


namespace TUTORIAL {

namespace UTILS {

/*!
\brief Apply the Dirichlet conditions to the global system and extract the reaction force

  Non-template function in header file must be declared inline to avoid duplicate definitions of functions.
  Define functions inline only when they are small, say, 10 lines or less.
  See LNM/Baci code style guide for further information.

\param stiff (out)    : pointer to system mytrix
\param RHS (out)      : pointer to right-hand side vector
\param inc (out)      : pointer to iteration increment
\param disp (out)     : pointer to displacment state vector
\param freact (out)   : pointer to reaction force vector
\param dirichcond (in): pointer to dirichlet condition

\return void

\note All input arguments are pointers (call by reference), since changes made in this method are supposed to be visible outside.
      Furthermore, we want to AVOID to make copies of data structures as often as possible.

*/
inline
void DoDirichletCond(LINALG::Matrix<7,7>* stiff, LINALG::Matrix<7,1>* RHS, LINALG::Matrix<7,1>* inc, LINALG::Matrix<7,1>* disp, LINALG::Matrix<7,1>* freact, const LINALG::Matrix<7,1>* dirichcond)
{
  // do dirichlet
  for(int dirich=0;dirich<(int)dirichcond->Rows();dirich++)
  {
    if((*dirichcond)(dirich)>-1234.0)
    {
      (*freact)(dirich)=(*RHS)(dirich);
      (*RHS)(dirich)=0.0;
      (*inc)(dirich)=0.0;
      (*disp)(dirich)=(*dirichcond)(dirich);


      for(int col=0;col<(int)stiff->Cols();++col)
      {
        (*stiff)(dirich,col) = 0.0;
      }
      for(int row=0;row<(int)stiff->Rows();++row)
      {
        (*stiff)(row,dirich) = 0.0;
      }
      (*stiff)(dirich,dirich)=1.0;
    }//
  }// do dirichlet
}; // DoDirichletCond


/*!
\brief Apply the Neumann condition to the right-hand side vector

  Non-template function in header file must be declared inline to avoid duplicate definitions of functions.
  Define functions inline only when they are small, say, 10 lines or less.
  See LNM/Baci code style guide for further information.

\param RHS (out): pointer to right-hand side vector
\param neumanncond (in): pointer to neumann condition

\return void

*/
inline
void DoNeumannCond(LINALG::Matrix<7,1>* RHS, const LINALG::Matrix<7,1>* neumanncond)
{
  // do neumann
  for(int neum=0;neum<(int)RHS->Rows();neum++)
  {
    if((*neumanncond)(neum)>-1234.0)
    {
      (*RHS)(neum)-=(*neumanncond)(neum);
    }
  } // do neumann

}; // DoNeumannCond


/*!
\brief Apply the Neumann condition to the right-hand side vector

  Non-template function in header file must be declared inline to avoid duplicate definitions of functions.
  Define functions inline only when they are small, say, 10 lines or less.
  See LNM/Baci code style guide for further information.

\param numberofelements (in): number of elements to be built
\param E (in): vector containing Young's moduli of elements to be built
\param X (in): vector containing the initial coordinates of the nodes of the elements to be built

\return discretization map associating the element id to the pointer to the element object

\note all input arguments are const, since they only provide information to this method and are not supposed to be changed

*/
inline
std::map<int, ELEMENTS::TutorialElement* >* BuildDiscretization(const int numberofelements, const LINALG::Matrix<6,1> E, const LINALG::Matrix<7,1> X)
{
  // build map containing element id and pointer to element
  std::map<int, ELEMENTS::TutorialElement* >* discret = new std::map<int, ELEMENTS::TutorialElement* >();

  // loop over desired amount of elements
  for (int ele=0; ele<numberofelements; ele++)
  {
    // add element to discretization map
    discret->insert( std::pair<int, ELEMENTS::TutorialElement* >(ele, (new ELEMENTS::TutorialElement(ele,E(ele),X(ele),X(ele+1))) ) );
  }

  return discret;

}; // BuildDiscretization


} // namespace UTILS
} // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef TUTORIAL_UTILS_H*/
