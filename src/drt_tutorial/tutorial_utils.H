/*!----------------------------------------------------------------------
\file tutorial_utils.H

\brief student's c++/baci tutorial utility methods

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef TUTORIAL_UTILS_H
#define TUTORIAL_UTILS_H

#include <map>
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_tutorial/tutorial_ele/tutorial_ele.H"

namespace TUTORIAL {

namespace UTILS {

void DoDirichletCond(LINALG::Matrix<7,7>* stiff, LINALG::Matrix<7,1>* RHS, LINALG::Matrix<7,1>* inc, LINALG::Matrix<7,1>* disp, LINALG::Matrix<7,1>* freact, LINALG::Matrix<7,1>* dirichcond)
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


  void DoNeumannCond(LINALG::Matrix<7,1>* RHS, LINALG::Matrix<7,1>* neumanncond)
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


  std::map<int, ELEMENTS::TutorialElement* >* BuildDiscretization(const int numberofelements, LINALG::Matrix<6,1> E, LINALG::Matrix<7,1> X)
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
