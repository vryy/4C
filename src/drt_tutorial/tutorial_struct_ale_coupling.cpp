/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci nonlinear truss tutorial

\maintainer  Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/

#include "tutorial_struct_ale_coupling.H"
#include "tutorial_nln_truss.H"


/// ctor
TUTORIAL::StructAleCoupling::StructAleCoupling() : TUTORIAL::NonlinearTruss() {}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::StructAleCoupling::PrintTutorialType()
{
  std::cout << "\n YOU CHOSE THE STRUCT-ALE COUPLING TOY PROBLEM ! \n\n" << std::endl;

  return;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::StructAleCoupling::PostEvaluate(LINALG::Matrix<numele + 1, 1> rhs,
    LINALG::Matrix<numele + 1, numele + 1> stiff, bool eval_rhs, bool eval_stiff)
{
  return;
}
