/*!----------------------------------------------------------------------
\file tutorial_struct_ale_coupling.cpp

\brief student's c++/baci nonlinear truss tutorial


\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
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
