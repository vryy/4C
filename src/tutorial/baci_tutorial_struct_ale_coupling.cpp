/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci nonlinear truss tutorial


\level 2

*/
/*---------------------------------------------------------------------*/

#include "baci_tutorial_struct_ale_coupling.H"

#include "baci_tutorial_nln_truss.H"


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
void TUTORIAL::StructAleCoupling::PostEvaluate(CORE::LINALG::Matrix<numele + 1, 1> rhs,
    CORE::LINALG::Matrix<numele + 1, numele + 1> stiff, bool eval_rhs, bool eval_stiff)
{
  return;
}
