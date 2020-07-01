/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial material St.Venant-Kirchhoff


\level 2

*/
/*---------------------------------------------------------------------*/

#include "tutorial_mat.H"
#include "../tutorial_ele/tutorial_ele.H"


TUTORIAL::MATERIAL::TutorialMat::TutorialMat(double young) : Young_(young) {}

/*----------------------------------------------------------------------*/

TUTORIAL::MATERIAL::TutorialMatStVenantKirchhoff::TutorialMatStVenantKirchhoff(const double young)
    : TutorialMat(young)
{
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::MATERIAL::TutorialMatStVenantKirchhoff::EvaluateMaterial(
    TUTORIAL::ELEMENTS::TutorialElement* ele, double* stress, LINALG::Matrix<2, 2>* linearization,
    const double strain, const double defgrd, const double J, const LINALG::Matrix<2, 1> B)
{
  // see NiliFEM Script page 101 !
  // calculate stress
  *stress = strain * YoungModulus();

  // how many nodes does the element have ?
  int numnodeperele = ele->NumNode();

  // see NiliFEM Script page 101 !
  // elastic and initial displacement stiffness
  for (int node = 0; node < numnodeperele; ++node)
    for (int colnode = 0; colnode < numnodeperele; ++colnode)
      (*linearization)(node, colnode) +=
          B(node) * defgrd * YoungModulus() * defgrd * B(colnode) * J;

  return;
}
