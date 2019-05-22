/*---------------------------------------------------------------------*/
/*!

\brief student's c++/baci tutorial element evaluation

\maintainer  Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/

#include <cstddef>

#include "../tutorial_ele/tutorial_ele.H"
#include "../tutorial_ele/tutorial_ele_calc.H"
#include "../tutorial_material/tutorial_mat.H"


/*----------------------------------------------------------------------*
 * construct/return pointer to singleton                                |
 *----------------------------------------------------------------------*/
TUTORIAL::ELEMENTS::TutorialEleCalc* TUTORIAL::ELEMENTS::TutorialEleCalc::Instance(bool create)
{
  static TutorialEleCalc* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new TutorialEleCalc();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 * clean up function to destroy instance in the end                     |
 *----------------------------------------------------------------------*/
void TUTORIAL::ELEMENTS::TutorialEleCalc::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 * constructor (protected)                                              |
 *----------------------------------------------------------------------*/
TUTORIAL::ELEMENTS::TutorialEleCalc::TutorialEleCalc() {}


/*----------------------------------------------------------------------*
 * evaluate terms                                                       |
 *----------------------------------------------------------------------*/
void TUTORIAL::ELEMENTS::TutorialEleCalc::Evaluate(TUTORIAL::ELEMENTS::TutorialElement* element,
    LINALG::Matrix<7, 1>* state, LINALG::Matrix<2, 2>* elematrix, LINALG::Matrix<2, 1>* elevector,
    bool eval_mat, bool eval_vec)
{
  // get number of nodes from the element
  int numnodeperele = element->NumNode();

  // So far, the gauss points (or generally called integration points)
  // are hard coded. For now we have 2 int. points per element.
  // number of gp per ele
  int numgp = 2;
  // gp coordinates
  LINALG::Matrix<2, 1> gpcoord(true);
  gpcoord(0) = -1.0 / sqrt(3.0);
  gpcoord(1) = 1.0 / sqrt(3.0);

  // gaus point loop
  for (int gp = 0; gp < numgp; gp++)
  {
    // eval shape functions at gp
    LINALG::Matrix<2, 1> N(true);
    N(0) = 0.5 * (1 - gpcoord(gp));
    N(1) = 0.5 * (1 + gpcoord(gp));
    // eval shape function derivs w.r.t. xi at gp
    LINALG::Matrix<2, 1> N_xi(true);
    N_xi(0) = -0.5;
    N_xi(1) = 0.5;
    // dX_dxi
    double dX_dxi = N_xi(0) * (element->X(0)) + N_xi(1) * (element->X(1));
    // eval shape function derivs w.r.t. X at gp
    LINALG::Matrix<2, 1> N_X(true);
    N_X(0) = N_xi(0) / dX_dxi;
    N_X(1) = N_xi(1) / dX_dxi;
    // evaluate defgrd
    double defgrd = 1.0 + (((*state)((element->Id()) + 1) - (*state)(element->Id())) /
                              (element->X(1) - element->X(0)));
    // evaluate GL-Strain at gp
    double EGL = 0.5 * (pow(defgrd, 2) - 1.0);
    // integral transform
    double J = dX_dxi;
    // dF/dd
    LINALG::Matrix<2, 1> B(true);
    B(0) = N_xi(0) / J;
    B(1) = N_xi(1) / J;

    // MATERIAL LAW EVALUATION
    ///////////////////////////////////////////////////////////
    // evaluate Material Law

    // initialize passive stress
    double Spassive = -1234.0;

    // initialize linearization of material law
    static LINALG::Matrix<2, 2> cmat(true);
    cmat.Clear();

    // call evaluate
    element->Material()->EvaluateMaterial(element, &Spassive, &cmat, EGL, defgrd, J, B);
    ///////////////////////////////////////////////////////////

    if (eval_vec)
    {
      // evaluate residual

      for (int node = 0; node < numnodeperele; ++node)
      {
        // eval right hand side
        (*elevector)(node) += (Spassive)*defgrd * B(node) * J;
      }
    }  // eval_rhs


    if (eval_mat)
    {
      // evaluate stiffness matrix
      // see NiliFEM Script page 102 !
      // geometric stiffness
      for (int node = 0; node < numnodeperele; ++node)
        for (int colnode = 0; colnode < numnodeperele; ++colnode)
          (*elematrix)(node, colnode) += B(node) * B(colnode) * (Spassive)*J;


      // combine material stress linearization and geometric stiffness
      elematrix->Update(1.0, cmat, 1.0);

    }  // eval_mat

  }  // gp loop


}  // Evaluate
