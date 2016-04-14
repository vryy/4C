/*!----------------------------------------------------------------------
\file tutorial_nln_truss.H

\brief student's c++/baci nonlinear truss tutorial

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef TUTORIAL_NONLINEAR_TRUSS_H
#define TUTORIAL_NONLINEAR_TRUSS_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_tutorial/tutorial_ele/tutorial_ele.H"
#include "../drt_tutorial/tutorial_material/tutorial_mat.H"

/// define the number of elements
#define numele 6


namespace TUTORIAL {

class NonlinearTruss
{

public:

  /// constructor (ctor)
  explicit NonlinearTruss();

  /// here you define your geometry and boundary conditions directly in the source code
  std::vector<LINALG::Matrix<7,1> > ProblemDefinition();

  /// build discretization
  void SetupProblem();

  /// execute the time loop
  void TimeLoop();

  /// build solver object
  template < unsigned int numrow, unsigned int numcol >
  void SetupSolver();

  /// newton iteration (find equilibrium in every time/load step)
  void Newton();

  /// evaluate rhs and stiffness matrix
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff, bool eval_rhs, bool eval_stiff);

  /// evaluate rhs and stiffness matrix
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff)
  {Evaluate(rhs, stiff,true, true); return;};

  /// evaluate rhs and stiffness matrix
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff, bool eval_rhs)
  {Evaluate(rhs, stiff, eval_rhs, true); return;};

  /// post evaluate operator ( e.g. for manipulation of rhs and stiffness matrix after usual evaluate)
  /// can be overriden by specialized nonlinear truss problem, which inherits from this class.
  /// here the PostEvaluate() operator does nothing.
  virtual void PostEvaluate(LINALG::Matrix<numele+1,1> rhs, LINALG::Matrix<numele+1,numele+1> stiff, bool eval_rhs, bool eval_stiff)
  {return;};

  /// solve the linear system of equations stiff*disp=rhs
  void Solve(LINALG::Matrix<numele+1,1> rhs, LINALG::Matrix<numele+1,numele+1> stiff, LINALG::Matrix<numele+1,1> inc);

  /// check if newton is converged
  bool ConvergenceCheck(double tol);

  /// return pointer to the solver object
  Teuchos::RCP<LINALG::FixedSizeSerialDenseSolver<numele+1,numele+1,1> > Solver(){return nxnsolver_;};

  /// print the results to screen
  void PrintResults();

  /// apply dirichlet and neumann conditions
  void DoBoundaryConditions();

  /// clean up the tutorial; delete objects built with 'new' and singleton instances
  void TutorialDone();


protected :

  /// print the type of the tutorial to screen
  virtual void PrintTutorialType();

private :

  /// print geometry and boundary conditions
  void PrintProblem(std::vector<LINALG::Matrix<7,1> >& conditions);

  /// print used material law
  void PrintMaterial();

  /// print the elements and their members (output example for std::map)
  void PrintDiscretization();

  /// assemble an element stiffness matrix into the global stiffness matrix
  void AssembleMatrix(TUTORIAL::ELEMENTS::TutorialElement* ele, LINALG::Matrix<2,2>* ele_stiff_mat, LINALG::Matrix<numele+1,numele+1>* glob_stiff);

  /// assemble an element residual vector into the global right-hand-side
  void AssembleVector(TUTORIAL::ELEMENTS::TutorialElement* ele, LINALG::Matrix<2,1>* ele_rhs_vec, LINALG::Matrix<numele+1,1>* glob_rhs);

  /// discretization map ele id => TutorialElement
  std::map<int, ELEMENTS::TutorialElement* >* discretization_;

  // pointer to dirichlet condition
  LINALG::Matrix<numele+1,1>* dirichletvalues_;

  // pointer to neumann condition
  LINALG::Matrix<numele+1,1>* neumannvalues_;

  // the solver object
  Teuchos::RCP<LINALG::FixedSizeSerialDenseSolver<numele+1,numele+1,1> > nxnsolver_;

  // define and initialize vectors and matrices
  // RHS vector
  LINALG::Matrix<numele+1,1> rhs_;
  // stiffness matrix
  LINALG::Matrix<numele+1,numele+1> stiff_;
  // increment vector
  LINALG::Matrix<numele+1,1> inc_;
  // displacement vector
  LINALG::Matrix<numele+1,1> disp_;
  // reaction force vector
  LINALG::Matrix<numele+1,1> freact_;

};// class NonlinearTruss
} // namespace TUTORIAL

#endif
