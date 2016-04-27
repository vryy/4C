/*!----------------------------------------------------------------------
\file tutorial_nln_truss.H

\brief student's c++/baci nonlinear truss tutorial

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/

/// All header files should have #define guards to prevent multiple inclusion.
/// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef TUTORIAL_NONLINEAR_TRUSS_H
#define TUTORIAL_NONLINEAR_TRUSS_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_globalproblem.H"

//////////////////////////////////////////////////////
///  DEFINE THE NUMBER OF ELEMENTS
#define numele 6
//////////////////////////////////////////////////////


namespace TUTORIAL {

/// FORWARD DECLARATION
/// Here we declare the class 'TutorialElement'.
/// This is needed because input arguments of this type
/// are declared in this .H-file. We could also include
/// the tutorial_ele.H-file but it should be avoided
/// to include headers in headers. This increases
/// dependencies in the code and slows down compilation
/// and linking.
namespace ELEMENTS {
class TutorialElement;
}

class NonlinearTruss
{

public:

  /*!
  \brief constructor (ctor)
  */
  NonlinearTruss();

  /*!
  \brief destructor (dtor)

  \note This Destructor is declared virtual, since class NonlinearTruss contains virtual methods.
        This is an indicator, that this class is supposed to be part of some inheritance line.
        Check the internet on why it is then important to declare the destructor virtual.
  */
  virtual ~NonlinearTruss(){};

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

  /*!
  \brief evaluate rhs and stiffness matrix

  \param rhs (out)       : right-hand side vector to be evaluated
  \param stiff (out)     : system (stiffness-) matrix to be evaluated
  \param eval_rhs (in)   : bool indicating whether rhs should be evaluated
  \param eval_stiff (in) : bool indicating whether stiff should be evaluated

  */
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff, bool eval_rhs, bool eval_stiff);


  /*!
  \brief evaluate rhs and stiffness matrix

  \note evaluates rhs and system matrix by default

  \param rhs (out)       : right-hand side vector to be evaluated
  \param stiff (out)     : system (stiffness-) matrix to be evaluated

  */
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff)
  {Evaluate(rhs, stiff,true, true); return;};

  /*!
  \brief evaluate rhs and stiffness matrix

  \note evaluates system matrix by default

  \param rhs (out)       : right-hand side vector to be evaluated
  \param stiff (out)     : system (stiffness-) matrix to be evaluated
  \param eval_rhs (in)   : bool indicating whether rhs should be evaluated

  */
  void Evaluate(LINALG::Matrix<numele+1,1>* rhs, LINALG::Matrix<numele+1,numele+1>* stiff, bool eval_rhs)
  {Evaluate(rhs, stiff, eval_rhs, true); return;};

  /*!
  \brief Post evaluate operator ( e.g. for manipulation of rhs and stiffness matrix after usual evaluate)

  \note Supposed to be overriden by specialized nonlinear truss problem, which inherits from this class.
        Here, the PostEvaluate() operator does nothing.

  \param rhs (out)       : right-hand side vector to be evaluated
  \param stiff (out)     : system (stiffness-) matrix to be evaluated
  \param eval_rhs (in)   : bool indicating whether rhs should be evaluated
  \param eval_stiff (in) : bool indicating whether stiff should be evaluated

  */
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
  void PrintProblem(std::vector<LINALG::Matrix<numele+1,1> >& conditions);

  /// print used material law
  void PrintMaterial();

  /// print the elements and their members (output example for std::map)
  void PrintDiscretization();

  /*!
  \brief assemble an element stiffness matrix into the global stiffness matrix

  \param ele (in)           : element whose contribution is assembled
  \param ele_stiff_mat (in) : element matrix to be assembled
  \param glob_stiff (out)   : global matrix into which we assemble

  */
  void AssembleMatrix(TUTORIAL::ELEMENTS::TutorialElement* ele, LINALG::Matrix<2,2>* ele_stiff_mat, LINALG::Matrix<numele+1,numele+1>* glob_stiff);

  /*!
  \brief assemble an element residual vector into the global right-hand-side

  \param ele (in)         : element whose contribution is assembled
  \param ele_rhs_vec (in) : element vector to be assembled
  \param glob_rhs (out)   : global vector into which we assemble

  */
  void AssembleVector(TUTORIAL::ELEMENTS::TutorialElement* ele,
                      LINALG::Matrix<2,1>* ele_rhs_vec,
                      LINALG::Matrix<numele+1,1>* glob_rhs);

  /// discretization map ele id => TutorialElement
  std::map<int, ELEMENTS::TutorialElement* >* discretization_;

  /// pointer to dirichlet condition
  LINALG::Matrix<numele+1,1>* dirichletvalues_;

  /// pointer to neumann condition
  LINALG::Matrix<numele+1,1>* neumannvalues_;

  /// the solver object
  Teuchos::RCP<LINALG::FixedSizeSerialDenseSolver<numele+1,numele+1,1> > nxnsolver_;

  //! @name Vectors and Matrices
  //@{
  LINALG::Matrix<numele+1,1> rhs_; //!< RHS vector
  LINALG::Matrix<numele+1,numele+1> stiff_; //!< stiffness matrix
  LINALG::Matrix<numele+1,1> inc_; //!< increment vector
  LINALG::Matrix<numele+1,1> disp_; //!< displacement vector
  LINALG::Matrix<numele+1,1> freact_; //!< reaction force vector
  //@}

};// class NonlinearTruss
} // namespace TUTORIAL

#endif
