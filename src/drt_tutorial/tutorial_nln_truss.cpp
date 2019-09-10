/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci nonlinear truss tutorial

\maintainer  Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/

#include "tutorial_utils.H"
#include "tutorial_nln_truss.H"
#include "tutorial_ele/tutorial_ele_calc.H"
#include "../drt_tutorial/tutorial_ele/tutorial_ele.H"
#include "../drt_tutorial/tutorial_material/tutorial_mat.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_fixedsizematrix.H"


/// ctor
TUTORIAL::NonlinearTruss::NonlinearTruss()
{
  /// initialize global vectors and matrices
  rhs_.PutScalar(0.0);
  stiff_.PutScalar(0.0);
  inc_.PutScalar(0.0);
  freact_.PutScalar(0.0);
  disp_.PutScalar(0.0);

  /// initialize pointers with NULL
  discretization_ = NULL;
  neumannvalues_ = NULL;
  dirichletvalues_ = NULL;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::PrintTutorialType()
{
  std::cout << "\n YOU CHOSE THE NONLINEAR TRUSS TOY PROBLEM ! \n\n" << std::endl;

  return;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
std::vector<LINALG::Matrix<7, 1>> TUTORIAL::NonlinearTruss::ProblemDefinition()
{
  /// print tutorial type
  PrintTutorialType();

  /**************************************************************/
  /*              DEFINE YOUR VALUES                            */

  // GEOMETRY
  //
  //    o-------o-------o-------o-------o--~
  //

  // Define number of elements
  //
  // Do this in the header of this file. In line 20 you find
  // the global definition of numele which equals 6 as default value.
  // this needs to be done, since the size of LINALG::Matrix needs to
  // known at compile time. You cannot put a variable as template argument.
  // try for example to declare a LINALG::Matrix like:
  //
  // int numrows = 6;
  // LINALG::MATRIX<numrows,1> mytestmatrix(true);
  //
  // You will see, that this results in an compile error !
  //
  // A remedy could be to use matrix objects other than LINALG::Matrix,
  // which allow a dynamic adaption of their size.
  //

  // BOUNDARY CONDITIONS

  // initialize variables
  LINALG::Matrix<numele + 1, 1> dirichletcond(true);
  dirichletcond.PutScalar(-1234.0);  // initialized with -1234.0
  LINALG::Matrix<numele + 1, 1> neumanncond(true);
  neumanncond.PutScalar(-1234.0);  // initialized with -1234.0

  // DIRICHLET conditions
  dirichletcond(0) = 0.0;
  dirichletcond(6) = 0.0;

  // NEUMANN conditions
  neumanncond(1) = 1.5;
  neumanncond(5) = -1.5;

  // NODE LOCATIONS

  // material configuration
  LINALG::Matrix<numele + 1, 1> X(true);
  X(0) = 0.0;
  X(1) = 1.0;
  X(2) = 2.0;
  X(3) = 3.0;
  X(4) = 4.0;
  X(5) = 5.0;
  X(6) = 6.0;


  // MATERIAL DEFINITION

  // Young's Modulus
  LINALG::Matrix<numele, 1> E(true);
  E.PutScalar(900.0);

  /**************************************************************/
  // DO NOT TOUCH

  // save condition matices in vector and return this vector
  std::vector<LINALG::Matrix<numele + 1, 1>> vector_of_condition_matrices;
  vector_of_condition_matrices.push_back(dirichletcond);
  vector_of_condition_matrices.push_back(neumanncond);

  /**************************************************************/
  // DO NOT TOUCH

  // build the discretization consisting 'numele' elements,
  // with material parameters E,
  // and with initial coordinates X
  discretization_ = TUTORIAL::UTILS::BuildDiscretization(numele, E, X);

  /**************************************************************/

  // good bye! And take this!
  return vector_of_condition_matrices;

}  // GeometryDefinition


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::SetupProblem()
{
  // define geometry and bc's and save conditions in vector 'vector_of_condition_matrices'
  static std::vector<LINALG::Matrix<numele + 1, 1>> vector_of_condition_matrices =
      ProblemDefinition();


  // extract dirichlet (position 0) and neumann (position 1) matrices from vector
  //
  // The data are stored in 'vector_of_condition_matrices', and we just want to
  // extract the dirichlet and neumann parts for convenience.
  // Therefore we do not copy the data, we rather have pointers pointing to the data !
  dirichletvalues_ = &(vector_of_condition_matrices[0]);
  neumannvalues_ = &(vector_of_condition_matrices[1]);

  // build the solver
  SetupSolver<numele + 1, numele + 1>();

  // print important information on the problem
  PrintProblem(vector_of_condition_matrices);
  PrintDiscretization();
  PrintMaterial();
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::TimeLoop()
{
  double dt = 1.0;
  double endtime = 1.0;
  double time = dt;

  int step = 1;

  // loop over
  while (time <= endtime)
  {
    std::cout << "\nTIMESTEP " << step << "/" << endtime / dt << " time=" << time << std::endl;

    // update step
    step++;

    // equilibrium iterations for this time step
    Newton();

    // update time n->n+1
    time = time + dt;
  }

}  // TimeLoop


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
template <unsigned int numrow, unsigned int numcol>
void TUTORIAL::NonlinearTruss::SetupSolver()
{
  nxnsolver_ = Teuchos::rcp(new LINALG::FixedSizeSerialDenseSolver<numrow, numcol, 1>());
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::PrintProblem(std::vector<LINALG::Matrix<7, 1>>& conditions)
{
  std::cout << "This is your Geometry and your Boundary Conditions:\n" << std::endl;
  std::cout << "D := Dirichlet Node" << std::endl;
  std::cout << "N := Neumann Node" << std::endl;
  std::cout << "\n" << std::endl;

  for (int c = 0; c < numele + 1; c++)
  {
    if ((conditions[0])(c) > -1234.0)
      std::cout << "D"
                << "      ";
    else if (conditions[1](c) > -1234.0)
      std::cout << "N"
                << "      ";
    else
      std::cout << " "
                << "      ";
  }
  std::cout << "" << std::endl;

  for (int c = 0; c < numele; c++) std::cout << "o------";
  std::cout << "o" << std::endl;
  for (int c = 0; c < numele + 1; c++) std::cout << c << "      ";
  std::cout << "\n" << std::endl;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::PrintMaterial()
{
  std::cout << "Material St.Venant-Kirchhoff with nu=" << 0.0 << std::endl;
  std::cout << "Strain Measure: Green-Lagrange" << std::endl;
  std::cout << "\n" << std::endl;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::PrintDiscretization()
{
  std::map<int, TUTORIAL::ELEMENTS::TutorialElement*>::iterator it = discretization_->begin();
  std::cout << "The Discretization contains the following elements:\n";
  for (it = discretization_->begin(); it != discretization_->end(); ++it)
    std::cout << "ID " << (it->second)->Id() << " => "
              << " numnodes=" << (it->second)->NumNode()
              << " Young's Modulus=" << (it->second)->Material()->YoungModulus()
              << " Dofs=" << (it->second)->Dof(0) << " "
              << " Dofs=" << (it->second)->Dof(1) << '\n';
  std::cout << "\n" << std::endl;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::Newton()
{
  // initialize converged flag to be unconverged
  bool converged = false;
  // initialize iteration counter
  int iter = 0;
  // get max. number of iterations from input file
  int maxiter = 30;
  // get convergence tolerance from input file
  double tol = 1e-11;


  // newton equilibrium loop
  while (not converged and iter < maxiter)
  {
    // i -> i+1
    iter++;

    // evaluate stiffness matrix and rhs with displacements from previous solve d^{n+1}_{i}
    stiff_.Clear();
    rhs_.Clear();
    Evaluate(&rhs_, &stiff_);
    PostEvaluate(rhs_, stiff_, false, false);

    // do bcs
    DoBoundaryConditions();

    // solve for new displacmements d^{n+1}_{i+1}
    Solve(rhs_, stiff_, inc_);

    // update new displacements
    disp_.Update(1.0, inc_, 1.0);

    // evaluate stiffness rhs with new displacements d^{n+1}_{i} for convergence check
    rhs_.Clear();
    Evaluate(&rhs_, &stiff_, true, false);
    PostEvaluate(rhs_, stiff_, false, false);

    // do bcs
    DoBoundaryConditions();

    // perform convergence check
    converged = ConvergenceCheck(tol);

    // print convergence over iterations
    std::cout << "iter " << iter << ": res=" << std::setprecision(7) << rhs_.Norm2()
              << "  inc=" << std::setprecision(7) << inc_.Norm2() << " -> converged=" << converged
              << std::endl;
  }

  if (not converged) dserror("Newton did not converge in %d iterations.", maxiter);
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::Evaluate(LINALG::Matrix<numele + 1, 1>* rhs,
    LINALG::Matrix<numele + 1, numele + 1>* stiff, bool eval_rhs, bool eval_stiff)
{
  // initialize element stiffness matrix that is filled during element evaluation.
  // Static objects are built only once. The second time we arrive here, this line is skipped.
  // Thus, the memory for this matrix is allocated only once.
  static LINALG::Matrix<2, 2> ele_stiff(true);

  // initialize element residual vector that is filled during element evaluation.
  static LINALG::Matrix<2, 1> ele_rhs(true);

  // get iterator to first map entry of discretization
  std::map<int, TUTORIAL::ELEMENTS::TutorialElement*>::iterator ele = discretization_->begin();

  // loop over all elements in discretization and call Evaluate() on element level
  for (ele = discretization_->begin(); ele != discretization_->end(); ++ele)
  {
    // call the element evaluate
    ele->second->Evaluate(&disp_, &ele_stiff, &ele_rhs, eval_stiff, eval_rhs);

    // assemble the element matrix into the global matrix
    if (eval_stiff) AssembleMatrix(ele->second, &ele_stiff, stiff);

    // assemble the element vector into the global vector
    if (eval_rhs) AssembleVector(ele->second, &ele_rhs, rhs);

    // clear the element matrices to prepare them for the next element evaluation
    ele_stiff.Clear();
    ele_rhs.Clear();
  }  // loop over all elements in discretization


  return;
}  // TUTORIAL::NonlinearTruss::Evaluate



/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::Solve(LINALG::Matrix<numele + 1, 1> rhs,
    LINALG::Matrix<numele + 1, numele + 1> stiff, LINALG::Matrix<numele + 1, 1> inc)
{
  // solve

  // clear increment ( just feels good :-) )
  inc_.Clear();
  // scale rhs with -1.0 (this is what the Newton scheme demands if you do the correct math!)
  rhs_.Scale(-1.0);
  // set matrix and vectors in nxnsolver
  Solver()->SetMatrix(stiff_);
  Solver()->SetVectors(inc_, rhs_);

  // finally we solve
  int err = Solver()->Solve();

  // any errors returned from the solver?
  if (err != 0) dserror("Solver returned err=%d", err);
}

/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
bool TUTORIAL::NonlinearTruss::ConvergenceCheck(double tol)
{
  if (rhs_.Norm2() <= tol)
    return true;
  else
    return false;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::DoBoundaryConditions()
{
  // apply the neumann conditions
  TUTORIAL::UTILS::DoNeumannCond(&rhs_, neumannvalues_);
  // apply the dirichlet conditions
  TUTORIAL::UTILS::DoDirichletCond(&stiff_, &rhs_, &inc_, &disp_, &freact_, dirichletvalues_);
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::PrintResults()
{
  std::cout << "\nRESULTS :\n" << std::endl;
  std::cout << "disp=" << std::setprecision(10) << disp_ << std::endl;
  std::cout << "freact=" << std::setprecision(10) << freact_ << std::endl;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::TutorialDone()
{
  std::map<int, TUTORIAL::ELEMENTS::TutorialElement*>::iterator it = discretization_->begin();

  for (it = discretization_->begin(); it != discretization_->end(); ++it) delete it->second;

  // also delete the calc instance
  TUTORIAL::ELEMENTS::TutorialEleCalc::Instance(false);
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::AssembleMatrix(TUTORIAL::ELEMENTS::TutorialElement* ele,
    LINALG::Matrix<2, 2>* ele_stiff_mat, LINALG::Matrix<numele + 1, numele + 1>* glob_stiff)
{
  // get number of nodes from element
  int numnode = ele->NumNode();
  int numdof = ele->NumDof();

  if (numnode != numdof)
    dserror(
        "Number of dofs per element != number of nodes per element!\n "
        "This is not supported for Truss TutorialElements!");

  // write the element matrix entries to the corresponding position in the global matrix
  // loop over nodes
  for (int rowdof = 0; rowdof < numdof; rowdof++)
  {
    for (int coldof = 0; coldof < numdof; coldof++)
    {
      // get the dof gid at this nodes
      int globalrow = ele->Dof(rowdof);
      int globalcol = ele->Dof(coldof);
      (*glob_stiff)(globalrow, globalcol) += (*ele_stiff_mat)(rowdof, coldof);
    }
  }

  return;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::NonlinearTruss::AssembleVector(TUTORIAL::ELEMENTS::TutorialElement* ele,
    LINALG::Matrix<2, 1>* ele_rhs_vec, LINALG::Matrix<numele + 1, 1>* glob_rhs)
{
  // get number of nodes from element
  int numnode = ele->NumNode();
  int numdof = ele->NumDof();

  if (numnode != numdof)
    dserror(
        "Number of dofs per element != number of nodes per element!\n "
        "This is not supported for Truss TutorialElements!");

  // write the element vector entries to the corresponding position in the global matrix
  // loop over nodes
  for (int rowdof = 0; rowdof < numdof; rowdof++)
  {
    // get the dof gid at this nodes
    int globalrow = ele->Dof(rowdof);
    (*glob_rhs)(globalrow) += (*ele_rhs_vec)(rowdof);
  }

  return;
}
