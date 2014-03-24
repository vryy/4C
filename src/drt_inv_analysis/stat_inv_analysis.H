/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANALYSIS_H
#define STAT_INV_ANALYSIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

// forward declarations

namespace COMM_UTILS
{
  class NestedParGroup;
}

namespace STR
{
  class TimInt;
}

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
  class InputControl;
  class OutputControl;
}

namespace LINALG
{
  class Solver;
}

namespace DRT
{
  class Condition;
  class Discretization;
  class ResultTest;
  namespace UTILS
  {
    template <typename> class TimIntMStep;
  }
}

namespace STR
{

namespace INVANA
{
class ObjectiveFunct;

class MatParManager;

  /*!
  \brief A class for the inverse analysis

  This class provides all the general stuff needed for some inverse analysis procedure, which is:

  - access to the forward problem time integrator
  - access to the adjoint problem time integrator
  - methods to run the forward problem
  - methods to run the adjoint problem
  - access to some objective function
  - access to the gradient of some objective function
  - an objective function to be minimized or whatever
  - a MatParManager that has control of the material parameters
  - final testing of the optimized material parmaters

  What exactly is to be done with this information must be provided by derived classes. This can be
  implementations like steepest descent, limited memory bfgs, monto carlo sampling, etc.

  A simple Thikonov regularization is included, when specified in the input file. This should be moved
  some standalone class eventually.
 */

class StatInvAnalysis {

protected:

  /// discretisation
  Teuchos::RCP<DRT::Discretization> discret_;

  //! value of the objective function
  double  objval_;

  //! value of the objective function of the prevous step
  double objval_o_;

  //! increment of the objective function
  double error_incr_;

  //! write restart every so often
  int restartevry_;

  //! solve the primal problem
  void SolveForwardProblem();

  //! solve the adjoint problem
  void SolveAdjointProblem();

  //! Evaluate the objective functions gradient w.r.t the material parameters
  void EvaluateGradient();

  //! Upate Gradient: push current to old
  void UpdateGradient(){objgrad_o_->Update(1.0, *objgrad_, 0.0);};

  //! Evaluate the objective function
  void EvaluateError();

  //! Evaluate FD approximation of the objective functions gradient w.r.t the design parameters
  void EvaluateGradientFD();

  //! Reset discretization
  void ResetDiscretization();

  //! a reference to the material parameter manager
  const Teuchos::RCP<STR::INVANA::MatParManager> Matman(){return matman_; };

  //! a reference to the discretization writer
  const Teuchos::RCP<IO::DiscretizationWriter> Writer(){return output_;};

  //! a reference to the file to restart from
  const Teuchos::RCP<IO::InputControl> RestartFromFile(){return inputfile_;};

  //! get current gradient
  const Teuchos::RCP<Epetra_MultiVector> GetGradient(){return objgrad_;};

  //! get current gradient
  const Teuchos::RCP<Epetra_MultiVector> GetGradientOld(){return objgrad_o_;};

private:

  //! a dofrowmap of the discretization
  const Epetra_Map* dofrowmap_;

  //! primal variables (displacements)
  Teuchos::RCP<Epetra_MultiVector> dis_;

  //! dual variables (lagrange multiplier)
  Teuchos::RCP<Epetra_MultiVector> disdual_;

  //! time of steps of the forward problem
  Teuchos::RCP<std::vector<double> > time_;

  //! Objective function to be optimized
  Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct_;

  //! Handler for calls to material routines
  Teuchos::RCP<STR::INVANA::MatParManager> matman_;

  //! number of simulation steps in the primal problem
  int msteps_;

  //! gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_;

  //! previous gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_o_;

  //! outputer for the inverse analysis
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! the file to start from in case of restarts
  Teuchos::RCP<IO::InputControl> inputfile_;

  //! have regularization?
  bool havereg_;

  //! weight of the regularization when being added to the objective function
  double regweight_;

  //! MStep EpetraVector to EpetraMultiVector
  void MStepEpetraToEpetraMulti(Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > mstepvec,
                                Teuchos::RCP<Epetra_MultiVector> multivec);

  //! Mstep double to std::vector<double>
  void MStepDToStdVecD(Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > mstepvec,
                       Teuchos::RCP<std::vector<double> > stdvec);

  // DEBUG PRINT OUT
  void PrintDataToScreen(Epetra_MultiVector& vec);

public:

  //! standard constructor
  StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis);

  //! destructor
  virtual ~StatInvAnalysis() {};

  //! use some specific optimization routine
  virtual void Optimize() = 0;

  //! read restart information
  virtual void ReadRestart(int run);

  //! write restart information
  virtual void WriteRestart();

  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_;};

  Teuchos::RCP<STR::INVANA::MatParManager> MatParManager() {return matman_;};

  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  // return the value of the objective function
  double GetError(){return objval_;};

  // return the value of the gradient 2-norm
  double GetGrad2Norm();


}; // class StatInvAnalysis
}  // namespace INVANA
}  // namespace STR

#endif /*STAT_INV_ANALYSIS_H_*/
