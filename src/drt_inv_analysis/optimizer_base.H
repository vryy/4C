/*----------------------------------------------------------------------*/
/*!
 \file optimizer_base.H
 \brief optimization algorithm base class

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef OPTIMIZER_BASE_H
#define OPTIMIZER_BASE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

// forward declarations
namespace Teuchos
{
class ParameterList;
}

namespace IO
{
class InputControl;
class DiscretizationWriter;
}

namespace DRT
{
class ResultTest;
}

namespace INVANA
{

class InvanaBase;

  /*!
  \brief A base optimization algorithm

  This class hold the quantities common to all kind of optimization algorithms, i.e.:\n
  - a reference to an optimization problem
  - the value of some objective function\n
  - the gradient of some opbjective function\n
  - the current and old solution vector to the optimization problem\n

  It's main purpose is the access to the optimization problem. It offers the evaluation of the value and
  the gradient of the optimization problem to be used by the derived implementations. So these, e.g.
  OptimizerBase::OptimizerGradDesc, OptimizerBase::OptimizerLBFGS or OptimizerBase::OptimizerMC, are
  completely decoupled from the optimization problem to be solved.

  \author kehl
  */

class OptimizerBase {

public:

  //! standard constructor
  OptimizerBase(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~OptimizerBase() {}

  //! Set member variables after construction.
  /*!
    \param In
    layoutrowmap - the unique layout of the optimization parameters.
    \param In
    layoutcolmap - the possibly non-unique layout of the optimization parameters
    \param In
    numvecs - number of vectors of the optimization parameters layout
    \param In
    optprob - a reference to the optimization problem for call to evaluate value and gradient of some objective function
  */
  void Init(Teuchos::RCP<Epetra_Map> layoutrowmap, Teuchos::RCP<Epetra_Map> layoutcolmap, int numvecs, Teuchos::RCP<InvanaBase> optprob);

  /*----------------------------------------------------------------------*/
  //! functions to be provided by derived classes

  //! setup up algorithm specific stuff
  virtual void Setup() = 0;

  //! use some specific optimization routine
  virtual void Integrate() = 0;

  //! read restart information
  virtual void ReadRestart(int run) = 0;

  //! write restart information
  virtual void WriteRestart() = 0;

  /*----------------------------------------------------------------------*/
  //! write access functions

  //! get reference to soution vector for write access
  Teuchos::RCP<Epetra_MultiVector> GetSolution() const {return sol_;}

  Teuchos::RCP<Epetra_MultiVector> GetGradient() const {return objgrad_;}

  //! get write access to objective function value
  double* GetObjFunctVal() {return &objval_;}

  /*----------------------------------------------------------------------*/
  //! read only access functions

  //! get reference to solution vector
  const Epetra_MultiVector& GetSolutionView() const {return *sol_;}

  //! get reference to the old solution vector
  const Epetra_MultiVector& GetSolutionOldView() const {return *sol_o_;}

  //! get reference to the gradient vector
  const Epetra_MultiVector& GetGradientView() const {return *objgrad_;}

  //! get reference to the old gradient vector
  const Epetra_MultiVector& GetGradientOldView() const {return *objgrad_o_;}

  //! get "reference" to objective function value
  double GetObjFunctValView() const {return objval_;}

  //! get "reference" to old objective function value
  double GetObjFunctValOldView() const {return objval_o_;}

  //! return the map of the vector layout
  const Epetra_Map& SolLayoutMap() const {return *solcolmap_;}

  //! return the unique map of the vector layout
  const Epetra_Map& SolLayoutMapUnique() const {return *solrowmap_;}

  //! current run counter
  int Runc() const {return runc_;}

  /*----------------------------------------------------------------------*/
  //! general access functions

  //! a reference to the underlying optimization problem
  Teuchos::RCP<InvanaBase> OptProb() const {return optprob_;}

protected:

  /*----------------------------------------------------------------------*/
  //! perform some sort of optimization using these functions

  //! get inital guess from underlying optimization problem
  void SetInitialGuess();

  /*! \brieff Evaluate the objective function
   *
   * \param Out
   * val - to be filled with the value of the objective function.
   * \param Out
   * gradient - to be filled with the gradient of the objective function
   */
  void Evaluate(double* val, Teuchos::RCP<Epetra_MultiVector> gradient);

  //! upate Gradient: push current to old
  void UpdateGradient(){objgrad_o_->Update(1.0, *objgrad_, 0.0);}

  //! upate objective function value: push current to old
  void UpdateObjFunctValue(){objval_o_=objval_;}

  /*! \brief update the solution (and push current to old)
   *
   * \param In
   * toadd - (properly scaled) step to be added to the solution vector
   */
  void UpdateSolution(const Epetra_MultiVector& toadd);

  //! undo the last UpdateSolution()-call
  void UndoUpdateSolution();

  //! question if base class is initialized
  bool IsInit(){return isinit_;}


  //! @name access to input and output
  //@{

  //! a reference to the file to restart from
  const Teuchos::RCP<IO::InputControl> RestartFromFile(){return inputfile_;}

  //! a reference to the discretization writer
  const Teuchos::RCP<IO::DiscretizationWriter> Writer(){return output_;}

  //@}


  //! @name members directly accessible from derived classes
  //@{

  //! number of vectors in the multivector
  int numvecs_;

  //! increment of the objective function
  double error_incr_;

  //! write restart every so often
  const int restartevry_;

  //! maximum number of iterations
  const int maxiter_;

  //! step size
  const double stepsize_;

  //! iteration counter
  int runc_;

  //! tolerance for convergence criterion
  double convtol_;

  //! current convergence criterion value
  double convcritc_;

  //@}

private:

  /*----------------------------------------------------------------------*/
  //! members accessible via access functions

  //! current value of the objective function
  double objval_;

  //! value of the objective function of the previous step
  double objval_o_;

  //! current gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_;

  //! gradient of the objective function of the previous step
  Teuchos::RCP<Epetra_MultiVector> objgrad_o_;

  //! current solution
  Teuchos::RCP<Epetra_MultiVector> sol_;

  //! solution of the previous step
  Teuchos::RCP<Epetra_MultiVector> sol_o_;

  //! a reference to the underlying optimization problem to be solved
  Teuchos::RCP<InvanaBase> optprob_;

  //! unique vector layout map
  Teuchos::RCP<Epetra_Map> solrowmap_;

  //! possibly non-unique vector layout map
  Teuchos::RCP<Epetra_Map> solcolmap_;

  //! outputer for the inverse analysis
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! the file to start from in case of restarts
  Teuchos::RCP<IO::InputControl> inputfile_;

  //! initialization bool
  bool isinit_;

};
}  // namespace INVANA

#endif /*OPTIMIZER_BASE_H_*/
