/*----------------------------------------------------------------------*/
/*!
 \file optimizer_base.H
 \brief optimization algorithm base class

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef OPTIMIZER_BASE_H
#define OPTIMIZER_BASE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "invana_writer.H"

// forward declarations
namespace Teuchos
{
  class ParameterList;
}

namespace IO
{
  class InputControl;
  class DiscretizationWriter;
}  // namespace IO

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

  It's main purpose is the access to the optimization problem. It offers the evaluation of the value
  and the gradient of the optimization problem to be used by the derived implementations. So these,
  e.g. OptimizerBase::OptimizerGradDesc, OptimizerBase::OptimizerLBFGS or
  OptimizerBase::OptimizerMC, are completely decoupled from the optimization problem to be solved.

  \author kehl
  */

  class OptimizerBase
  {
   public:
    //! @name Construction/destruction
    //@{

    //! standard constructor
    OptimizerBase(const Teuchos::ParameterList& invp);

    //! destructor
    virtual ~OptimizerBase() {}

    //! Set member variables after construction.
    void Init(Teuchos::RCP<InvanaBase> optprob  // [in] a reference to an optimization problem
    );

    //@}


    //! @name Initialize the execution of the algorithm
    //@{

    //! setup up algorithm specific stuff
    virtual void Setup() = 0;

    //! use some specific optimization routine
    virtual void Integrate() = 0;

    //! read restart information
    virtual void ReadRestart(int run) = 0;

    //! write restart information
    virtual void WriteRestart() = 0;

    //@}


    //! @name Write access
    //@{

    //! Access to solution
    Teuchos::RCP<Epetra_MultiVector> GetSolution() { return sol_; }

    //! Access to gradient
    Teuchos::RCP<Epetra_MultiVector> GetGradient() { return objgrad_; }

    //! Access to objective function value
    double* GetObjFunctVal() { return &objval_; }

    //@}


    //! @name Read only access
    //@{

    //! get reference to solution vector
    const Epetra_MultiVector& GetSolutionView() { return *sol_; }

    //! get reference to the old solution vector
    const Epetra_MultiVector& GetSolutionOldView() { return *sol_o_; }

    //! get reference to the gradient vector
    const Epetra_MultiVector& GetGradientView() { return *objgrad_; }

    //! get reference to the old gradient vector
    const Epetra_MultiVector& GetGradientOldView() { return *objgrad_o_; }

    //! get the objective function value
    double GetObjFunctValView() { return objval_; }

    //! get the old objective function value
    double GetObjFunctValOldView() { return objval_o_; }

    //! return the unique map of the vector layout
    const Epetra_Map& SolLayoutMap() { return *solrowmap_; }

    //! current run counter
    int Runc() { return runc_; }

    //@}


   protected:
    //! @name General purpose functions for nonlinear solution
    //@{

    //! evaluate the objective function
    void Evaluate(double* val,                     /// [out] value of the objective function
        Teuchos::RCP<Epetra_MultiVector> gradient  /// [out] gradient of the objective function
    );

    //! upate gradient: push current to old
    void UpdateGradient() { objgrad_o_->Update(1.0, *objgrad_, 0.0); }

    //! upate objective function value: push current to old
    void UpdateObjFunctValue() { objval_o_ = objval_; }

    //! update the solution (and push current to old)
    void UpdateSolution(const Epetra_MultiVector&
            toadd  /// [in] (properly scaled) step to be added to the solution vector
    );

    //! undo the last UpdateSolution()-call
    void UndoUpdateSolution();

    //@}

    //! @name General query
    //@{

    //! query if base class is initialized
    bool IsInit() { return isinit_; }

    //! a reference to the underlying optimization problem
    Teuchos::RCP<InvanaBase> OptProb() { return optprob_; }

    //@}


    //! @name Access to input and output
    //@{

    //! the input file parameters
    const Teuchos::ParameterList& Inpar() { return inpar_; }

    //! a reference to the file to restart from
    const Teuchos::RCP<IO::InputControl> RestartFromFile() { return inputfile_; }

    //! a reference to the discretization writer
    const Teuchos::RCP<InvanaWriter> Writer() { return writer_; }

    //@}


    //! @name interesting members running an iterative optimization algorithm
    //@{

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
    //! @name current and old state of the algorithm iterations
    //@{

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

    //@}

    //! @name Problem to be solved
    //@{
    //! a reference to the underlying optimization problem to be solved
    Teuchos::RCP<InvanaBase> optprob_;

    //! unique vector layout map
    Teuchos::RCP<Epetra_Map> solrowmap_;

    //@}

    //! @name Input and output
    //@{

    //! Input file parameters
    const Teuchos::ParameterList& inpar_;

    //! output control
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! An inverse analysis specific wrapper for the output
    Teuchos::RCP<InvanaWriter> writer_;

    //! the file to start from in case of restarts
    Teuchos::RCP<IO::InputControl> inputfile_;

    //@}

    //! initialization bool
    bool isinit_;
  };
}  // namespace INVANA

#endif /*OPTIMIZER_BASE_H_*/
