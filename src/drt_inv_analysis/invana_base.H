/*----------------------------------------------------------------------*/
/*!
 * \file invana_base.H
 * \brief inverse analysis base class
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_BASE_H
#define INVANA_BASE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

// forward declarations
namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
  class InputControl;
  class OutputControl;
}

namespace DRT
{
  class ResultTest;
  class Discretization;
}

namespace STR
{

namespace INVANA
{
class ObjectiveFunct;
class MatParManager;
class RegularizationBase;
class OptimizerBase;

  /*!
  \brief A base class for the inverse analysis
  \author kehl

  This class is a base class to setup an optimzation problem. It's main purpose is to hold
  all base class pointers to all the functional classes used to make up an optimization
  problem, i.e.:\n
  - the discretization (more for convenience)\n
  - an objective function\n
  - a parameterization of the material parameters to be estimated\n
  - a regularization to smooth the objective function\n
  - an optimization algorithm\n

  It triggers the optimization algorithm to start the optimization and it provides and interface function
  Evaluate(sol, val, gradient) to be used by the optimization algorithm to query the optimization problem
  for the actual value "val" and gradient "gradient" of the objective function given the current solution "sol".
  Since this highly depends on the actual implementation this is pure virtual here.
 */

class InvanaBase {

public:

  //! standard constructor
  InvanaBase();

  //! destructor
  virtual ~InvanaBase() {};

  //! Set member variables after construction.
  /*!
    \param In
    discret - a reference to the discretization.
    \param In
    objfunct - a reference to an objective function.
    \param In
    matman - a reference to the material parameterization.
    \param In
    regman - a reference to a regularization strategy for the objective function.
    \param In
    optimizer - a reference to an optimization algorithm.
    \param In
    optprob - a reference to itself to be passed to the optimization algorithm.
  */
  void Init(Teuchos::RCP<DRT::Discretization> discret,
            Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct,
            Teuchos::RCP<STR::INVANA::MatParManager> matman,
            Teuchos::RCP<STR::INVANA::RegularizationBase> regman,
            Teuchos::RCP<STR::INVANA::OptimizerBase> optimizer,
            Teuchos::RCP<STR::INVANA::InvanaBase> optprob);

  //! Setup of the derived classes.
  virtual void Setup() = 0;

  //! Solve the optimization problem using the optimization algorithm kept as member
  void Solve(int restart);

  //! Evaluate value and gradient of the objective function defined in the derived class
  /*!
    \param In
    sol - the current solution of the optimization problem
    \param Out
    val - the value of the objective function
    \param Out
    optprob - the gradient of the objective function wrt to the optimization parameters
  */
  virtual void Evaluate(const Epetra_MultiVector& sol, double* val, Teuchos::RCP<Epetra_MultiVector> gradient) = 0;

  //! a reference to the discretization
  Teuchos::RCP<DRT::Discretization> Discret() {return discret_;}

  //! a reference to the objective function
  Teuchos::RCP<STR::INVANA::ObjectiveFunct> ObjectiveFunct() {return objfunct_;}

  //! a reference to the material parameter manager
  Teuchos::RCP<STR::INVANA::MatParManager> Matman(){return matman_; }

  //! a reference to the regularization of the objective function
  Teuchos::RCP<STR::INVANA::RegularizationBase> Regman(){return regman_; }

  //! a reference to the optimization algorithm
  Teuchos::RCP<STR::INVANA::OptimizerBase> Optimizer(){return optimizer_; }

  //! initial guess to be used by the optimizer
  const Epetra_MultiVector& InitialGuess();

  //! a reference to the programs communicator
  const Epetra_Comm& Comm();

  //! the unique layout of the optimization parameters
  Teuchos::RCP<Epetra_Map> VectorRowLayout();

  //! the possibly non-unique layout of the optimization parameters
  Teuchos::RCP<Epetra_Map> VectorColLayout();

  //! number of vectors in the optimization parameter Epetra_MultiVector
  double NumVectors();

  //! return invana specific result test to the global problem field tests
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

private:

  //! a reference to the discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! objective function to be optimized
  Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct_;

  //! handler for calls to the parametrization
  Teuchos::RCP<STR::INVANA::MatParManager> matman_;

  //! handler for calls to the regularization
  Teuchos::RCP<STR::INVANA::RegularizationBase> regman_;

  //! optimization algorithm to solve the optimization problem
  Teuchos::RCP<STR::INVANA::OptimizerBase> optimizer_;

  //! initialization bool
  bool isinit_;



};
}  // namespace INVANA
}  // namespace STR

#endif /*INVANA_BASE_H_*/
