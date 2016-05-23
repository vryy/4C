/*----------------------------------------------------------------------*/
/*!
\file invana_base.H

\brief inverse analysis base class

<pre>
\level 3
\maintainer Sebastian Kehl
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
namespace DRT
{
  class ResultTest;
  class Discretization;
}

namespace INVANA
{
class ObjectiveFunct;
class MatParManager;
class RegularizationBase;

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
 */

class InvanaBase {

public:

  //! standard constructor
  InvanaBase();

  //! destructor
  virtual ~InvanaBase() {};

  /*! Initialize class members.
   *
   *  Setup of the neccessary elements to be found in every
   *  standard inverse analysis procedure
   */
  void Init(
      Teuchos::RCP<DRT::Discretization> discret, ///< the underlying discretization
      Teuchos::RCP<INVANA::ObjectiveFunct> objfunct, ///< an objective function
      Teuchos::RCP<INVANA::MatParManager> matman, ///< the material parameterization
      Teuchos::RCP<INVANA::RegularizationBase> regman ///< a regularization strategy
      );

  //! Setup of the derived classes.
  virtual void Setup() = 0;

  /*! \brief Value and gradient of the functional
    \param In
    sol - the current solution of the optimization problem
    \param Out
    val - the value of the objective function
    \param Out
    optprob - the gradient of the objective function wrt to the optimization parameters
  */
  virtual void Evaluate(
      const Epetra_MultiVector& sol, ///< input to the problem
      double* val, ///> value of the functional to be filled
      Teuchos::RCP<Epetra_MultiVector> gradient ///> gradient of the functional
      ) = 0;

  //! a reference to the discretization
  Teuchos::RCP<DRT::Discretization> Discret() {return discret_;}

  //! a reference to the objective function
  Teuchos::RCP<INVANA::ObjectiveFunct> ObjectiveFunct() {return objfunct_;}

  //! a reference to the material parameter manager
  Teuchos::RCP<INVANA::MatParManager> Matman(){return matman_; }

  //! a reference to the regularization of the objective function
  Teuchos::RCP<INVANA::RegularizationBase> Regman(){return regman_; }

  //! initial guess to be used by the optimizer
  const Epetra_MultiVector& InitialGuess();

  //! a reference to the programs communicator
  const Epetra_Comm& Comm();

  //! the unique layout of the optimization parameters
  Teuchos::RCP<Epetra_Map> VectorRowLayout();

private:

  //! a reference to the discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! objective function to be optimized
  Teuchos::RCP<INVANA::ObjectiveFunct> objfunct_;

  //! handler for calls to the parametrization
  Teuchos::RCP<INVANA::MatParManager> matman_;

  //! handler for calls to the regularization
  Teuchos::RCP<INVANA::RegularizationBase> regman_;

  //! initialization bool
  bool isinit_;



};
}  // namespace INVANA

#endif /*INVANA_BASE_H_*/
