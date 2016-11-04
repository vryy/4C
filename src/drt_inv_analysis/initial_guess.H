/*----------------------------------------------------------------------*/
/*!
 * \file initial_guess.H
 * \brief Initial guess for optimization problems
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef INITIAL_GUESS_H_
#define INITIAL_GUESS_H_

#include "Teuchos_RCP.hpp"

namespace Teuchos{
class ParameterList;
}

class Epetra_Map;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsMatrix;

namespace IO{
class DiscretizationReader;
}

namespace TIMINT
{
  template <typename> class TimIntMStep;
}

namespace DRT
{
  class Discretization;
}

namespace INVANA
{

class MatParManager;
class ObjectiveFunct;
class DcsMatrix;

class InvanaBase;

/*! \class InitialGuess
 * \brief Initial guess for optimization problems
 *
 * \author kehl \date 08/2016
*/
class InitialGuess
{

 public:

  //! @name Construction and setup
  //@{

  //! Constructor
  InitialGuess(const Teuchos::ParameterList& invp);

  //! Destructor
  virtual ~InitialGuess() {};
  //@}

  //! @ name Compute initial guess for the optimization problem
  //@{

  //! Compute initial guess
  void Compute(Teuchos::RCP<DRT::Discretization> discret, /// [in] discret to build a reader
      Teuchos::RCP<MatParManager> matman, /// [in] parameterization
      Teuchos::RCP<ObjectiveFunct> objfunct /// [in] objective function
      );
  //@}


  //! @name Provide initial guess after computation
  //@{

  /*! \brief Initial solution
   *
   * Refering to the possibly variation Bayesian interpretation
   * of some solution methods it is called mean.
   *
   * \author kehl \date 08/2016
   */
  Teuchos::RCP<Epetra_MultiVector> Mean();

  /*! \brief Initial inverse Hessian approximation
   *
   * Refering to the possibly variation Bayesian interpretation
   * of some solution methods this is called covariance
   *
   * \author kehl \date 08/2016
   */
  Teuchos::RCP<Epetra_CrsMatrix> Covariance();

  //@}

private:

  //! Parameters with input
  const Teuchos::ParameterList& params_;

  //! has been computed
  bool iscomputed_;

  //! @name Initial guess data
  //@{

  //! mean
  Teuchos::RCP<Epetra_MultiVector> mean_;

  //! precision
  Teuchos::RCP<Epetra_CrsMatrix> covariance_;
  //@}

  //! @name internal compute initial guess
  //@{

  //! init from dat file info
  void InitfromDat(Teuchos::RCP<DRT::Discretization> discret, /// [in] discret to build a reader
      Teuchos::RCP<MatParManager> matman, /// [in] parameterization
      Teuchos::RCP<ObjectiveFunct> objfunct /// [in] objective function
      );

  //! init from map approximation
  void InitfromMap(Teuchos::RCP<DRT::Discretization> discret, /// [in] discret to build a reader
      Teuchos::RCP<MatParManager> matman, /// [in] parameterization
      Teuchos::RCP<ObjectiveFunct> objfunct /// [in] objective function
      );

  //! read storage from lbfgs
  void ReadLBFGSStorage(const Epetra_Map& readmap, /// [in] use this map for reading
      IO::DiscretizationReader& reader, /// [in] reader to use
      Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore, /// [out] sstore
      Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore, /// [out] ystore
      Teuchos::RCP<Epetra_MultiVector> solution ///[out] map solution
      );

  /* \brief apply projection
   *
   * project matrix toproject using the projector and
   * return result.
   */
  void ProjecttoBasis(
      const Epetra_CrsMatrix& projector, ///[in] projection matrix
      DcsMatrix& toproject ///[in] matrix to be projected
      );
  //@}



};

}
/*----------------------------------------------------------------------*/
#endif /* INITIAL_GUESS_H_ */
