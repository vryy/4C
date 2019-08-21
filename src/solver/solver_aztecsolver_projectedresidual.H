/*----------------------------------------------------------------------*/
/*! \file


\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

Created on: Mar 4, 2013
*/
/*---------------------------------------------------------------------------*/
#ifndef SOLVER_AZTECSOLVER_PROJECTEDRESIDUAL_H_
#define SOLVER_AZTECSOLVER_PROJECTEDRESIDUAL_H_

#include "AztecOO_StatusTestResNorm.h"
#include "Teuchos_RCP.hpp"

namespace LINALG
{
  class KrylovProjector;
}

//! AztecOO_StatusTestProjResNorm: An implementation of AztecOO_StatusTestResNorm
//  applying a projection on the residual prior to convergence check

/*!
 *  AztecOO_StatusTestProjResNorm is an implementation of
 *  AztecOO_StatusTestResNorm that allows a user to project the residual prior
 *  to the convergence check. This must be done in case that certain
 *  undetermined zero energy modes were projected with the help of a Krylov
 *  space projector to make the problem analytically and numerically
 *  non-singular.
 */

class AztecOO_StatusTestProjResNorm : public AztecOO_StatusTestResNorm
{
 public:
  //@{ \name Constructors/destructors.
  //! Constructor
  /*!
   * This class works exactly the same way AztecOO_StatusTestResNorm, except
   * that residual is projected prior to the convergence check. Thus the
   * Krylov space projector has to be available in this class as well.
   *
   * \param Operator  (In) The original linear operator that was passed in to
   *                       the AztecOO solver object.
   * \param LHS       (In) The original left hand side vector that was passed
   *                       in to the AztecOO solver object.  NOTE: AztecOO
   *                       accepts multivector objects, but
   *                       AztecOO_StatusTestResNorm does not handle residual
   *                       tests for multiple vectors.  Most AztecOO users tend
   *                       to have Epetra_Vectors, even though they are passing
   *                       them in as Epetra_MultiVectors, so this should not
   *                       be an issue.  If you are truly using
   *                       Epetra_MultiVector objects, remember that for a
   *                       multivector object mv, the ith vector can be
   *                       extracted as mv(i).
   * \param RHS       (In) The original right hand side vector that was passed
   *                       in to the AztecOO solver object. See note for LHS.
   * \param projector (In) Krylov space projector of the calling solver
   * \param Tolerance (In) A double value that is used to test the convergence
   *                       criterion.  Can be reset using ResetTolerance().
   */
  AztecOO_StatusTestProjResNorm(const Epetra_Operator& Operator, const Epetra_Vector& LHS,
      const Epetra_Vector& RHS, const Teuchos::RCP<LINALG::KrylovProjector>& projector,
      double Tolerance);

  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met. This is
    done after the residual has been projected apprpriately. Depending on how
    the residual test is constructed this method will return the appropriate
    status type.

    \param CurrentIter (In) Current iteration of iterative method.

    \param CurrentResVector (In) The current residuals of the iterative process.

    \param CurrentResNormEst (In) Estimate of the two-norm of the
    residual. Not available, leads to dserror.

    \param SolutionUpdated (In) If this argument is true, then the solution vector that is part
    of the Epetra_LinearProblem
    object being solved is consistent with the residual.

    \return AztecOO_StatusType: Unconverged, Converged or Failed.
  */
  AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector* CurrentResVector,
      double CurrentResNormEst, bool SolutionUpdated);

 private:
  //! Krylov space projector to project the residual
  Teuchos::RCP<LINALG::KrylovProjector> projector_;
};

#endif  // SOLVER_AZTECSOLVER_PROJECTEDRESIDUAL_h_
