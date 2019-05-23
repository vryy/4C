/*----------------------------------------------------------------------------*/
/*!

\brief Nonlinear GMRES

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_NGMRES_H
#define NLN_OPERATOR_NGMRES_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;

namespace NLNSOL
{
  class LineSearchBase;
  class NlnProblem;
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorNGmres
   *  \brief Nonlinear GMRES
   *
   *  We do a nonlinear Krylov acceleration, also referred to as nonlinear
   *  GMRES (NGMRES). By recombining previous iterates/solutions, we hope to
   *  accelerate the convergence of the nonlinear solver. The algorithm consists
   *  of 3 steps:
   *
   *  Step 1: Generate a new tentative iterate/solution
   *
   *  Step 2: Build a linear combination of the previous iterates and the most
   *          recent one. The weights of each iterate are determined by solving
   *          a least squares problem with the goal of reducing the L2-norm of
   *          the residual.
   *
   *  Step 3: Line search between tentative solution (cf. step 1) and the one
   *          obtained by linear combination (cf. step 2).
   *
   *  Basically, a linear system is not required. Evaluation of the residual is
   *  sufficient.
   *
   *  <h3> References </h3>
   *  <ul>
   *  <li>[1] Washio T and Oosterlee CW: Krylov subspace acceleration for
   *          nonlinear multigrid schemes, Electronic Transactions on Numerical
   *          Analysis, Vol. 6, pp. 271 - 290, 1997 </li>
   *  <li>[2] Sterck, Hans De: A Nonlinear GMRES Optimization Algorithm for
   *          Canonical Tensor Decomposition, SIAM J. Sci. Comput., Vol. 34,
   *          pp. A1351-1379, 2012 </li>
   *  </ul>
   *
   *  \note Both references, [Washio (1997)] and [Sterck (2012)], describe the
   *        same algorithm. However, some signs are different at certain places.
   *        The present implementation follows the notation in [Sterck (2012)].
   *
   *  \sa NlnOperatorBase
   *
   *  \author mayr.mt \date 03/2014
   */
  class NlnOperatorNGmres : public NlnOperatorBase
  {
   public:
    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorNGmres();

    //! Destructor
    virtual ~NlnOperatorNGmres() {}

    //@}

    //! Setup of nonlinear GMRES algorithm
    virtual void Setup();

    /*! \brief Rebuild preconditioner
     *
     *  No action required.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrec();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Do the Krylov acceleration procedure.
     *
     *  The algorithm consists of three steps:
     *    Step 1: Generate a new tentative iterate/solution
     *
     *    Step 2: Build a linear combination of the previous iterates and the
     *            most recent one. The weights of each iterate are determined by
     *            solving a least squares problem with the goal of reducing the
     *            L2-norm of the residual.
     *
     *    Step 3: Line search between tentative solution (cf. step 1) and the
     *            one obtained by linear combination (cf. step 2).
     *
     *  The implementation follows [Sterck (2012)].
     *
     *  \return Integer error code, set to 0 if successful.
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f,  ///< residual
        Epetra_MultiVector& x                              ///< solution
        ) const;

    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNGmres"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "NGmres"; }

    //@}

   protected:
   private:
    //! Compute line search parameter (=step length)
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const;

    //! @name Setup
    //@{

    //! Setup line search
    virtual void SetupLineSearch();

    //! Setup nonlinear preconditioner
    virtual void SetupPreconditioner();

    //@}

    //! @name Contributions to the algorithm
    //@{

    /*! \brief Compute a new tentative iterate
     *
     *  Apply the preconditioner once and store result in \p xbar.
     *
     *  \return Integer error code, set to 0 if successful
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual int ComputeTentativeIterate(const Epetra_MultiVector& fbar,  ///< current residual
        Epetra_MultiVector& xbar                                         ///< current solution
        ) const;

    /*! \brief Compute the accelerated iterate
     *
     *  Generate an accelerated iterate by nonlinear GMRES step. This requires
     *  solving a least squares problem for the coefficients alpha_j that are
     *  needed to compute the accelerated iterate as linear combination of
     *  previous iterates.
     *
     *  The least squares problem is solved using normal equations where the
     *  corresponding system is build purely on information from previous and
     *  current residuals. The system size equals the window size and, thus, is
     *  rather small, such that it can be solved with a direct solver
     *  efficiently.
     *
     *  The implementation follows [Sterck (2012)].
     *
     *  \return Boolean error code (set to true if everything is ok, false
     *          otherwise)
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual bool ComputeAcceleratedIterate(const Teuchos::RCP<const Epetra_MultiVector>
                                               xbar,  ///< tentative iterate / current solution
        const Teuchos::RCP<const Epetra_MultiVector>
            fbar,  ///< residual based on tentative iterate / current solution
        const std::vector<Teuchos::RCP<const Epetra_MultiVector>>&
            sol,  ///< collection of previous solution iterates
        const std::vector<Teuchos::RCP<const Epetra_MultiVector>>&
            res,                               ///< collection of previous residuals
        Teuchos::RCP<Epetra_MultiVector> xhat  ///< accelerated iterate (to be determined)
        ) const;

    /*! \brief Perform line search
     *
     *  We perform a line search along the search direction betweeen \p xbar and
     *  \p xhat. For storage and algorithmic reasons, the result is written in
     *  \p xbar on output.
     *
     *  The implementation follows [Sterck (2012)].
     *
     *  \return Boolean error code (set to true if everything is ok, false
     *          otherwise)
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void PerformLineSearchStep(
        Teuchos::RCP<Epetra_MultiVector> xbar,  ///< previous iterate (input), next iterate (output)
        Teuchos::RCP<const Epetra_MultiVector> xhat,  ///< accelerated solution
        Teuchos::RCP<Epetra_MultiVector> fbar         ///< residual at previous iterate
        ) const;

    //@}

    //! @name Utilities
    //@{

    /*! \brief Add iterate and residual vectors to window history
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void AddToWindow(
        Teuchos::RCP<const Epetra_MultiVector> vec,  ///< vector to add to history
        std::vector<Teuchos::RCP<const Epetra_MultiVector>>&
            history  ///< collection of history vectors
        ) const;

    //@}

    //! @name Access routines
    //@{

    //! Get the maximally allowed window size
    virtual unsigned int GetMaxWindowSize() const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! Line search algorithm
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //! Nonlinear preconditioner
    mutable Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnprec_;

    //@}

  };  // class NlnOperatorNGmres
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_NGMRES_H */
