/*----------------------------------------------------------------------*/
/*!
\file nln_operator_ngmres.H

\brief Nonlinear GMRES

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_NGMRES_H
#define NLN_OPERATOR_NGMRES_H

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;

namespace LINALG
{
  class Solver;
} // namespace LINALG

namespace NLNSOL
{
  class LineSearchBase;
  class NlnProblem;
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Nonlinear GMRES
   *
   *  We do a nonlinear Krylov acceleration, also referred to as nonlinear GMRES (NGMRES).
   *  By recombining previous iterates/solutions, we hope to accelerate the
   *  convergence of the nonlinear solver. The algorithm consists of 3 steps:
   *
   *  Step 1: Generate a new tentative iterate/solution
   *
   *  Step 2: Build a linear combination of the previous iterates and the most
   *          recent one. The weights of each iterate are determined by solving a
   *          least squares problem with the goal of reducing the L2-norm of the
   *          residual.
   *
   *  Step 3: Line search between tentative solution (cf. step 1) and the one
   *          obtained by linear combination (cf. step 2).
   *
   *  Basically, a linear system is not required. Evaluation of the residual is
   *  sufficient.
   *
   *
   *  <h3> References </h3>
   *
   *  - Washio, T. and Oosterlee, C. W.: Krylov subspace acceleration for
   *  nonlinear multigrid schemes, Electronic Transactions on Numerical Analysis,
   *  Vol. 6, pp. 271 - 290, 1997
   *
   *  - Sterck, Hans De: A Nonlinear GMRES Optimization Algorithm for
   *  Canonical Tensor Decomposition, SIAM J. Sci. Comput., Vol. 34, pp. A1351-
   *  1379, 2012
   *
   * \sa NLNSOL::NlnOperatorBase()
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
    virtual ~NlnOperatorNGmres(){}

    //@}

    //! Setup the nonlinear GMRES algorithm
    virtual void Setup();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Do the Krylov acceleration procedure.
     *
     *  The algorithm consists of three steps:
     *    Step 1: Generate a new tentative iterate/solution
     *
     *    Step 2: Build a linear combination of the previous iterates and the most
     *            recent one. The weights of each iterate are determined by solving a
     *           least squares problem with the goal of reducing the L2-norm of the
     *           residual.
     *
     *    Step 3: Line search between tentative solution (cf. step 1) and the one
     *           obtained by linear combination (cf. step 2).
     *
     *  Implementation follows this publication: Sterck, Hans De: A Nonlinear
     *  GMRES Optimization Algorithm for Canonical Tensor Decomposition,
     *  SIAM J. Sci. Comput., Vol. 34, pp. A1351-1379, 2012
     *
     * \return Integer error code, set to 0 if successful.
     *
     * \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual
        Epetra_MultiVector& x ///< solution
        )  const;

    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNGmres"; }

    //@}

  protected:

    //! @name Attribute set methods
    //@{

    //@}

  private:

    //! @name Setup
    //@{

    //! Setup nonlinear preconditioner
    virtual void SetupPreconditioner();

    //! Setup linear solver
    virtual void SetupLinearSolver();

    //! Setup line search
    virtual void SetupLineSearch();

    //@}

    //! @name Contributions to the algorithm
    //@{

    /*! \brief Compute a new tentative iterate
     *
     *  Apply the preconditioner once and store result in #xbar. Afterwards,
     *  evaluate the residual and store it in #fbar.
     *
     *  \return Integer error code, set to 0 if successful
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual int ComputeTentativeIterate
    (
      const Epetra_MultiVector& fbar,  ///< current residual
      Epetra_MultiVector& xbar   ///< current solution
    ) const;

    /*! \brief Compute the accelerated iterate
     *
     *  Generate tan accelerated iterate by nonlinear GMRES step. This requires
     *  solving a least squares problem for the coefficients alpha_j that are
     *  needed to compute the accelerated iterate as linear combination of
     *  previous iterates.
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual Teuchos::RCP<Epetra_MultiVector> ComputeAcceleratedIterate
    (
      const Teuchos::RCP<Epetra_MultiVector>& xbar,               ///< tentative iterate / current solution
      const Teuchos::RCP<Epetra_MultiVector>& fbar,               ///< residual based on tentative iterate / current solution
      const std::vector<Teuchos::RCP<Epetra_MultiVector> >& sol,  ///< collection of previous solution iterates
      const std::vector<Teuchos::RCP<Epetra_MultiVector> >& res,  ///< collection of previous residuals
      const int wsize                                             ///< current window size
    ) const;

    //@}

    //! @name Utilities
    //@{

    //! Print summary of current iteration
    virtual void PrintIterSummary(const int iter, ///< iteration counter
        const double resnorm  ///< residual norm
        ) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! linear solver from baci
    Teuchos::RCP<LINALG::Solver> linsolver_;

    //! Line search algorithm
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //! Nonlinear preconditioner
    Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnprec_;

    //@}

    //! max number of iterations
    int itermax_;

    //! max window size;
    int windowsize_;

  }; // class NlnOperatorNGmres
} // namespace NLNSOL

#endif /* NLN_OPERATOR_NGMRES_H */
