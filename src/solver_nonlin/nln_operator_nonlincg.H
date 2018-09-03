/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_nonlincg.H

\brief Nonlinear Conjugate Gradient

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_NONLIN_CG_H
#define NLN_OPERATOR_NONLIN_CG_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"

#include "../drt_inpar/inpar_solver_nonlin.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace LINALG
{
  class Solver;
  class SparseOperator;
}  // namespace LINALG

namespace NLNSOL
{
  class LineSearchBase;
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorNonlinCG
   *  \brief Nonlinear Conjugate Gradient as nonlinear preconditioner
   *
   *  <h3>References:</h3>
   *  <ul>
   *  <li>[1] Shewchuk, J. R.: An Introduction to the Conjugate Gradient
   *          Method Without the Agonizing Pain (1994) </li>
   *  <li>[2] Hager, W. W. and Zhang, H.: A survey of nonlinear conjugate
   *          gradient methods, Pacific Journal of Optimization (2), issue 1,
   *          pp. 35-58, 2006</li>
   *  <li>[3] Gee, M. W and Tuminaro, R. S.: Nonlinear algebraic multigrid for
   *          constrained solid mechanics problems, Sandia National
   *          Laboratories, Technical Report Nr. SAND2006-2256, 2006 </li>
   *  </ul>
   *
   *  \sa NlnOperatorBase
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorNonlinCG : public NlnOperatorBase
  {
   public:
    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorNonlinCG();

    //! Destructor
    virtual ~NlnOperatorNonlinCG() {}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
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

    /*! \brief Apply the preconditioner
     *
     *  Implementation follows reference [1].
     *
     *  \return Integer error code, set to 0 if successful.
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(
        const Epetra_MultiVector&
            f,  ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
        Epetra_MultiVector& x  ///< solution
        ) const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNonlinCG"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "NonlinCG"; }

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

    //! Setup of line search
    virtual void SetupLineSearch();

    //! Setup nonlinear preconditioner
    virtual void SetupPreconditioner();

    //@}

    /*! \brief Apply the preconditioner once
     *
     *  \return integer error code, set to 0 if successful
     */
    virtual int ApplyPreconditioner(const Epetra_MultiVector& f,  ///< residual
        Epetra_MultiVector& x                                     ///< solution
        ) const;

    //! @name Computation of beta
    //@{

    /*! \brief Compute the parameter beta
     *
     *  There are several methods of computing beta:
     *  - Polak-Ribiere
     *  - Fletcher-Reeves
     *  - Hestenes-Stiefel
     *
     *  Decide internally, which one is used.
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBeta(double& beta,            ///< parameter to be computed
        Teuchos::RCP<const Epetra_MultiVector> fnew,  ///< most recent residual
        Teuchos::RCP<const Epetra_MultiVector> fold,  ///< residual from previous iteration
        Teuchos::RCP<const Epetra_MultiVector> s,     ///< most recent preconditioned residual
        Teuchos::RCP<const Epetra_MultiVector>
            sold  ///< preconditioned residual from previous iteration
        ) const;

    /*! \brief Compute the parameter beta according to Fletcher-Reeves formula
     *
     *  \f[ \beta_{k+1}^{FR} = \frac{f_{k+1}^TM^{-1}f_{k+1}}{f_{k}^TM^{-1}f_{k}} \f]
     *
     *  \note This is the same as in linear CG.
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaFletcherReeves(double& beta,  ///< parameter to be computed
        Teuchos::RCP<const Epetra_MultiVector> fnew,      ///< most recent residual
        Teuchos::RCP<const Epetra_MultiVector> fold,      ///< residual from previous iteration
        Teuchos::RCP<const Epetra_MultiVector> s,         ///< most recent preconditioned residual
        Teuchos::RCP<const Epetra_MultiVector>
            sold  ///< preconditioned residual from previous iteration
        ) const;

    /*! \brief Compute the parameter beta according to Polak-Ribiere formula
     *
     *  \f[ \beta_{k+1}^{PR} = \frac{f_{k+1}^TM^{-1}\left(f_{k+1}-f_{k}\right)}{f_{k}^TM^{-1}f_{k}}
     * \f]
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaPolakRibiere(double& beta,  ///< parameter to be computed
        Teuchos::RCP<const Epetra_MultiVector> fnew,    ///< most recent residual
        Teuchos::RCP<const Epetra_MultiVector> fold,    ///< residual from previous iteration
        Teuchos::RCP<const Epetra_MultiVector> s,       ///< most recent preconditioned residual
        Teuchos::RCP<const Epetra_MultiVector>
            sold  ///< preconditioned residual from previous iteration
        ) const;

    /*! \brief Compute the parameter beta according to Hestenes-Stiefel formula
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaHestenesStiefel(double& beta) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //! Nonlinear preconditioner
    mutable Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnprec_;

    //@}

    //! Formula that is used to compute the parameter \f$\beta\f$
    INPAR::NLNSOL::NONLINCG::BetaType betatype_;

    //! Restart the iteration every #restartevery_ iterations
    int restartevery_;

  };  // class NlnOperatorNonlinCG
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_NONLIN_CG_H */
