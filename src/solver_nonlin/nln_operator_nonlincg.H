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

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace LINALG
{
  class Solver;
  class SparseOperator;
} // namespace LINALG

namespace NLNSOL
{
  class LineSearchBase;
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Nonlinear Conjugate Gradient as nonlinear preconditioner
   *
   *  \sa NLNSOL::NlnOperatorBase()
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
    virtual ~NlnOperatorNonlinCG(){}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner, i.e. do one step of Quasi-Newton
     *
     * \return Integer error code, set to 0 if successful.
     *
     * \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
        Epetra_MultiVector& x ///< solution
        )  const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNonlinCG"; }

    //@}

  protected:

  private:

    //! @name Setup
    //@{

    //! Setup of linear solver
    virtual void SetupLinearSolver();

    //! Setup of line search
    virtual void SetupLineSearch();

    //! Setup nonlinear preconditioner
    virtual void SetupPreconditioner();

    //@}

    //! Apply the preconditioner once
    virtual void ApplyPreconditioner
    (
        const Epetra_MultiVector& f, ///< residual
        Epetra_MultiVector& x ///< solution
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
    virtual void ComputeBeta(double& beta, ///< parameter to be computed
        const Epetra_MultiVector& fnew,  ///< most recent residual
        const Epetra_MultiVector& fold   ///< residual from previous iteration
        ) const;

    /*! \brief Compute the parameter beta according to Fletcher-Reeves formula
     *
     *  \f$ \beta_{k+1}^{FR} = \frac{f_{k+1}^Tf_{k+1}}{f_{k}^Tf_{k}} \f$
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaFletcherReeves(double& beta, ///< parameter to be computed
        const Epetra_MultiVector& fnew, ///< most recent residual
        const Epetra_MultiVector& fold ///< residual from previous iteration
        ) const;

    /*! \brief Compute the parameter beta according to Polak-Ribiere formula
     *
     *  \f$ \beta_{k+1}^{PR} = \frac{f_{k+1}^T\left(f_{k+1}-f_{k}\right)}{\vert f_{k}\vert^2} \f$
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaPolakRibiere(double& beta) const;

    /*! \brief Compute the parameter beta according to Hestenes-Stiefel formula
     *
     *  \author mayr \date 11/2013
     */
    virtual void ComputeBetaHestenesStiefel(double& beta) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! linear solver from baci
    Teuchos::RCP<LINALG::Solver> linsolver_;

    //! line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //! Nonlinear preconditioner
    Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnprec_;

    //@}

  }; // class NlnOperatorNonlinCG
} // namespace NLNSOL

#endif /* NLN_OPERATOR_NONLIN_CG_H */
