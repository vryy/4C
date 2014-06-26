/*----------------------------------------------------------------------*/
/*!
\file nln_operator_quasinewton.H

\brief Quasi-Newton

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_QUASI_NEWTON_H
#define NLN_OPERATOR_QUASI_NEWTON_H

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace LINALG
{
  class Solver;
  class SparseOperator;
}

namespace NLNSOL
{
  class LineSearchBase;
}

namespace Teuchos
{
  class ParameterList;
}

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Quasi-Newton as nonlinear preconditioner
   *
   *  \sa NLNSOL::NlnOperator()
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorQuasiNewton : public NlnOperator
  {
  public:

    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorQuasiNewton();

    //! Destructor
    virtual ~NlnOperatorQuasiNewton(){}

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
    virtual const char* Label() const { return "NlnPrecQuasiNewton"; }

    //@}

  protected:

  private:

    //! @name Setup
    //@{

    //! Setup of linear solver
    virtual void SetupLinearSolver();

    //! Setup of line search
    virtual void SetupLineSearch();

    //@}

    //! @name Access routines
    //@{

    //! Get max allowed number of iterations
    virtual const int GetMaxIter() const { return maxiter_; }

    //@}

    //! linear solver from baci
    Teuchos::RCP<LINALG::Solver> linsolver_;

    //! line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //! max number of iterations
    int maxiter_;

  };
}

#endif /* NLN_OPERATOR_QUASI_NEWTON_H */
