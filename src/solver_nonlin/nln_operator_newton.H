/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_newton.H

\brief Newton-type algorithm

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_NEWTON_H
#define NLN_OPERATOR_NEWTON_H

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
class Epetra_Operator;

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
  /*! \brief Newton-type algorithm as nonlinear preconditioner
   *
   *  \sa NLNSOL::NlnOperatorBase()
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorNewton : public NlnOperatorBase
  {
  public:

    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorNewton();

    //! Destructor
    virtual ~NlnOperatorNewton(){}

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
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNewton"; }

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

    //! Compute the search direction
    virtual const int ComputeSearchDirection(
        Teuchos::RCP<Epetra_MultiVector>& rhs, ///< residual
        Teuchos::RCP<Epetra_MultiVector>& inc, ///< increment (to be filled)
        const int iter ///< current Newton iteration
        ) const;

    //! Use a fixed Jacobian (true) or update it in every iteration (false)?
    virtual const bool FixedJacobian() const { return Params().get<bool>("Newton: Fixed Jacobian"); }

    //! @name Necessary sub-algorithms
    //@{

    //! Linear solver from baci
    Teuchos::RCP<LINALG::Solver> linsolver_;

    //! Line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //@}

    //! Fixed Jacobian matrix
    Teuchos::RCP<Epetra_Operator> jac_;

  }; // class NlnOperatorNewton
} // namespace NLNSOL

#endif /* NLN_OPERATOR_NEWTON_H */
