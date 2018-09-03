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
  /*! \class NlnOperatorNewton
   *  \brief Newton-type algorithm as nonlinear preconditioner
   *
   *  This Newton-type algorithm enables several versions that differ in the
   *  frequency of updating the Jacobian. This update frequency is controlled
   *  via #jacevery_. The following methods can be represented
   *  (c.f. [Kelley (1995), chaps. 5.3 and 5.4]):
   *  <ul>
   *  <li> Full Newton with exact linearization </li>
   *  <li> Chord's method (keeps the Jacobian fixed) </li>
   *  </ul>
   *
   *  <h3> References </h3>
   *  <ul>
   *  <li>[1] Kelley CT: Iterative Methods for Linear and Nonlinear Equations,
   *          SIAM, Philadelphia (1995)
   *  </ul>
   *
   *  \sa NlnOperatorBase
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
    virtual ~NlnOperatorNewton() {}

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
    virtual void RebuildPrec(){/* no action required */};

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner
     *
     *  Implementation of Newton-type algorithm follows
     *  [Kelley (1995), p. 86, Algorithm 5.6.1].
     *
     *  The algorithmic parameter #jacevery_ controls, whether the Jacobian
     *  needs to be updated or not. We allow for
     *  <ul>
     *  <li> Chord's method: #jacevery_ = 0 </li>
     *  <li> Full Newton: #jacevery_ = 1 </li>
     *  </ul>
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
    virtual const char* Label() const { return "NLNSOL::NlnOperatorNewton"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "Newton"; }

    //@}

   protected:
   private:
    //! @name Setup of internal objects/sub-algorithms
    //@{

    //! Setup of linear solver
    virtual void SetupLinearSolver();

    //! Setup of line search
    virtual void SetupLineSearch();

    //@}

    //! @name Search direction and step length
    //@{

    /*! \brief Compute the search direction
     *
     *  Compute search direction \f$p\f$ by solving \f$Jp=f\f$ with the Jacobian
     *  matrix \f$J\f$ and the residual \f$f\f$.
     *
     *  \note We test for success of the linear solver only, if Newton is
     *  applied as the outer nonlinear solver. If Newton is used as
     *  preconditioner, we allow incomplete linear solves (e.g. with limited
     *  Krylov subspace) to save computational time.
     *
     *  \sa ComputeStepLength
     *
     *  \author mayr \date 02/2015
     */
    virtual int ComputeSearchDirection(Teuchos::RCP<Epetra_MultiVector> rhs,  ///< residual
        Teuchos::RCP<Epetra_MultiVector> inc,  ///< increment (to be filled)
        const bool refactor                    ///< re-factor Jacobian?
        ) const;

    //! Compute line search parameter / step length
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! @name Jacobian update strategy
    //@{

    /*! \brief Update Jacobian matrix
     *
     *  Updating the Jacobian matrix if desired by the user. The update strategy
     *  is controlled via #jacevery_.
     *
     *  \author mayr.mt \date 02/2015
     */
    virtual void UpdateJacobian(bool& refactor  ///< Do we need to re-factor the Jacobian?
        ) const;

    /*! \brief Recompute Jacobian every #jacevery_ iterations
     *
     *  This number specifies the exact Newton-type algorithm. There are several
     *  options (c.f. [Kelley (1995), chaps. 5.3 and 5.4]):
     *  <ul>
     *  <li> Full Newton with exact linearization: set to 1 </li>
     *  <li> Chord's method keeps the Jacobian fixed: set to 0 </li>
     *  </ul>
     */
    int jacevery_;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! Linear solver from baci
    Teuchos::RCP<LINALG::Solver> linsolver_;

    //! Line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //@}

  };  // class NlnOperatorNewton
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_NEWTON_H */
