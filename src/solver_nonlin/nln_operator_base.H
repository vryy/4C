/*----------------------------------------------------------------------*/
/*!
\file nln_operator_base.H

\brief General nonlinear operator base class

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_BASE_H
#define NLN_OPERATOR_BASE_H

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace NLNSOL
{
  class NlnProblem;
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definition of classes */
//! NLNSOL: a nonlinear solver tool box
namespace NLNSOL
{
  /*! \brief General nonlinear operator base class
   *
   *  This is a base class for all kinds of nonlinear solution and preconditioning
   *  algorithms. The central method is ApplyInverse(), where the actual work is
   *  done.
   *
   *  The operator is fed with:
   *  - an Epetra_Comm
   *  - a Teuchos::ParameterList for configuration
   *  - a Teuchos::RCP to the nonlinear problem interface class (connection to time integration)
   *
   *  Setup() constructs the actual algorithm and has to be implemented by every
   *  specific algorithm itself.
   *
   *  \author mayr.mt \date 05/2014
   */
  class NlnOperatorBase
  {
    public:

      //! @name Construction / Destruction
      //@{

      //! Constructor (empty)
      NlnOperatorBase();

      //! Destructor
      virtual ~NlnOperatorBase(){};

      //@}

      //! @name Setup
      //@{

      // fill base class member variables
      virtual void Init(const Epetra_Comm& comm, ///< communicator
          const Teuchos::ParameterList& params, ///< parameter list
          Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem ///< nonlinear problem to solve
          );

      // setup the specific nonlinear operator
      virtual void Setup() = 0;

      //@}

      //! @name Mathematical functions
      //@{

      /*! \brief Returns the result of a NlnOperatir inverse applied to an Epetra_MultiVector X in Y.
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual int ApplyInverse(const Epetra_MultiVector& f, ///< current residual
          Epetra_MultiVector& x ///< current solution (to be filled)
          ) const = 0;

      //@}

      //! @name Attribute access functions
      //@{

      //! Returns a character string describing the operator
      virtual const char* Label() const = 0;

      //@}

    protected:

      //! @name Set attributes
      //@{

      //! Init() has been called
      virtual void SetIsInit() { isinit_ = true; return; }

      //! Setup() has been called
      virtual void SetIsSetup() { issetup_ = true; return; }

      //@}

      //! @name Get attributes
      //@{

      //! Has Init() already been called?
      virtual const bool IsInit() const { return isinit_; }

      //! Has Setup() already been called?
      virtual const bool IsSetup() const { return issetup_; }

      //! Access to communicator
      virtual const Epetra_Comm& Comm() const;

      //! Access to parameter list
      virtual const Teuchos::ParameterList& Params() const;

      //! Access to the nonlinear problem
      virtual Teuchos::RCP<NLNSOL::NlnProblem> NlnProblem() const;

      //! Is this nonlinear operator a solver?
      virtual bool IsSolver() const { return issolver_; }

      //@}

    private:

      //! Init() has been called
      bool isinit_;

      //! Setup() has been called
      bool issetup_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! parameter list
      Teuchos::RCP<const Teuchos::ParameterList> params_;

      //! Nonlinear problem to be evaluated and solved
      Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem_;

      /*! \brief Is this nonlinear operator a solver or a smoother?
       *
       *  If this nonlinear operator is a solver, we perform a convergence check
       *  and throw an dserror() if convergence fails. If this nonlinear
       *  operator is a smoother, we do not care about convergence and just
       *  advance with the smoothed solution.
       */
      bool issolver_;

  }; // class NlnOperatorBase
} // namespace NLNSOL

#endif /* NLN_OPERATOR_BASE_H */
