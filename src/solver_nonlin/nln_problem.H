/*----------------------------------------------------------------------*/
/*!
\file nln_problem.H

\brief Base class for interface of nonlinear solver to BACI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_H
#define NLN_PROBLEM_H

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;
class Epetra_Operator;

namespace NOX
{
  namespace Abstract
  {
    class Group;
  } // namespace Abstract
} // namespace NOX

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Base class for interface of nonlinear solver to BACI
   *
   *  Interface class between the nonlinear solver and BACI. We utilize the
   *  NOX::Abstract::Group() to access the time integration within BACI for
   *  Evaluate() / computeF().
   *
   *  Whenever a nonlinear operator wants to call Evaluate() or ConvergenceCheck(),
   *  the NLNSOL::NlnProblem provides access to the time integration scheme in
   *  BACI.
   *
   *  NLNSOL::NlnProblem acts on the DofRowMap() of the (fine) level, only. For
   *  coarse level representations, we use NLNSOL::NlnProblemCoarseLevel().
   *
   *  \sa NLNSOL::NlnProblemCoarseLevel()
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnProblem
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblem();

    //! Destructor
    virtual ~NlnProblem(){};

    //@}

    //! @name Setup
    //@{

    //! Initialize member variables
    virtual void Init(const Epetra_Comm& comm, ///< communicator
        const Teuchos::ParameterList& params, ///< parameter list
        NOX::Abstract::Group& noxgrp, ///< NOX group for call back routines to time integration
        Teuchos::RCP<Epetra_Operator> jac  /// Jacobian matrix of underlying problem
        );

    //! Setup of nonlinear problem
    virtual void Setup();

    //! Set coarse level residual corrections for FAS
    virtual void SetFHatFBar(Teuchos::RCP<Epetra_MultiVector> fhat, ///< \f$\hat{F}\f$
        Teuchos::RCP<Epetra_MultiVector> fbar ///< \f$\bar{F}\f$
        )
    {
      dserror("Not implemented in base class. This routine is needed only in "
          "case of a coarse level problem within FAS-AMG.");
    }

    //@}

    //! @name Actions on the nonlinear problem at hand
    //@{

    /*! \brief Evaluate the residual
     *
     *  Set the most recent solution 'x' in the NOX group #noxgrp_. Afterwards,
     *  the marching time integrator is asked to evaluate the residual via the
     *  computeF()-call.
     *
     *  Updating the solution variables in the marching time integrator is done
     *  during computeF().
     *
     *  Return the most recent residual in 'f'.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void Evaluate(const Epetra_MultiVector& x, ///< current solution
        Epetra_MultiVector& f ///< current residual (to be filled)
        ) const;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f ///< current residual
        ) const;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f, ///< current residual
        double& fnorm2 ///< residual L2 norm (to be filled)
        ) const;

    //@}

    //! @name Get attributes
    //@{

    //! communicator
    virtual const Epetra_Comm& Comm() const;

    //! parameter list
    virtual const Teuchos::ParameterList& Params() const;

    //! NOX group
    virtual NOX::Abstract::Group& NOXGroup() const;

    //@}

    //! @name Access to underlying NOX group elements
    //@{

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator();

    //@}

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblem"; }

    //@}

  protected:

    //! @name Attribute set methods
    //@{

    //! Init() has been called
    virtual void SetIsInit() { isinit_ = true; return; }

    //! Setup() has been called
    virtual void SetIsSetup() { issetup_ = true; return; }

    //@}

    //! @name Get attributes
    //@{

    //! Has Init() already been called?
    virtual bool IsInit() const { return isinit_; }

    //! Has Setup() already been called?
    virtual bool IsSetup() const { return issetup_; }

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

    //! NOX group for access to time integrator via call-back routines
    Teuchos::RCP<NOX::Abstract::Group> noxgrp_;

    /*! \brief Pointer to Jacobian matrix of underlying problem
     *
     *  \warning Access to time integration WIHTOUT using the NOX group #noxgrp_.
     *  Seems to be inconsistent since all other accesses to time integration
     *  are done via the NOX group #noxgrp_. However, usage of NOX group might
     *  go away at some point.
     */
    Teuchos::RCP<Epetra_Operator> jac_;

    //! @name Tolerances for convergence check
    //@{

    //! residual L2-norm tolerance
    double tolresl2_;

    //@}

    //! @name Configuration flags for printing
    //@{

    //! Print the convergence check? (yes / no)
    bool printconvcheck_;

    //@}

  }; // class NlnProblem
} // namespace NLNSOL

#endif /* NLN_PROBLEM_H */
