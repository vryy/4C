/*----------------------------------------------------------------------------*/
/*!
\file linesearch_base.H

\brief Line search base class

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_BASE_H
#define LINESEARCH_BASE_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_MultiVector;

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

namespace NLNSOL
{
  class NlnOperator;
  class NlnProblem;
} // namespace NLNSOL

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*!
   * \brief Line search base class
   *
   * \author mayr.mt
   * \date 05/2014
   */
  class LineSearchBase
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchBase();

    //! Destructor
    virtual ~LineSearchBase(){};

    //@}

    //! @name Setup
    //@{

    //! Initialize some member variables
    virtual void Init(Teuchos::RCP<const NLNSOL::NlnProblem> nlnproblem, ///< nonlinear problem to solve
        const Teuchos::ParameterList& params, ///< parameter list for configuration of line search algorithm
        const Epetra_MultiVector& xold, ///< last iterate
        const Epetra_MultiVector& inc, ///< the Newton direction (full step)
        const double resnormold ///< residual norm of last iterate
        );

    //! Setup the line search object
    virtual void Setup() = 0;

    //@}

    //! @name Computation
    //@{

    //! Compute the step length / line search parameter
    virtual const double ComputeLSParam() const = 0;

    //@}

    //! @name Access to member variables
    //@{

    //! Has Init() already been called?
    virtual const bool IsInit() const { return isinit_; }

    //! Has SetupLineSearch() already been called?
    virtual const bool IsSetup() const { return issetup_; }

    //@}

  protected:

    /*! \brief Check sufficient decrease condition
     *
     *  <h3> References </h3>
     *
     *  - Kelley CT: Solving Nonlinear Euqations with Newton's
     *  method, SIAM, p. 12 (2003)
     *
     *  \author mayr.mt \date 11/2013
     */
    virtual bool IsSufficientDecrease(const double fnorm2, ///< current L2-norm of residual
        const double lsparam ///< current line search parameter
        ) const;

    /*! \brief Safeguard strategy for line search parameter
     *
     *  Set upper and lower bound for ration of two subsequent line search
     *  parameters based on heuristic formulas.
     *
     *  <h3> References </h3>
     *
     *  - Kelley CT: Solving Nonlinear Euqations with Newton's
     *  method, SIAM, p. 12 (2003)
     *
     *  \author mayr.mt \date 11/2013
     */
    virtual void Safeguard(double& lsparam, ///< current line search parameter
        const double lsparamold ///< previous line search parameter
        ) const;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  Use the governing nonlinear solver object to evaluate the current
     *  residual.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void ComputeF(const Epetra_MultiVector& x, ///< current solution
        Epetra_MultiVector& f  ///< residual (to be filled)
        ) const;

    /*! \brief Evaluate residual
     *
     *  Use the governing nonlinear solver object to evaluate the current
     *  residual.
     *
     *  \return Bool indicating whether convergence has been achieved or not
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f, ///< current residual
        double& fnorm2  ///< L2-norm of residual (to be filled)
        ) const;

    //! @name Attribute access functions
    //@{

    //! Access to last iterate
    virtual const Epetra_MultiVector& GetXOld() const { return *xold_; }

    //! Access to Newton direction
    virtual const Epetra_MultiVector& GetXInc() const { return *inc_; }

    //! Access to residual norm of last iterate
    virtual const double GetFNormOld() const { return resnormold_; }

    //! Access to parameter list
    virtual const Teuchos::ParameterList& Params() const;

    //@}

    //! @name Set attributes
    //@{

    //! Init() has been called.
    virtual void SetIsInit() { isinit_ = true; return; }

    //! SetupLineSearch() has been called.
    virtual void SetIsSetup() { issetup_ = true; return; }

    //@}

  private:

    //! nonlinear operator that has initiated this line search
    Teuchos::RCP<const NLNSOL::NlnProblem> nlnproblem_;

    //! parameter list
    Teuchos::RCP<const Teuchos::ParameterList> params_;

    //!< the last iterate
    Teuchos::RCP<const Epetra_MultiVector> xold_;

    //!< the Newton direction (full step)
    Teuchos::RCP<const Epetra_MultiVector> inc_;

    //!< residual norm of last iterate
    double resnormold_;

    //! @name Setup / Initialization
    //@{

    //! Has Init() already been called?
    bool isinit_;

    //! Has SetupLineSearch() already been called?
    bool issetup_;

    //@}

  }; // class LineSearchBase
} // namespace NLNSOL

#endif /* LINESEARCH_BASE_H */
