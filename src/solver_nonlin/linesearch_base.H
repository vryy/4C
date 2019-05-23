/*----------------------------------------------------------------------------*/
/*!

\brief Line search base class

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_BASE_H
#define LINESEARCH_BASE_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "nln_utils.H"
#include "../drt_inpar/inpar_solver_nonlin.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_MultiVector;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace NLNSOL
{
  class NlnOperator;
  class NlnProblemBase;
}  // namespace NLNSOL

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class LineSearchBase
   *  \brief Line search base class
   *
   *  This base class for different line search strategies provides some common
   *  routines, that are part of every line search strategy. The core routine
   *  ComputeLSParam() has to be implemented by derived classes since there are
   *  the main differences between all line search strategies.
   *
   *  \author mayr.mt \date 05/2014
   */
  class LineSearchBase : public Teuchos::VerboseObject<LineSearchBase>
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
    virtual void Init(
        Teuchos::RCP<const NLNSOL::NlnProblemBase> nlnproblem,  ///< nonlinear problem to solve
        Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,    ///< configuration
        const std::string listname,      ///< name of sublist with configuration of this operator
        const Epetra_MultiVector& xold,  ///< last iterate
        const Epetra_MultiVector& fold,  ///< residual at last iterate
        const Epetra_MultiVector& inc,   ///< the Newton direction (full step)
        const double resnormold          ///< residual norm of last iterate
    );

    //! Setup the line search object
    virtual void Setup() = 0;

    //@}

    //! @name Computation
    //@{

    //! Compute the step length / line search parameter
    virtual void ComputeLSParam(double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr  ///< achieved sufficient decreas? (to be determined)
        ) const = 0;

    //! Compute the step length / line search parameter
    virtual void ComputeLSParam(double& lsparam  ///< line search parameter (to be determined)
        ) const;

    //@}

    //! @name Access to member variables
    //@{

    //! Has Init() already been called?
    virtual bool IsInit() const { return isinit_; }

    //! Has SetupLineSearch() already been called?
    virtual bool IsSetup() const { return issetup_; }

    //@}

    //! @name Description of this class
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const = 0;

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const = 0;

    //@}

   protected:
    /*! \brief Check sufficient decrease condition
     *
     *  <h3> References </h3>
     *  <ul>
     *  <li>[1] Kelley CT: Solving Nonlinear Equations with Newton's method,
     *          SIAM, p. 12 (2003) </li>
     *  </ul>
     *
     *  \author mayr.mt \date 11/2013
     */
    virtual bool IsSufficientDecrease(const double fnorm2,  ///< current L2-norm of residual
        const double lsparam                                ///< current line search parameter
        ) const;

    /*! \brief Safeguard strategy for line search parameter
     *
     *  Set upper and lower bound for ration of two subsequent line search
     *  parameters based on heuristic formulas.
     *
     *  <h3> References </h3>
     *  <ul>
     *  <li>[1] Kelley CT: Solving Nonlinear Equations with Newton's
     *          method, SIAM, p. 12 (2003) </li>
     *  </ul>
     *
     *  \author mayr.mt \date 11/2013
     */
    virtual void Safeguard(double& lsparam,  ///< current line search parameter
        const double lsparamold              ///< previous line search parameter
        ) const;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  Use the governing nonlinear solver object to evaluate the current
     *  residual.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void ComputeF(const Epetra_MultiVector& x,  ///< current solution
        Epetra_MultiVector& f                           ///< residual (to be filled)
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
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f,  ///< current residual
        double& fnorm2  ///< L2-norm of residual (to be filled)
        ) const;

    //! @name Attribute access functions
    //@{

    //! Access to last iterate
    virtual const Epetra_MultiVector& GetXOld() const { return *xold_; }

    //! Access to residual at last iterate
    virtual const Epetra_MultiVector& GetFOld() const { return *fold_; }

    //! Access to Newton direction
    virtual const Epetra_MultiVector& GetXInc() const { return *inc_; }

    //! Access to residual norm of last iterate
    virtual double GetFNormOld() const { return resnormold_; }

    //! Access to global configuration
    virtual Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> Configuration() const;

    //! Access to a parameter of my configuration
    template <typename ParamType>
    ParamType MyGetParameter(const std::string param) const
    {
      return Configuration()->GetParameter<ParamType>(MyListName(), param);
    }

    //! Get name of this operator's parameter list in the global configuration
    virtual const std::string& MyListName() const { return listname_; }

    //@}

    //! @name Set attributes
    //@{

    //! Init() has been called.
    virtual void SetIsInit()
    {
      isinit_ = true;
      return;
    }

    //! SetupLineSearch() has been called.
    virtual void SetIsSetup()
    {
      issetup_ = true;
      return;
    }

    //@}

   private:
    /*! \brief Test for sufficient decrease using Armijo's rule
     *
     *  <h3> References </h3>
     *  <ul>
     *  <li>[1] Kelley CT: Iterative Methods for Linear and Nonlinear
     *          Equations, SIAM, Philadelphia (1995) </li>
     *  <li>[2] Kelley CT: Solving Nonlinear Equations with Newton's method,
     *          SIAM, Philadelphia, p. 12 (2003)
     *  </ul>
     *
     *  \note This is also known as Goldstein conditions
     *
     *  \author mayr.mt \date 02/2015
     */
    virtual bool SufficientDecreaseArmijo(
        const double fnorm2,  ///< L2-norm of residual at trial solution
        const double lsparam  ///< current line search parameter
        ) const;

    /*! \brief Test for sufficient decrease using actual and predicted reductions
     *
     *  ToDo (mayr) include forcing term parameter
     *
     *  <h3> References </h3>
     *  <ul>
     *  <li>[1] Shadid JN, Tuminaro RS and Walter HF: An Inexact Newton Method
     *          for Fully Coupled Solution of the Navier-Stokes Equations with
     *          Heat and Mass Transport, Journal of Computational Physics,
     *          Vol. 137(1), pp. 155--185 (1997)</li>
     *  <li>[2] Eisenstat S and Walker H: Globally Convergent Inexact Newton
     *          Methods, SIAM J. Optim., Vol. 4(2), pp. 393--422 (1994)</li>
     *  <li>[3] Eisenstat S and Walker H: Choosing the Forcing Terms in an
     *          Inexact Newton Method, SIAM J. Sci. Comput., Vol. 17(1),
     *          pp. 16--32 (1996)</li>
     *  </ul>
     *
     *  \author mayr.mt \date 02/2015
     */
    virtual bool SufficientDecreaseARedPRed(
        const double fnorm2  ///< L2-norm of residual at trial solution
        ) const;

    /*! \brief Test for sufficient decrease very loosely
     *
     *  \author mayr.mt \date 02/2015
     */
    virtual bool SufficientDecreaseLoose(
        const double fnorm2  ///< L2-norm of residual at trial solution
        ) const;

    //! nonlinear operator that has initiated this line search
    Teuchos::RCP<const NLNSOL::NlnProblemBase> nlnproblem_;

    //! @name Configuration
    //@{

    //! Configuration of this line search operator
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

    //! Name of sublist in #config_ for this operator
    std::string listname_;

    //@}

    //! last iterate
    Teuchos::RCP<const Epetra_MultiVector> xold_;

    //! residual at the last iterate
    Teuchos::RCP<const Epetra_MultiVector> fold_;

    //! search direction
    Teuchos::RCP<const Epetra_MultiVector> inc_;

    //! residual norm of last iterate
    double resnormold_;

    //! type of sufficient decrease condition
    INPAR::NLNSOL::LINESEARCH::SufficientDeacreaseType suffdecrtype_;

    //! @name Setup / Initialization
    //@{

    //! Has Init() already been called?
    bool isinit_;

    //! Has SetupLineSearch() already been called?
    bool issetup_;

    //@}

  };  // class LineSearchBase
}  // namespace NLNSOL

#endif /* LINESEARCH_BASE_H */
