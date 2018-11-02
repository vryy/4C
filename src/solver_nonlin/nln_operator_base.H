/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_base.H

\brief General nonlinear operator base class

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_BASE_H
#define NLN_OPERATOR_BASE_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "nln_utils.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace LINALG
{
  class Solver;
}  // namespace LINALG

namespace NLNSOL
{
  class NlnProblemBase;
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* definition of classes */

namespace NLNSOL
{
  /*! \class NlnOperatorBase
   *  \brief General nonlinear operator base class
   *
   *  This is a base class for all kinds of nonlinear solution and
   *  preconditioning algorithms. The central method is ApplyInverse(), where
   *  the actual work is done.
   *
   *  Init() is the main initialization routine, that handles all base class
   *  initializations. The operator is fed with:
   *  <ul>
   *  <li> an Epetra_Comm </li>
   *  <li> a Teuchos::RCP<const Teuchos::NLNSOL::UTILS::NlnConfig> object that
   *       controls access to the global configuration of the NLNSOL framework
   *       </li>
   *  <li> a const std::string that describes the name of the sublist of this
   *       operator in the global configuration object. </li>
   *  <li> a Teuchos::RCP to the nonlinear problem interface class (connection
   *       to time integration) </li>
   *  <li> a parameter that reflects how deep the operator is nested in the
   *       overall algorithm </li>
   *  </ul>
   *
   *  Setup() constructs the actual algorithm and has to be implemented by every
   *  specific algorithm itself.
   *
   *  \note The flag #isinit_ has to be checked at the entry of Setup() to make
   *  sure, Init() has already been called.
   *
   *  \note The flags #isinit_ and #issetup_ indicate, whether Init() and
   *  Setup() have already been called. We need to check them when entering
   *  every public method to make sure, that the operator is in a valid state
   *  and can perform its actions properly. This is not necessary in protected
   *  and private routines.
   *
   *  ApplyInverse() is the central routine. It needs to implement the actual
   *  solution algorithm and iterative process.
   *
   *  After finishing iterations, output information is provided via an output
   *  parameter list #outparams_ which can be accessed via GetOutParams().
   *
   *  \author mayr.mt \date 05/2014
   */
  class NlnOperatorBase : public Teuchos::VerboseObject<NlnOperatorBase>
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

    //! Initialization: Fill base class member variables
    virtual void Init(const Epetra_Comm& comm,                ///< communicator
        Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,  ///< configuration
        const std::string listname,                           ///< name of sublist for this operator
        Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem,      ///< nonlinear problem to solve
        Teuchos::RCP<LINALG::Solver> bacisolver = Teuchos::null,  ///< linear solver from Baci
        const int nested = 0                                      ///< How nested is this algorithm?
    );

    //! Setup the specific nonlinear operator
    virtual void Setup() = 0;

    /*! \brief Rebuild preconditioner
     *
     *  Rebuild the nonlinear operator in case it is used as a nonlinear
     *  preconditioner. Some Operator do not require any action, some do.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrec() = 0;

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Returns the result of a NlnOperator inverse applied to an
     *  Epetra_MultiVector X in Y.
     *
     *  \return Integer error code, set to 0 if successful.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f,  ///< current residual
        Epetra_MultiVector& x                              ///< current solution (to be filled)
        ) const = 0;

    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const = 0;

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const = 0;

    //! Access to output parameter list
    virtual Teuchos::RCP<const Teuchos::ParameterList> GetOutParams() const;

    //! Has Init() already been called?
    virtual bool IsInit() const { return isinit_; }

    //! Has Setup() already been called?
    virtual bool IsSetup() const { return issetup_; }

    //@}

   protected:
    //! @name Set attributes
    //@{

    //! Init() has been called
    virtual void SetIsInit()
    {
      isinit_ = true;
      return;
    }

    /*! \brief Setup() has been called
     *
     *  We allow for the input argument 'issetup = false' to reset the setup
     *  state of the operator and to enforce a re-setup with a new state.
     */
    virtual void SetIsSetup(const bool issetup = true  ///< flag indicating validity of setup state
    )
    {
      issetup_ = issetup;
      return;
    }

    //@}

    //! @name Get attributes
    //@{

    //! Access to communicator
    virtual const Epetra_Comm& Comm() const;

    //! Access to global configuration
    virtual Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> Configuration() const;

    //! Access to a parameter of this operator's configuration
    template <typename ParamType>
    ParamType MyGetParameter(const std::string param  ///< name of parameter
        ) const
    {
      return Configuration()->GetParameter<ParamType>(MyListName(), param);
    }

    //! Get name of this operator's parameter list in the global configuration
    virtual const std::string& MyListName() const { return listname_; }

    //! Access to the nonlinear problem
    virtual Teuchos::RCP<NLNSOL::NlnProblemBase> NlnProblem() const;

    //! Access to Baci linear solver
    virtual Teuchos::RCP<LINALG::Solver> BaciLinearSolver() const;

    /*! \brief Is this nonlinear operator a solver or a smoother?
     *
     *  If this nonlinear operator is a solver, we perform a convergence check
     *  and throw an dserror() if convergence fails. If this nonlinear
     *  operator is a smoother, we do not care about convergence and just
     *  advance with the smoothed solution.
     */
    virtual bool IsSolver() const;

    //! Get maximally allowed number of iterations
    virtual int GetMaxIter() const;

    //! Get nested level
    int Nested() const { return nested_; }

    //! Get level of indentation for formatted screen output
    int Indentation() const { return 4 * nested_; }

    //@}

    /*! @name Iteration loop Utilities
     *
     *  To enable efficient implementations of nonlinear solvers, the
     *  iteration loops have to be implemented by each nonlinear operator
     *  itself. This is important especially in the case when history
     *  quantities are needed (as in NLNSOL::NlnOperatorNGmres or
     *  NLNSOL::NlnOperatorNonlinCG). To have a common base of these iteration
     *  loops, these helper routines provide actions that are needed in every
     *  iteration loop.
     *
     *  For future implementation of additional nonlinear operators, please
     *  make use of these generalized routines to make the code as unified
     *  as possible
     *
     *  \author mayr.mt \date 09/2014
     */
    //@{

    /*! \brief Check if we still continue to iterate
     *
     *  Use this in the while-statement of the iteration loop.
     *
     *  \return Bool that indicates whether iterations have to be continued or
     *  stopped
     *
     *  \author mayr.mt \date 09/2014
     */
    virtual bool ContinueIterations(const int iter,  ///< current iteration counter
        const bool converged                         ///< convergence flag
        ) const;

    //! Print summary of current iteration
    virtual void PrintIterSummary(const int iter,  ///< iteration counter
        const double fnorm2                        ///< L2-norm of residual
        ) const;

    /*! \brief Determine error code for ApplyInverse()
     *
     *  \note Here, we care for stagnation of iterative process.
     *
     *  \return error code of type NLNSOL::UTILS::OperatorStatus
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual NLNSOL::UTILS::OperatorStatus ErrorCode(const int iter,  ///< current iteration counter
        const bool converged,                                        ///< convergence status
        const bool error = false,                                    ///< any unspecified error
        const bool stagnation = false                                ///< stagnation status
        ) const;

    //! Write number of iterations to output parameter list
    virtual void SetOutParameterIter(const int iter) const;

    //! Write status of convergence to output parameter list
    virtual void SetOutParameterConverged(const bool converged) const;

    //! Write residual norm to output parameter list
    virtual void SetOutParameterResidualNorm(const double norm) const;

    //! Write operator status / error code to output parameter list
    virtual void SetOutParameterErrorCode(
        const NLNSOL::UTILS::OperatorStatus errorcode  ///< error code
        ) const;

    //! Write stagnation status to output parameter list
    virtual void SetOutParameterStagnation(Teuchos::RCP<const Teuchos::ParameterList>
            stagparams  ///< parameter list from stagnation detection
        ) const;

    //@}

   private:
    //! Compute step length, i.e. line search parameter
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const = 0;

    /*! \brief Check successful convergence within max number of iterations
     *
     *  Use this after exiting the iteration loop to check whether it
     *  converged successfully or not. The routine takes care of the decision
     *  whether the current nonlinear operator is used as a solver or just as
     *  a smoother.
     *
     *  \return Bool that indicates whether iteration loop converged
     *  successfully within the given number of maximally allowed iterations
     *
     *  \author mayr.mt \date 09/2014
     */
    virtual bool CheckSuccessfulConvergence(const int iter,  ///< number of performed iterations
        const bool converged                                 ///< convergence flag
        ) const;

    //! Set entry in output parameter list (templated on data type of entry)
    template <typename T>
    void SetOutParameter(const std::string param,  ///< name of parameter entry
        T value                                    ///< value of parameter entry
        ) const
    {
      if (outparams_.is_null())
        dserror(
            "Output parameter list has not been created, yet.\n"
            "Create and initialize it before writing to it!");
      outparams_->set<T>(param, value);
      return;
    }

    //! Init() has been called
    bool isinit_;

    //! Setup() has been called
    bool issetup_;

    //! Communicator
    Teuchos::RCP<const Epetra_Comm> comm_;

    //! @name Configuration
    //@{

    //! Configuration of this nonlinear operator
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

    //! Name of sublist in #config_ for this operator
    std::string listname_;

    //@}

    //! Nonlinear problem to be evaluated and solved
    Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem_;

    /*! \brief Linear solver from Baci
     *
     *  A linear solver from Baci can be passed in. These solvers provide
     *  additional features like AMG using ML. To use the linear solver
     *  created by Baci, choose '0' as linear solver number in xml input and
     *  set the desired linear solver number in the single field input block
     *  via the LINEAR_SOVLER input flag.
     *
     *  \author mayr.mt \date 03/2015
     */
    Teuchos::RCP<LINALG::Solver> bacisolver_;

    /*! \brief Parameter list for output purposes
     *
     *  \note Needs to be mutable to allow for writing access during
     *  Epeta_Operator::ApplyInverse()
     *
     *  \author mayr.mt \date 03/2015
     */
    mutable Teuchos::RCP<Teuchos::ParameterList> outparams_;

    //! How nested is this algorithm? 0 = outer algorithm
    int nested_;

  };  // class NlnOperatorBase
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_BASE_H */
