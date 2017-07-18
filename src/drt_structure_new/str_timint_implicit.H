/*-----------------------------------------------------------*/
/*!
\file str_timint_implicit.H

\brief Implicit structural time integration strategy.

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_TIMINT_IMPLICIT_H_
#define STR_TIMINT_IMPLICIT_H_

#include "str_timint_implicitbase.H"    // base class

// forward declarations ...
namespace NOX {
namespace NLN {
namespace INNER {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
} // namespace INNER
} // namespace NLN
} // namespace NOX
namespace STR {
namespace IMPLICIT {
  class Generic;
} // namespace IMPLICIT
namespace PREDICT{
  class Generic;
} // namespace PREDICT
namespace NLN {
namespace SOLVER {
  class Generic;
namespace INTERFACE {
  class Required;
} // namespace INTERFACE
} // namespace SOLVER
} // namespace NLN

namespace TIMINT {
  /** \brief Implicit time integration strategy
   *
   * \author Michael Hiermeier */
  class Implicit : public ImplicitBase
  {
  public:
    //! constructor
    Implicit();

    //! destructor
    virtual ~Implicit() {};

    virtual void Setup();

    virtual int Integrate();

    virtual int IntegrateStep();

    /// set the state of the nox group and the global state data container (implicit only)
    virtual void SetState(const Teuchos::RCP<Epetra_Vector> & x);

    /*! \brief nonlinear solve
     *
     *  Do the nonlinear solve, i.e. (multiple) corrector,
     *  for the time step. All boundary conditions have
     *  been set. */
    virtual INPAR::STR::ConvergenceStatus Solve();

    /** \brief Identify residual
     *
     *  This method does not predict the target solution but
     *  evaluates the residual and the stiffness matrix.
     *  In partitioned solution schemes, it is better to keep the current
     *  solution instead of evaluating the initial guess (as the predictor)
     *  does. */
    virtual void PreparePartitionStep();

    //! Prepare time step
    virtual void PrepareTimeStep();

    //! @name Accessors
    //! @{
    //! return the predictor
    const STR::PREDICT::Generic& Predictor() const
    {
      CheckInitSetup();
      return *predictor_ptr_;
    }

    //! @}
    Teuchos::RCP<const STR::NLN::SOLVER::Generic> GetNlnSolverPtr() const
    {
      return nlnsolver_ptr_;
    };

    //! do something in case nonlinear solution does not converge for some reason
     INPAR::STR::ConvergenceStatus PerformErrorAction(INPAR::STR::ConvergenceStatus nonlinsoldiv);


    //! check, if according to divercont flag time step size can be increased
    void CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status);

    //! returns pointer to generic implicit object
    const Teuchos::RCP<STR::IMPLICIT::Generic> ImplIntPtr()
    {
      CheckInitSetup();
      return implint_ptr_;
    };

    ///  Evaluate routine for coupled problems with monolithic approach
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc);
    virtual void Evaluate();

    /** \brief Print structural jacobian matrix into a text file for later use
     *  in MATLAB
     *
     *  This routine can be activated via the input parameter
     *  %STRUCT_JACOBIAN_MATLAB. See the corresponding inpar section for more
     *  details.
     *
     *  The text file can be found in the user-provided output folder using the
     *  following file name extension:
     *
     *  [OUTPUT-FOLDER]/[OUTPUT FILE NAME]_str_jacobian_step-[STEP]_nlniter-[NEWTON-ITERATION].mtl
     *
     *  \author hiermeier \date 06/17 */
    void PrintJacobianInMatlabFormat( const NOX::NLN::Group& curr_grp ) const;

  protected:
    //! returns the current solution group
    virtual const NOX::Abstract::Group& GetSolutionGroup() const;

    //! returns the current solution group ptr
    virtual Teuchos::RCP<NOX::Abstract::Group> SolutionGroupPtr();

    STR::IMPLICIT::Generic& ImplInt()
    {
      CheckInitSetup();
      return *implint_ptr_;
    };

    STR::PREDICT::Generic& Predictor()
    {
      CheckInitSetup();
      return *predictor_ptr_;
    };

    const Teuchos::RCP<STR::PREDICT::Generic> PredictorPtr()
    {
      CheckInitSetup();
      return predictor_ptr_;
    };

    const STR::NLN::SOLVER::Generic& NlnSolver() const
    {
      CheckInitSetup();
      return *nlnsolver_ptr_;
    };

    STR::NLN::SOLVER::Generic& NlnSolver()
    {
      CheckInitSetup();
      return *nlnsolver_ptr_;
    };

    const Teuchos::RCP<STR::NLN::SOLVER::Generic> NlnSolverPtr()
    {
      CheckInitSetup();
      return nlnsolver_ptr_;
    };



  private:
    //! ptr to the implicit time integrator object
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

    //! ptr to the predictor object
    Teuchos::RCP<STR::PREDICT::Generic> predictor_ptr_;

    //! ptr to the non-linear solver object
    Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_ptr_;

    //! ptr to the nox group object
    Teuchos::RCP<NOX::Abstract::Group> grp_ptr_;
  };
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_IMPLICIT_H_ */
