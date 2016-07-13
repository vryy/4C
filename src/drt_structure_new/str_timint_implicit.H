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

    /*! \brief linear structure solve with just a interface load
     *
     * The very special solve done in steepest descent relaxation
     * calculation (and matrix free Newton Krylov).
     *
     * \note Can only be called after a valid structural solve. */
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear();

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
      CheckInitSetup();
      return nlnsolver_ptr_;
    };

    //! do something in case nonlinear solution does not converge for some reason
     INPAR::STR::ConvergenceStatus PerformErrorAction(INPAR::STR::ConvergenceStatus nonlinsoldiv);

    //! check, if according to divercont flag time step size can be increased
    void CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status);

  protected:
    // returns the current solution group
    virtual const NOX::Abstract::Group& GetSolutionGroup() const;

    STR::IMPLICIT::Generic& ImplInt()
    {
      CheckInitSetup();
      return *implint_ptr_;
    };

    const Teuchos::RCP<STR::IMPLICIT::Generic> ImplIntPtr()
    {
      CheckInitSetup();
      return implint_ptr_;
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
