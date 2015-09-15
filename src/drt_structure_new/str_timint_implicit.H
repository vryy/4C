/*
 * str_timint_implicit.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TIMINT_IMPLICIT_H_
#define STR_TIMINT_IMPLICIT_H_

#include "str_timint_base.H"
#include "str_nln_solver_interface_required.H"

// forward declarations ...
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
  class Implicit : public Base,
                   public STR::NLN::SOLVER::INTERFACE::Required
  {
  public:
    //! constructor
    Implicit();

    //! destructor
    virtual ~Implicit() {};

    virtual int IntegrateStep();

    virtual INPAR::STR::ConvergenceStatus Solve();

    //! Return the rhs-vector (negative sign for Newton is already included.)
    Teuchos::RCP<const Epetra_Vector> RHS()
    {
      return rhs_;
    }

  protected:
    virtual void Setup();

    STR::IMPLICIT::Generic& ImplInt()
    {
      CheckInitSetup();
      return *implint_;
    };

    const Teuchos::RCP<STR::IMPLICIT::Generic> ImplIntPtr()
    {
      CheckInitSetup();
      return implint_;
    };

    STR::PREDICT::Generic& Predictor()
    {
      CheckInitSetup();
      return *predictor_;
    };

    const Teuchos::RCP<STR::PREDICT::Generic> PredictorPtr()
    {
      CheckInitSetup();
      return predictor_;
    };

    STR::NLN::SOLVER::Generic& NlnSolver()
    {
      CheckInitSetup();
      return *nlnsolver_;
    };

    const Teuchos::RCP<STR::NLN::SOLVER::Generic> NlnSolverPtr()
    {
      CheckInitSetup();
      return nlnsolver_;
    };

    //! Prepare time step
    virtual void PrepareTimeStep()
    {
      // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
      DataGlobalStatePtr()->SetTimenp(DataGlobalStatePtr()->GetUpdatedTime());

      // things that need to be done before Predict
      PrePredict();

      // TODO prepare contact for new time step
      // PrepareStepContact();
      Predictor().Predict();

      // things that need to be done after Predict
      PostPredict();

      return;
    }

  private:
    //! ptr to the implicit time integrator object
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_;

    //! ptr to the predictor object
    Teuchos::RCP<STR::PREDICT::Generic> predictor_;

    //! ptr to the non-linear solver object
    Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_;
  };
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_IMPLICIT_H_ */
