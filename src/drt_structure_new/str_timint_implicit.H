/*
 * str_timint_implicit.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TIMINT_IMPLICIT_H_
#define STR_TIMINT_IMPLICIT_H_

#include "str_timint_base.H"

// forward declarations ...
namespace STR {
namespace IMPLICIT {
  class Generic;
} // namespace IMPLICIT
namespace ModelEvaluator {
  class Generic;
} // namespace  ModelEvaluator
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

    std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >& ModelEvaluators()
    {
      CheckInitSetup();
      return *modelevaluators_;
    }

    const Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > ModelEvaluatorsPtr()
    {
      CheckInitSetup();
      return modelevaluators_;
    }

    STR::PREDICT::Generic& Predictor()
    {
      CheckInitSetup();
      return *predictor_;
    };

    const Teuchos::RCP<STR::NLN::SOLVER::Generic> PredictorPtr()
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

  private:
    //! ptr to the implicit time integrator object
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_;

    //! ptr to the model evaluator map
    Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > modelevaluators_;

    //! ptr to the predictor object
    Teuchos::RCP<STR::PREDICT::Generic> predictor_;

    //! ptr to the non-linear solver object
    Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_;
  };
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_IMPLICIT_H_ */
