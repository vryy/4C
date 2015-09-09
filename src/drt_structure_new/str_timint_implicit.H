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
} // namespace SOLVER
} // namespace NLN

namespace TIMINT {
  class Implicit : public Base
  {
  public:
    //! constructor
    Implicit();

    //! destructor
    virtual ~Implicit() {};

    virtual void IntegrateStep();

    void CreateNlnSolver();

    void CreatePredictor();

  protected:
    virtual void Setup();

    STR::IMPLICIT::Generic& ImplInt()
    {
      CheckInitSetup();
      return *implint_;
    };

    Teuchos::RCP<STR::IMPLICIT::Generic> ImplIntPtr()
    {
      CheckInitSetup();
      return implint_;
    };

    std::map<INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >& ModelEvaluators()
    {
      CheckInitSetup();
      return *modelevaluators_;
    }

    Teuchos::RCP<std::map<INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > ModelEvaluatorsPtr()
    {
      CheckInitSetup();
      return modelevaluators_;
    }

    STR::NLN::SOLVER::Generic& NlnSolver()
    {
      CheckInitSetup();
      return *nlnsolver_;
    };

    Teuchos::RCP<STR::NLN::SOLVER::Generic> NlnSolverPtr()
    {
      CheckInitSetup();
      return nlnsolver_;
    };

  private:
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_;
    Teuchos::RCP<std::map<INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > modelevaluators_;
    Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolver_;
  };
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_IMPLICIT_H_ */
