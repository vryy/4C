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

  };
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_IMPLICIT_H_ */
