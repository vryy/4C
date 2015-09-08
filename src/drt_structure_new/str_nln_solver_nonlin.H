/*
 * str_nln_solver_nonlin.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_NONLIN_H_
#define STR_NLN_SOLVER_NONLIN_H_

#include "str_nln_solver_generic.H"

namespace STR {
namespace NLN {
namespace SOLVER {
  class Nonlin : public Generic
  {
  public:
    Nonlin() {};

    virtual ~Nonlin() {};
  };
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_NONLIN_H_ */
