/*
 * str_nln_solver_FullNewton.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_FULLNEWTON_H_
#define STR_NLN_SOLVER_FULLNEWTON_H_

#include "str_nln_solver_generic.H"

namespace STR {
namespace NLN {
namespace SOLVER {
  class FullNewton : public Generic
  {
  public:
    FullNewton() {};

    virtual ~FullNewton() {};
  };
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_FULLNEWTON_H_ */
