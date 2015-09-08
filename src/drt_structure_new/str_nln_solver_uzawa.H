/*
 * str_nln_solver_Uzawa.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_UZAWA_H_
#define STR_NLN_SOLVER_UZAWA_H_

#include "str_nln_solver_generic.H"

namespace STR {
namespace NLN {
namespace SOLVER {
  class Uzawa : public Generic
  {
  public:
    Uzawa() {};

    virtual ~Uzawa() {};
  };
} // namespace SOLVER
} // namespace NLN
} // namespace STR

#endif /* STR_NLN_SOLVER_UZAWA_H_ */
