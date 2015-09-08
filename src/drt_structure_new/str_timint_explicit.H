/*
 * str_timint_explicit.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TIMINT_EXPLICIT_H_
#define STR_TIMINT_EXPLICIT_H_


#include "str_timint_base.H"

namespace STR {
namespace TIMINT {
  class Explicit : public Base
  {
  public:
    //! constructor
    Explicit();

    //! destructor
    virtual ~Explicit() {};

    virtual void IntegrateStep(){};
  };
} // namespace TimInt
} // namespace STR

#endif /* STR_TIMINT_EXPLICIT_H_ */
