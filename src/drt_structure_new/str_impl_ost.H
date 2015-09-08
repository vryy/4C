/*
 * str_impl_ost.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_IMPL_OST_H_
#define STR_IMPL_OST_H_

#include "str_impl_generic.H"

namespace STR {
namespace IMPLICIT {
  class OneStepTheta : public Generic
  {
  public:
    //! constructor
    OneStepTheta();

    //! destructor
    virtual ~OneStepTheta() {};
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_OST_H_ */
