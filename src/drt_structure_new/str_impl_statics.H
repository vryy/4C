/*
 * str_impl_statics.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_IMPL_STATICS_H_
#define STR_IMPL_STATICS_H_

#include "str_impl_generic.H"

namespace STR {
namespace IMPLICIT {
  class Statics : public Generic
  {
  public:
    //! constructor
    Statics();

    //! destructor
    virtual ~Statics() {};
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_STATICS_H_ */
