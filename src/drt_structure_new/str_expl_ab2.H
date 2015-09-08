/*
 * str_expl_ab2.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_EXPL_AB2_H_
#define STR_EXPL_AB2_H_

#include "str_expl_generic.H"

namespace STR {
namespace EXPLICIT {
  class AdamsBashforth2 : public Generic
  {
  public:
    //! constructor
    AdamsBashforth2();

    //! destructor
    virtual ~AdamsBashforth2(){};
  };
} // namespace EXPLICIT
} // namespace STR



#endif /* STR_EXPL_AB2_H_ */
