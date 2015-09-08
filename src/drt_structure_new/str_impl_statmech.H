/*
 * str_impl_statmech.H
 *
 *  Created on: Sep 2, 2015
 *      Author: hiermeier
 */

#ifndef STR_IMPL_STATMECH_H_
#define STR_IMPL_STATMECH_H_

#include "str_impl_ost.H"

namespace STR {
namespace IMPLICIT {
  class StatMech : public OneStepTheta
  {
  public:
    //! constructor
    StatMech();

    //! destructor
    virtual ~StatMech() {};
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_STATMECH_H_ */
