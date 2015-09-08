/*
 * str_expl_centrdiff.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_EXPL_CENTRDIFF_H_
#define STR_EXPL_CENTRDIFF_H_

#include "str_expl_generic.H"

namespace STR {
namespace EXPLICIT {
  class CentrDiff : public Generic
  {
  public:
    //! constructor
    CentrDiff();

    //! destructor
    virtual ~CentrDiff() {};
  };
} // namespace EXPLICIT
} // namespace STR


#endif /* STR_EXPL_CENTRDIFF_H_ */
