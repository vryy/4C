/*
 * str_expl_forwardeuler.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_EXPL_FORWARDEULER_H_
#define STR_EXPL_FORWARDEULER_H_

#include "str_expl_generic.H"

namespace STR {
namespace EXPLICIT {
  class ForwardEuler : public Generic
  {
  public:
    //! constructor
    ForwardEuler();

    //! destructor
    virtual ~ForwardEuler() {};
  };
} // namespace EXPLICIT
} // namespace STR


#endif /* STR_EXPL_FORWARDEULER_H_ */
