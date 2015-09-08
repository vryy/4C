/*
 * str_impl_gemm.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_IMPL_GEMM_H_
#define STR_IMPL_GEMM_H_

#include "str_impl_generic.H"

namespace STR {
namespace IMPLICIT {
  class GEMM : public Generic
  {
  public:
    //! constructor
    GEMM() {};

    //! destructor
    virtual ~GEMM() {};
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_GEMM_H_ */
