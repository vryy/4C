/*-----------------------------------------------------------*/
/*!
\file str_impl_prestress.H

\maintainer Michael Hiermeier

\date Sep 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_PRESTRESS_H_
#define STR_IMPL_PRESTRESS_H_

#include "str_impl_statics.H"

namespace STR
{
  namespace IMPLICIT
  {
    class PreStress : public Statics
    {
     public:
      //! constructor
      PreStress(){};

      //! destructor
      virtual ~PreStress(){};
    };
  }  // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_PRESTRESS_H_ */
