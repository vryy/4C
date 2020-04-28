/*-----------------------------------------------------------*/
/*! \file

\brief This file contains the declaration of a time integrator for prestressing

\level 3

\maintainer Amadeus Gebauer

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_PRESTRESS_H_
#define STR_IMPL_PRESTRESS_H_

#include "str_impl_statics.H"

namespace STR
{
  namespace IMPLICIT
  {
    /*!
     * \brief The time integrator
     *
     */
    class PreStress : public Statics
    {
     public:
      //! constructor
      PreStress() = default;

      //! destructor
      ~PreStress() override = default;

      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;
    };
  }  // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_PRESTRESS_H_ */
