/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_CONSTDISVELACCPRESS_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_CONSTDISVELACCPRESS_HPP

#include "4C_config.hpp"

#include "4C_structure_new_predict_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace Predict
  {
    class ConstDisVelAccPress : public Generic
    {
     public:
      //! constructor
      ConstDisVelAccPress();

      //! setup class specific stuff
      void Setup() override;

      //! do the class specific predictor step
      void Compute(::NOX::Abstract::Group& grp) override;

     private:
      Teuchos::RCP<STR::Predict::Generic> tangdis_ptr_;
    };  // class ConstDisVelAccPress
  }     // namespace Predict
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
