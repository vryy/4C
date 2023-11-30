/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_PREDICT_CONSTDISVELACCPRESS_H
#define BACI_STRUCTURE_NEW_PREDICT_CONSTDISVELACCPRESS_H

#include "baci_structure_new_predict_generic.H"

namespace STR
{
  namespace PREDICT
  {
    class ConstDisVelAccPress : public Generic
    {
     public:
      //! constructor
      ConstDisVelAccPress();


      //! setup class specific stuff
      void Setup() override;

      //! do the class specific predictor step
      void Compute(NOX::Abstract::Group& grp) override;

     private:
      Teuchos::RCP<STR::PREDICT::Generic> tangdis_ptr_;
    };  // class ConstDisVelAccPress
  }     // namespace PREDICT
}  // namespace STR


#endif  // STRUCTURE_NEW_PREDICT_CONSTDISVELACCPRESS_H
