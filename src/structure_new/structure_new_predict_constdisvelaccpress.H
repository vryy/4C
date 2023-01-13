/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_PREDICT_CONSTDIS_H_
#define STR_PREDICT_CONSTDIS_H_

#include "structure_new_predict_generic.H"

namespace STR
{
  namespace PREDICT
  {
    class ConstDisVelAccPress : public Generic
    {
     public:
      //! constructor
      ConstDisVelAccPress();

      //! destructor
      virtual ~ConstDisVelAccPress(){};

      //! setup class specific stuff
      void Setup() override;

      //! do the class specific predictor step
      void Compute(NOX::Abstract::Group& grp) override;

     private:
      Teuchos::RCP<STR::PREDICT::Generic> tangdis_ptr_;
    };  // class ConstDisVelAccPress
  }     // namespace PREDICT
}  // namespace STR


#endif /* STR_PREDICT_CONSTDIS_H_ */
