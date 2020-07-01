/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_PREDICT_CONSTDIS_H_
#define STR_PREDICT_CONSTDIS_H_

#include "str_predict_generic.H"

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
      virtual void Setup();

      //! do the class specific predictor step
      virtual void Compute(NOX::Abstract::Group& grp);

     private:
      Teuchos::RCP<STR::PREDICT::Generic> tangdis_ptr_;
    };  // class ConstDisVelAccPress
  }     // namespace PREDICT
}  // namespace STR


#endif /* STR_PREDICT_CONSTDIS_H_ */
