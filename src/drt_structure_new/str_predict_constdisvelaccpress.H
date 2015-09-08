/*
 * str_predict_constdisvelaccpress.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_PREDICT_CONSTDIS_H_
#define STR_PREDICT_CONSTDIS_H_

#include "str_predict_generic.H"

namespace STR {
namespace PREDICT {
  class ConstDisVelAccPress : public Generic
  {
  public:
    //! constructor
    ConstDisVelAccPress();

    //! destructor
    virtual ~ConstDisVelAccPress() {};

    //! do the class specific predictor step
    void PredictStep();

  protected:
    //! setup class specific stuff
    virtual void Setup();
  };
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_CONSTDIS_H_ */
