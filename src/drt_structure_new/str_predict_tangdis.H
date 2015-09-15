/*
 * str_predict_tangdis.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_PREDICT_TANGDIS_H_
#define STR_PREDICT_TANGDIS_H_

#include "str_predict_generic.H"

namespace STR {
namespace PREDICT {
  class TangDis : public Generic
  {
  public:
    TangDis() {};

    virtual ~TangDis() {};

    //! do the class specific predictor step
    virtual void Predict();

  protected:
    //! setup class specific stuff
    virtual void Setup();
  };
} // namespace PREDICT
} // namespace STR

#endif /* STR_PREDICT_TANGDIS_H_ */
