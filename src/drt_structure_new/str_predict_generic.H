/*
 * str_predict_generic.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_PREDICT_GENERIC_H_
#define STR_PREDICT_GENERIC_H_

#include "str_predict_generic.H"

namespace STR {
namespace PREDICT {
  class Generic
  {
  public:
      Generic();

      virtual ~Generic() {};

      //! Get the predictor type enum
      const INPAR::STR::PredEnum& GetType() const { return type_; };

      //! returns the name of the used predictor
      virtual const std::string& Name() const;

      //! Calculate the specific predictor step
      void PredictStep() = 0;

  protected:
      //! initialize the base class variables
      virtual void Init(const INPAR::STR::PredEnum& type);

      //! setup of the specific predictor
      virtual void Setup() = 0;

      //! returns init state
      const bool& IsInit() const { return isInit_; };

      //! returns setup state
      const bool& IsSetup() const { return isSetup_; };

  protected:
      //! indicates if the Init() function has been called
      bool isInit_;

      //! indicates if the Setup() function has been called
      bool isSetup_;
  private:
      //! predictor type
      INPAR::STR::PredEnum type_;
  };
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_GENERIC_H_ */
