/*
 * str_predict_generic.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_PREDICT_GENERIC_H_
#define STR_PREDICT_GENERIC_H_

#include "../drt_inpar/inpar_structure.H"

namespace STR {
namespace PREDICT {
  class Generic
  {
  public:
      Generic();

      virtual ~Generic() {};

      //! initialize the base class variables
      virtual void Init(const enum INPAR::STR::PredEnum& type);

      //! setup of the specific predictor
      virtual void Setup() = 0;

      //! Get the predictor type enum
      const INPAR::STR::PredEnum& GetType() const { return type_; };

      //! returns the name of the used predictor
      virtual const std::string& Name() const;

      //! Calculate the specific predictor step
      virtual void Predict() = 0;

  protected:
      //! returns init state
      const bool& IsInit() const { return isinit_; };

      //! returns setup state
      const bool& IsSetup() const { return issetup_; };

  protected:
      //! indicates if the Init() function has been called
      bool isinit_;

      //! indicates if the Setup() function has been called
      bool issetup_;

  private:
      //! predictor type
      enum INPAR::STR::PredEnum type_;
  };
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_GENERIC_H_ */
