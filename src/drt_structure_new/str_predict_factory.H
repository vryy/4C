/*
 * str_predict_factory.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_PREDICT_FACTORY_H_
#define STR_PREDICT_FACTORY_H_

#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace STR {
namespace PREDICT {
  class Generic;

  class Factory
  {
  public:
    //! constructor
    Factory() {};

    //! destructor
    virtual ~Factory() {};

    //! build the desired predictor
    Teuchos::RCP<STR::PREDICT::Generic>  BuildPredictor(const enum INPAR::STR::PredEnum& predType) const;
  };

  //! non-member function
  Teuchos::RCP<STR::PREDICT::Generic> BuildPredictor(const enum INPAR::STR::PredEnum& predType);
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_FACTORY_H_ */
