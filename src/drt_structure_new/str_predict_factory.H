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
    Factory() {};

    virtual ~Factory() {};

    Teuchos::RCP<STR::PREDICT::Generic>  BuildPredictor(const INPAR::STR::PredEnum& predType) const;
  };

  Teuchos::RCP<STR::PREDICT::Generic> BuildPredictor(const INPAR::STR::PredEnum& predType);
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_FACTORY_H_ */
