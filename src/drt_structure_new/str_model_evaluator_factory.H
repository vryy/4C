/*
 * str_model_evaluator_factory.H
 *
 *  Created on: Sep 10, 2015
 *      Author: hiermeier
 */

#ifndef STR_MODEL_EVALUATOR_FACTORY_H_
#define STR_MODEL_EVALUATOR_FACTORY_H_


#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_structure.H"

namespace STR {
namespace MODELEVALUATOR {
class Generic;

  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > >
    BuildModelEvaluators(const std::vector<const enum INPAR::STR::ModelType>& modeltypes) const;

  private:
  };  // class Factory

  //! non-member function
  Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > >
  BuildModelEvaluators(const std::vector<const enum INPAR::STR::ModelType>& modeltypes);

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_FACTORY_H_ */
