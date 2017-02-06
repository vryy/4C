/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_factory.H

\brief Factory to create the desired model evaluators.

\maintainer Michael Hiermeier

\date Sep 10, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_MODEL_EVALUATOR_FACTORY_H_
#define STR_MODEL_EVALUATOR_FACTORY_H_


#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_structure.H"
#include "str_model_evaluator.H"  // typedef

namespace STR {
namespace MODELEVALUATOR {
  class Generic;

  /*! Factory to build the desired model evaluator std::map
   *
   *  \author Michael Hiermeier */
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};


    Teuchos::RCP<STR::ModelEvaluator::Map>
    BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes,
        const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr) const;

  private:
    //! return the proper type for the contact model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> BuildContactModelEvaluator() const;

    //! return the proper type for the standard structural model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::Generic> BuildStructureModelEvaluator() const;

  };  // class Factory

  //! non-member function, which relates to the STR::MODELEVALUATOR::Factory
  Teuchos::RCP<STR::ModelEvaluator::Map>
  BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes,
      const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr);

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_FACTORY_H_ */
