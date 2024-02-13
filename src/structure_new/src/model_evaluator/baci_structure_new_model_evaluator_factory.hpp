/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired model evaluators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP
#define BACI_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"
#include "baci_structure_new_model_evaluator.hpp"  // typedef

#include <Teuchos_RCP.hpp>

#include <set>

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace MODELEVALUATOR
  {
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
      virtual ~Factory() = default;


      Teuchos::RCP<STR::ModelEvaluator::Map> BuildModelEvaluators(
          const std::set<enum INPAR::STR::ModelType>& modeltypes,
          const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr) const;

     private:
      //! return the proper type for the contact model evaluator
      Teuchos::RCP<STR::MODELEVALUATOR::Generic> BuildContactModelEvaluator() const;

      //! return the proper type for the standard structural model evaluator
      Teuchos::RCP<STR::MODELEVALUATOR::Generic> BuildStructureModelEvaluator() const;

    };  // class Factory

    //! non-member function, which relates to the STR::MODELEVALUATOR::Factory
    Teuchos::RCP<STR::ModelEvaluator::Map> BuildModelEvaluators(
        const std::set<enum INPAR::STR::ModelType>& modeltypes,
        const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr);

  }  // namespace MODELEVALUATOR
}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_H
