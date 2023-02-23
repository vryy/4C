/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired model evaluators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_H
#define STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_H

#include <set>

#include <Teuchos_RCP.hpp>
#include "inpar_structure.H"
#include "structure_new_model_evaluator.H"  // typedef

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
      virtual ~Factory(){};


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

#endif  // STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_H
