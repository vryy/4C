/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired model evaluators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_structure_new_model_evaluator.hpp"  // typedef

#include <Teuchos_RCP.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Solid
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


      Teuchos::RCP<Solid::ModelEvaluator::Map> build_model_evaluators(
          const std::set<enum Inpar::Solid::ModelType>& modeltypes,
          const Teuchos::RCP<Solid::MODELEVALUATOR::Generic>& coupling_model_ptr) const;

     private:
      //! return the proper type for the contact model evaluator
      Teuchos::RCP<Solid::MODELEVALUATOR::Generic> build_contact_model_evaluator() const;

      //! return the proper type for the standard structural model evaluator
      Teuchos::RCP<Solid::MODELEVALUATOR::Generic> build_structure_model_evaluator() const;

    };  // class Factory

    //! non-member function, which relates to the Solid::MODELEVALUATOR::Factory
    Teuchos::RCP<Solid::ModelEvaluator::Map> build_model_evaluators(
        const std::set<enum Inpar::Solid::ModelType>& modeltypes,
        const Teuchos::RCP<Solid::MODELEVALUATOR::Generic>& coupling_model_ptr);

  }  // namespace MODELEVALUATOR
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
