/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired submodel evaluators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_HPP


#include "baci_config.hpp"

#include "baci_beaminteraction_str_model_evaluator.hpp"  // typedef

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  namespace SUBMODELEVALUATOR
  {
    class Generic;

    /*! Factory to build the desired model evaluator std::map
     *
     *  \author Jonas Eichinger */
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory() = default;

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> BuildModelEvaluators(
          const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes) const;

     private:
    };

    //! non-member function, which relates to the STR::MODELEVALUATOR::Factory
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> BuildModelEvaluators(
        const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes);

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
