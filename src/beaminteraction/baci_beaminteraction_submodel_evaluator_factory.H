/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired submodel evaluators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_H
#define BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_H


#include "baci_beaminteraction_str_model_evaluator.H"  // typedef

#include <Teuchos_RCP.hpp>

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
      virtual ~Factory(){};

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> BuildModelEvaluators(
          const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes) const;

     private:
    };

    //! non-member function, which relates to the STR::MODELEVALUATOR::Factory
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> BuildModelEvaluators(
        const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes);

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
