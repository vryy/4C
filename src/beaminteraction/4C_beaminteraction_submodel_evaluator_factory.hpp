// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_FACTORY_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_str_model_evaluator.hpp"  // typedef

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

      Teuchos::RCP<Solid::ModelEvaluator::BeamInteraction::Map> build_model_evaluators(
          const std::set<enum Inpar::BEAMINTERACTION::SubModelType>& submodeltypes) const;

     private:
    };

    //! non-member function, which relates to the Solid::ModelEvaluator::Factory
    Teuchos::RCP<Solid::ModelEvaluator::BeamInteraction::Map> build_model_evaluators(
        const std::set<enum Inpar::BEAMINTERACTION::SubModelType>& submodeltypes);

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
