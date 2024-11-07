// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_SSIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_SSIWRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"
#include "4C_ssi_str_model_evaluator_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Solid
{
  namespace Aux
  {
    class MapExtractor;
  }
}  // namespace Solid

namespace Solid
{
  namespace ModelEvaluator
  {
    class PartitionedSSI;
  }
}  // namespace Solid


namespace Adapter
{
  class SSIStructureWrapper : public StructureWrapper
  {
   public:
    /// constructor
    explicit SSIStructureWrapper(std::shared_ptr<Structure> structure);

    /// set pointer to model evaluator
    void set_model_evaluator_ptr(std::shared_ptr<Solid::ModelEvaluator::PartitionedSSI> me)
    {
      ssi_model_evaluator_ = me;
      return;
    }

   protected:
    /// The structural model evaluator object.
    /// Your SSI algorithm calls methods in this adapter.
    /// If this method is related to the structural field,
    /// a corresponding method in the model evaluator may be
    /// called, if necessary.
    std::shared_ptr<Solid::ModelEvaluator::PartitionedSSI> ssi_model_evaluator_;

    /// access the fsi model evaluator
    std::shared_ptr<Solid::ModelEvaluator::PartitionedSSI> ssi_model_evaluator()
    {
      return ssi_model_evaluator_;
    };

  };  // class SSIStructureWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
