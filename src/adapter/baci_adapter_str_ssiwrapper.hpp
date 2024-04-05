/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for SSI problems


\level 1
*/

#ifndef FOUR_C_ADAPTER_STR_SSIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_SSIWRAPPER_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_wrapper.hpp"
#include "baci_ssi_str_model_evaluator_partitioned.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}  // namespace STR

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedSSI;
  }
}  // namespace STR


namespace ADAPTER
{
  class SSIStructureWrapper : public StructureWrapper
  {
   public:
    /// constructor
    explicit SSIStructureWrapper(Teuchos::RCP<Structure> structure);

    /// set pointer to model evaluator
    void SetModelEvaluatorPtr(Teuchos::RCP<STR::MODELEVALUATOR::PartitionedSSI> me)
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
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedSSI> ssi_model_evaluator_;

    /// access the fsi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedSSI> SSIModelEvaluator()
    {
      return ssi_model_evaluator_;
    };

  };  // class SSIStructureWrapper
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_STR_SSIWRAPPER_H
