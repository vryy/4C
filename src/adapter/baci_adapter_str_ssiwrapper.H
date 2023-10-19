/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for SSI problems


\level 1
*/

#ifndef BACI_ADAPTER_STR_SSIWRAPPER_H
#define BACI_ADAPTER_STR_SSIWRAPPER_H

#include "baci_adapter_str_wrapper.H"
#include "baci_ssi_str_model_evaluator_partitioned.H"

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

    /// Destructor
    ~SSIStructureWrapper() override { ; }

    /// set pointer to model evaluator
    void SetModelEvaluatorPtr(Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedSSI> me)
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
    Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedSSI> ssi_model_evaluator_;

    /// access the fsi model evaluator
    Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedSSI> SSIModelEvaluator()
    {
      return ssi_model_evaluator_;
    };

  };  // class SSIStructureWrapper
}  // namespace ADAPTER

#endif  // ADAPTER_STR_SSIWRAPPER_H
