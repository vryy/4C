/*----------------------------------------------------------------------*/
/*!

\brief Structural adapter for SSI problems

\maintainer Christoph Schmidt

\level 1
*/

#ifndef SSI_STRUCTURE_WRAPPER_H
#define SSI_STRUCTURE_WRAPPER_H

#include "ad_str_wrapper.H"
#include "../drt_ssi/ssi_str_model_evaluator_partitioned.H"

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
    virtual ~SSIStructureWrapper() { ; }

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

#endif /* SSI_STRUCTURE_WRAPPER_H */
