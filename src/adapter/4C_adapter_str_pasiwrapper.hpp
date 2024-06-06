/*----------------------------------------------------------------------*/
/*! \file

\brief structural adapter for PASI problems

\level 3



*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_STR_PASIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_PASIWRAPPER_HPP

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"
#include "4C_pasi_str_model_evaluator_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace STR
{
  class MapExtractor;

  namespace MODELEVALUATOR
  {
    class PartitionedPASI;
  }  // namespace MODELEVALUATOR
}  // namespace STR

/*----------------------------------------------------------------------*
 | class declarations                                                   |
 *----------------------------------------------------------------------*/
namespace Adapter
{
  class PASIStructureWrapper : public StructureWrapper
  {
   public:
    //! constructor
    explicit PASIStructureWrapper(Teuchos::RCP<Structure> structure);

    //! set pointer to model evaluator
    void set_model_evaluator_ptr(Teuchos::RCP<STR::MODELEVALUATOR::PartitionedPASI> me)
    {
      pasi_model_evaluator_ = me;
      return;
    }

    //! communication object at the interface
    virtual Teuchos::RCP<const STR::MapExtractor> Interface() const { return interface_; }

    //! apply interface force to structure interface
    void ApplyInterfaceForce(Teuchos::RCP<const Epetra_Vector> intfforce);

   protected:
    //! The structural model evaluator object.
    //! Your PASI algorithm calls methods in this adapter.
    //! If this method is related to the structural field,
    //! a corresponding method in the model evaluator may be
    //! called, if necessary.
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedPASI> pasi_model_evaluator_;

    //! access the pasi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedPASI> pasi_model_evaluator()
    {
      return pasi_model_evaluator_;
    };

    //! the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::MapExtractor> interface_;
  };

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
