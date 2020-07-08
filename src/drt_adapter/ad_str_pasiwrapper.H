/*----------------------------------------------------------------------*/
/*! \file

\brief structural adapter for PASI problems

\level 3



*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/
#ifndef PASI_STRUCTURE_WRAPPER_H
#define PASI_STRUCTURE_WRAPPER_H

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "ad_str_wrapper.H"
#include "../drt_pasi/pasi_str_model_evaluator_partitioned.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }  // namespace AUX
  namespace MODELEVALUATOR
  {
    class PartitionedPASI;
  }  // namespace MODELEVALUATOR
}  // namespace STR

/*----------------------------------------------------------------------*
 | class declarations                                                   |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class PASIStructureWrapper : public StructureWrapper
  {
   public:
    //! constructor
    explicit PASIStructureWrapper(Teuchos::RCP<Structure> structure);

    //! Destructor
    virtual ~PASIStructureWrapper(){};

    //! set pointer to model evaluator
    void SetModelEvaluatorPtr(Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedPASI> me)
    {
      pasi_model_evaluator_ = me;
      return;
    }

    //! communication object at the interface
    virtual Teuchos::RCP<const ::STR::AUX::MapExtractor> Interface() const { return interface_; }

    //! apply interface force to structure interface
    void ApplyInterfaceForce(Teuchos::RCP<const Epetra_Vector> intfforce);

   protected:
    //! The structural model evaluator object.
    //! Your PASI algorithm calls methods in this adapter.
    //! If this method is related to the structural field,
    //! a corresponding method in the model evaluator may be
    //! called, if necessary.
    Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedPASI> pasi_model_evaluator_;

    //! access the pasi model evaluator
    Teuchos::RCP<::STR::MODELEVALUATOR::PartitionedPASI> PASIModelEvaluator()
    {
      return pasi_model_evaluator_;
    };

    //! the interface map setup for interface <-> full translation
    Teuchos::RCP<::STR::AUX::MapExtractor> interface_;
  };

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
