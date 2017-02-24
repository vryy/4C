/*!----------------------------------------------------------------------
\file ad_str_pasiwrapper.H

\brief structural adapter for PASI problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_STRUCTURE_WRAPPER_H
#define PASI_STRUCTURE_WRAPPER_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "ad_str_wrapper.H"
#include "../drt_pasi/pasi_str_model_evaluator_partitioned.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace STR
{
namespace MODELEVALUATOR
{
class PartitionedPASI;
}
}

/*----------------------------------------------------------------------*
 | pasi adapter                                          sfuchs 01/2017 |
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
  void SetModelEvaluatorPtr(Teuchos::RCP< ::STR::MODELEVALUATOR::PartitionedPASI> me)
  {
    pasi_model_evaluator_=me;
    return;
  }

protected:

  //! The structural model evaluator object.
  //! Your PASI algorithm calls methods in this adapter.
  //! If this method is related to the structural field,
  //! a corresponding method in the model evaluator may be
  //! called, if necessary.
  Teuchos::RCP< ::STR::MODELEVALUATOR::PartitionedPASI> pasi_model_evaluator_;

  //! access the pasi model evaluator
  Teuchos::RCP< ::STR::MODELEVALUATOR::PartitionedPASI> PASIModelEvaluator(){return pasi_model_evaluator_;};

}; // class PASIStructureWrapper : public StructureWrapper

} // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif // PASI_STRUCTURE_WRAPPER_H
