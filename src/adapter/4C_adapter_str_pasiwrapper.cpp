/*----------------------------------------------------------------------*/
/*! \file

\brief structural adapter for PASI problems

\level 3



*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_adapter_str_pasiwrapper.hpp"

#include "4C_lib_discret.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/
ADAPTER::PASIStructureWrapper::PASIStructureWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
  // set-up PASI interface
  interface_ = Teuchos::rcp(new STR::MapExtractor);

  interface_->Setup(*Discretization(), *Discretization()->DofRowMap());
}

void ADAPTER::PASIStructureWrapper::ApplyInterfaceForce(Teuchos::RCP<const Epetra_Vector> intfforce)
{
  PASIModelEvaluator()->GetInterfaceForceNpPtr()->PutScalar(0.0);

  if (intfforce != Teuchos::null)
    interface_->AddPASICondVector(intfforce, PASIModelEvaluator()->GetInterfaceForceNpPtr());
}

FOUR_C_NAMESPACE_CLOSE
