/*----------------------------------------------------------------------*/
/*! \file

\brief structural adapter for PASI problems

\level 3

\maintainer  Sebastian Fuchs


*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "ad_str_pasiwrapper.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 | pasi adapter                                          sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
ADAPTER::PASIStructureWrapper::PASIStructureWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
  // set-up PASI interface
  interface_ = Teuchos::rcp(new STR::AUX::MapExtractor);

  interface_->Setup(*Discretization(), *Discretization()->DofRowMap());

}  // ADAPTER::PASIStructureWrapper::PASIStructureWrapper()

/*----------------------------------------------------------------------*
 | apply particle wall force to structure interface      sfuchs 03/2017 |
 *----------------------------------------------------------------------*/
void ADAPTER::PASIStructureWrapper::ApplyInterfaceForce(Teuchos::RCP<Epetra_Vector> wallforce)
{
  PASIModelEvaluator()->GetInterfaceForceNpPtr()->Scale(0.0);

  if (wallforce != Teuchos::null)
    interface_->AddPASICondVector(wallforce, PASIModelEvaluator()->GetInterfaceForceNpPtr());

}  // ADAPTER::PASIStructureWrapper::ApplyInterfaceForce()
