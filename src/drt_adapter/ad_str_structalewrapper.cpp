/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Structure-ALE problems.



\level 3
*/
/*----------------------------------------------------------------------*/


#include "ad_str_structalewrapper.H"
#include "ad_str_structure_new.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructAleWrapper::StructAleWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure),
      structure_(Teuchos::rcp_dynamic_cast<ADAPTER::StructureNew>(structure, true)),
      struct_ale_model_evaluator_(
          Teuchos::rcpFromRef(structure_->ModelEvaluator(INPAR::STR::model_structure)))
{
  // empty
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector>& ADAPTER::StructAleWrapper::GetMaterialDisplacementNpPtr()
{
  return GetStructAleModelEvaluatorPtr()->GetMaterialDisplacementNpPtr();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructAleWrapper::UpdateMaterialDisplacements(Teuchos::RCP<Epetra_Vector> dispmat)
{
  GetStructAleModelEvaluatorPtr()->GetMaterialDisplacementNpPtr()->Update(1.0, *dispmat, 0.0);
}
