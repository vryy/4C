/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for Structure-ALE problems.


\level 3
*/
/*-----------------------------------------------------------*/


#include "struct_ale_str_model_evaluator.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::StructAle::StructAle() : material_displacements_np_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::StructAle::Setup()
{
  // call Setup() in base class
  STR::MODELEVALUATOR::Structure::Setup();

  // construct material displacements at \f$t_{n+1}\f$
  material_displacements_np_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView(), true));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::StructAle::PreEvaluateInternal()
{
  // set state
  DiscretPtr()->SetState(0, "material_displacement", material_displacements_np_);
  return;
}
