/*
 * str_model_evaluator_structure.cpp
 *
 *  Created on: Aug 11, 2015
 *      Author: farah
 */

#include "str_model_evaluator_structure.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Structure::Structure()
    : STR::MODELEVALUATOR::Generic(),
      discret_(Teuchos::null),
      fintn_(Teuchos::null),
      fintnp_(Teuchos::null),
      fextn_(Teuchos::null),
      fextnp_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Init()
{
  isinit_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  issetup_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ApplyForce(
    Teuchos::ParameterList& p,
    const Teuchos::RCP<Epetra_Vector>& fresnp,
    Teuchos::RCP<const Epetra_Vector> disnp,
    Teuchos::RCP<const Epetra_Vector> disn)
{
  // ---------------------------------------
  // (1) EXTERNAL FORCES
  // ---------------------------------------
  // initialize external forces
  fextnp_->PutScalar(0.0);
  ApplyForceExternal(p, disnp, disn);

  // ---------------------------------------
  // (2) INTERNAL FORCES
  // ---------------------------------------
  // initialize internal forces
  fintnp_->PutScalar(0.0);

  // ordinary internal force
  ApplyForceInternal(p, disnp, disn);

  // ******************** Finally, put everything together ********************

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fresnp->Update(-1.0, *fextnp_, 0.0);
  fresnp->Update(1.0, *fintnp_, 1.0);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ApplyForceInternal(
    Teuchos::ParameterList& p,
    Teuchos::RCP<const Epetra_Vector> disnp,
    Teuchos::RCP<const Epetra_Vector> disn)
{
  // action for elements
  std::string action = "calc_struct_internalforce";
  p.set("action", action);
//  // other parameters that might be needed by the elements
//  p.set("total time", time);
//  p.set("delta time", dt);
//  p.set<int>("young_temp", young_temp_);

  // set vector values needed by elements
  discret_->ClearState();
//  discret_->SetState("residual displacement", disi);  // these are incremental
  discret_->SetState(0,"displacement", disnp);

  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     fintnp_, Teuchos::null, Teuchos::null);

  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ApplyForceExternal(
    Teuchos::ParameterList& p,
    Teuchos::RCP<const Epetra_Vector> disnp,
    Teuchos::RCP<const Epetra_Vector> disn)
{
//  // set target time t_n+1
//  Teuchos::ParameterList p;
//  p.set("total time", time);
  // Set to default value, because it is unnecessary for the
  // EvaluateNeumann routine.
  std::string action = "calc_none";
  p.set("action", action);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0, "displacement old", disn);
  discret_->SetState(0, "displacement", disnp);
  discret_->EvaluateNeumann(p, *fextnp_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ApplyForceAndJacobian(
    const Teuchos::RCP<LINALG::SparseOperator>& stiff,
    const Teuchos::RCP<Epetra_Vector>& fresnp,
    Teuchos::RCP<const Epetra_Vector> disnp,
    Teuchos::RCP<const Epetra_Vector> disn)
{

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetInternalForceN() const
{
  if (fintn_.is_null())
    dserror("The member variable fintn_ is not initialized!");

  return fintn_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetInternalForceNp() const
{
  if (fintnp_.is_null())
    dserror("The member variable fintnp_ is not initialized!");

  return fintnp_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetExternalForceN() const
{
  if (fextn_.is_null())
    dserror("The member variable fextn_ is not initialized!");

  return fextn_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetExternalForceNp() const
{
  if (fextnp_.is_null())
    dserror("The member variable fextnp_ is not initialized!");

  return fextnp_;
}
