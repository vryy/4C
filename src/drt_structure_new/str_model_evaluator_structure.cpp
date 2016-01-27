/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_structure.cpp

\maintainer Michael Hiermeier

\date Aug 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_structure.H"
#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Structure::Structure()
    : fintnp_ptr_(Teuchos::null),
      fextnp_ptr_(Teuchos::null),
      disnp_ptr_(Teuchos::null),
      stiff_ptr_(Teuchos::null),
      dt_ele_ptr_(NULL)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // setup the internal forces and the external force pointers
  fintnp_ptr_ = gstate_ptr_->GetMutableFintNp();
  fextnp_ptr_ = gstate_ptr_->GetMutableFextNp();
  // setup the displacement pointer
  disnp_ptr_ = gstate_ptr_->GetMutableDisNp();
  // structural element evaluation time
  dt_ele_ptr_ =
      &(gstate_ptr_->GetMutableElementEvaluationTime());

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Reset(const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  stiff_ptr_ = gstate_ptr_->ExtractDisplBlock(jac);

  Reset(x);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  Teuchos::RCP<const Epetra_Vector> disnp =
      gstate_ptr_->ExportDisplEntries(x);
  disnp_ptr_->Scale(1.0,*disnp);

  // reset external forces
  fextnp_ptr_->PutScalar(0.0);

  // reset internal forces
  fintnp_ptr_->PutScalar(0.0);

  // set evaluation time back to zero
  *dt_ele_ptr_ = 0.0;
}

// TODO split this for all other implicit time integrators!!!
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForce(
    const Epetra_Vector& x,
    Epetra_Vector& f)
{
  CheckInitSetup();
  Reset(x);
  bool ok = true;
  // ---------------------------------------
  // (1) EXTERNAL FORCES
  // ---------------------------------------
  ok = ApplyForceExternal();

  // ---------------------------------------
  // (2) INTERNAL FORCES
  // ---------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceInternal() : false);

  // ******************** Finally, put everything together ********************

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  STR::AssembleForce(f,*fextnp_ptr_,-1.0);
  STR::AssembleForce(f,*fintnp_ptr_, 1.0);
  return ok;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceInternal()
{
  CheckInitSetup();
  // action for elements
  Teuchos::ParameterList p;
  p.set<bool>("tolerate_errors",true);
  p.set("action", "calc_struct_internalforce");
  // other parameters that might be needed by the elements
  p.set<double>("total time", gstate_ptr_->GetTimeNp());
  p.set("delta time", gstate_ptr_->GetDeltaTime());
//  p.set<int>("young_temp", young_temp_);

  // set vector values needed by elements
  discret_ptr_->ClearState();
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  // FixMe do we really need these residual displacements?
  Teuchos::RCP<Epetra_Vector> delta_dis_ptr =
      Teuchos::rcp(new Epetra_Vector(disnp_ptr_->Map()));
  delta_dis_ptr->PutScalar(1.0e+50);
  discret_ptr_->SetState(0,"residual displacement", delta_dis_ptr);
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  discret_ptr_->SetState(0,"displacement", disnp_ptr_);

  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_ptr_->Evaluate(p, Teuchos::null, Teuchos::null,
      fintnp_ptr_, Teuchos::null, Teuchos::null);

  discret_ptr_->ClearState();

  return EvalErrorCheck(p);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceExternal()
{
  CheckInitSetup();
  // set target time t_n+1
  Teuchos::ParameterList p;
  p.set<double>("total time", gstate_ptr_->GetTimeNp());
  // Set to default value, because it is unnecessary for the
  // EvaluateNeumann routine.
  p.set("action", "calc_none");

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0, "displacement", gstate_ptr_->GetDisN());
  discret_ptr_->SetState(0, "displacement new", disnp_ptr_);
  discret_ptr_->EvaluateNeumann(p, *fextnp_ptr_);

  return EvalErrorCheck(p);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x,jac);
  bool ok = true;

  /* We use the same routines as for the ApplyForceStiff case, but we
   * do not update the global force vector, which is used for the
   * solution process in the NOX library.
   * This is meaningful, since the computational overhead, which is
   * generated by evaluating the right hand side is negligible */
  // *********** time measurement ***********
  double dtcpu = gstate_ptr_->GetTimer()->WallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = ApplyForceStiffExternal();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceStiffInternal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ +=
      gstate_ptr_->GetTimer()->WallTime() - dtcpu;
  // *********** time measurement ***********

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x,jac);
  bool ok = true;

  // *********** time measurement ***********
  double dtcpu = gstate_ptr_->GetTimer()->WallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = ApplyForceStiffExternal();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceStiffInternal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ +=
      gstate_ptr_->GetTimer()->WallTime() - dtcpu;
  // *********** time measurement ***********

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  STR::AssembleForce(f,*fextnp_ptr_,-1.0);
  STR::AssembleForce(f,*fintnp_ptr_, 1.0);

  // that's it
  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffExternal()
{
  CheckInitSetup();
  const enum INPAR::STR::DampKind& damping_type =
      timint_ptr_->GetDataSDyn().GetDampingType();

  Teuchos::ParameterList p;
  // other parameters needed by the elements
  p.set("total time",gstate_ptr_->GetTimeNp());

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"displacement",gstate_ptr_->GetDisN());

  if (damping_type == INPAR::STR::damp_material)
    discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelN());

  // get load vector
  if (!timint_ptr_->GetDataSDyn().GetLoadLin())
    discret_ptr_->EvaluateNeumann(p, *fextnp_ptr_);
  else
  {
    discret_ptr_->SetState(0,"displacement new", disnp_ptr_);
    /* Add the linearization of the external force to the stiffness
     * matrix. */
    discret_ptr_->EvaluateNeumann(p, fextnp_ptr_, stiff_ptr_);
  }

  return EvalErrorCheck(p);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffInternal()
{
  CheckInitSetup();

  const enum INPAR::STR::DampKind& damping_type =
      timint_ptr_->GetDataSDyn().GetDampingType();

  // action for elements
  Teuchos::ParameterList p;
  p.set("action", "calc_struct_nlnstiff");

  // other parameters that might be needed by the elements
  p.set("total time", gstate_ptr_->GetTimeNp());
  p.set("delta time", gstate_ptr_->GetDeltaTime());
  p.set("damping", damping_type);

  // set vector values needed by elements
  discret_ptr_->ClearState();
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  // FixMe do we really need these residual displacements?
  Teuchos::RCP<Epetra_Vector> delta_dis_ptr =
      Teuchos::rcp(new Epetra_Vector(disnp_ptr_->Map()));
  delta_dis_ptr->PutScalar(1.0e+50);
  discret_ptr_->SetState(0,"residual displacement", delta_dis_ptr);
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  discret_ptr_->SetState(0,"displacement", disnp_ptr_);

  if (damping_type == INPAR::STR::damp_material)
    discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());

  discret_ptr_->Evaluate(
      p,
      stiff_ptr_,
      gstate_ptr_->GetMutableDampMatrix(),
      fintnp_ptr_,
      Teuchos::null,
      Teuchos::null);
  discret_ptr_->ClearState();

  return EvalErrorCheck(p);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepState()
{
  CheckInitSetup();
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  gstate_ptr_->GetMutableMultiDis()->UpdateSteps(*disnp_ptr_);

  // new velocities at t_{n+1} -> t_{n}
  //    V_{n} := V_{n+1}
  gstate_ptr_->GetMutableMultiVel()->UpdateSteps(*gstate_ptr_->GetVelNp());

  // new at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  gstate_ptr_->GetMutableMultiAcc()->UpdateSteps(*gstate_ptr_->GetAccNp());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepElement()
{
  CheckInitSetup();
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", gstate_ptr_->GetTimeNp());
  p.set("delta time", (*gstate_ptr_->GetDeltaTime())[0]);
  //p.set("alpha f", theta_);
  // action for elements
  p.set("action", "calc_struct_update_istep");
  // go to elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState("displacement",gstate_ptr_->GetDisN());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::OutputStepState()
{
  CheckInitSetup();

  // write now
  if (gio_ptr_->IsOutputEveryIter())
  {
    dserror("Not yet implemented!");
    gio_ptr_->GetMutableOutputPtr()->NewStep(gio_ptr_->GetOEI_OutputCounter(), (double) gio_ptr_->GetOEI_OutputCounter());
    gio_ptr_->GetMutableOutputPtr()->WriteVector("displacement",Teuchos::rcp_static_cast<Epetra_MultiVector>(gstate_ptr_->GetMutableDisNp()));
  }
  else
  {
    gio_ptr_->GetMutableOutputPtr()->NewStep(gstate_ptr_->GetStepN(), gstate_ptr_->GetTimeNp());
    gio_ptr_->GetMutableOutputPtr()->WriteVector("displacement", gstate_ptr_->GetDisN());
  }
}
