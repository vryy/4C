/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_structure.cpp

\maintainer Michael Hiermeier

\date Aug 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"
#include "str_utils.H"
#include "str_integrator.H"

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
      finertialnp_ptr_(Teuchos::null),
      fvisconp_ptr_(Teuchos::null),
      disnp_ptr_(Teuchos::null),
      stiff_ptr_(Teuchos::null),
      mass_ptr_(Teuchos::null),
      damp_ptr_(Teuchos::null),
      dt_ele_ptr_(NULL),
      masslin_type_(INPAR::STR::ml_none),
      dis_incr_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // get the global state content
  {
    // setup the internal forces and the external force pointers
    fintnp_ptr_ = gstate_ptr_->GetMutableFintNp();
    fextnp_ptr_ = gstate_ptr_->GetMutableFextNp();
    // setup the inertial force vector
    finertialnp_ptr_ = gstate_ptr_->GetMutableFinertialNp();
    // setup the viscous force vector
    fvisconp_ptr_ = gstate_ptr_->GetMutableFviscoNp();
    // setup the displacement pointer
    disnp_ptr_ = gstate_ptr_->GetMutableDisNp();
    // setup the dynamic matrix pointers
    mass_ptr_ = gstate_ptr_->GetMutableMassMatrix();
    damp_ptr_ = gstate_ptr_->GetMutableDampMatrix();
    // structural element evaluation time
    dt_ele_ptr_ =
        &(gstate_ptr_->GetMutableElementEvaluationTime());
  }
  // get the structural dynamic content
  {
    // setup important evaluation booleans
    masslin_type_ = timint_ptr_->GetDataSDyn().GetMassLinType();
  }
  // setup new variables
  {
    dis_incr_ptr_ = Teuchos::rcp(new Epetra_Vector(disnp_ptr_->Map(),true));
  }
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

  // update the state variables of the current time integrator
  int_ptr_->SetState(x);

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
  STR::AssembleVector(1.0,f,-1.0,*fextnp_ptr_);
  STR::AssembleVector(1.0,f,1.0,*fintnp_ptr_);
  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceInternal()
{
  CheckInitSetup();

  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec [3] =
      {Teuchos::null,Teuchos::null,Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] =
      {Teuchos::null,Teuchos::null};

  // set default matrices and vectors
  eval_vec[0] = fintnp_ptr_;

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"residual displacement", dis_incr_ptr_);
  discret_ptr_->SetState(0,"displacement", disnp_ptr_);

  switch (EvalData().GetDampingType())
  {
    case INPAR::STR::damp_material:
    {
      // -----------------------------------------------------------------
      /* evaluate the current (static) state and material damping effects
       * (Note: material damping and non-linear mass effects cannot be
       * considered at the same time at the moment) */
      // -----------------------------------------------------------------
      // action for elements
      EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
      // set the discretization state
      discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
      // reset stiffness matrix
      stiff_ptr_->Zero();
      // reset damping matrix
      damp_ptr_->Zero();
      // set stiffness matrix
      eval_mat[0] = stiff_ptr_;
      // set damping matrix
      eval_mat[1] = damp_ptr_;
      // evaluate ...
      EvaluateInternal(&eval_mat[0],&eval_vec[0]);
      break;
    }
    case INPAR::STR::damp_rayleigh:
    {
      // reset stiffness matrix
      stiff_ptr_->Zero();
      // set stiffness matrix
      eval_mat[0] = stiff_ptr_;

      switch (masslin_type_)
      {
        case INPAR::STR::ml_none:
        {
          // -------------------------------------------------------------
          // evaluate the current (static) state
          // -------------------------------------------------------------
          // action for elements
          EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
          // evaluate ...
          EvaluateInternal(&eval_mat[0],&eval_vec[0]);

          break;
        }
        case INPAR::STR::ml_standard:
        case INPAR::STR::ml_rotations:
        {
          // -------------------------------------------------------------
          /* evaluate the current (static) state and non-linear inertia
           * effects */
          // -------------------------------------------------------------
          // action for elements
          EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);
          // reset the inertial stuff...
          mass_ptr_->Zero();
          finertialnp_ptr_->PutScalar(0.0);
          // set the discretization state
          discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
          discret_ptr_->SetState(0,"acceleration", gstate_ptr_->GetAccNp());
          // set mass matrix
          eval_mat[1] = mass_ptr_;
          // set inertial force
          eval_vec[1] = finertialnp_ptr_;
          // evaluate ...
          EvaluateInternal(&eval_mat[0],&eval_vec[0]);
          // complete the mass matrix
          mass_ptr_->Complete();

          break;
        }
        default:
          dserror("Unknown mass linearization type!");
          break;
      }
      break;
    }
    case INPAR::STR::damp_none:
    {
      switch (masslin_type_)
      {
        case INPAR::STR::ml_none:
        {
          // -------------------------------------------------------------
          // evaluate the current (static) state
          // -------------------------------------------------------------
          // action for elements
          EvalData().SetActionType(DRT::ELEMENTS::struct_calc_internalforce);
          // evaluate ...
          EvaluateInternal(&eval_mat[0],&eval_vec[0]);

          break;
        }
        case INPAR::STR::ml_standard:
        case INPAR::STR::ml_rotations:
        {
          // -------------------------------------------------------------
          /* evaluate the current (static) state and non-linear inertia
           * effects */
          // -------------------------------------------------------------
          // action for elements
          EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);
          // reset the inertial stuff...
          finertialnp_ptr_->PutScalar(0.0);
          // set the discretization state
          discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
          discret_ptr_->SetState(0,"acceleration", gstate_ptr_->GetAccNp());
          // set inertial force
          eval_vec[1] = finertialnp_ptr_;
          // evaluate ...
          EvaluateInternal(&eval_mat[0],&eval_vec[0]);

          break;
        }
        default:
          dserror("Unknown mass linearization type!");
          break;
      }

      break;
    }
    default:
      dserror("Unsupported damping type!");
      break;
  }

  /* We have to do it here, since we assemble directly into the
   * stiffness block during the call of the remaining model
   * evaluators. */
  if (EvalData().GetDampingType()==INPAR::STR::damp_rayleigh)
  {
    stiff_ptr_->Complete();
    const double& dampk =
        timint_ptr_->GetDataSDyn().GetDampingStiffnessFactor();
    const double& dampm =
        timint_ptr_->GetDataSDyn().GetDampingMassFactor();

    damp_ptr_->Add(*stiff_ptr_,false,dampk,0.0);
    damp_ptr_->Add(*mass_ptr_,false,dampm,1.0);
    damp_ptr_->Complete();
  }

  // calculate the inertial force at t_{n+1}
  mass_ptr_->Multiply(false,
      *gstate_ptr_->GetAccNp(),*finertialnp_ptr_);

  // calculate the viscous/damping force at t_{n+1}
  if (EvalData().GetDampingType()!=INPAR::STR::damp_none)
  {
    damp_ptr_->Multiply(false,
        *gstate_ptr_->GetVelNp(),*fvisconp_ptr_);
  }

  return EvalErrorCheck();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceExternal()
{
  CheckInitSetup();

  // Set to default value, because it is unnecessary for the
  // EvaluateNeumann routine.
  EvalData().SetActionType(DRT::ELEMENTS::none);
  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0, "displacement", gstate_ptr_->GetDisN());
  if (EvalData().GetDampingType() == INPAR::STR::damp_material)
    discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelN());
  discret_ptr_->SetState(0, "displacement new", disnp_ptr_);
  EvaluateNeumann(fextnp_ptr_,Teuchos::null);

  return EvalErrorCheck();
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
  STR::AssembleVector(1.0,f,-1.0,*fextnp_ptr_);
  STR::AssembleVector(1.0,f,1.0,*fintnp_ptr_);

  // that's it
  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffExternal()
{
  CheckInitSetup();

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"displacement",gstate_ptr_->GetDisN());

  if (EvalData().GetDampingType() == INPAR::STR::damp_material)
    discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelN());

  // get load vector
  if (!timint_ptr_->GetDataSDyn().GetLoadLin())
    EvaluateNeumann(fextnp_ptr_,Teuchos::null);
  else
  {
    discret_ptr_->SetState(0,"displacement new", disnp_ptr_);
    /* Add the linearization of the external force to the stiffness
     * matrix. */
    EvaluateNeumann(fextnp_ptr_,stiff_ptr_);
  }

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffInternal()
{
  CheckInitSetup();

  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec [3] =
      {Teuchos::null,Teuchos::null,Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] =
      {Teuchos::null,Teuchos::null};

  // set default matrices and vectors
  eval_mat[0] = stiff_ptr_;
  eval_vec[0] = fintnp_ptr_;

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"residual displacement", dis_incr_ptr_);
  discret_ptr_->SetState(0,"displacement", disnp_ptr_);

  // -------------------------------------------------------------
  // evaluate initial dynamic state
  // -------------------------------------------------------------
  if (int_ptr_->IsEquilibriateInitialState())
  {
    // overwrite the standard element action, if lumping is desired
    if (timint_ptr_->GetDataSDyn().IsMassLumping())
      EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstifflmass);
    else
      EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);

    // reset the mass matrix
    mass_ptr_->Zero();
    // set the discretization state
    discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
    discret_ptr_->SetState(0,"acceleration", gstate_ptr_->GetAccNp());
    // set mass matrix
    eval_mat[1] = mass_ptr_;
    // evaluate ...
    EvaluateInternal(&eval_mat[0],&eval_vec[0]);
    mass_ptr_->Complete();
  }
  else
  {
    switch (EvalData().GetDampingType())
    {
      case INPAR::STR::damp_material:
      {
        // -----------------------------------------------------------------
        /* evaluate the current (static) state and material damping effects
         * (Note: material damping and non-linear mass effects cannot be
         * considered at the same time at the moment) */
        // -----------------------------------------------------------------
        // action for elements
        EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
        // reset the damping matrix
        damp_ptr_->Zero();
        // set the discretization state
        discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
        // set damping matrix
        eval_mat[1] = damp_ptr_;
        // evaluate ...
        EvaluateInternal(&eval_mat[0],&eval_vec[0]);
        // complete the damping matrix
        damp_ptr_->Complete();
        break;
      }
      case INPAR::STR::damp_none:
      case INPAR::STR::damp_rayleigh:
      {
        switch (masslin_type_)
        {
          case INPAR::STR::ml_none:
          {
            // -------------------------------------------------------------
            // evaluate the current (static) state
            // -------------------------------------------------------------
            // action for elements
            EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
            // evaluate ...
            EvaluateInternal(&eval_mat[0],&eval_vec[0]);

            break;
          }
          case INPAR::STR::ml_standard:
          case INPAR::STR::ml_rotations:
          {
            // -------------------------------------------------------------
            /* evaluate the current (static) state and non-linear inertia
             * effects */
            // -------------------------------------------------------------
            // action for elements
            EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);
            // reset the inertial stuff...
            mass_ptr_->Zero();
            finertialnp_ptr_->PutScalar(0.0);
            // set the discretization state
            discret_ptr_->SetState(0,"velocity", gstate_ptr_->GetVelNp());
            discret_ptr_->SetState(0,"acceleration", gstate_ptr_->GetAccNp());
            // set mass matrix
            eval_mat[1] = mass_ptr_;
            // set inertial force
            eval_vec[1] = finertialnp_ptr_;
            // evaluate ...
            EvaluateInternal(&eval_mat[0],&eval_vec[0]);
            // complete the mass matrix
            mass_ptr_->Complete();

            break;
          }
          default:
            dserror("Unknown mass linearization type!");
            break;
        }
        break;
      }
      default:
        dserror("Unsupported damping type!");
        break;
    }
  }
  /* We have to do it here, since we assemble directly into the
   * stiffness block during the call of the remaining model
   * evaluators. */
  if (EvalData().GetDampingType()==INPAR::STR::damp_rayleigh)
  {
    stiff_ptr_->Complete();
    const double& dampk =
        timint_ptr_->GetDataSDyn().GetDampingStiffnessFactor();
    const double& dampm =
        timint_ptr_->GetDataSDyn().GetDampingMassFactor();

    damp_ptr_->Add(*stiff_ptr_,false,dampk,0.0);
    damp_ptr_->Add(*mass_ptr_,false,dampm,1.0);
    damp_ptr_->Complete();
  }

  // calculate the inertial force at t_{n+1}
  mass_ptr_->Multiply(false,
      *gstate_ptr_->GetAccNp(),*finertialnp_ptr_);

  // calculate the viscous/damping force at t_{n+1}
  if (EvalData().GetDampingType()!=INPAR::STR::damp_none)
  {
    damp_ptr_->Multiply(false,
        *gstate_ptr_->GetVelNp(),*fvisconp_ptr_);
  }

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
    Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> >("interface",
      eval_data_ptr_);
  EvaluateInternal(p,eval_mat,eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(
    Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
    Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  if (p.numParams()>1)
    dserror("Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  discret_ptr_->Evaluate(p, eval_mat[0], eval_mat[1],
      eval_vec[0], eval_vec[1], eval_vec[2]);
  discret_ptr_->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(
    Teuchos::RCP<Epetra_Vector> eval_vec,
    Teuchos::RCP<LINALG::SparseOperator> eval_mat)
{
  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> >("interface",
      eval_data_ptr_);
  EvaluateNeumann(p,eval_vec,eval_mat);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(
    Teuchos::ParameterList& p,
    Teuchos::RCP<Epetra_Vector> eval_vec,
    Teuchos::RCP<LINALG::SparseOperator> eval_mat)
{
  if (p.numParams()>1)
    dserror("Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  discret_ptr_->EvaluateNeumann(p,eval_vec,eval_mat);
  discret_ptr_->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  CheckInitSetup();
  Reset(xnew);
  /* set the class internal displacement increment vector. Check if it is
   * meaningful in some cases, like incremental strains etc. */
  dis_incr_ptr_ = gstate_ptr_->ExportDisplEntries(dir);
  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"residual displacement",dis_incr_ptr_);
  discret_ptr_->SetState(0,"displacement",disnp_ptr_);
  // set the element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_recover);
  // set the matrix and vector pointers to Teuchos::null
  Teuchos::RCP<Epetra_Vector> eval_vec [3] =
      {Teuchos::null,Teuchos::null,Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] =
      {Teuchos::null,Teuchos::null};

  EvaluateInternal(eval_mat,eval_vec);
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

  // new at t_{n+1} -> t_n
  //    F^{struct}_{n} := F^{struct}_{n+1}
  Teuchos::RCP<Epetra_Vector>& fstructn_ptr =
      gstate_ptr_->GetMutableFstructureN();
  fstructn_ptr->Update(1.0,*fintnp_ptr_,1.0);
  fstructn_ptr->Update(-1.0,*fextnp_ptr_,1.0);

  // set the displacement increment back to zero
  dis_incr_ptr_->Scale(0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepElement()
{
  CheckInitSetup();
  // other parameters that might be needed by the elements
  EvalData().SetTotalTime(gstate_ptr_->GetTimeNp());
  EvalData().SetDeltaTime((*gstate_ptr_->GetDeltaTime())[0]);
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_update_istep);
  // go to elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState("displacement",gstate_ptr_->GetDisN());

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec [3] =
      {Teuchos::null,Teuchos::null,Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] =
      {Teuchos::null,Teuchos::null};
  EvaluateInternal(eval_mat,eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineStressStrain()
{
  CheckInitSetup();
  // set all parameters in the evaluation data container
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_stress);
  EvalData().SetTotalTime(gstate_ptr_->GetTimeNp());
  EvalData().SetDeltaTime((*gstate_ptr_->GetDeltaTime())[0]);
  EvalData().SetStressData(Teuchos::rcp(new std::vector<char>()));
  EvalData().SetStrainData(Teuchos::rcp(new std::vector<char>()));
  EvalData().SetPlasticStrainData(Teuchos::rcp(new std::vector<char>()));

  // set vector values needed by elements
  discret_ptr_->ClearState();
  discret_ptr_->SetState(0,"displacement",gstate_ptr_->GetDisN());
  discret_ptr_->SetState(0,"residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec [3] =
      {Teuchos::null,Teuchos::null,Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] =
      {Teuchos::null,Teuchos::null};

  EvaluateInternal(eval_mat,eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineEnergy()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::OutputStepState()
{
  CheckInitSetup();
  IO::DiscretizationWriter& owriter =
      *(gio_ptr_->GetMutableOutputPtr());

  // write output every iteration for debug purposes
  if (gio_ptr_->IsOutputEveryIter())
  {
    owriter.WriteVector("displacement",
        Teuchos::rcp_static_cast<Epetra_MultiVector>(gstate_ptr_->GetMutableDisNp()));
    /* for visualization of vel and acc do not forget to comment in
     * corresponding lines in StructureEnsightWriter */
    if (gio_ptr_->IsWriteVelAcc())
    {
      owriter.WriteVector("velocity", gstate_ptr_->GetVelNp());
      owriter.WriteVector("acceleration", gstate_ptr_->GetAccNp());
    }
    return;
  }
  else
  {
    // write default output...
    owriter.WriteVector("displacement", gstate_ptr_->GetDisN());

    /* for visualization of vel and acc do not forget to comment in
     * corresponding lines in StructureEnsightWriter */
    if (gio_ptr_->IsWriteVelAcc())
    {
      owriter.WriteVector("velocity", gstate_ptr_->GetVelN());
      owriter.WriteVector("acceleration", gstate_ptr_->GetAccN());
    }
  }
  return;
}
