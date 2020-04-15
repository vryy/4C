/*----------------------------------------------------------------------*/
/*! \file
\brief Implicit time integration for spatial discretised
       thermal dynamics
\level 1
\maintainer Sebastian Proell
*/

/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_impl.H"
#include "thr_aux.H"
#include "thermo_ele_action.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*
 | constructor                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
THR::TimIntImpl::TimIntImpl(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : TimInt(ioparams, tdynparams, xparams, actdis, solver, output),
      pred_(DRT::INPUT::IntegralValue<INPAR::THR::PredEnum>(tdynparams, "PREDICT")),
      itertype_(DRT::INPUT::IntegralValue<INPAR::THR::NonlinSolTech>(tdynparams, "NLNSOL")),
      normtypetempi_(DRT::INPUT::IntegralValue<INPAR::THR::ConvNorm>(tdynparams, "NORM_TEMP")),
      normtypefres_(DRT::INPUT::IntegralValue<INPAR::THR::ConvNorm>(tdynparams, "NORM_RESF")),
      combtempifres_(
          DRT::INPUT::IntegralValue<INPAR::THR::BinaryOp>(tdynparams, "NORMCOMBI_RESFTEMP")),
      iternorm_(DRT::INPUT::IntegralValue<INPAR::THR::VectorNorm>(tdynparams, "ITERNORM")),
      itermax_(tdynparams.get<int>("MAXITER")),
      itermin_(tdynparams.get<int>("MINITER")),
      divcontype_(DRT::INPUT::IntegralValue<INPAR::THR::DivContAct>(tdynparams, "DIVERCONT")),
      divcontrefinelevel_(0),
      divcontfinesteps_(0),
      toltempi_(tdynparams.get<double>("TOLTEMP")),
      tolfres_(tdynparams.get<double>("TOLRES")),
      iter_(-1),
      resetiter_(0),
      normcharforce_(0.0),
      normchartemp_(0.0),
      normfres_(0.0),
      normtempi_(0.0),
      tempi_(Teuchos::null),
      tempinc_(Teuchos::null),
      timer_(actdis->Comm()),
      fres_(Teuchos::null),
      freact_(Teuchos::null),
      fmelt_(Teuchos::null),
      heatint_(DRT::INPUT::IntegralValue<int>(tdynparams, "HEATINTEGRATION") == 1),
      melttol_(tdynparams.get<double>("TOLMELT"))
{
  if (heatint_)
  {
    if (!lumpcapa_) dserror("Heat integration is only possible with a lumped capacity.");
  }
  // create empty residual force vector
  fres_ = LINALG::CreateVector(*discret_->DofRowMap(), false);

  // create empty reaction force vector of full length
  freact_ = LINALG::CreateVector(*discret_->DofRowMap(), false);

  // iterative temperature increments IncT_{n+1}
  // also known as residual temperatures
  tempi_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // incremental temperature increments IncT_{n+1}
  tempinc_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // artificial melting force accountig for latent heat
  fmelt_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  if (heatint_)
  {
    // apply dirichlet BCs first to initialize latent heat history with them
    ApplyDirichletBC((*time_)[0], (*temp_)(0), (*rate_)(0), false);
    // trigger evaluation and storage of element-wise available latent heat on element level
    discret_->ClearState();
    Teuchos::ParameterList params;
    params.set<double>("total time", 0);
    params.set<int>("action", calc_thermo_totallatentheat);
    params.set<double>("delta time", (*dt_)[0]);
    discret_->SetState(0, "temperature", (*temp_)(0));
    // this call is required for correct evaluation of location array
    discret_->Evaluate(params, Teuchos::null, Teuchos::null);
  }


  // setup mortar coupling
  if (DRT::Problem::Instance()->GetProblemType() == prb_thermo)
  {
    DRT::Condition* mrtrcond = actdis->GetCondition("Mortar");
    if (mrtrcond != NULL)
    {
      adaptermeshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());

      std::vector<int> coupleddof(1, 1);
      adaptermeshtying_->Setup(
          actdis, actdis, Teuchos::null, coupleddof, "Mortar", actdis->Comm(), false, false, 0, 0);
      adaptermeshtying_->Evaluate();
    }
  }

  // done so far
  return;
}

/*----------------------------------------------------------------------*
 | integrate step                                           bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::IntegrateStep()
{
  Predict();
  Solve();
  return;
}

/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual   bborn 08/09 |
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::Evaluate(Teuchos::RCP<const Epetra_Vector> tempi)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  //  if (temp != Teuchos::null)
  //  {
  //    // residual temperatures (or iteration increments or iteratively
  //    // incremental temperatures)
  //    Teuchos::RCP<Epetra_Vector> tempi = Teuchos::rcp(new Epetra_Vector(*temp));
  //    tempi->Update(-1.0, *tempinc_, 1.0);
  //
  //    // update incremental temperature member to provided step increments
  //    // shortly: tempinc_^<i> := temp^<i+1>
  //    tempinc_->Update(1.0, *temp, 0.0);
  //
  //    // do thermal update with provided residual temperatures
  //    // recent increment: tempi == tempi_ = \f$\Delta{T}^{<k>}_{n+1}\f$
  //    thermo_->UpdateIterIncrementally(tempi);
  //  }
  //  else
  //  {
  //    thermo_->UpdateIterIncrementally(Teuchos::null);
  //  }

  // TSI does not use NOX --> the Newton increment is passed to the field solver
  UpdateIterIncrementally(tempi);

  // builds tangent, residual and applies DBC
  Evaluate();
}

/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual    dano 02/11 |
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::Evaluate()
{
  // builds tangent, residual and applies DBC
  EvaluateRhsTangResidual();
  PrepareSystemForNewtonSolve();
}

/*----------------------------------------------------------------------*
 | predict solution                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::Predict()
{
  // choose predictor
  if (pred_ == INPAR::THR::pred_consttemp)
  {
    PredictConstTempConsistRate();
    normtempi_ = 1.0e6;
  }
  else if (pred_ == INPAR::THR::pred_consttemprate)
  {
    PredictConstTempRate();
    normtempi_ = 1.0e6;
  }
  else if (pred_ == INPAR::THR::pred_tangtemp)
  {
    PredictTangTempConsistRate();
    // normtempi_ has been set
  }
  else
  {
    dserror("Trouble in determining predictor %i", pred_);
  }

  // apply Dirichlet BCs
  //  ApplyDirichletBC(timen_, temon_, raten_, Teuchos::null, false);
  ApplyDirichletBC(timen_, tempn_, raten_, false);

  // compute residual forces fres_ and tangent tang_
  EvaluateRhsTangResidual();

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // build residual force norm
  normfres_ = THR::AUX::CalculateVectorNorm(iternorm_, fres_);

  // determine characteristic norms
  // we set the minimum of CalcRefNormForce() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchartemp_ = CalcRefNormTemperature();
  if (normchartemp_ == 0.0) normchartemp_ = toltempi_;

  // output
  PrintPredictor();

  // enjoy your meal
  return;
}

/*----------------------------------------------------------------------*
 | prepare partiton step                                     dano 12/10 |
 | like Predict() but without predict the unknown variables T,R         |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PreparePartitionStep()
{
  // set iteration step to 0
  iter_ = 0;

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, tempn_, raten_, false);

  // compute residual forces fres_ and stiffness tang_
  EvaluateRhsTangResidual();

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // split norms
  // build residual force norm
  normfres_ = THR::AUX::CalculateVectorNorm(iternorm_, fres_);

  // determine characteristic norms
  // we set the minumum of CalcRefNormForce() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchartemp_ = CalcRefNormTemperature();
  if (normchartemp_ == 0.0) normchartemp_ = toltempi_;

  // output
  PrintPredictor();

  // enjoy your meal
  return;
}


/*----------------------------------------------------------------------*
 | predict solution as constant temperatures,               bborn 08/09 |
 | temperature rates                                                    |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PredictConstTempRate()
{
  // constant predictor
  tempn_->Update(1.0, *(*temp_)(0), 0.0);
  raten_->Update(1.0, *(*rate_)(0), 0.0);

  // see you next time step
  return;
}

/*----------------------------------------------------------------------*
 | Predict solution as constant temperatures,               bborn 08/09 |
 | temperature rates and tangent                                        |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PredictTangTempConsistRate()
{
  // initialise
  tempn_->Update(1.0, *(*temp_)(0), 0.0);
  raten_->Update(1.0, *(*rate_)(0), 0.0);
  tempi_->PutScalar(0.0);

  // for temperature increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // copy last converged temperatures
  dbcinc->Update(1.0, *(*temp_)(0), 0.0);

  // get Dirichlet values at t_{n+1}
  ApplyDirichletBC(timen_, dbcinc, Teuchos::null, false);

  // subtract the temperatures of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *(*temp_)(0), 1.0);

  // compute residual forces fres_ and tangent tang_
  // at tempn_, etc which are unchanged
  EvaluateRhsTangResidual();

  // add linear reaction forces to residual
  {
    // linear reactions
    Teuchos::RCP<Epetra_Vector> freact = LINALG::CreateVector(*discret_->DofRowMap(), true);
    tang_->Multiply(false, *dbcinc, *freact);

    // add linear reaction forces due to prescribed Dirichlet BCs
    fres_->Update(1.0, *freact, 1.0);
  }

  // extract reaction forces
  freact_->Update(-1.0, *fres_, 0.0);  // reactions are negative
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // make negative residual
  // K . DT = -fres = -(fint - fext)
  fres_->Scale(-1.0);

  // apply Dirichlet BCs to system of equations
  tempi_->PutScalar(0.0);
  tang_->Complete();
  LINALG::ApplyDirichlettoSystem(
      tang_, tempi_, fres_, Teuchos::null, zeros_, *(dbcmaps_->CondMap()));

  // solve for tempi_
  // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
  solver_->Reset();
  solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, true, true);
  solver_->Reset();

  // build residual temperature norm
  normtempi_ = THR::AUX::CalculateVectorNorm(iternorm_, tempi_);

  // set Dirichlet increments in temperature increments
  tempi_->Update(1.0, *dbcinc, 1.0);

  // update end-point temperatures etc
  UpdateIterIncrementally();
  // tempn_->Update(1.0, *tempi_, 1.0);

  // MARK:
  // temperature rates unset on Dirichlet boundary

  // reset to zero
  tempi_->PutScalar(0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    p.set<int>("action", THR::calc_thermo_reset_istep);
    // set the total time
    p.set("total time", (*time_)[0]);
    // go to elements
    discret_->Evaluate(
        p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // shalom
  return;
}

/*----------------------------------------------------------------------*
 | prepare time step                                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  // predict
  Predict();

  // initialise incremental temperatures
  tempinc_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | converged                                                bborn 08/09 |
 *----------------------------------------------------------------------*/
bool THR::TimIntImpl::Converged()
{
  // verify: #normcharforce_ has been delivered strictly larger than zero
  if (normcharforce_ <= 0.0)
  {
    dserror("Characteristic force norm %g must be strictly larger than 0", normcharforce_);
  }
  // verify: #normchartemp_ has been delivered strictly larger than zero
  if (normchartemp_ <= 0.0)
  {
    dserror("Characteristic temperature norm %g must be strictly larger than 0", normchartemp_);
  }

  // check for single norms
  bool convtemp = false;
  bool convfres = false;

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::THR::convnorm_abs:
      convfres = normfres_ < tolfres_;
      break;
    case INPAR::THR::convnorm_rel:
      convfres = normfres_ < std::max(normcharforce_ * tolfres_, 1e-15);
      break;
    case INPAR::THR::convnorm_mix:
      convfres =
          ((normfres_ < tolfres_) or (normfres_ < std::max(normcharforce_ * tolfres_, 1e-15)));
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }

  // residual temperature
  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_abs:
      convtemp = normtempi_ < toltempi_;
      break;
    case INPAR::THR::convnorm_rel:
      convtemp = normtempi_ < std::max(normchartemp_ * toltempi_, 1e-15);
      break;
    case INPAR::THR::convnorm_mix:
      convtemp =
          ((normtempi_ < toltempi_) or (normtempi_ < std::max(normchartemp_ * toltempi_, 1e-15)));
      break;
    default:
      dserror("Cannot check for convergence of residual temperatures!");
      break;
  }

  // combine temperature-like and force-like residuals
  bool conv = false;
  if (combtempifres_ == INPAR::THR::bop_and)
    conv = convtemp and convfres;
  else if (combtempifres_ == INPAR::THR::bop_or)
    conv = convtemp or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;
}

/*----------------------------------------------------------------------*
 | solve equilibrium                                        bborn 08/09 |
 *----------------------------------------------------------------------*/
INPAR::THR::ConvergenceStatus THR::TimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  switch (itertype_)
  {
    case INPAR::THR::soltech_newtonfull:
      return NewtonFull();
    // catch problems
    default:
      dserror("Solution technique \"%s\" is not implemented",
          INPAR::THR::NonlinSolTechString(itertype_).c_str());
      return INPAR::THR::conv_nonlin_fail;  // compiler happiness
  }
}

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration              bborn 08/09 |
 *----------------------------------------------------------------------*/
INPAR::THR::ConvergenceStatus THR::TimIntImpl::NewtonFull()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #tang_ is the effective dynamic tangent matrix

  // check whether we have a sanely filled tangent matrix
  if (not tang_->Filled())
  {
    dserror("Effective tangent matrix must be filled here");
  }

  fmelt_->PutScalar(0);

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normtempi_ was already set in predictor; this is strictly >0
  timer_.ResetStartTime();

  // Do mortar condensation
  if (adaptermeshtying_ != Teuchos::null) adaptermeshtying_->MortarCondensation(tang_, fres_);

  // equilibrium iteration loop
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    // make negative residual
    fres_->Scale(-1.0);

#ifdef THRASOUTPUT
    // finite difference check
    FDCheck();
#endif

    // apply Dirichlet BCs to system of equations
    tempi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(
        tang_, tempi_, fres_, Teuchos::null, zeros_, *(dbcmaps_->CondMap()));

    // Solve for tempi_
    // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, true, iter_ == 1);
    solver_->ResetTolerance();

    // recover condensed variables
    if (adaptermeshtying_ != Teuchos::null) adaptermeshtying_->MortarRecover(tang_, tempi_);

    // update end-point temperatures etc
    UpdateIter(iter_);

    if (heatint_)
    {
#ifdef LATENTHEAT_FIXPOINT_ITER
      EvaluateRhsTangResidual();
      BlankDirichletAndCalcNorms();

      if (Converged())
      {
        std::cout << "inner Newton converged (res=" << normfres_ << ") -> applying latent heat"
                  << std::endl;
        ApplyLatentHeatIntegration();
        // reevlaute residuals and norms
        EvaluateRhsTangResidual();
        BlankDirichletAndCalcNorms();
      }
#else
      ApplyLatentHeatIntegration();
#endif
    }

    // compute residual forces #fres_ and tangent #tang_
    // whose components are globally oriented
    EvaluateRhsTangResidual();

    BlankDirichletAndCalcNorms();

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  return NewtonFullErrorCheck();
}


void THR::TimIntImpl::ApplyLatentHeatIntegration()
{
  Teuchos::ParameterList p;
  const THR::Action action = THR::calc_thermo_phasechangeinc;
  p.set<int>("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  p.set("melt tolerance", melttol_);
  // apply the source term for melting
  discret_->ClearState();
  // SetState(0,...) in case of multiple dofsets (e.g. TSI)
  discret_->SetState(0, "temperature", tempn_);
  // required for linearization of T-dependent capacity
  discret_->SetState(0, "last temperature", (*temp_)(0));
  // TODO create these vectors only once
  Teuchos::RCP<Epetra_Vector> fmeltinc = LINALG::CreateVector(*discret_->DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> tempnmod = LINALG::CreateVector(*discret_->DofRowMap(), true);
  discret_->Evaluate(p, Teuchos::null, Teuchos::null, fmeltinc, tempnmod, Teuchos::null);
  discret_->ClearState();
  fmelt_->Update(1.0, *fmeltinc, 1.0);
  tempn_->Update(1.0, *tempnmod, 0.0);
  // reapply Dirichlet to fix possibly modified nodes
  ApplyDirichletBC(timen_, tempn_, raten_, false);
}


void THR::TimIntImpl::BlankDirichletAndCalcNorms()
{
  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  // copie the dbc onto freact_,
  // everything that is not DBC node ("OtherVector") is blanked
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  // DBC node do not enter the residual, because values are known at the nodes
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // do mortar condensation
  if (adaptermeshtying_ != Teuchos::null) adaptermeshtying_->MortarCondensation(tang_, fres_);

  // build residual force norm
  normfres_ = THR::AUX::CalculateVectorNorm(iternorm_, fres_);
  // build residual temperature norm
  normtempi_ = THR::AUX::CalculateVectorNorm(iternorm_, tempi_);
}



INPAR::THR::ConvergenceStatus THR::TimIntImpl::NewtonFullErrorCheck()
{
  // do some error checks
  if ((iter_ >= itermax_) and (divcontype_ == INPAR::THR::divcont_stop))
  {
    // write restart output of last converged step before stopping
    Output(true);

    dserror("Newton unconverged in %d iterations", iter_);
    return INPAR::THR::conv_nonlin_fail;
  }
  else if ((iter_ >= itermax_) and (divcontype_ == INPAR::THR::divcont_continue))
  {
    if (myrank_ == 0)
      IO::cout << "Newton unconverged in " << iter_ << " iterations, continuing" << IO::endl;
    return INPAR::THR::conv_success;
  }
  else if ((iter_ >= itermax_) and divcontype_ == INPAR::THR::divcont_halve_step)
  {
    HalveTimeStep();
    return INPAR::THR::conv_fail_repeat;
  }
  else if (divcontype_ == INPAR::THR::divcont_repeat_step or
           divcontype_ == INPAR::THR::divcont_repeat_simulation)
  {
    if (myrank_ == 0)
      dserror(
          "Fatal failure in NewtonFullErrorCheck()! divcont_repeat_step and "
          "divcont_repeat_simulation not implemented for THR");
    return INPAR::THR::conv_nonlin_fail;
  }
  // if everything is fine print to screen and return
  if (Converged())
  {
    CheckForTimeStepIncrease();
    return INPAR::THR::conv_success;
  }
  else
    return INPAR::THR::conv_nonlin_fail;

}  // NewtonFull()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::HalveTimeStep()
{
  const double old_dt = Dt();
  const double new_dt = old_dt * 0.5;
  const int endstep = NumStep() + (NumStep() - Step()) + 1;
  SetDt(new_dt);
  SetTimen(TimeOld() + new_dt);
  SetNumStep(endstep);
  ResetStep();
  // TODO limit the maximum number of refinement levels?
  // go down one refinement level
  divcontrefinelevel_++;
  divcontfinesteps_ = 0;

  // remember number of iterations
  resetiter_ += iter_;
  if (Comm().MyPID() == 0)
    IO::cout << "Nonlinear solver failed to converge in step " << Step()
             << ". Divide timestep in half. "
             << "Old time step: " << old_dt << IO::endl
             << "New time step: " << new_dt << IO::endl
             << IO::endl;
}

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                            proell 09/18
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void THR::TimIntImpl::CheckForTimeStepIncrease()
{
  const int maxnumfinestep = 4;

  if (divcontype_ != INPAR::THR::divcont_halve_step)
    return;
  else if (divcontrefinelevel_ != 0)
  {
    // increment for the current, converged step
    divcontfinesteps_++;
    if (divcontfinesteps_ >= maxnumfinestep)
    {
      // increase the step size if the remaining number of steps is a even number
      if (((NumStep() - Step()) % 2) == 0 and NumStep() != Step())
      {
        if (Comm().MyPID() == 0)
          IO::cout << "Nonlinear solver successful. Double timestep size!" << IO::endl;

        // step up one refinement level
        divcontrefinelevel_--;
        divcontfinesteps_ = 0;
        // update total number of steps and next time step
        const int endstep = NumStep() - (NumStep() - Step()) / 2;
        SetNumStep(endstep);
        SetDt(Dt() * 2.0);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | Prepare system for solving with Newton's method          bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrepareSystemForNewtonSolve()
{
  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // make the residual negative
  fres_->Scale(-1.0);
  // blank residual at DOFs on Dirichlet BCs, fres_=0 at nodes with DBC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // apply Dirichlet BCs to system of equations
  tempi_->PutScalar(0.0);  // Useful? depends on solver and more
  // at dofs with DBC change tang_:
  // blank all off-diagonal terms and put 1s at diagonal terms of tang_
  LINALG::ApplyDirichlettoSystem(tang_, tempi_, fres_, zeros_, *(dbcmaps_->CondMap()));

  // final sip
  return;
}  // PrepareSystemForNewtonSolve()


/*----------------------------------------------------------------------*
 | Update iteration                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::UpdateIter(const int iter  //!< iteration counter
)
{
  // we need to do an incremental update (expensive)
  // in the very first iteration (i.e. predictor) of a Newton loop
  // to protect the Dirichlet BCs and to achieve consistent
  // behaviour across all predictors
  // HINT: Sorry, this comment was added delayed and might be inaccurate.
  if (iter <= 1)
  {
    UpdateIterIncrementally();
  }
  else
  {
    UpdateIterIteratively();
  }

  // morning is broken
  return;
}  // UpdateIter()


/*----------------------------------------------------------------------*
 | Update iteration incrementally with prescribed           bborn 08/09 |
 | residual temperatures                                                |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::UpdateIterIncrementally(
    const Teuchos::RCP<const Epetra_Vector> tempi  //!< input residual temperatures
)
{
  // select residual temperatures
  if (tempi != Teuchos::null)
    // tempi_ = \f$\Delta{T}^{<k>}_{n+1}\f$
    tempi_->Update(1.0, *tempi, 0.0);  // set the new solution we just got
  else
    tempi_->PutScalar(0.0);

  // Update using #tempi_
  UpdateIterIncrementally();

  // leave this place
  return;
}  // UpdateIterIncrementally()


/*----------------------------------------------------------------------*
 | update time step                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::Update()
{
  // update temperature and temperature rate
  // after this call we will have tempn_ == temp_ (temp_{n+1} == temp_n), etc.
  UpdateStepState();
  // update everything on the element level
  UpdateStepElement();
  // update time and step
  UpdateStepTime();
  // correct iteration counter by adding all reset iterations
  iter_ += resetiter_;
  resetiter_ = 0;
  return;

}  // Update()


/*----------------------------------------------------------------------*
 | update Newton step                                        dano 02/11 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::UpdateNewton(Teuchos::RCP<const Epetra_Vector> tempi)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  UpdateIterIncrementally(tempi);
  return;

}  // UpdateNewton()


/*----------------------------------------------------------------------*
 | print to screen                                          bborn 08/09 |
 | originally by lw 12/07                                               |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintPredictor()
{
  // only master processor
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    // relative check of force residual
    if (normtypefres_ == INPAR::THR::convnorm_rel)
    {
      std::cout << "Predictor thermo scaled res-norm " << normfres_ / normcharforce_ << std::endl;
    }
    // absolute check of force residual
    else if (normtypefres_ == INPAR::THR::convnorm_abs)
    {
      std::cout << "Predictor thermo absolute res-norm " << normfres_ << std::endl;
    }
    // mixed absolute-relative check of force residual
    else if (normtypefres_ == INPAR::THR::convnorm_mix)
    {
      std::cout << "Predictor thermo mixed res-norm "
                << std::min(normfres_, normfres_ / normcharforce_) << std::endl;
    }
    // default
    else
    {
      dserror("You should not turn up here.");
    }
    // print it, now
    fflush(stdout);
  }

  // leave your hat on
  return;

}  // PrintPredictor()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  bborn 08/09 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIter()
{
  // print to standard out
  if ((myrank_ == 0) and printscreen_ and printiter_ and (StepOld() % printscreen_ == 0))
  {
    if (iter_ == 1) PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  // see you
  return;

}  // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print header                                             bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

  // different style due relative or absolute error checking
  // temperature
  switch (normtypefres_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << "rel-res-norm";
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << "abs-res-norm";
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << "mix-res-norm";
      break;
    default:
      dserror("Unknown type of convergence check for residual forces.");
      break;
  }

  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << "rel-temp-norm";
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << "abs-temp-norm";
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << "mix-temp-norm";
      break;
    default:
      dserror("Unknown type of convergence check for residual temperatures.");
      break;
  }

  // add solution time
  oss << std::setw(14) << "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                 bborn 08/09 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking
  // temperature
  switch (normtypefres_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_ / normcharforce_;
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_;
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normfres_, normfres_ / normcharforce_);
      break;
    default:
      dserror("Unknown type of convergence check for residual forces.");
      break;
  }

  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_ / normchartemp_;
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_;
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normtempi_, normtempi_ / normchartemp_);
      break;
    default:
      dserror("Unknown type of convergence check for residual temperatures.");
      break;
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // PrintNewtonIterText()


/*----------------------------------------------------------------------*
 | print statistics of converged NRI                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonConv()
{
  // somebody did the door
  return;
}  // PrintNewtonConv()


/*----------------------------------------------------------------------*
 | print step summary                                       bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintStep()
{
  // print out (only on master CPU)
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    PrintStepText(stdout);
  }

  // print to error file (on every CPU involved)
  if (printerrfile_)
  {
    PrintStepText(errfile_);
  }

  // fall asleep
  return;
}  // PrintStep()


/*----------------------------------------------------------------------*
 | print step summary                                       bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintStepText(FILE* ofile)
{
  // the text
  fprintf(ofile,
      "Finalised: step %6d"
      " | nstep %6d"
      " | time %-14.8E"
      " | dt %-14.8E"
      " | numiter %3d\n",
      step_, stepmax_, (*time_)[0], (*dt_)[0], iter_ + resetiter_);
  // print a beautiful line made exactly of 80 dashes
  fprintf(ofile,
      "--------------------------------------------------------------"
      "------------------\n");
  // do it, print now!
  fflush(ofile);

  // fall asleep
  return;
}  // PrintStepText()


/*----------------------------------------------------------------------*
 | finite difference check of thermal tangent                dano 09/13 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::FDCheck()
{
  // value of disturbance
  const double delta = 1.0e-8;
  // disturb the current temperature increment

  // ------------------------------------------ initialise matrices and vectors

  // initialise discurbed increment vector
  Teuchos::RCP<Epetra_Vector> disturbtempi = LINALG::CreateVector(*DofRowMap(), true);
  const int dofs = disturbtempi->GlobalLength();
  disturbtempi->PutScalar(0.0);
  disturbtempi->ReplaceGlobalValue(0, 0, delta);

  // initialise rhs
  Teuchos::RCP<Epetra_Vector> rhs_old =
      Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap(), true));
  rhs_old->Update(1.0, *fres_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy =
      Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap(), true));

  // initialise approximation of tangent
  Teuchos::RCP<Epetra_CrsMatrix> tang_approx = LINALG::CreateMatrix((tang_->RowMap()), 81);

  Teuchos::RCP<LINALG::SparseMatrix> tang_copy =
      Teuchos::rcp(new LINALG::SparseMatrix(tang_->EpetraMatrix(), LINALG::Copy));
  std::cout << "\n****************** THR finite difference check ******************" << std::endl;
  std::cout << "thermo field has " << dofs << " DOFs" << std::endl;

  // loop over columns
  // in case of pure thermal problem, start at 0,
  // BUT in case of TSI vector is filled first with STR DOFs followed by THR
  // i.e. insert maximal value of i=STR_DOFs+dofs
  for (int i = 0; i < dofs; ++i)  // TSI: j=STR_DOFs+dofs
  {
    // DOFs that have DBC are not disturbed, i.e. set to zero
    if (dbcmaps_->CondMap()->MyGID(i))
    {
      disturbtempi->ReplaceGlobalValue(i, 0, 0.0);
    }
    // evaluate the element with disturb temperature increment
    Evaluate(disturbtempi);
    rhs_copy->Update(1.0, *fres_, 0.0);
    tempi_->PutScalar(0.0);
    LINALG::ApplyDirichlettoSystem(
        tang_copy, disturbtempi, rhs_copy, Teuchos::null, zeros_, *(dbcmaps_->CondMap()));
    // finite difference approximation of partial derivative
    // rhs_copy = ( rhs_disturb - rhs_old ) . (-1)/delta with rhs_copy==rhs_disturb
    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    int* index = &i;
    // loop over rows
    for (int j = 0; j < dofs; ++j)  // TSI: j=STR_DOFs+dofs
    {
      // insert approximate values using FD into tang_approx
      double value = (*rhs_copy)[j];
      tang_approx->InsertGlobalValues(j, 1, &value, index);
    }  // loop over rows

    // free DOFs (no DBC) get the value (-delta)
    if (not dbcmaps_->CondMap()->MyGID(i))
      disturbtempi->ReplaceGlobalValue(i, 0, -delta);  // row: i, vector index: 0, value: -delta

    // TODO 2013-09-18 was machen diese drei Zeilen??
    disturbtempi->ReplaceGlobalValue(i - 1, 0, 0.0);
    if (i != dofs - 1) disturbtempi->ReplaceGlobalValue(i + 1, 0, delta);
  }  // loop over columns

  // evaluate the element with changed disturbed incremental vector
  Evaluate(disturbtempi);
  tang_approx->FillComplete();
  // copy tang_approx
  Teuchos::RCP<LINALG::SparseMatrix> tang_approx_sparse =
      Teuchos::rcp(new LINALG::SparseMatrix(tang_approx, LINALG::Copy));
  // tang_approx_sparse = tang_approx_sparse - tang_copy
  tang_approx_sparse->Add(*tang_copy, false, -1.0, 1.0);

  // initialise CRSMatrices for the two tangents
  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = tang_copy->EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> error_crs = tang_approx_sparse->EpetraMatrix();
  error_crs->FillComplete();
  sparse_crs->FillComplete();

  // ------------------------------------- initialise values for actual FDCheck
  bool success = true;
  double error_max = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    // only do the check for DOFs which have NO Dirichlet boundary condition
    if (not dbcmaps_->CondMap()->MyGID(i))
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not dbcmaps_->CondMap()->MyGID(j))
        {
          double tang_approx_ij = 0.0;
          double sparse_ij = 0.0;
          double error_ij = 0.0;

          // --------------------------------- get errors of tangent difference
          int errornumentries = 0;
          int errorlength = error_crs->NumMyEntries(i);
          std::vector<double> errorvalues(errorlength);
          std::vector<int> errorindices(errorlength);
          error_crs->ExtractGlobalRowCopy(
              i, errorlength, errornumentries, &errorvalues[0], &errorindices[0]);
          for (int k = 0; k < errorlength; ++k)
          {
            if (errorindices[k] == j)
            {
              error_ij = errorvalues[k];
              break;
            }
          }

          // -------------------------------------- get exact values of tangent
          // get errors of exact tangent
          int sparsenumentries = 0;
          int sparselength = sparse_crs->NumMyEntries(i);
          std::vector<double> sparsevalues(sparselength);
          std::vector<int> sparseindices(sparselength);
          sparse_crs->ExtractGlobalRowCopy(
              i, sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);
          for (int k = 0; k < sparselength; ++k)
          {
            if (sparseindices[k] == j)
            {
              sparse_ij = sparsevalues[k];
              break;
            }
            // else sparse_ij = 0.0;
          }

          // ---------------------------- get approximate values of tang_approx
          int approxnumentries = 0;
          int approxlength = tang_approx->NumMyEntries(i);
          std::vector<double> approxvalues(approxlength);
          std::vector<int> approxindices(approxlength);
          tang_approx->ExtractGlobalRowCopy(
              i, approxlength, approxnumentries, &approxvalues[0], &approxindices[0]);
          for (int k = 0; k < approxlength; ++k)
          {
            if (approxindices[k] == j)
            {
              tang_approx_ij = approxvalues[k];
              break;
            }
            // else tang_approx_ij = 0.0;
          }

          // check value of
          double error = 0.0;
          if (abs(tang_approx_ij) > 1e-7)
            error = error_ij / tang_approx_ij;
          else if (abs(sparse_ij) > 1e-7)
            error = error_ij / sparse_ij;
          // in case current error is higher than maximal, permitted one
          // --> set error_max to current error
          if (abs(error) > abs(error_max)) error_max = abs(error);

          // ---------------------------------------- control values of FDCheck
          if ((abs(error) > 1e-6) and (abs(error_ij) > 1e-7))
          {
            // FDCheck of tangent was NOT successful
            success = false;

            std::cout << "finite difference check failed!\n"
                      << "entry (" << i << "," << j << ") of tang = " << sparse_ij
                      << " and of approx. tang = " << tang_approx_ij
                      << ".\nAbsolute error = " << error_ij << ", relative error = " << error
                      << std::endl;
          }  // control the error values
        }
      }
    }  // FDCheck only for DOFs which have NO DBC
  }    // loop over dofs of successful FD check

  // --------------------------------------------------- FDCheck was successful
  // i.e. tang and its approxiamation are equal w.r.t. given tolerance
  if (success == true)
  {
    std::cout.precision(12);
    std::cout << "finite difference check successful! Maximal relative error = " << error_max
              << std::endl;
    std::cout << "****************** finite difference check done ***************\n\n" << std::endl;
  }
  else
    dserror("FDCheck of thermal tangent failed!");

  return;

}  // FDCheck()


/*----------------------------------------------------------------------*/
