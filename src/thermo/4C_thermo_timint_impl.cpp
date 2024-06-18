/*----------------------------------------------------------------------*/
/*! \file
\brief Implicit time integration for spatial discretised
       thermal dynamics
\level 1
*/

/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "4C_thermo_timint_impl.hpp"

#include "4C_contact_nitsche_strategy_tsi.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_io_pstream.hpp"
#include "4C_thermo_aux.hpp"
#include "4C_thermo_ele_action.hpp"
#include "4C_thermo_timint.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
THR::TimIntImpl::TimIntImpl(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : TimInt(ioparams, tdynparams, xparams, actdis, solver, output),
      pred_(Core::UTILS::IntegralValue<Inpar::THR::PredEnum>(tdynparams, "PREDICT")),
      itertype_(Core::UTILS::IntegralValue<Inpar::THR::NonlinSolTech>(tdynparams, "NLNSOL")),
      normtypetempi_(Core::UTILS::IntegralValue<Inpar::THR::ConvNorm>(tdynparams, "NORM_TEMP")),
      normtypefres_(Core::UTILS::IntegralValue<Inpar::THR::ConvNorm>(tdynparams, "NORM_RESF")),
      combtempifres_(
          Core::UTILS::IntegralValue<Inpar::THR::BinaryOp>(tdynparams, "NORMCOMBI_RESFTEMP")),
      iternorm_(Core::UTILS::IntegralValue<Inpar::THR::VectorNorm>(tdynparams, "ITERNORM")),
      itermax_(tdynparams.get<int>("MAXITER")),
      itermin_(tdynparams.get<int>("MINITER")),
      divcontype_(Core::UTILS::IntegralValue<Inpar::THR::DivContAct>(tdynparams, "DIVERCONT")),
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
      timer_("", true),
      fres_(Teuchos::null),
      freact_(Teuchos::null)
{
  // create empty residual force vector
  fres_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), false);

  // create empty reaction force vector of full length
  freact_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), false);

  // iterative temperature increments IncT_{n+1}
  // also known as residual temperatures
  tempi_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);

  // incremental temperature increments IncT_{n+1}
  tempinc_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);

  // setup mortar coupling
  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::thermo)
  {
    Core::Conditions::Condition* mrtrcond = actdis->GetCondition("Mortar");
    if (mrtrcond != nullptr)
    {
      adaptermeshtying_ =
          Teuchos::rcp(new Core::Adapter::CouplingMortar(Global::Problem::Instance()->NDim(),
              Global::Problem::Instance()->mortar_coupling_params(),
              Global::Problem::Instance()->contact_dynamic_params(),
              Global::Problem::Instance()->spatial_approximation_type()));

      std::vector<int> coupleddof(1, 1);
      adaptermeshtying_->setup(actdis, actdis, Teuchos::null, coupleddof, "Mortar", actdis->Comm(),
          Global::Problem::Instance()->FunctionManager(), false, false, 0, 0);
      adaptermeshtying_->evaluate();
    }
  }
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
void THR::TimIntImpl::evaluate(Teuchos::RCP<const Epetra_Vector> tempi)
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
  //    thermo_->update_iter_incrementally(tempi);
  //  }
  //  else
  //  {
  //    thermo_->update_iter_incrementally(Teuchos::null);
  //  }

  // TSI does not use NOX --> the Newton increment is passed to the field solver
  update_iter_incrementally(tempi);

  // builds tangent, residual and applies DBC
  evaluate();
}

/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual    dano 02/11 |
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::evaluate()
{
  // builds tangent, residual and applies DBC
  evaluate_rhs_tang_residual();
  prepare_system_for_newton_solve();
}

/*----------------------------------------------------------------------*
 | predict solution                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::Predict()
{
  // choose predictor
  if (pred_ == Inpar::THR::pred_consttemp)
  {
    predict_const_temp_consist_rate();
    normtempi_ = 1.0e6;
  }
  else if (pred_ == Inpar::THR::pred_consttemprate)
  {
    predict_const_temp_rate();
    normtempi_ = 1.0e6;
  }
  else if (pred_ == Inpar::THR::pred_tangtemp)
  {
    predict_tang_temp_consist_rate();
    // normtempi_ has been set
  }
  else
  {
    FOUR_C_THROW("Trouble in determining predictor %i", pred_);
  }

  // integrate the contact to obtain the contact-induced boundary condition
  if (contact_strategy_nitsche_ != Teuchos::null)
  {
    // re-evaluate the contact
    contact_strategy_nitsche_->Integrate(*contact_params_interface_);
  }

  // apply Dirichlet BCs
  //  apply_dirichlet_bc(timen_, temon_, raten_, Teuchos::null, false);
  apply_dirichlet_bc(timen_, tempn_, raten_, false);

  // compute residual forces fres_ and tangent tang_
  evaluate_rhs_tang_residual();

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // build residual force norm
  normfres_ = THR::Aux::calculate_vector_norm(iternorm_, fres_);

  // determine characteristic norms
  // we set the minimum of calc_ref_norm_force() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchartemp_ = calc_ref_norm_temperature();
  if (normchartemp_ == 0.0) normchartemp_ = toltempi_;

  // output
  print_predictor();

  // enjoy your meal
  return;
}

/*----------------------------------------------------------------------*
 | prepare partiton step                                     dano 12/10 |
 | like Predict() but without predict the unknown variables T,R         |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::prepare_partition_step()
{
  // set iteration step to 0
  iter_ = 0;

  // apply Dirichlet BCs
  apply_dirichlet_bc(timen_, tempn_, raten_, false);

  // compute residual forces fres_ and stiffness tang_
  evaluate_rhs_tang_residual();

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // split norms
  // build residual force norm
  normfres_ = THR::Aux::calculate_vector_norm(iternorm_, fres_);

  // determine characteristic norms
  // we set the minumum of calc_ref_norm_force() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchartemp_ = calc_ref_norm_temperature();
  if (normchartemp_ == 0.0) normchartemp_ = toltempi_;

  // output
  print_predictor();

  // enjoy your meal
  return;
}


/*----------------------------------------------------------------------*
 | predict solution as constant temperatures,               bborn 08/09 |
 | temperature rates                                                    |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::predict_const_temp_rate()
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
void THR::TimIntImpl::predict_tang_temp_consist_rate()
{
  // initialise
  tempn_->Update(1.0, *(*temp_)(0), 0.0);
  raten_->Update(1.0, *(*rate_)(0), 0.0);
  tempi_->PutScalar(0.0);

  // for temperature increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);

  // copy last converged temperatures
  dbcinc->Update(1.0, *(*temp_)(0), 0.0);

  // get Dirichlet values at t_{n+1}
  apply_dirichlet_bc(timen_, dbcinc, Teuchos::null, false);

  // subtract the temperatures of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *(*temp_)(0), 1.0);

  // compute residual forces fres_ and tangent tang_
  // at tempn_, etc which are unchanged
  evaluate_rhs_tang_residual();

  // add linear reaction forces to residual
  {
    // linear reactions
    Teuchos::RCP<Epetra_Vector> freact = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
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
  Core::LinAlg::apply_dirichlet_to_system(*tang_, *tempi_, *fres_, *zeros_, *(dbcmaps_->CondMap()));

  // solve for tempi_
  // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
  solver_->Reset();
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, solver_params);
  solver_->Reset();

  // build residual temperature norm
  normtempi_ = THR::Aux::calculate_vector_norm(iternorm_, tempi_);

  // set Dirichlet increments in temperature increments
  tempi_->Update(1.0, *dbcinc, 1.0);

  // update end-point temperatures etc
  update_iter_incrementally();
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
    discret_->evaluate(
        p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // shalom
  return;
}

/*----------------------------------------------------------------------*
 | prepare time step                                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::prepare_time_step()
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
    FOUR_C_THROW("Characteristic force norm %g must be strictly larger than 0", normcharforce_);
  }
  // verify: #normchartemp_ has been delivered strictly larger than zero
  if (normchartemp_ <= 0.0)
  {
    FOUR_C_THROW(
        "Characteristic temperature norm %g must be strictly larger than 0", normchartemp_);
  }

  // check for single norms
  bool convtemp = false;
  bool convfres = false;

  // residual forces
  switch (normtypefres_)
  {
    case Inpar::THR::convnorm_abs:
      convfres = normfres_ < tolfres_;
      break;
    case Inpar::THR::convnorm_rel:
      convfres = normfres_ < std::max(normcharforce_ * tolfres_, 1e-15);
      break;
    case Inpar::THR::convnorm_mix:
      convfres =
          ((normfres_ < tolfres_) or (normfres_ < std::max(normcharforce_ * tolfres_, 1e-15)));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }

  // residual temperature
  switch (normtypetempi_)
  {
    case Inpar::THR::convnorm_abs:
      convtemp = normtempi_ < toltempi_;
      break;
    case Inpar::THR::convnorm_rel:
      convtemp = normtempi_ < std::max(normchartemp_ * toltempi_, 1e-15);
      break;
    case Inpar::THR::convnorm_mix:
      convtemp =
          ((normtempi_ < toltempi_) or (normtempi_ < std::max(normchartemp_ * toltempi_, 1e-15)));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual temperatures!");
      break;
  }

  // combine temperature-like and force-like residuals
  bool conv = false;
  if (combtempifres_ == Inpar::THR::bop_and)
    conv = convtemp and convfres;
  else if (combtempifres_ == Inpar::THR::bop_or)
    conv = convtemp or convfres;
  else
    FOUR_C_THROW("Something went terribly wrong with binary operator!");

  // return things
  return conv;
}

/*----------------------------------------------------------------------*
 | solve equilibrium                                        bborn 08/09 |
 *----------------------------------------------------------------------*/
Inpar::THR::ConvergenceStatus THR::TimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  switch (itertype_)
  {
    case Inpar::THR::soltech_newtonfull:
      return NewtonFull();
    // catch problems
    default:
      FOUR_C_THROW("Solution technique \"%s\" is not implemented",
          Inpar::THR::NonlinSolTechString(itertype_).c_str());
      return Inpar::THR::conv_nonlin_fail;  // compiler happiness
  }
}

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration              bborn 08/09 |
 *----------------------------------------------------------------------*/
Inpar::THR::ConvergenceStatus THR::TimIntImpl::NewtonFull()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #tang_ is the effective dynamic tangent matrix

  // check whether we have a sanely filled tangent matrix
  if (not tang_->Filled())
  {
    FOUR_C_THROW("Effective tangent matrix must be filled here");
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normtempi_ was already set in predictor; this is strictly >0
  timer_.reset();

  // Do mortar condensation
  if (adaptermeshtying_ != Teuchos::null) adaptermeshtying_->MortarCondensation(tang_, fres_);

  // equilibrium iteration loop
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    // make negative residual
    fres_->Scale(-1.0);

#ifdef THRASOUTPUT
    // finite difference check
    fd_check();
#endif

    // apply Dirichlet BCs to system of equations
    tempi_->PutScalar(0.0);  // Useful? depends on solver and more
    Core::LinAlg::apply_dirichlet_to_system(
        *tang_, *tempi_, *fres_, *zeros_, *(dbcmaps_->CondMap()));

    // Solve for tempi_
    // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
    Core::LinAlg::SolverParams solver_params;
    if (solveradapttol_ and (iter_ > 1))
    {
      solver_params.nonlin_tolerance = tolfres_;
      solver_params.nonlin_residual = normfres_;
      solver_params.lin_tol_better = solveradaptolbetter_;
    }
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, solver_params);
    solver_->ResetTolerance();

    // recover condensed variables
    if (adaptermeshtying_ != Teuchos::null) adaptermeshtying_->MortarRecover(tang_, tempi_);

    // update end-point temperatures etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and tangent #tang_
    // whose components are globally oriented
    evaluate_rhs_tang_residual();

    blank_dirichlet_and_calc_norms();

    // print stuff
    print_newton_iter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  return newton_full_error_check();
}


void THR::TimIntImpl::blank_dirichlet_and_calc_norms()
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
  normfres_ = THR::Aux::calculate_vector_norm(iternorm_, fres_);
  // build residual temperature norm
  normtempi_ = THR::Aux::calculate_vector_norm(iternorm_, tempi_);
}



Inpar::THR::ConvergenceStatus THR::TimIntImpl::newton_full_error_check()
{
  // do some error checks
  if ((iter_ >= itermax_) and (divcontype_ == Inpar::THR::divcont_stop))
  {
    // write restart output of last converged step before stopping
    output(true);

    FOUR_C_THROW("Newton unconverged in %d iterations", iter_);
    return Inpar::THR::conv_nonlin_fail;
  }
  else if ((iter_ >= itermax_) and (divcontype_ == Inpar::THR::divcont_continue))
  {
    if (myrank_ == 0)
      Core::IO::cout << "Newton unconverged in " << iter_ << " iterations, continuing"
                     << Core::IO::endl;
    return Inpar::THR::conv_success;
  }
  else if ((iter_ >= itermax_) and divcontype_ == Inpar::THR::divcont_halve_step)
  {
    halve_time_step();
    return Inpar::THR::conv_fail_repeat;
  }
  else if (divcontype_ == Inpar::THR::divcont_repeat_step or
           divcontype_ == Inpar::THR::divcont_repeat_simulation)
  {
    if (myrank_ == 0)
      FOUR_C_THROW(
          "Fatal failure in newton_full_error_check()! divcont_repeat_step and "
          "divcont_repeat_simulation not implemented for THR");
    return Inpar::THR::conv_nonlin_fail;
  }
  // if everything is fine print to screen and return
  if (Converged())
  {
    check_for_time_step_increase();
    return Inpar::THR::conv_success;
  }
  else
    return Inpar::THR::conv_nonlin_fail;

}  // NewtonFull()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::halve_time_step()
{
  const double old_dt = Dt();
  const double new_dt = old_dt * 0.5;
  const int endstep = NumStep() + (NumStep() - Step()) + 1;
  set_dt(new_dt);
  SetTimen(TimeOld() + new_dt);
  SetNumStep(endstep);
  reset_step();
  // TODO limit the maximum number of refinement levels?
  // go down one refinement level
  divcontrefinelevel_++;
  divcontfinesteps_ = 0;

  // remember number of iterations
  resetiter_ += iter_;
  if (Comm().MyPID() == 0)
    Core::IO::cout << "Nonlinear solver failed to converge in step " << Step()
                   << ". Divide timestep in half. "
                   << "Old time step: " << old_dt << Core::IO::endl
                   << "New time step: " << new_dt << Core::IO::endl
                   << Core::IO::endl;
}

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                            proell 09/18
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void THR::TimIntImpl::check_for_time_step_increase()
{
  const int maxnumfinestep = 4;

  if (divcontype_ != Inpar::THR::divcont_halve_step)
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
          Core::IO::cout << "Nonlinear solver successful. Double timestep size!" << Core::IO::endl;

        // step up one refinement level
        divcontrefinelevel_--;
        divcontfinesteps_ = 0;
        // update total number of steps and next time step
        const int endstep = NumStep() - (NumStep() - Step()) / 2;
        SetNumStep(endstep);
        set_dt(Dt() * 2.0);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | Prepare system for solving with Newton's method          bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::prepare_system_for_newton_solve()
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
  Core::LinAlg::apply_dirichlet_to_system(*tang_, *tempi_, *fres_, *zeros_, *(dbcmaps_->CondMap()));

  // final sip
  return;
}  // prepare_system_for_newton_solve()


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
    update_iter_incrementally();
  }
  else
  {
    update_iter_iteratively();
  }

  // morning is broken
  return;
}  // UpdateIter()


/*----------------------------------------------------------------------*
 | Update iteration incrementally with prescribed           bborn 08/09 |
 | residual temperatures                                                |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::update_iter_incrementally(
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
  update_iter_incrementally();

  // leave this place
  return;
}  // update_iter_incrementally()


/*----------------------------------------------------------------------*
 | update time step                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::update()
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

}  // update()


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
  update_iter_incrementally(tempi);
  return;

}  // UpdateNewton()


/*----------------------------------------------------------------------*
 | print to screen                                          bborn 08/09 |
 | originally by lw 12/07                                               |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_predictor()
{
  // only master processor
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    // relative check of force residual
    if (normtypefres_ == Inpar::THR::convnorm_rel)
    {
      std::cout << "Predictor thermo scaled res-norm " << normfres_ / normcharforce_ << std::endl;
    }
    // absolute check of force residual
    else if (normtypefres_ == Inpar::THR::convnorm_abs)
    {
      std::cout << "Predictor thermo absolute res-norm " << normfres_ << std::endl;
    }
    // mixed absolute-relative check of force residual
    else if (normtypefres_ == Inpar::THR::convnorm_mix)
    {
      std::cout << "Predictor thermo mixed res-norm "
                << std::min(normfres_, normfres_ / normcharforce_) << std::endl;
    }
    // default
    else
    {
      FOUR_C_THROW("You should not turn up here.");
    }
    // print it, now
    fflush(stdout);
  }

  // leave your hat on
  return;

}  // print_predictor()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  bborn 08/09 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_newton_iter()
{
  // print to standard out
  if ((myrank_ == 0) and printscreen_ and printiter_ and (StepOld() % printscreen_ == 0))
  {
    if (iter_ == 1) print_newton_iter_header(stdout);
    print_newton_iter_text(stdout);
  }

}  // print_newton_iter()


/*----------------------------------------------------------------------*
 | print header                                             bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_newton_iter_header(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

  // different style due relative or absolute error checking
  // temperature
  switch (normtypefres_)
  {
    case Inpar::THR::convnorm_rel:
      oss << std::setw(18) << "rel-res-norm";
      break;
    case Inpar::THR::convnorm_abs:
      oss << std::setw(18) << "abs-res-norm";
      break;
    case Inpar::THR::convnorm_mix:
      oss << std::setw(18) << "mix-res-norm";
      break;
    default:
      FOUR_C_THROW("Unknown type of convergence check for residual forces.");
      break;
  }

  switch (normtypetempi_)
  {
    case Inpar::THR::convnorm_rel:
      oss << std::setw(18) << "rel-temp-norm";
      break;
    case Inpar::THR::convnorm_abs:
      oss << std::setw(18) << "abs-temp-norm";
      break;
    case Inpar::THR::convnorm_mix:
      oss << std::setw(18) << "mix-temp-norm";
      break;
    default:
      FOUR_C_THROW("Unknown type of convergence check for residual temperatures.");
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
}  // print_newton_iter_header()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                 bborn 08/09 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_newton_iter_text(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking
  // temperature
  switch (normtypefres_)
  {
    case Inpar::THR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_ / normcharforce_;
      break;
    case Inpar::THR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_;
      break;
    case Inpar::THR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normfres_, normfres_ / normcharforce_);
      break;
    default:
      FOUR_C_THROW("Unknown type of convergence check for residual forces.");
      break;
  }

  switch (normtypetempi_)
  {
    case Inpar::THR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_ / normchartemp_;
      break;
    case Inpar::THR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_;
      break;
    case Inpar::THR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normtempi_, normtempi_ / normchartemp_);
      break;
    default:
      FOUR_C_THROW("Unknown type of convergence check for residual temperatures.");
      break;
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.totalElapsedTime(true);

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // print_newton_iter_text()


/*----------------------------------------------------------------------*
 | print statistics of converged NRI                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_newton_conv()
{
  // somebody did the door
  return;
}  // print_newton_conv()


/*----------------------------------------------------------------------*
 | print step summary                                       bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintStep()
{
  // print out (only on master CPU)
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    print_step_text(stdout);
  }
}  // print_step()


/*----------------------------------------------------------------------*
 | print step summary                                       bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::print_step_text(FILE* ofile)
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
}  // print_step_text()


/*----------------------------------------------------------------------*
 | finite difference check of thermal tangent                dano 09/13 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::fd_check()
{
  // value of disturbance
  const double delta = 1.0e-8;
  // disturb the current temperature increment

  // ------------------------------------------ initialise matrices and vectors

  // initialise discurbed increment vector
  Teuchos::RCP<Epetra_Vector> disturbtempi = Core::LinAlg::CreateVector(*dof_row_map(), true);
  const int dofs = disturbtempi->GlobalLength();
  disturbtempi->PutScalar(0.0);
  disturbtempi->ReplaceGlobalValue(0, 0, delta);

  // initialise rhs
  Teuchos::RCP<Epetra_Vector> rhs_old =
      Teuchos::rcp(new Epetra_Vector(*discret_->dof_row_map(), true));
  rhs_old->Update(1.0, *fres_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy =
      Teuchos::rcp(new Epetra_Vector(*discret_->dof_row_map(), true));

  // initialise approximation of tangent
  Teuchos::RCP<Epetra_CrsMatrix> tang_approx = Core::LinAlg::CreateMatrix((tang_->RowMap()), 81);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tang_copy =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(tang_->EpetraMatrix(), Core::LinAlg::Copy));
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
    evaluate(disturbtempi);
    rhs_copy->Update(1.0, *fres_, 0.0);
    tempi_->PutScalar(0.0);
    Core::LinAlg::apply_dirichlet_to_system(
        *tang_copy, *disturbtempi, *rhs_copy, *zeros_, *(dbcmaps_->CondMap()));
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
  evaluate(disturbtempi);
  tang_approx->FillComplete();
  // copy tang_approx
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tang_approx_sparse =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(tang_approx, Core::LinAlg::Copy));
  // tang_approx_sparse = tang_approx_sparse - tang_copy
  tang_approx_sparse->Add(*tang_copy, false, -1.0, 1.0);

  // initialise CRSMatrices for the two tangents
  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = tang_copy->EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> error_crs = tang_approx_sparse->EpetraMatrix();
  error_crs->FillComplete();
  sparse_crs->FillComplete();

  // ------------------------------------- initialise values for actual fd_check
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
              i, errorlength, errornumentries, errorvalues.data(), errorindices.data());
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
              i, sparselength, sparsenumentries, sparsevalues.data(), sparseindices.data());
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
              i, approxlength, approxnumentries, approxvalues.data(), approxindices.data());
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

          // ---------------------------------------- control values of fd_check
          if ((abs(error) > 1e-6) and (abs(error_ij) > 1e-7))
          {
            // fd_check of tangent was NOT successful
            success = false;

            std::cout << "finite difference check failed!\n"
                      << "entry (" << i << "," << j << ") of tang = " << sparse_ij
                      << " and of approx. tang = " << tang_approx_ij
                      << ".\nAbsolute error = " << error_ij << ", relative error = " << error
                      << std::endl;
          }  // control the error values
        }
      }
    }  // fd_check only for DOFs which have NO DBC
  }    // loop over dofs of successful FD check

  // --------------------------------------------------- fd_check was successful
  // i.e. tang and its approxiamation are equal w.r.t. given tolerance
  if (success == true)
  {
    std::cout.precision(12);
    std::cout << "finite difference check successful! Maximal relative error = " << error_max
              << std::endl;
    std::cout << "****************** finite difference check done ***************\n\n" << std::endl;
  }
  else
    FOUR_C_THROW("fd_check of thermal tangent failed!");

  return;

}  // fd_check()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
