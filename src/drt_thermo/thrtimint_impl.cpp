/*----------------------------------------------------------------------*/
/*!
\file thrtimint_impl.cpp
\brief Implicit time integration for spatial discretised
       thermal dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_impl.H"
#include "thr_aux.H"

/*----------------------------------------------------------------------*
 | constructor                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
THR::TimIntImpl::TimIntImpl
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& tdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimInt
  (
    ioparams,
    tdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
  pred_(Teuchos::getIntegralValue<INPAR::THR::PredEnum>(tdynparams,"PREDICT")),
  itertype_(Teuchos::getIntegralValue<INPAR::THR::NonlinSolTech>(tdynparams,"NLNSOL")),
  normtypetempi_(Teuchos::getIntegralValue<INPAR::THR::ConvNorm>(tdynparams,"NORM_TEMP")),
  normtypefres_(Teuchos::getIntegralValue<INPAR::THR::ConvNorm>(tdynparams,"NORM_RESF")),
  combtempifres_(Teuchos::getIntegralValue<INPAR::THR::BinaryOp>(tdynparams,"NORMCOMBI_RESFTEMP")),
  iternorm_(Teuchos::getIntegralValue<INPAR::THR::VectorNorm>(tdynparams,"ITERNORM")),
  itermax_(tdynparams.get<int>("MAXITER")),
  itermin_(tdynparams.get<int>("MINITER")),
  iterdivercont_(Teuchos::getIntegralValue<int>(tdynparams,"DIVERCONT")==1),
  toltempi_(tdynparams.get<double>("TOLTEMP")),
  tolfres_(tdynparams.get<double>("TOLRES")),
  iter_(-1),
  normcharforce_(0.0),
  normchartemp_(0.0),
  normfres_(0.0),
  normtempi_(0.0),
  tempi_(Teuchos::null),
  timer_(actdis->Comm()),
  fres_(Teuchos::null),
  freact_(Teuchos::null)
{

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*dofrowmap_, false);

  // create empty reaction force vector of full length
  freact_ = LINALG::CreateVector(*dofrowmap_, false);

  // iterative temperature increments IncT_{n+1}
  // also known as residual temperatures
  tempi_ = LINALG::CreateVector(*dofrowmap_, true);

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
  Teuchos::RCP<Epetra_Vector> dbcinc
    = LINALG::CreateVector(*dofrowmap_, true);

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
    Teuchos::RCP<Epetra_Vector> freact
      = LINALG::CreateVector(*dofrowmap_, true);
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
  fres_->Scale(-1.0);

  // apply Dirichlet BCs to system of equations
  tempi_->PutScalar(0.0);
  tang_->Complete();
  LINALG::ApplyDirichlettoSystem(tang_, tempi_, fres_,
                                 Teuchos::null, zeros_, *(dbcmaps_->CondMap()));

  // solve for tempi_
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  solver_->Reset();
  solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, true, true);
  solver_->Reset();

  // build residual temperature norm
  normtempi_ = THR::AUX::CalculateVectorNorm(iternorm_, tempi_);

  // set Dirichlet increments in temperature increments
  tempi_->Update(1.0, *dbcinc, 1.0);

  // update end-point temperatures etc
  UpdateIterIncrementally();
  //tempn_->Update(1.0, *tempi_, 1.0);

  // MARK:
  // temperature rates unset on Dirichlet boundary

  // reset to zero
  tempi_->PutScalar(0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_thermo_reset_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // shalom
  return;
}

/*----------------------------------------------------------------------*
 | converged                                                bborn 08/09 |
 *----------------------------------------------------------------------*/
bool THR::TimIntImpl::Converged()
{
  // verify: #normcharforce_ has been delivered strictly larger than zero
  if (normcharforce_ <= 0.0)
  {
    dserror("Characteristic force norm %g must be strictly larger than 0",
            normcharforce_);
  }
  // verify: #normchartemp_ has been delivered strictly larger than zero
  if (normchartemp_ <= 0.0)
  {
    dserror("Characteristic temperature norm %g must be strictly larger than 0",
            normchartemp_);
  }

  // check for single norms
  bool convtemp = false;
  bool convfres = false;

  // residual temperature
  switch (normtypetempi_)
  {
  case INPAR::THR::convnorm_abs:
    convtemp = normtempi_ < toltempi_;
    break;
  case INPAR::THR::convnorm_rel:
    convtemp = normtempi_/normchartemp_ < toltempi_;
    break;
  case INPAR::THR::convnorm_mix:
    convtemp = ( (normtempi_ < toltempi_) or (normtempi_/normchartemp_ < toltempi_) );
    break;
  default:
    dserror("Cannot check for convergence of residual temperatures!");
  }

  // residual forces
  switch (normtypefres_)
  {
  case INPAR::THR::convnorm_abs:
    convfres = normfres_ < tolfres_;
    break;
  case INPAR::THR::convnorm_rel:
    convfres = normfres_/normcharforce_ < tolfres_;
    break;
  case INPAR::THR::convnorm_mix:
    convfres = ( (normfres_ < tolfres_) or (normfres_/normcharforce_ < tolfres_) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
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
void THR::TimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  switch (itertype_)
  {
  case INPAR::THR::soltech_newtonfull :
    NewtonFull();
    break;
  // catch problems
  default :
    dserror("Solution technique \"%s\" is not implemented",
            INPAR::THR::NonlinSolTechString(itertype_).c_str());
    break;
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration              bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::NewtonFull()
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

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normtempi_ was already set in predictor; this is strictly >0
  timer_.ResetStartTime();

  // equilibrium iteration loop
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {

    // make negative residual
    fres_->Scale(-1.0);

    // apply Dirichlet BCs to system of equations
    tempi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(tang_, tempi_, fres_,
                                   Teuchos::null, zeros_, *(dbcmaps_->CondMap()));

    // Solve for tempi_
    // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    solver_->Solve(tang_->EpetraMatrix(), tempi_, fres_, true, iter_==1);
    solver_->ResetTolerance();

    // update end-point temperatures etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and tangent #tang_
    // whose components are globally oriented
    EvaluateRhsTangResidual();

    // extract reaction forces
    // reactions are negative to balance residual on DBC
    freact_->Update(-1.0, *fres_, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

    // build residual force norm
    normfres_ = THR::AUX::CalculateVectorNorm(iternorm_, fres_);
    // build residual temperature norm
    normtempi_ = THR::AUX::CalculateVectorNorm(iternorm_, tempi_);

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (iter_ >= itermax_) and (not iterdivercont_) )
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }
  else if ( (iter_ >= itermax_) and (iterdivercont_) and (myrank_ == 0) )
  {
    printf("Newton unconverged in %d iterations ... continuing\n", iter_);
  }
  else if ( (Converged()) and (myrank_ == 0) )
  {
    PrintNewtonConv();
  }

  // get out of here
  return;
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
  // blank residual at DOFs on Dirichlet BCs
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // apply Dirichlet BCs to system of equations
  tempi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(tang_, tempi_, fres_,
                                 zeros_, *(dbcmaps_->CondMap()));

  // final sip
  return;
}

/*----------------------------------------------------------------------*
 | Update iteration                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::UpdateIter
(
  const int iter  //!< iteration counter
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
}

/*----------------------------------------------------------------------*
 |  Update iteration incrementally with prescribed          bborn 08/09 |
 |  residual temperatures                                               |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::UpdateIterIncrementally
(
  const Teuchos::RCP<const Epetra_Vector> tempi  //!< input residual temperatures
)
{
  // select residual temperatures
  if (tempi != Teuchos::null)
    tempi_->Update(1.0, *tempi, 0.0);  // set the new solution we just got
  else
    tempi_->PutScalar(0.0);

  // Update using #tempi_
  UpdateIterIncrementally();

  // leave this place
  return;
}

/*----------------------------------------------------------------------*
 |  print to screen                                         bborn 08/09 |
 |  originally by lw 12/07                                              |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintPredictor()
{
  // only master processor
  if ( (myrank_ == 0) and printscreen_ )
  {
    // relative check of force residual
    if ( normtypefres_ == INPAR::THR::convnorm_rel )
    {
      std::cout << "Predictor scaled res-norm "
                << normfres_/normcharforce_
                << std::endl;
    }
    // absolute check of force residual
    else if ( normtypefres_ == INPAR::THR::convnorm_abs )
    {
      std::cout << "Predictor absolute res-norm "
                << normfres_
                << std::endl;
    }
    // mixed absolute-relative check of force residual
    else if ( normtypefres_ == INPAR::THR::convnorm_mix )
    {
      std::cout << "Predictor mixed res-norm "
                << min(normfres_, normfres_/normcharforce_)
                << std::endl;
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
}

/*----------------------------------------------------------------------*
 |  print Newton-Raphson iteration to screen and error file bborn 08/09 |
 |  originally by lw 12/07, tk 01/08                                    |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIter()
{
  // print to standard out
  if ( (myrank_ == 0) and printscreen_ and printiter_ )
  {
    if (iter_== 1)
      PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if ( printerrfile_ and printiter_ )
  {
    if (iter_== 1)
      PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*
 |  print header                                            bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIterHeader
(
  FILE* ofile
)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter";

  // different style due relative or absolute error checking
  // temperature
  switch ( normtypefres_ )
  {
  case INPAR::THR::convnorm_rel:
    oss <<std::setw(18)<< "rel-res-norm";
    break;
  case INPAR::THR::convnorm_abs :
    oss <<std::setw(18)<< "abs-res-norm";
    break;
  case INPAR::THR::convnorm_mix :
    oss <<std::setw(18)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypetempi_ )
  {
  case INPAR::THR::convnorm_rel:
    oss <<std::setw(18)<< "rel-temp-norm";
    break;
  case INPAR::THR::convnorm_abs :
    oss <<std::setw(18)<< "abs-temp-norm";
    break;
  case INPAR::THR::convnorm_mix :
    oss <<std::setw(18)<< "mix-temp-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }


  // add solution time
  oss << std::setw(14)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}

/*----------------------------------------------------------------------*
 |  print Newton-Raphson iteration to screen                bborn 08/09 |
 |  originally by lw 12/07, tk 01/08                                    |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonIterText
(
  FILE* ofile
)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7)<< iter_;

  // different style due relative or absolute error checking
  // temperature
  switch ( normtypefres_ )
  {
  case INPAR::THR::convnorm_rel:
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_/normcharforce_;
    break;
  case INPAR::THR::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_;
    break;
  case INPAR::THR::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << min(normfres_, normfres_/normcharforce_);
    break;
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypetempi_ )
  {
  case INPAR::THR::convnorm_rel:
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_/normchartemp_;
    break;
  case INPAR::THR::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normtempi_;
    break;
  case INPAR::THR::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << min(normtempi_, normtempi_/normchartemp_);
    break;
  default:
    dserror("You should not turn up here.");
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  // fprintf(ofile, "%s\n", oss.str().c_str());
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}

/*----------------------------------------------------------------------*
 |  print statistics of converged NRI                       bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintNewtonConv()
{
  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*
 |  print step summary                                      bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintStep()
{

  // print out (only on master CPU)
  if ( (myrank_ == 0) and printscreen_ )
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
}

/*----------------------------------------------------------------------*
 |  print step summary                                      bborn 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::PrintStepText
(
  FILE* ofile
)
{
  // the text
  fprintf(ofile,
          "Finalised: step %6d"
          " | nstep %6d"
          " | time %-14.8E"
          " | dt %-14.8E"
          " | numiter %3d\n",
          step_, stepmax_, (*time_)[0], (*dt_)[0], iter_);
  // print a beautiful line made exactly of 80 dashes
  fprintf(ofile,
          "--------------------------------------------------------------"
          "------------------\n");
  // do it, print now!
  fflush(ofile);

  // fall asleep
  return;
}


#endif  // #ifdef CCADISCRET
