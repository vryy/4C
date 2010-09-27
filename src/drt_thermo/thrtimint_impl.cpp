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

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"

/*----------------------------------------------------------------------*
 | constructor                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
THR::TimIntImpl::TimIntImpl(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& tdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
  )
: TimInt(
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
    // set the total time
    p.set("total time",(*time_)[0]);
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
    //    04.03.10 like Georg did in scatra_timint_implicit.cpp
    /*
          // matrix printing options (DEBUGGING!)
          RCP<LINALG::SparseMatrix> A = SystemMatrix();
          if (A != Teuchos::null)
          {
            // print to file in matlab format
            const std::string fname = "sparsematrix.mtl";
            LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
            // print to screen
            (A->EpetraMatrix())->Print(cout);
            // print sparsity pattern to file
            LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
          }
          else
          {
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = BlockSystemMatrix();
            const std::string fname = "sparsematrix.mtl";
            LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
          }
          */
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
void THR::TimIntImpl::UpdateIter(
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
void THR::TimIntImpl::UpdateIterIncrementally(
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
void THR::TimIntImpl::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

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
void THR::TimIntImpl::PrintNewtonIterText(FILE* ofile)
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
void THR::TimIntImpl::PrintStepText(FILE* ofile)
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


/*----------------------------------------------------------------------*
 | recalculate thermal matrices for coupled problems         dano 05/10 |
 | (like STR::TimIntImpl::UseBlockMatrix)                               |
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::TSIMatrix()
{
  // recalculate thermal matrices

  Teuchos::RCP<Epetra_Vector> fint
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  // initialise matrix (in thermo only one global matrix required)
  tang_->Zero();
//  capa_->Zero();

  {
    // create the parameters for the discretization
    ParameterList p;

    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState(0,"residual temperature", zeros_);
    discret_->SetState(0,"temperature", (*temp_)(0));
    if (disn_!=Teuchos::null)
    {
      discret_->SetState(1,"displacement",disn_);
    }
    if (veln_!=Teuchos::null)
    {
      discret_->SetState(1,"velocity",veln_);
    }
    // set action for elements depending on time integrator
//    p.set("action", "calc_thermo_fintcond");
    // type of calling time integrator
    p.set("time integrator", MethodName());
    if (MethodName()== INPAR::THR::dyna_onesteptheta)
    {
      p.set("action", "calc_thermo_fintcapa");
      // calculate the dynamic tangent (save on place of capacity matrix)
      discret_->Evaluate(p, Teuchos::null, tang_, fint, Teuchos::null, Teuchos::null);
    }
    else if (MethodName()== INPAR::THR::dyna_statics)
    {
      p.set("action", "calc_thermo_fintcond");
      // calculate the tangent
      discret_->Evaluate(p, tang_, Teuchos::null, fint, Teuchos::null,
                         Teuchos::null);
    }
    discret_->ClearState();
  }

//  // finish mass matrix
//  mass_->Complete();

  // close stiffness matrix
  tang_->Complete();

  // We need to reset the stiffness matrix because its graph (topology)
  // is not finished yet in case of constraints and posssibly other side
  // effects (basically managers).
  tang_->Reset();
} // TSIMatrix

/*----------------------------------------------------------------------*
 |  Modify thermal system of equation towards thermal contact mgit 09/10|
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::ApplyThermoContact(Teuchos::RCP<LINALG::SparseMatrix>& tang,
                                         Teuchos::RCP<Epetra_Vector>& feff,
                                         Teuchos::RCP<Epetra_Vector>& temp)
{
  // only in the case of contact
  if(cmtman_==Teuchos::null)
    return;

   //contact / meshtying modifications need -fres
   feff->Scale(-1.0);

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  tang->Complete();

  // convert maps (from structure discretization to thermo discretization)
  // slave-, active-, inactive-, master-, activemaster-, n- smdofs
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  ConvertMaps (sdofs,adofs,mdofs);

  // map of active and master dofs
  amdofs = LINALG::MergeMap(adofs,mdofs,false);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);

  // row map of thermal problem
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(discret_->DofRowMap())));

  // split problemrowmap in n+am
  ndofs = LINALG::SplitMap(*problemrowmap,*smdofs);

  // modifications only for active nodes
  if (adofs->NumGlobalElements()==0)
  {
    feff->Scale(-1.0);
    return;
  }

  // assemble Mortar Matrices D and M in thermo dofs for active nodes
  RCP<LINALG::SparseMatrix> dmatrix = rcp(new LINALG::SparseMatrix(*sdofs,10));
  RCP<LINALG::SparseMatrix> mmatrix = rcp(new LINALG::SparseMatrix(*sdofs,100));

  AssembleDM(*dmatrix,*mmatrix);

  // FillComplete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*mdofs,*sdofs);

  // active part of dmatrix and mmatrix
  RCP<Epetra_Map> tmp;
  RCP<LINALG::SparseMatrix> dmatrixa,mmatrixa,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
  LINALG::SplitMatrix2x2(dmatrix,adofs,idofs,adofs,idofs,dmatrixa,tmp1,tmp2,tmp3);
  LINALG::SplitMatrix2x2(mmatrix,adofs,idofs,mdofs,tmp,mmatrixa,tmp4,tmp5,tmp6);

  // assemble mechanical dissipation
  RCP<Epetra_Vector> mechdissrate = LINALG::CreateVector(*mdofs,true);
  AssembleMechDissRate(*mechdissrate);

  // matrices from linearized thermal contact condition
  RCP<LINALG::SparseMatrix>  thermcontLM = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<LINALG::SparseMatrix>  thermcontTEMP = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<Epetra_Vector>         thermcontRHS = LINALG::CreateVector(*adofs,true);

  // assemble thermal contact contition
  AssembleThermContCondition(*thermcontLM,*thermcontTEMP,*thermcontRHS,*dmatrixa,*mmatrixa,adofs,mdofs);

  // complete the matrices
  thermcontLM->Complete(*adofs,*adofs);
  thermcontTEMP->Complete(*smdofs,*adofs);

  /**********************************************************************/
  /* Modification of the stiff matrix and rhs towards thermo contact    */
  /**********************************************************************/

  /**********************************************************************/
  /* Create inv(D)                                                      */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*sdofs,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  RCP<LINALG::SparseMatrix> mhatmatrix;
  mhatmatrix = LINALG::MLMultiply(*invd,false,*mmatrix,false,false,false,true);

  /**********************************************************************/
  /* Split tang into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  RCP<LINALG::SparseMatrix> tangmatrix = rcp(new LINALG::SparseMatrix(*tang));
  LINALG::SplitMatrix2x2(tangmatrix,smdofs,ndofs,smdofs,ndofs,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,sdofs,mdofs,sdofs,mdofs,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,sdofs,mdofs,ndofs,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,ndofs,tempmap,sdofs,mdofs,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*problemrowmap,*feff,smdofs,fsm,ndofs,fn);
  LINALG::SplitVector(*smdofs,*fsm,sdofs,fs,mdofs,fm);

  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,adofs,idofs,adofs,idofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,adofs,idofs,ndofs,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,adofs,idofs,mdofs,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,mdofs,tempmap,adofs,idofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*adofs));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*idofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*sdofs,*fs,adofs,fa,idofs,fi);

  // abbreviations for active and inactive set
  int aset = adofs->NumGlobalElements();
  int iset = idofs->NumGlobalElements();

  // active part of invd and mhatmatrix
  RCP<Epetra_Map> tmpmap;
  RCP<LINALG::SparseMatrix> invda,mhata;
  LINALG::SplitMatrix2x2(invd,adofs,idofs,adofs,idofs,invda,tmp1,tmp2,tmp3);
  LINALG::SplitMatrix2x2(mhatmatrix,adofs,idofs,mdofs,tmpmap,mhata,tmp1,tmp2,tmp3);

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmnmod->Add(*kmn,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhata,true,*kan,false,false,false,true);
  kmnmod->Add(*kmnadd,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmmmod->Add(*kmm,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhata,true,*kam,false,false,false,true);
  kmmmod->Add(*kmmadd,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmimod->Add(*kmi,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhata,true,*kai,false,false,false,true);
    kmimod->Add(*kmiadd,false,1.0,1.0);
    kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  }

  // kmi: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmamod->Add(*kma,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhata,true,*kaa,false,false,false,true);
    kmamod->Add(*kmaadd,false,1.0,1.0);
    kmamod->Complete(kma->DomainMap(),kma->RowMap());
  }

  // kan: thermcontlm*invd*kan
  RCP<LINALG::SparseMatrix> kanmod;
  if (aset)
  {
    kanmod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
    kanmod->Complete(kan->DomainMap(),kan->RowMap());
  }

  // kam: thermcontlm*invd*kam
  RCP<LINALG::SparseMatrix> kammod;
  if (aset)
  {
    kammod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
    kammod->Complete(kam->DomainMap(),kam->RowMap());
  }

  // kai: thermcontlm*invd*kai
  RCP<LINALG::SparseMatrix> kaimod;
  if (aset && iset)
  {
    kaimod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
    kaimod->Complete(kai->DomainMap(),kai->RowMap());
  }

  // kaa: thermcontlm*invd*kaa
  RCP<LINALG::SparseMatrix> kaamod;
  if (aset)
  {
    kaamod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
    kaamod->Complete(kaa->DomainMap(),kaa->RowMap());
  }

  // Modifications towards rhs
  // FIXGIT: pay attention to genalpha
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*mdofs));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  // fa: mutliply with thermcontlm
  RCP<Epetra_Vector> famod;
  {
    famod = rcp(new Epetra_Vector(*adofs));
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    temp->Multiply(false,*fa,*famod);
  }

  /**********************************************************************/
  /* Global setup of tangnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> tangnew = rcp(new LINALG::SparseMatrix(*problemrowmap,81,true,false,tangmatrix->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap);

  // add n submatrices to tangnew
  tangnew->Add(*knn,false,1.0,1.0);
  tangnew->Add(*knm,false,1.0,1.0);
  tangnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to tangnew
  tangnew->Add(*kmnmod,false,1.0,1.0);
  tangnew->Add(*kmmmod,false,1.0,1.0);
  if (iset) tangnew->Add(*kmimod,false,1.0,1.0);
  if (aset) tangnew->Add(*kmamod,false,1.0,1.0);

  // add i submatrices to tangnew
  if (iset) tangnew->Add(*kin,false,1.0,1.0);
  if (iset) tangnew->Add(*kim,false,1.0,1.0);
  if (iset) tangnew->Add(*kii,false,1.0,1.0);
  if (iset) tangnew->Add(*kia,false,1.0,1.0);

  // add a submatrices to tangnew
  if (aset) tangnew->Add(*kanmod,false,1.0,1.0);
  if (aset) tangnew->Add(*kammod,false,1.0,1.0);
  if (aset && iset) tangnew->Add(*kaimod,false,1.0,1.0);
  if (aset) tangnew->Add(*kaamod,false,1.0,1.0);

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add mechanical dissipation to feffnew
  RCP<Epetra_Vector> mechdissrateexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*mechdissrate,*mechdissrateexp);
  feffnew->Update(1.0,*mechdissrateexp,1.0);

  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp;
  if (iset)
  {
    fiexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*fi,*fiexp);
    feffnew->Update(1.0,*fiexp,1.0);
  }

  // add a subvector to feffnew
  RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*famod,*famodexp);
    feffnew->Update(1.0,*famodexp,+1.0);
  }

  // add linearized thermo contact condition
  tangnew->Add(*thermcontTEMP,false,-1.0,+1.0);

  // add rhs of thermal contact condition to feffnew
  RCP<Epetra_Vector> thermcontRHSexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*thermcontRHS,*thermcontRHSexp);
  feffnew->Update(-1.0,*thermcontRHSexp,1.0);

  // FillComplete tangnew (square)
  tangnew->Complete();

  /**********************************************************************/
  /* Replace tang and feff by tangnew and feffnew                     */
  /**********************************************************************/
  tang = tangnew;
  feff = feffnew;

  feff->Scale(-1.0);

  // leave this place
  return;
}

/*----------------------------------------------------------------------*
 | convert maps form structure dofs to thermo dofs            mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::ConvertMaps(RCP<Epetra_Map>& slavedofs,
                                 RCP<Epetra_Map>& activedofs,
                                 RCP<Epetra_Map>& masterdofs)
{

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes/dofs
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // define local variables
    int slavecountnodes = 0;
    vector<int> myslavegids(slavenodes->NumMyElements());

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      myslavegids[slavecountnodes] = (discretstruct_->Dof(1,node))[0];
      ++slavecountnodes;
    }

    // resize the temporary vectors
    myslavegids.resize(slavecountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gslavecountnodes;
    Comm().SumAll(&slavecountnodes,&gslavecountnodes,1);

    // create active node map and active dof map
    slavedofs = rcp(new Epetra_Map(gslavecountnodes,slavecountnodes,&myslavegids[0],0,Comm()));

    // active nodes/dofs
    const RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // define local variables
    int countnodes = 0;
    vector<int> mynodegids(activenodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mynodegids[countnodes] = (discretstruct_->Dof(1,node))[0];
      ++countnodes;
    }

    // resize the temporary vectors
    mynodegids.resize(countnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gcountnodes;
    Comm().SumAll(&countnodes,&gcountnodes,1);

    // create active node map and active dof map
    activedofs = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));

    // master nodes/dofs
    const RCP<Epetra_Map> masternodes = interface[m]->MasterRowNodes();

    // define local variables
    int mastercountnodes = 0;
    vector<int> mymastergids(masternodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mymastergids[mastercountnodes] = (discretstruct_->Dof(1,node))[0];
      ++mastercountnodes;
    }

    // resize the temporary vectors
    mymastergids.resize(mastercountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gmastercountnodes;
    Comm().SumAll(&mastercountnodes,&gmastercountnodes,1);

    // create active node map and active dof map
    masterdofs = rcp(new Epetra_Map(gmastercountnodes,mastercountnodes,&mymastergids[0],0,Comm()));
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mortar matrices in thermo dofs (active nodes)     mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleDM(LINALG::SparseMatrix& dmatrix,
                                LINALG::SparseMatrix& mmatrix)

{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // This is a little bit complicated and a lot of parallel stuff has to
  // be done here. The point is that, when assembling the mortar matrix
  // M, we need the temperature dof from the master node which can lie on
  // a complete different proc. For this reason, we have tho keep all procs
  // around

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveFullNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      /************************************************** D-matrix ******/
      if (Comm().MyPID()==cnode->Owner())
      {
        if ((cnode->MoData().GetD()).size()>0)
        {
          vector<map<int,double> > dmap = cnode->MoData().GetD();
          int rowdisp = cnode->Dofs()[0];
          double val = (dmap[0])[rowdisp];
          dmatrix.Assemble(val, rowtemp, rowtemp);
        }
      }

      /************************************************** M-matrix ******/
      set<int> mnodes;
      int mastergid=0;
      set<int>::iterator mcurr;
      int mastersize = 0;
      vector<map<int,double> > mmap;

      if (Comm().MyPID()==cnode->Owner())
      {
        mmap = cnode->MoData().GetM();
        mnodes = cnode->FriData().GetMNodes();
        mastersize = mnodes.size();
        mcurr = mnodes.begin();
      }

      // commiunicate number of master nodes
      Comm().Broadcast(&mastersize,1,cnode->Owner());

      // loop over all according master nodes
      for (int l=0;l<mastersize;++l)
      {
        if (Comm().MyPID()==cnode->Owner())
          mastergid=*mcurr;

        // communicate GID of masternode
        Comm().Broadcast(&mastergid,1,cnode->Owner());

        DRT::Node* mnode = (interface[m]->Discret()).gNode(mastergid);
        DRT::Node* mnodeges = discretstruct_->gNode(mastergid);

        // temperature and displacement dofs
        int coltemp = 0;
        int coldis = 0;
        if(Comm().MyPID()==mnode->Owner())
        {
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          coltemp = discretstruct_->Dof(1,mnodeges)[0];
          coldis = (cmnode->Dofs())[0];
        }

        // communicate temperature and displacement dof
        Comm().Broadcast(&coltemp,1,mnode->Owner());
        Comm().Broadcast(&coldis,1,mnode->Owner());

        // do the assembly
        if (Comm().MyPID()==cnode->Owner())
        {
          double val = mmap[0][coldis];
          if (abs(val)>1e-12) mmatrix.Assemble(val, rowtemp, coltemp);
          ++mcurr;
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mechanical dissipation for master nodes            mgit 08/10|
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleMechDissRate(Epetra_Vector& mechdissrate)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // time step size
  double dt = GetTimeStepSize();

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> masternodes = interface[m]->MasterRowNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      Epetra_SerialDenseVector mechdissiprate(1);
      vector<int> dof(1);
      vector<int> owner(1);

      mechdissiprate(0) = 1/dt*cnode->MechDiss();
      dof[0] = rowtemp;
      owner[0] = cnode->Owner();

      // do assembly
      if(abs(mechdissiprate(0))>1e-12)
        LINALG::Assemble(mechdissrate, mechdissiprate, dof, owner);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble the thermal contact conditions for slave nodes    mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                                 LINALG::SparseMatrix& thermcontTEMP,
                                                 Epetra_Vector& thermcontRHS,
                                                 LINALG::SparseMatrix& dmatrix,
                                                 LINALG::SparseMatrix& mmatrix,
                                                 RCP<Epetra_Map> activedofs,
                                                 RCP<Epetra_Map> masterdofs)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet and for one heat
  // transfer coefficient
  // FIXGIT: The heat transfer coefficient should be a condition on
  // the single interfaces!!
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::AssembleThermContCondition: Only for one interface yet.");

  // heat transfer coefficient for slave and master surface
  double heattranss = interface[0]->IParams().get<double>("HEATTRANSSLAVE");
  double heattransm = interface[0]->IParams().get<double>("HEATTRANSMASTER");

  if (heattranss <= 0 or heattransm <= 0)
   dserror("Error: Choose realistic heat transfer parameter");

  // time step size
  double dt = GetTimeStepSize();

  double beta = heattranss*heattransm/(heattranss+heattransm);
  double delta = heattranss/(heattranss+heattransm);

  // with respect to Lagrange multipliers
  thermcontLM.Add(dmatrix,false,1.0,1.0);

  // with respect to temperature
  thermcontTEMP.Add(dmatrix,false,-beta,1.0);
  thermcontTEMP.Add(mmatrix,false,+beta,1.0);

  RCP<Epetra_Vector> fa, fm, rest1, rest2;

  // row map of thermal problem
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(discret_->DofRowMap())));

  LINALG::SplitVector(*problemrowmap,*tempn_,activedofs,fa,masterdofs,fm);

  RCP <Epetra_Vector> DdotTemp = rcp(new Epetra_Vector(*activedofs));
  dmatrix.Multiply(false,*fa,*DdotTemp);
  thermcontRHS.Update(beta,*DdotTemp,1.0);

  RCP <Epetra_Vector> MdotTemp = rcp(new Epetra_Vector(*activedofs));
  mmatrix.Multiply(false,*fm,*MdotTemp);
  thermcontRHS.Update(-beta,*MdotTemp,1.0);

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      Epetra_SerialDenseVector mechdissiprate(1);
      vector<int> dof(1);
      vector<int> owner(1);

      mechdissiprate(0) = delta/dt*cnode->MechDiss();
      dof[0] = rowtemp;
      owner[0] = cnode->Owner();

      // do assembly
      if(abs(mechdissiprate(0)>1e-12 and cnode->Active()))
        LINALG::Assemble(thermcontRHS, mechdissiprate, dof, owner);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
