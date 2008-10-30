/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.cpp
\brief Implicit time integration for spatial discretised
       structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>

#include "strtimint.H"
#include "strtimint_impl.H"

/*----------------------------------------------------------------------*/
/* map enum term to input string */
enum STR::TimIntImpl::PredEnum STR::TimIntImpl::MapPredictorStringToEnum
(
  const std::string name  //!< identifier
)
{
  // for explanations please look at #PredEnum
  if (name == "Vague")
  {
    return pred_vague;
  }
  else if (name == "ConstDis")
  {
    return pred_constdis;
  }
  else if (name == "ConstDisVelAcc")
  {
    return pred_constdisvelacc;
  }
  else
  {
    dserror("Unknown kind of predictor %s", name.c_str());
    return pred_vague;
  }
}

/*----------------------------------------------------------------------*/
/* map solution technique identification string to enum term */
enum STR::TimIntImpl::SolTechEnum STR::TimIntImpl::MapSolTechStringToEnum
(
  const std::string name  //!< name identification string
)
{
  if (name == "vague")
  {
    return soltech_vague;
  }
  else if (name == "fullnewton")
  {
    return soltech_newtonfull;
  }
  else if (name == "modnewton")
  {
    return soltech_newtonmod;
  }
  else if (name == "newtonlinuzawa")
  {
    return soltech_uzawalinnewton;
  }
  else if (name == "augmentedlagrange")
  {
    return soltech_uzawanonlinnewton;
  }
  else
  {
    dserror("Unknown type of solution technique %s", name.c_str());
    return soltech_vague;
  }
}

/*----------------------------------------------------------------------*/
/* map solution technique identification string to enum term */
std::string  STR::TimIntImpl::MapSolTechEnumToString
(
  const enum SolTechEnum name  //!< identifying enum term
)
{
  switch (name)
  {
  case soltech_vague :
    return "vague";
    break;
  case soltech_newtonfull :
    return "fullnewton";
    break;
  case soltech_newtonmod :
    return "modnewton";
    break;
  case soltech_uzawalinnewton :
    return "newtonlinuzawa";
    break;
  case soltech_uzawanonlinnewton :
    return "augmentedlagrange";
    break;
  default :
    dserror("Unknown type of solution technique %d", name);
    return "vague";
  }
}

/*----------------------------------------------------------------------*/
/* map convergence check to enum term */
enum STR::TimIntImpl::ConvCheckEnum STR::TimIntImpl::MapConvCheckStringToEnum
(
  const std::string name  //!< name identification string
)
{
  if (name == "AbsRes_Or_AbsDis")
  {
    return convcheck_absres_or_absdis;
  }
  else if (name == "AbsRes_And_AbsDis")
  {
    return convcheck_absres_and_absdis;
  }
  else if (name == "RelRes_Or_AbsDis")
  {
    return convcheck_relres_or_absdis;
  }
  else if (name == "RelRes_And_AbsDis")
  {
    return convcheck_relres_and_absdis;
  }
  else if (name == "RelRes_Or_RelDis")
  {
    return convcheck_relres_or_reldis;
  }
  else if (name == "RelRes_And_RelDis")
  {
    return convcheck_relres_and_reldis;
  }
  else
  {
    dserror("Unknown type %s of iterative convergence check", name.c_str());
    return convcheck_vague;
  }
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntImpl::TimIntImpl
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
  fsisurface_(NULL),
  pred_(MapPredictorStringToEnum(sdynparams.get<string>("PREDICT"))),
  itertype_(MapSolTechStringToEnum(sdynparams.get<string>("NLNSOL"))),
  itercnvchk_(MapConvCheckStringToEnum(sdynparams.get<string>("CONV_CHECK"))),
  iternorm_(TimIntVector::MapNormStringToEnum("L2")),  // ADD INPUT FEATURE
  itermax_(sdynparams.get<int>("MAXITER")),
  toldisi_(sdynparams.get<double>("TOLDISP")),
  tolfres_(sdynparams.get<double>("TOLRES")),
  uzawaparam_(sdynparams.get<double>("UZAWAPARAM")),
  uzawaitermax_(sdynparams.get<int>("UZAWAMAXITER")),
  tolcon_(sdynparams.get<double>("TOLCONSTR")),
  iter_(-1),
  normcharforce_(0.0),
  normchardis_(0.0),
  normfres_(0.0),
  normdisi_(0.0),
  disi_(Teuchos::null),
  timer_(actdis->Comm()),
  fres_(Teuchos::null),
  fifc_(Teuchos::null)
{
  // verify: if system has constraints, then Uzawa-type solver is used
  if ( conman_->HaveConstraint()
       and ( (itertype_ != soltech_uzawalinnewton)
             and (itertype_ != soltech_uzawanonlinnewton) ) )
  {
    dserror("Chosen solution technique %s does not work constrained.",
            MapSolTechEnumToString(itertype_).c_str());
  }

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*dofrowmap_, false);

  // create empty interface force vector
  fifc_ = LINALG::CreateVector(*dofrowmap_, true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap_, true);

  // done so far
  return;
}

/*----------------------------------------------------------------------*/
/* integrate step */
void STR::TimIntImpl::IntegrateStep()
{
  Predict();
  Solve();
  return;
}

/*----------------------------------------------------------------------*/
/* predict solution */
void STR::TimIntImpl::Predict()
{

  // choose predictor
  if (pred_ == pred_constdis)
  {
    PredictConstDisConsistVelAcc();
  }
  else if (pred_ == pred_constdisvelacc)
  {
    PredictConstDisVelAcc();
  }
  else
  {
    dserror("Trouble in determing predictor %i", pred_);
  }

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, accn_, false);

  // possibly initialise Lagrange multiplicators to zero
  //  if ( (conman_->HaveConstraint())
  //       and (itertype_ == soltech_uzawalinnewton) )
  //  {
  //    conman_->ScaleLagrMult(0.0);
  //  }

  // compute residual forces fres_ and stiffness stiff_
  EvaluateForceStiffResidual();

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // determine residual norm of predictor
  normfres_ = TimIntVector::CalculateNorm(iternorm_, fres_);

  // determine characteristic norms
  // we set the minumum of CalcRefNormForce() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchardis_ = CalcRefNormDisplacement();
  if (normchardis_ == 0.0) normchardis_ = toldisi_;

  // output
  PrintPredictor();

  // enjoy your meal
  return;
}

/*----------------------------------------------------------------------*/
/* predict solution as constant displacements, velocities
 * and accelerations */
void STR::TimIntImpl::PredictConstDisVelAcc()
{
  // constant predictor
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  accn_->Update(1.0, *(*acc_)(0), 0.0);

  // see you next time step
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ surface stresses and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::ApplyForceStiffSurfstress
(
  const Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff
)
{
  // surface stress loads (but on internal force vector side)
  if (surfstressman_ != Teuchos::null)
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("surfstr_man", surfstressman_);

    surfstressman_->EvaluateSurfStress(p, dis, disn_, fint, stiff);
  }

  // bye bye
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ potential forces and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::ApplyForceStiffPotential
(
  const Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff
)
{
  // potential force loads (but on internal force vector side)
  if (potman_ != Teuchos::null)
  {
    ParameterList p; // create the parameters for manager
    p.set("pot_man", potman_);
    potman_->EvaluatePotential(p, dis, fint, stiff);
  }

  // bye bye
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate forces due to constraints */
void STR::TimIntImpl::ApplyForceStiffConstraint
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis,
  const Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff,
  Teuchos::ParameterList pcon
)
{
  if (conman_->HaveConstraint())
  {
    conman_->StiffnessAndInternalForces(time, dis, disn, fint, stiff, pcon);
  }

  // wotcha
  return;
}

/*----------------------------------------------------------------------*/
/* check convergence
 * originally by lw 12/07 and tk 01/08 */
bool STR::TimIntImpl::Converged()
{
  // verify: #normcharforce_ has been delivered strictly larger than zero
  if (normcharforce_ <= 0.0)
  {
    dserror("Characteristic force norm %g must be strictly larger than 0",
            normcharforce_);
  }
  // verify: #normchardis_ has been delivered strictly larger than zero
  if (normchardis_ <= 0.0)
  {
    dserror("Characteristic displacement norm %g must be strictly larger than 0",
            normchardis_);
  }

  // check force and displacement residuals
  bool fdc = false;
  switch (itercnvchk_)
  {
  case convcheck_absres_or_absdis:
    fdc = ( (normdisi_ < toldisi_) or (normfres_ < tolfres_) );
    break;
  case convcheck_absres_and_absdis:
    fdc = ( (normdisi_ < toldisi_) and (normfres_ < tolfres_) );
    break;
  case convcheck_relres_or_absdis:
    fdc = ( (normdisi_ < toldisi_)
            or (normfres_/normcharforce_ < tolfres_) );
    break;
  case convcheck_relres_and_absdis:
    fdc = ( (normdisi_ < toldisi_)
            and (normfres_/normcharforce_ < tolfres_) );
    break;
  case convcheck_relres_or_reldis:
    fdc = ( (normdisi_/normchardis_ < toldisi_)
            or (normfres_/normcharforce_ < tolfres_) );
    break;
  case convcheck_relres_and_reldis:
    fdc = ( (normdisi_/normchardis_ < toldisi_)
            and (normfres_/normcharforce_ < tolfres_) );
    break;
  default:
    dserror("Requested convergence check %i is not (yet) implemented",
            itercnvchk_);
    fdc = true;
  }

  // check constraint
  bool cc = true;
  if (conman_->HaveConstraint())
  {
    cc = normcon_ < tolcon_;
  }

  // return things
  return (fdc and cc);
}

/*----------------------------------------------------------------------*/
/* solve equilibrium */
void STR::TimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  switch (itertype_)
  {
  case soltech_newtonfull :
    NewtonFull();
    break;
  case soltech_uzawanonlinnewton :
    UzawaNonLinearNewtonFull();
    break;
  case soltech_uzawalinnewton :
    UzawaLinearNewtonFull();
    break;
  // catch problems
  default :
    dserror("Solution technique %s is not implemented",
            MapSolTechEnumToString(itertype_).c_str());
    break;
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/* solution with full Newton-Raphson iteration */
void STR::TimIntImpl::NewtonFull()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #stiff_ is the effective dynamic stiffness matrix

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled())
  {
    dserror("Effective stiffness matrix must be filled here");
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  normdisi_ = 1.0e6;  // this is strictly >0,toldisi_
  timer_.ResetStartTime();
  //bool print_unconv = true;

  // equilibrium iteration loop
  while ( (not Converged()) and (iter_ <= itermax_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    //LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
    //                               zeros_, dirichtoggle_);
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   zeros_, *(dbcmaps_->CondMap()));

    // solve for disi_
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    solver_->Solve(stiff_->EpetraMatrix(), disi_, fres_, true, iter_==1);
    solver_->ResetTolerance();

    // update end-point displacements etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and stiffness #stiff_
    EvaluateForceStiffResidual();

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

    // build residual force norm
    normfres_ = TimIntVector::CalculateNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = TimIntVector::CalculateNorm(iternorm_, disi_);

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if (iter_ >= itermax_)
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }
  else if ( (Converged()) and (myrank_ == 0) )
  {
    PrintNewtonConv();
  }

  // get out of here
  return;
}

/*----------------------------------------------------------------------*/
/* do non-linear Uzawa iteration within a full NRI is called,
 * originally by tk */
void STR::TimIntImpl::UzawaNonLinearNewtonFull()
{
  // now or never, break it
  dserror("Sorry dude, non-linear Uzawa with full Newton-Raphson"
          " iteration is available in source, but it has not been"
          " tested in silico and should not be used overcredulously."
          " Feel free to remove this dserror but be careful and check"
          " if things run as expected.");

  // do Newton-Raphson iteration, which contains here effects of
  // constraint forces and stiffness
  // this call ends up with new displacements etc on \f$D_{n+1}\f$ etc
  NewtonFull();

  // compute constraint error ...
  conman_->ComputeError(timen_, disn_);
  // ... and its norm
  normcon_ = conman_->GetErrorNorm();
  // talk to user
  if (myrank_ == 0)
  {
    std::cout << "Constraint error for Newton solution: " << normcon_
              << std::endl;
  }

  // Uzawa iteration loop
  int uziter = 0;
  while ( (normcon_ > tolcon_) and (uziter <= uzawaitermax_) )
  {
    // Lagrange multiplier is increased by #uzawaparam_ times ConstrError
    conman_->UpdateLagrMult(uzawaparam_);

    // Keep new Lagrange multiplier fixed and solve for new displacements

    // REALLY NECESSARY, OR EVEN COUNTERPRODUCTIVE ???
    Predict();

    // do Newton-Raphson iteration, which contains here effects of
    // constraint forces and stiffness
    // this call ends up with new displacements etc on \f$D_{n+1}\f$ etc
    NewtonFull();

    // compute constraint error ...
    conman_->ComputeError(timen_, disn_);
    // ... and its norm
    normcon_ = conman_->GetErrorNorm();
    // talk to user
    if (myrank_ == 0)
    {
      std::cout << "Constraint error for computed displacement: " << normcon_
                << std::endl;
    }

    // increment loop counter
    uziter += 1;
  }

  // for output
  iter_ = uziter + 1;
}

/*----------------------------------------------------------------------*/
/* do linearised Uzawa iterations with full NRI
 * originally by tk 11/07 */
void STR::TimIntImpl::UzawaLinearNewtonFull()
{
  // allocate additional vectors and matrices
  Teuchos::RCP<LINALG::SparseMatrix> conmatrix
    = Teuchos::rcp(new LINALG::SparseMatrix(*(conman_->GetConstrMatrix())));
  Teuchos::RCP<Epetra_Vector> conrhs
    = Teuchos::rcp(new Epetra_Vector(*(conman_->GetError())));

  Teuchos::RCP<Epetra_Vector> lagrincr
    = Teuchos::rcp(new Epetra_Vector(*(conman_->GetConstraintMap())));

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled())
  {
    dserror("Effective stiffness matrix must be filled here");
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  normdisi_ = 1.0e6;
  normcon_ = conman_->GetErrorNorm();
  timer_.ResetStartTime();
  
  // equilibrium iteration loop
  while ( (not Converged()) and (iter_ <= itermax_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

//    // uncomplete stiffness matrix, so stuff can be inserted again
//    stiff_->UnComplete();

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   zeros_, *(dbcmaps_->CondMap()));
    // prepare residual Lagrange multiplier
    lagrincr->PutScalar(0.0);
    // Call constraint solver to solve system with zeros on diagonal
    consolv_->Solve(stiff_, conmatrix,
                    disi_, lagrincr,
                    fres_, conrhs);

    
    // update Lagrange multiplier
    conman_->UpdateLagrMult(lagrincr);
    // update end-point displacements etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and stiffness #stiff_
    // which contain forces and stiffness of constraints
    EvaluateForceStiffResidual();
    // compute residual and stiffness of constraint equations
    conmatrix = conman_->GetConstrMatrix();
    conrhs = Teuchos::rcp(new Epetra_Vector(*(conman_->GetError())));

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

    // build residual force norm
    normfres_ = TimIntVector::CalculateNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = TimIntVector::CalculateNorm(iternorm_, disi_);
    // build residual Lagrange multiplier norm
    normcon_ = conman_->GetErrorNorm();

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if (iter_ >= itermax_)
  {
     dserror("Newton unconverged in %d iterations", iter_);
  }
  else
  {
    // monitor values
    if (conman_->HaveMonitor())
    {
       conman_->ComputeMonitorValues(disn_);
    }

    // print message
    if ( Converged() and (myrank_ == 0) )
    {
      PrintNewtonConv();
    }
  }

  // good evening
  return;
}

/*----------------------------------------------------------------------*/
/* Update iteration */
void STR::TimIntImpl::UpdateIter
(
  const int iter  //!< iteration counter
)
{
  if (iter <= 1)
  {
    UpdateIterIncrementally();
  }
  else
  {
    UpdateIterIteratively();
  }
}

/*----------------------------------------------------------------------*/
/* print to screen
 * lw 12/07 */
void STR::TimIntImpl::PrintPredictor()
{
  // only master processor
  if ( (myrank_ == 0) and printscreen_ )
  {
    // relative check of force residual
    if ( (itercnvchk_ != convcheck_absres_or_absdis)
         and (itercnvchk_ != convcheck_absres_and_absdis) )
    {
      std::cout << "Predictor scaled res-norm "
                << normfres_/normcharforce_
                << std::endl;
    }
    // absolute check of force residual
    else
    {
      std::cout << "Predictor absolute res-norm "
                << normfres_
                << std::endl;
    }
    // print it, now
    fflush(stdout);
  }

  // leave your hat on
  return;
}

/*----------------------------------------------------------------------*/
/* print Newton-Raphson iteration to screen and error file
 * originally by lw 12/07, tk 01/08 */
void STR::TimIntImpl::PrintNewtonIter()
{
  // print to standard out
  if ( (myrank_ == 0) and printscreen_ and printiter_ )
  {
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if ( printerrfile_ and printiter_ )
  {
    PrintNewtonIterText(errfile_);
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/* print Newton-Raphson iteration to screen
 * originally by lw 12/07, tk 01/08 */
void STR::TimIntImpl::PrintNewtonIterText
(
  FILE* ofile
)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << " numiter " << std::setw(2) << iter_;

  // different style due relative or absolute error checking
  switch (itercnvchk_)
  {
  // relative residual forces AND displacements
  case convcheck_relres_and_reldis:
  case convcheck_relres_or_reldis:
    oss << " rel-res-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normfres_/normcharforce_
        << " rel-dis-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normdisi_/normchardis_;
    break;
  // relative residual forces
  case convcheck_relres_and_absdis:
  case convcheck_relres_or_absdis:
    oss << " rel-res-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normfres_/normcharforce_
        << " abs-dis-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normdisi_;
    break;
  // absolute forces and displacements
  case convcheck_absres_and_absdis:
  case convcheck_absres_or_absdis:
    oss << " abs-res-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normfres_
        << " abs-dis-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normdisi_;
    break;
  default:
    dserror("Cannot handle requested convergence check %i", itercnvchk_);
    break;
  }  // end switch

  // add constraint norm
  if (conman_->HaveConstraint())
  {
    oss << " abs-constr-norm "
        << std::setw(10) << std::setprecision(5) << std::scientific
        << normcon_;
  }

  // add solution time
  oss << " time " <<  timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);
}

/*----------------------------------------------------------------------*/
/* print statistics of converged NRI */
void STR::TimIntImpl::PrintNewtonConv()
{
  // print constraint manager
  if (conman_->HaveMonitor())
  {
    conman_->PrintMonitorValues();
  }

  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << "Newton iteration converged:"
      << " numiter " << iter_;

  switch (itercnvchk_)
  {
  // relative residual forces AND displacements
  case convcheck_relres_and_reldis:
  case convcheck_relres_or_reldis:
    oss << " scaled res-norm " << normfres_/normcharforce_
        << " scaled dis-norm " << normdisi_/normchardis_;
    break;
  // relative residual forces
  case convcheck_relres_and_absdis:
  case convcheck_relres_or_absdis:
    oss << " scaled res-norm " << normfres_/normcharforce_
        << " absolute dis-norm " << normdisi_;
    break;
  // absolute forces and displacements
  case convcheck_absres_and_absdis:
  case convcheck_absres_or_absdis:
    oss << " absolute res-norm " << normfres_
        << " absolute dis-norm " << normdisi_;
    break;
  default:
    dserror("Cannot handle requested convergence check %i", itercnvchk_);
    break;
  }  // end switch

  // add constraint norm
  if (conman_->HaveConstraint())
  {
    oss << " absolute constr_norm " << normcon_;
  }

  // add solution time
  oss << " time " <<  timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  printf("%s\n", oss.str().c_str());

  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntImpl::PrintStep()
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

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntImpl::PrintStepText
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

/*----------------------------------------------------------------------*/
/* introduce (robin) fsi surface extractor object */
void STR::TimIntImpl::SetSurfaceFSI
(
  const LINALG::MapExtractor* fsisurface  //!< the FSI surface
)
{
  fsisurface_ = fsisurface;
}

/*----------------------------------------------------------------------*/
/* Set forces due to interface with fluid */
void STR::TimIntImpl::SetForceInterface
(
  const LINALG::MapExtractor& extractor,
  Teuchos::RCP<Epetra_Vector> iforce  ///< the force on interface
)
{
  fifc_->PutScalar(0.0);
  extractor.AddCondVector(iforce, fifc_);
}

/*----------------------------------------------------------------------*/
/* Set forces due to interface with fluid,
 * the force is expected external-force-like */
void STR::TimIntImpl::SetForceInterface
(
  Teuchos::RCP<Epetra_Vector> iforce  ///< the force on interface
)
{
  fifc_->Update(1.0, *iforce, 0.0);
}

/*----------------------------------------------------------------------*/
/* Linear structure solve with just an interface load */
Teuchos::RCP<Epetra_Vector> STR::TimIntImpl::SolveRelaxationLinear()
{
  // Evaluate/define the residual force vector #fres_ for
  // relaxation solution with SolveRelaxationLinear
  EvaluateForceStiffResidualRelax();

  // negative residual
  fres_->Scale(-1.0);

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 zeros_, *(dbcmaps_->CondMap()));

  // solve for #disi_
  solver_->Solve(stiff_->EpetraMatrix(), disi_, fres_, true, true);

  // go back
  return disi_;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
