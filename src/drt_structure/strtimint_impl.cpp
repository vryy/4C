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
  fres_(Teuchos::null)
{
  // verify: if system has constraints, then Uzawa-type solver is used
  if ( conman_->HaveConstraint()
       and ( (itertype_ != soltech_uzawalinnewton)
             and (itertype_ != soltech_uzawanonlinnewton) ) )
  {
    dserror("Chosen solution technique %s does not work constraints.",
            MapSolTechEnumToString(itertype_).c_str());
  }

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*dofrowmap_, false);

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
  ApplyDirichletBC(timen_, disn_, veln_, accn_);

  // initialise Lagrange multiplicators to zero
  if ( (conman_->HaveConstraint())
       and (itertype_ == soltech_uzawalinnewton) )
  {
    conman_->ScaleLagrMult(0.0);
  }

  // compute residual forces fres_ and stiffness stiff_
  EvaluateForceStiffResidual();

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector frescopy(*fres_);
    fres_->Multiply(1.0, *invtoggle_, frescopy, 0.0);
  }

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

    surfstressman_->EvaluateSurfStress(p, dis, fint, stiff);
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
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff
)
{
  if (conman_->HaveConstraint())
  {
    conman_->StiffnessAndInternalForces(time, dis, fint, stiff);
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
  // --> On #fres_ is the negative force residuum
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
    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_, 
                                   zeros_, dirichtoggle_);

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
    {
      Epetra_Vector frescopy(*fres_);
      fres_->Multiply(1.0, *invtoggle_, frescopy, 0.0);
    }

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
  else if (Converged())
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
          " tested in silico and should not be used overcredulously");

  // do Newton-Raphson iteration, which contains here effects of
  // constraint forces and stiffness
  // this call ends up with new displacements etc on \f$D_{n+1}\f$ etc
  NewtonFull();

  // compute constraint error ...
  conman_->ComputeError(timen_, disn_);
  // ... and its norm
  normcon_ = conman_->GetErrorNorm();
  // talk to user
  std::cout << "Constraint error for Newton solution: " << normcon_
            << std::endl;

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
    std::cout << "Constraint error for computed displacement: " << normcon_
              << std::endl;
    
    // increment loop counter
    uziter += 1;
  }

  // SENSIBLE??? FOR OUTPUT???
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
  //bool print_unconv = true;

  // equilibrium iteration loop
  while ( (not Converged()) and (iter_ <= itermax_) )
  {
    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   zeros_, dirichtoggle_);
    // prepare residual Lagrange multiplier
    lagrincr->PutScalar(0.0);
    // Call Uzawa algorithm to solve system with zeros on diagonal
    uzawasolv_->Solve(stiff_, conmatrix,
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
    {
      Epetra_Vector frescopy(*fres_);
      fres_->Multiply(1.0, *invtoggle_, frescopy, 0.0);
    }

    // build residual force norm
    normfres_ = TimIntVector::CalculateNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = TimIntVector::CalculateNorm(iternorm_, disi_);
    // build residual Lagrange multipilcator norm
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
      // WARNING: THIS WAS dism_, BUT WE DO NOT WANT THIS HERE!!!
      //          NEED TO TALK TO THOMAS KLOEPPEL, IF disn_ WORKS AS WELL
      conman_->ComputeMonitorValues(disn_);
    }

    // print message
    if (Converged())
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
         and (itercnvchk_ != convcheck_absres_or_absdis) )
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
    if (conman_->HaveMonitor()) conman_->PrintMonitorValues();
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
/* output to file
 * originally by mwgee 03/07 */
void STR::TimIntImpl::OutputStep()
{
  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  OutputRestart(datawritten);

  // output results (not necessary if restart in same step)
  OutputState(datawritten);

  // output stress & strain
  OutputStressStrain(datawritten);

  // what's next?
  return;
}

/*----------------------------------------------------------------------*/
/* write restart
 * originally mwgee 03/07 */
void STR::TimIntImpl::OutputRestart
(
  bool& datawritten
)
{
  // write restart step
  if (writerestartevery_ and (step_%writerestartevery_ == 0) )
  {
    // Yes, we are going to write...
    datawritten = true;

    // write restart output, please
    output_->WriteMesh(step_, (*time_)[0]);
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
    //output_->WriteVector("fexternal", fext_);  // CURRENTLY NOT AVAILABLE THINK OF SCENARIO

    // surface stress
    if (surfstressman_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Map> surfrowmap 
        = surfstressman_->GetSurfRowmap();
      Teuchos::RCP<Epetra_Vector> A 
        = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
      Teuchos::RCP<Epetra_Vector> con 
        = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
      surfstressman_->GetHistory(A,con);
      output_->WriteVector("Aold", A);
      output_->WriteVector("conquot", con);
    }
    
    // potential forces
    if (potman_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Map> surfrowmap = potman_->GetSurfRowmap();
      Teuchos::RCP<Epetra_Vector> A 
        = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
      potman_->GetHistory(A);
      output_->WriteVector("Aold", A);
    }

    // constraints
    if (conman_->HaveConstraint())
    {
      output_->WriteDouble("uzawaparameter",
                          uzawasolv_->GetUzawaParameter());
    }

    // info dedicated to user's eyes staring at standard out
    if ( (myrank_ == 0) and printscreen_)
    { 
      printf("====== Restart written in step %d\n", step_);
      fflush(stdout);
    }

    // info dedicated to processor error file
    if (printerrfile_)
    {
      fprintf(errfile_, "====== Restart written in step %d\n", step_);
      fflush(errfile_);
    }
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations
 * originally mwgee 03/07 */
void STR::TimIntImpl::OutputState
(
  bool& datawritten
)
{
  if ( writestate_ 
       and writestateevery_ and (step_%writestateevery_ == 0)
       and (not datawritten) )
  {
    // Yes, we are going to write...
    datawritten = true;

    // write now
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
    //output_->WriteVector("fexternal",fext_);  // CURRENTLY NOT AVAILABLE
    output_->WriteElementData();
  }

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* stress output
 * originally by lw */
void STR::TimIntImpl::OutputStressStrain
(
  bool& datawritten
)
{
  // do stress calculation and output
  if ( writestrevery_
       and ( (writestress_ != stress_none)
             or (writestrain_ != strain_none) )
       and (step_%writestrevery_ == 0) )
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    
    // stress
    if (writestress_ == stress_cauchy)
    {
      // output of Cauchy stresses instead of 2PK stresses
      p.set("cauchy", true);
    }
    else
    {
      // this will produce 2nd PK stress ????
      p.set("cauchy", false);
    }
    Teuchos::RCP<std::vector<char> > stressdata
      = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stressdata);

    // strain
    if (writestrain_ == strain_ea)
    {
      p.set("iostrain", "euler_almansi");
    }
    else if (writestrain_ == strain_gl)
    {
      // WILL THIS CAUSE TROUBLE ????
      // THIS STRING DOES NOT EXIST IN SO3
      p.set("iostrain", "green_lagrange");
    }
    else
    {
      p.set("iostrain", "none");
    }
    Teuchos::RCP<std::vector<char> > straindata
      = Teuchos::rcp(new std::vector<char>());
    p.set("strain", straindata);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement", zeros_);
    discret_->SetState("displacement", (*dis_)(0));
    discret_->Evaluate(p, null, null, null, null, null);
    discret_->ClearState();

    if (not datawritten)
    {
      output_->NewStep(step_, (*time_)[0]);
    }
    datawritten = true;

    // write stress
    if (writestress_ != stress_none)
    {
      std::string stresstext = "";
      if (writestress_ == stress_cauchy)
      {
        stresstext = "gauss_cauchy_stresses_xyz";
      }
      else if (writestress_ == stress_pk2)
      {
        stresstext = "gauss_2PK_stresses_xyz";
      }
      output_->WriteVector(stresstext, *stressdata, 
                           *(discret_->ElementColMap()));
    }

    // write strain
    if (writestrain_ != strain_none)
    {
      std::string straintext = "";
      if (writestrain_ == strain_ea)
      {
        straintext = "gauss_EA_strains_xyz";
      }
      else
      {
        straintext = "gauss_GL_strains_xyz";
      }
      output_->WriteVector(straintext, *straindata,
                           *(discret_->ElementColMap()));
    }
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
