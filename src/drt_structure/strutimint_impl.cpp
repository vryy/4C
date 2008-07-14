/*----------------------------------------------------------------------*/
/*!
\file strutimint_impl.cpp
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

#include "strutimint.H"
#include "strutimint_impl.H"

/*----------------------------------------------------------------------*/
/* map enum term to input string */
enum StruTimIntImpl::PredEnum StruTimIntImpl::MapPredictorStringToEnum
(
  const std::string name  //!< identifier
)
{
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
enum StruTimIntImpl::SolTechEnum StruTimIntImpl::MapSolTechStringToEnum
(
  const std::string name  //!< name identification string
)
{
  if (name == "Vague")
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
    return soltech_newtonuzawalin;
  }
  else if (name == "augmentedlagrange")
  {
    return soltech_newtonuzawanonlin;
  }
  else
  {
    dserror("Unknown type of solution technique %s", name.c_str());
    return soltech_vague;
  }
}

/*----------------------------------------------------------------------*/
/* map convergence check to enum term */
enum StruTimIntImpl::ConvCheckEnum StruTimIntImpl::MapConvCheckStringToEnum
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
StruTimIntImpl::StruTimIntImpl
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  DRT::Discretization& actdis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: StruTimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
  pred_(MapPredictorStringToEnum(sdynparams.get<string>("PREDICT"))),
  constrman_(Teuchos::null),
  uzawasolv_(Teuchos::null),
  surfstressman_(Teuchos::null),
  potman_(Teuchos::null),
  itertype_(MapSolTechStringToEnum(sdynparams.get<string>("NLNSOL"))),
  itercnvchk_(MapConvCheckStringToEnum(sdynparams.get<string>("CONV_CHECK"))),
  iternorm_(StruTimIntVector::MapNormStringToEnum("L2")),  // ADD INPUT FEATURE
  itermax_(sdynparams.get<int>("MAXITER")),
  toldisi_(sdynparams.get<double>("TOLDISP")),
  tolfres_(sdynparams.get<double>("TOLRES")),
  tolcon_(sdynparams.get<double>("TOLCONSTR")),
  iter_(-1),
  normcharforce_(0.0),
  normchardis_(0.0),
  normfres_(0.0),
  normdisi_(0.0),
  disi_(Teuchos::null),
  timer_(actdis.Comm()),
  fres_(Teuchos::null)
{

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*dofrowmap_, false);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap_, true);

  // initialize constraint manager
  constrman_ = rcp(new ConstrManager(Discretization(), dis_, sdynparams));
  // initialize Uzawa solver
  uzawasolv_ = rcp(new UzawaSolver(Discretization(), solver_, 
                                   dirichtoggle_, invtoggle_, 
                                   sdynparams));
  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_.DofRowMap();

  // Check for surface stress conditions due to interfacial phenomena
  {
    vector<DRT::Condition*> surfstresscond(0);
    discret_.GetCondition("SurfaceStress",surfstresscond);
    if (surfstresscond.size())
    {
      surfstressman_ = rcp(new DRT::SurfStressManager(discret_));
    }
  }
  
  // Check for potential conditions 
  {
    vector<DRT::Condition*> potentialcond(0);
    discret_.GetCondition("Potential",potentialcond);
    if (potentialcond.size())
    {
      potman_ = rcp(new DRT::PotentialManager(discret_));
    }
  }

  // done so far
  return;
}

/*----------------------------------------------------------------------*/
/* integrate step */
void StruTimIntImpl::IntegrateStep()
{
  Predict();
  Solve();
  return;
}

/*----------------------------------------------------------------------*/
/* predict solution */
void StruTimIntImpl::Predict()
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

  // compute residual forces fres_ and stiffness stiff_
  EvaluateForceStiffResidual();

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector frescopy(*fres_);
    fres_->Multiply(1.0, *invtoggle_, frescopy, 0.0);
  }

  // determine residual norm of predictor
  normfres_ = StruTimIntVector::CalculateNorm(iternorm_, fres_);

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
void StruTimIntImpl::PredictConstDisVelAcc()
{
  // constant predictor
  disn_->Update(1.0, *dis_, 0.0);
  veln_->Update(1.0, *vel_, 0.0);
  accn_->Update(1.0, *acc_, 0.0);  
  
  // see you next time step
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ surface stresses and stiffness
 * evaluation happens internal-force like */
void StruTimIntImpl::ApplyForceStiffSurfstress
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
void StruTimIntImpl::ApplyForceStiffPotential
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
void StruTimIntImpl::ApplyForceStiffConstraint
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff
)
{
  if (constrman_->HaveConstraint())
  {
    constrman_->StiffnessAndInternalForces(time, dis, fint, stiff);
  }

  // wotcha
  return;
}

/*----------------------------------------------------------------------*/
/* check convergence
 * originally by lw 12/07 and tk 01/08 */
bool StruTimIntImpl::Converged()
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
  if (constrman_->HaveConstraint())
  {
    cc = normcon_ < tolcon_;
  }

  // return things
  return (fdc and cc);
}

/*----------------------------------------------------------------------*/
/* solve equilibrium */
void StruTimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  if (itertype_ == soltech_newtonfull)
  {
    NewtonFull();
  }
  // catch problems
  else
  {
    dserror("Solution technique %i is not available", itertype_);
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/* solution with full Newton-Raphson iteration */
void StruTimIntImpl::NewtonFull()
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
      solver_.AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    solver_.Solve(stiff_->EpetraMatrix(), disi_, fres_, true, iter_==1);
    solver_.ResetTolerance();

    // update end-point displacements etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and stiffness #stiff_
    EvaluateForceStiffResidual();

    // build residual force norm
    normfres_ = StruTimIntVector::CalculateNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = StruTimIntVector::CalculateNorm(iternorm_, disi_);

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
/* Update iteration */
void StruTimIntImpl::UpdateIter
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
void StruTimIntImpl::PrintPredictor()
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
void StruTimIntImpl::PrintNewtonIter()
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
void StruTimIntImpl::PrintNewtonIterText
(
  FILE* ofile
)
{
  // different style due relative or absolute error checking
  switch (itercnvchk_)
  {
  // relative residual forces AND displacements
  case convcheck_relres_and_reldis:
  case convcheck_relres_or_reldis:
    fprintf(ofile,
            "numiter %2d"
            " scaled res-norm %10.5e"
            " scaled dis-norm %20.15E\n",
            iter_, normfres_/normcharforce_, normdisi_/normchardis_);
    break;
  // relative residual forces
  case convcheck_relres_and_absdis:
  case convcheck_relres_or_absdis:
    fprintf(ofile, 
            "numiter %2d"
            " scaled res-norm %10.5e"
            " absolute dis-norm %20.15E\n",
            iter_, normfres_/normcharforce_, normdisi_);
    break;
  // absolute forces and displacements
  case convcheck_absres_and_absdis:
  case convcheck_absres_or_absdis:
    fprintf(ofile,
            "numiter %2d"
            " absolute res-norm %10.5e"
            " absolute dis-norm %20.15E\n",
            iter_, normfres_, normdisi_);
    break;
  default:
    dserror("Cannot handle requested convergence check %i", itercnvchk_);
    break;
  }  // end switch
  // print it, now
  fflush(ofile);
}

/*----------------------------------------------------------------------*/
/* print statistics of converged NRI */
void StruTimIntImpl::PrintNewtonConv()
{
  // print constraint manager
  if (constrman_->HaveMonitor())
  {
    constrman_->PrintMonitorValues();
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
  if (constrman_->HaveConstraint())
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
void StruTimIntImpl::PrintStep()
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
void StruTimIntImpl::PrintStepText
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
          step_, stepmax_, time_, dt_, iter_);
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
void StruTimIntImpl::OutputStep()
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
void StruTimIntImpl::OutputRestart
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
    output_.WriteMesh(step_, time_);
    output_.NewStep(step_, time_);
    output_.WriteVector("displacement", dis_);
    output_.WriteVector("velocity", vel_);
    output_.WriteVector("acceleration", acc_);
    //output_.WriteVector("fexternal", fext_);  // CURRENTLY NOT AVAILABLE THINK OF SCENARIO

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
      output_.WriteVector("Aold", A);
      output_.WriteVector("conquot", con);
    }
    
    // potential forces
    if (potman_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Map> surfrowmap = potman_->GetSurfRowmap();
      Teuchos::RCP<Epetra_Vector> A 
        = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
      potman_->GetHistory(A);
      output_.WriteVector("Aold", A);
    }

    // constraints
    if (constrman_->HaveConstraint())
    {
      output_.WriteDouble("uzawaparameter",
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
void StruTimIntImpl::OutputState
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
    output_.NewStep(step_, time_);
    output_.WriteVector("displacement", dis_);
    output_.WriteVector("velocity", vel_);
    output_.WriteVector("acceleration", acc_);
    //output_.WriteVector("fexternal",fext_);  // CURRENTLY NOT AVAILABLE
    output_.WriteElementData();
  }

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* stress output
 * originally by lw */
void StruTimIntImpl::OutputStressStrain
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
    p.set("total time", time_);
    p.set("delta time", dt_);
    
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
    discret_.ClearState();
    discret_.SetState("residual displacement", zeros_);
    discret_.SetState("displacement", dis_);
    discret_.Evaluate(p, null, null, null, null, null);
    discret_.ClearState();

    if (not datawritten)
    {
      output_.NewStep(step_, time_);
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
      output_.WriteVector(stresstext, *stressdata, 
                          *discret_.ElementColMap());
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
      output_.WriteVector(straintext, *straindata,
                          *discret_.ElementColMap());
    }
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
