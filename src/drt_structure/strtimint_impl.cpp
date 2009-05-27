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
#include "str_aux.H"


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
  pred_(Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdynparams,"PREDICT")),
  itertype_(Teuchos::getIntegralValue<INPAR::STR::NonlinSolTech>(sdynparams,"NLNSOL")),
  normtypedisi_(Teuchos::getIntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_DISP")),
  normtypefres_(Teuchos::getIntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_RESF")),
  normtypepres_(Teuchos::getIntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_PRES")),
  normtypepfres_(Teuchos::getIntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_INCO")),
  combdispre_(Teuchos::getIntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_DISPPRES")),
  combfrespfres_(Teuchos::getIntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_RESFINCO")),
  combdisifres_(Teuchos::getIntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_RESFDISP")),
  iternorm_(Teuchos::getIntegralValue<INPAR::STR::VectorNorm>(sdynparams,"ITERNORM")),
  itermax_(sdynparams.get<int>("MAXITER")),
  itermin_(sdynparams.get<int>("MINITER")),
  iterdivercont_(Teuchos::getIntegralValue<int>(sdynparams,"DIVERCONT")==1),
  toldisi_(sdynparams.get<double>("TOLDISP")),
  tolfres_(sdynparams.get<double>("TOLRES")),
  tolpfres_(sdynparams.get<double>("TOLINCO")),
  tolpres_(sdynparams.get<double>("TOLPRE")),
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
  freact_(Teuchos::null),
  fifc_(Teuchos::null)
{

  // verify: Old-style convergence check has to be 'vague' to
  if(Teuchos::getIntegralValue<INPAR::STR::ConvCheck>(sdynparams,"CONV_CHECK")!=INPAR::STR::convcheck_vague)
  {
    if (pressure_!=Teuchos::null)
      dserror("For new structural time integration and pressure formulation, please choose CONV_CHECK = None");
    else
      ConvertConvCheck(Teuchos::getIntegralValue<INPAR::STR::ConvCheck>(sdynparams,"CONV_CHECK"));      
  }

  // verify: if system has constraints, then Uzawa-type solver is used
  if ( conman_->HaveConstraint())
  {
    if ( (itertype_ != INPAR::STR::soltech_newtonuzawalin)
               and (itertype_ != INPAR::STR::soltech_newtonuzawanonlin) )
    {
      dserror("Chosen solution technique %s does not work constrained.",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
    }
  }
  else if ( (itertype_ == INPAR::STR::soltech_newtonuzawalin)
            or (itertype_ == INPAR::STR::soltech_newtonuzawanonlin) )
  {
      dserror("Chosen solution technique %s does only work constrained.",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
  }

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*dofrowmap_, false);

  // create empty reaction force vector of full length
  freact_ = LINALG::CreateVector(*dofrowmap_, false);

  // create empty interface force vector
  fifc_ = LINALG::CreateVector(*dofrowmap_, true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap_, true);

  // setup NOX parameter lists
  if (itertype_ == INPAR::STR::soltech_noxnewtonlinesearch)
    NoxSetup();
  else if (itertype_ == INPAR::STR::soltech_noxgeneral)
    NoxSetup(xparams.sublist("NOX"));

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
  if (pred_ == INPAR::STR::pred_constdis)
  {
    PredictConstDisConsistVelAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if (pred_ == INPAR::STR::pred_constdisvelacc)
  {
    PredictConstDisVelAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if (pred_ == INPAR::STR::pred_tangdis)
  {
    PredictTangDisConsistVelAcc();
    // normdisi_ has been set
  }
  else
  {
    dserror("Trouble in determining predictor %i", pred_);
  }

  // zerofy pressure DOFs
  if (pressure_ != Teuchos::null)
    pressure_->InsertCondVector(pressure_->ExtractCondVector(zeros_), disn_);

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

  // rotate to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // rotate back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(fres_);

  // split norms
  if (pressure_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> fres = pressure_->ExtractOtherVector(fres_);
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres);
    Teuchos::RCP<Epetra_Vector> fpres = pressure_->ExtractCondVector(fres_);
    normpfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fpres);

  }
  else
  {
    // build residual force norm
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
  }
 
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
void STR::TimIntImpl::PredictTangDisConsistVelAcc()
{
  // initialise
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  accn_->Update(1.0, *(*acc_)(0), 0.0);
  disi_->PutScalar(0.0);

  // for displacement increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc
    = LINALG::CreateVector(*dofrowmap_, true);

  // copy last converged displacements
  dbcinc->Update(1.0, *(*dis_)(0), 0.0);

  // get Dirichlet values at t_{n+1}
  ApplyDirichletBC(timen_, dbcinc, Teuchos::null, Teuchos::null, false);

  // subtract the displacements of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *(*dis_)(0), 1.0);

  // compute residual forces fres_ and stiffness stiff_
  // at disn_, etc which are unchanged
  EvaluateForceStiffResidual();

  // add linear reaction forces to residual
  {
    // linear reactions
    Teuchos::RCP<Epetra_Vector> freact
      = LINALG::CreateVector(*dofrowmap_, true);
    stiff_->Multiply(false, *dbcinc, *freact);

    // add linear reaction forces due to prescribed Dirichlet BCs
    fres_->Update(1.0, *freact, 1.0);
  }

  // rotate to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  freact_->Update(-1.0, *fres_, 0.0);  // reactions are negative
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // rotate back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(fres_);

  // make negative residual
  fres_->Scale(-1.0);

  // transform to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(stiff_, fres_);

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);
  stiff_->Complete();
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

  // solve for disi_
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  solver_->Reset();
  solver_->Solve(stiff_->EpetraMatrix(), disi_, fres_, true, true);
  solver_->Reset();

  // extract norm of disi_
  if (pressure_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(disi_);
    Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(disi_);
    normpres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);
  }
  else
  {
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
  }

  // set Dirichlet increments in displacement increments
  disi_->Update(1.0, *dbcinc, 1.0);

  // update end-point displacements etc
  UpdateIterIncrementally();
  //disn_->Update(1.0, *disi_, 1.0);

  // MARK:
  // velocities and accelerations unset on Dirichlet boundary

  // reset to zero
  disi_->PutScalar(0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_struct_reset_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // shalom
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ surface stresses and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::ApplyForceStiffSurfstress
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dism,
  const Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseMatrix>& stiff
)
{
  // surface stress loads (but on internal force vector side)
  if (surfstressman_->HaveSurfStress())
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("surfstr_man", surfstressman_);
    p.set("total time", time);
    p.set("delta time", dt);
    surfstressman_->EvaluateSurfStress(p, dism, disn, fint, stiff);
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

  // wooop
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


  // check for single norms
  bool convdis = false;
  bool convfres = false;
    
  // residual displacement
  switch (normtypedisi_)
  {
  case INPAR::STR::convnorm_abs:
    convdis = normdisi_ < toldisi_;
    break;
  case INPAR::STR::convnorm_rel:
    convdis = normdisi_/normchardis_ < toldisi_;
    break;
  case INPAR::STR::convnorm_mix:
    convdis = ( (normdisi_ < toldisi_) or (normdisi_/normchardis_ < toldisi_) );
    break;
  default:
    dserror("Cannot check for convergence of residual displacements!");
  }
  
  // residual forces
  switch (normtypefres_)
  {
  case INPAR::STR::convnorm_abs:
    convfres = normfres_ < tolfres_;
    break;
  case INPAR::STR::convnorm_rel:
    convfres = normfres_/normcharforce_ < tolfres_;
    break;
  case INPAR::STR::convnorm_mix:
    convfres = ( (normfres_ < tolfres_) or (normfres_/normcharforce_ < tolfres_) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
  }
  
  // check constraint
  bool cc = true;
  if (conman_->HaveConstraint())
  {
    cc = normcon_ < tolcon_;
  }
  
  //pressure related stuff
  if (pressure_ != Teuchos::null)
  {
    bool convpre = false;
    bool convfpre = false;
    
    //pressure
    switch (normtypepres_)
    {
    case INPAR::STR::convnorm_abs:
      convpre = normpres_ < tolpres_;
      break;
     default:
      dserror("Cannot check for convergence of residual pressures! Only for absolute residuals implemeted so far!");
    }
    
    // incompressible residual
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs:
      convfpre = normpfres_ < tolpfres_;
      break;
    default:
      dserror("Cannot check for convergence of incompressible force residuals!");
    }
    
    
    // combine fields
    if (combdispre_==INPAR::STR::bop_and)
      convdis = convdis and convpre;
    else if (combdispre_==INPAR::STR::bop_or)
      convdis = convdis or convpre;
    else
      dserror("Something went terribly wrong with binary operator!");
    
    if (combfrespfres_==INPAR::STR::bop_and)
      convfres = convfres and convfpre;
    else if (combfrespfres_==INPAR::STR::bop_or)
      convfres = convfres or convfpre;
    else
      dserror("Something went terribly wrong with binary operator!");
    
  }
  
  // combine displacement-like and force-like residuals
  bool conv = false;
  if (combdisifres_ == INPAR::STR::bop_and)
     conv = convdis and convfres;
   else if (combdisifres_ == INPAR::STR::bop_or)
     conv = convdis or convfres;
   else
     dserror("Something went terribly wrong with binary operator!");

  
  // return things
  return (conv and cc);
}

/*----------------------------------------------------------------------*/
/* solve equilibrium */
void STR::TimIntImpl::Solve()
{
  // choose solution technique in accordance with user's will
  switch (itertype_)
  {
  case INPAR::STR::soltech_newtonfull :
    NewtonFull();
    break;
  case INPAR::STR::soltech_newtonuzawanonlin :
    UzawaNonLinearNewtonFull();
    break;
  case INPAR::STR::soltech_newtonuzawalin :
    UzawaLinearNewtonFull();
    break;
  case INPAR::STR::soltech_noxnewtonlinesearch :
  case INPAR::STR::soltech_noxgeneral :
    NoxSolve();
    break;
  // catch problems
  default :
    dserror("Solution technique \"%s\" is not implemented",
            INPAR::STR::NonlinSolTechString(itertype_).c_str());
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
  // normdisi_ was already set in predictor; this is strictly >0
  timer_.ResetStartTime();

  // equilibrium iteration loop
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {

    // make negative residual
    fres_->Scale(-1.0);

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(stiff_, fres_);

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

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

#if 0
    if (pressure_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(disn_);
      Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(disn_);
      double normpres = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      double normdisp = STR::AUX::CalculateVectorNorm(iternorm_, disp);
      if (myrank_ == 0)
        cout << std::scientific << "DispCurr=" << normdisp << " PresCurr=" << normpres;
    }
#endif

    // compute residual forces #fres_ and stiffness #stiff_
    // whose components are globally oriented
    EvaluateForceStiffResidual();

    // blank residual at (locally oriented) Dirichlet DOFs
    // rotate to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(fres_);

    // extract reaction forces
    // reactions are negative to balance residual on DBC
    freact_->Update(-1.0, *fres_, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
    // rotate reaction forces back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(freact_);

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
    // rotate back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(fres_);


    if (pressure_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(fres_);
      Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(fres_);
      normpfres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);

      pres = pressure_->ExtractCondVector(disi_);
      disp = pressure_->ExtractOtherVector(disi_);
      normpres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);
    }
    else
    {
      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    }

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // call monitor
  if (conman_->HaveMonitor())
  {
    conman_->ComputeMonitorValues(disn_);
  }

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
  // normdisi_ was already set in predictor; this is strictly >0
  normcon_ = conman_->GetErrorNorm();
  timer_.ResetStartTime();

  // equilibrium iteration loop
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

//    // uncomplete stiffness matrix, so stuff can be inserted again
//    stiff_->UnComplete();

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(stiff_, fres_);

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

    // prepare residual Lagrange multiplier
    lagrincr->PutScalar(0.0);
    // Call constraint solver to solve system with zeros on diagonal
    consolv_->Solve(stiff_, conmatrix,
                    disi_, lagrincr,
                    fres_, conrhs);
    
    // transform back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(disi_);

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

    // blank residual at (locally oriented) Dirichlet DOFs
    // rotate to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(fres_);

    // extract reaction forces
    // reactions are negative to balance residual on DBC
    freact_->Update(-1.0, *fres_, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
    // rotate reaction forces back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(freact_);

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
    // rotate back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(fres_);

    // build residual force norm
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    // build residual Lagrange multiplier norm
    normcon_ = conman_->GetErrorNorm();
    
    if (pressure_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(fres_);
      Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(fres_);
      normpfres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);

      pres = pressure_->ExtractCondVector(disi_);
      disp = pressure_->ExtractOtherVector(disi_);
      normpres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);
    }
    else
    {
      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    }

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
  else if ( (iter_ >= itermax_) and (iterdivercont_) )
  {
    // print newton message on proc 0
    if (myrank_ == 0)
      printf("Newton unconverged in %d iterations ... continuing\n", iter_);
    
    // compute and print monitor values
    if (conman_->HaveMonitor())
      conman_->ComputeMonitorValues(disn_);
    
  }
  else if ( Converged() )
  {
    // compute and print monitor values
    if (conman_->HaveMonitor())
    {
      conman_->ComputeMonitorValues(disn_);
    }
    
    
    // print newton message on proc 0
    if (myrank_ == 0)
      conman_->PrintMonitorValues();
    
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

/*----------------------------------------------------------------------*/
/* Update iteration incrementally with prescribed residual displacements */
void STR::TimIntImpl::UpdateIterIncrementally
(
  const Teuchos::RCP<const Epetra_Vector> disi  //!< input residual displacements
)
{
  // select residual displacements
  if (disi != Teuchos::null)
    disi_->Update(1.0, *disi, 0.0);  // set the new solution we just got
  else
    disi_->PutScalar(0.0);

  // Update using #disi_
  UpdateIterIncrementally();

  // leave this place
  return;
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
    if ( normtypefres_ == INPAR::STR::convnorm_rel )
    {
      std::cout << "Predictor scaled res-norm "
                << normfres_/normcharforce_
                << std::endl;
    }
    // absolute check of force residual
    else if ( normtypefres_ == INPAR::STR::convnorm_abs )
    {
      std::cout << "Predictor absolute res-norm "
                << normfres_
                << std::endl;
    }
    // mixed absolute-relative check of force residual
    else if ( normtypefres_ == INPAR::STR::convnorm_mix )
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

/*----------------------------------------------------------------------*/
/* print Newton-Raphson iteration to screen and error file
 * originally by lw 12/07, tk 01/08 */
void STR::TimIntImpl::PrintNewtonIter()
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

void STR::TimIntImpl::PrintNewtonIterHeader
(
  FILE* ofile
)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter";
  
  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::STR::convnorm_rel:
    oss <<std::setw(18)<< "rel-res-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss <<std::setw(18)<< "abs-res-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss <<std::setw(18)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }
  
  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss <<std::setw(18)<< "abs-pre-res-norm";
      break;
    default:
      dserror("You should not turn up here.");
    }
  }
  
  switch ( normtypedisi_ )
  {
  case INPAR::STR::convnorm_rel:
    oss <<std::setw(18)<< "rel-dis-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss <<std::setw(18)<< "abs-dis-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss <<std::setw(18)<< "mix-dis-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }
  
  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss <<std::setw(18)<< "abs-pre-norm";
      break;
    default:
      dserror("You should not turn up here.");
    }
  }
  
  // add constraint norm
  if (conman_->HaveConstraint())
  {
    oss << std::setw(18)<< "abs-constr-norm";
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
  oss << std::setw(7)<< iter_;
  
  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::STR::convnorm_rel:
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_/normcharforce_;
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normfres_;
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << min(normfres_, normfres_/normcharforce_);
    break;
  default:
    dserror("You should not turn up here.");
  }
  
  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normpfres_;
      break;
    default:
      dserror("You should not turn up here.");
    }
  }
  
  switch ( normtypedisi_ )
  {
  case INPAR::STR::convnorm_rel:
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normdisi_/normchardis_;
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normdisi_;
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << min(normdisi_, normdisi_/normchardis_);
    break;
  default:
    dserror("You should not turn up here.");
  }
  
  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss << std::setw(18) << std::scientific << normpres_;
      break;
    default:
      dserror("You should not turn up here.");
    }
  }
  
  // add constraint norm
  if (conman_->HaveConstraint())
  {
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normcon_;
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


/*----------------------------------------------------------------------*/
/* print statistics of converged NRI */
void STR::TimIntImpl::PrintNewtonConv()
{
  // print constraint manager's lore
  if (conman_->HaveMonitor())
  {
    conman_->PrintMonitorValues();
  }

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
/* Prepare system for solving with Newton's method */
void STR::TimIntImpl::PrepareSystemForNewtonSolve()
{

  // rotate residual to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);

  // make the residual negative
  fres_->Scale(-1.0);
  // blank residual at DOFs on Dirichlet BCs
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

  // final sip
  return;
}

void STR::TimIntImpl::ConvertConvCheck
(
    INPAR::STR::ConvCheck convcheck 
)
{
  if (myrank_ == 0)
  {
    // print a beautiful line made exactly of 80 dashes
    std::cout << endl << "--------------------------------------------------------------------------------" << endl;
    std::cout << "WARNING! Duplicated input!"<<endl;
    std::cout << "For new structural time integration input variable 'CONV_CHECK' should be set to 'None',"<<endl;
    std::cout << "and convergence check should be defined via 'NORM_DISP', 'NORM_FRES' and 'NORMCOMBI_RESFDISP'. "<<endl;
    std::cout << "These values will now be overwritten as you can see below! "<<endl<<endl;
  }
  
  // conversion
  switch (convcheck)  
  {
  case INPAR::STR::convcheck_absres_or_absdis:
    normtypefres_ = INPAR::STR::convnorm_abs;
    normtypedisi_ = INPAR::STR::convnorm_abs;
    combdisifres_ = INPAR::STR::bop_or;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Abs"<<endl;
      std::cout << "NORM_DISP                Abs"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       Or"<<endl;
    }
    break;
  case INPAR::STR::convcheck_absres_and_absdis:
    normtypefres_ = INPAR::STR::convnorm_abs;
    normtypedisi_ = INPAR::STR::convnorm_abs;
    combdisifres_ = INPAR::STR::bop_and;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Abs"<<endl;
      std::cout << "NORM_DISP                Abs"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       And"<<endl;
    }
    break;
  case INPAR::STR::convcheck_relres_or_absdis:
    normtypefres_ = INPAR::STR::convnorm_rel;
    normtypedisi_ = INPAR::STR::convnorm_abs;
    combdisifres_ = INPAR::STR::bop_or;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Rel"<<endl;
      std::cout << "NORM_DISP                Abs"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       Or"<<endl;
    }
    break;
  case INPAR::STR::convcheck_relres_and_absdis:
    normtypefres_ = INPAR::STR::convnorm_rel;
    normtypedisi_ = INPAR::STR::convnorm_abs;
    combdisifres_ = INPAR::STR::bop_and;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Rel"<<endl;
      std::cout << "NORM_DISP                Abs"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       And"<<endl;
    }
    break;
  case INPAR::STR::convcheck_relres_or_reldis:
    normtypefres_ = INPAR::STR::convnorm_rel;
    normtypedisi_ = INPAR::STR::convnorm_rel;
    combdisifres_ = INPAR::STR::bop_or;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Rel"<<endl;
      std::cout << "NORM_DISP                Rel"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       Or"<<endl;
    }
    break;
  case INPAR::STR::convcheck_relres_and_reldis:
    normtypefres_ = INPAR::STR::convnorm_rel;
    normtypedisi_ = INPAR::STR::convnorm_rel;
    combdisifres_ = INPAR::STR::bop_and;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Rel"<<endl;
      std::cout << "NORM_DISP                Rel"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       And"<<endl;
    }
    break;
  case INPAR::STR::convcheck_mixres_or_mixdis:
    normtypefres_ = INPAR::STR::convnorm_mix;
    normtypedisi_ = INPAR::STR::convnorm_mix;
    combdisifres_ = INPAR::STR::bop_or;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Mix"<<endl;
      std::cout << "NORM_DISP                Mix"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       Or"<<endl;
    }
    break;
  case INPAR::STR::convcheck_mixres_and_mixdis:
    normtypefres_ = INPAR::STR::convnorm_mix;
    normtypedisi_ = INPAR::STR::convnorm_mix;
    combdisifres_ = INPAR::STR::bop_and;
    if (myrank_ == 0)
    {
      std::cout << "NORM_RESF                Mix"<<endl;
      std::cout << "NORM_DISP                Mix"<<endl;
      std::cout << "NORMCOMBI_RESFDISP       And"<<endl;
    }
    break;
  default:
    dserror("Requested convergence check %i cannot (yet) be converted",
            convcheck);
  }
  std::cout << "--------------------------------------------------------------------------------" << endl << endl;
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
