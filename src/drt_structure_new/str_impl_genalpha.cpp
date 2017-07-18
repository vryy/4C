/*-----------------------------------------------------------*/
/*!
\file str_impl_genalpha.cpp

\brief Generalized Alpha time integrator.

\maintainer Michael Hiermeier

\date Dec 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_genalpha.H"
#include "str_dbc.H"
#include "str_utils.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::GenAlpha::GenAlpha()
    : beta_(-1.0),
      gamma_(-1.0),
      alphaf_(-1.0),
      alpham_(-1.0),
      rhoinf_(-1.0),
      const_vel_acc_update_ptr_(Teuchos::null),
      fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null)

{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::Setup()
{
  CheckInit();
  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  const STR::TIMINT::GenAlphaDataSDyn& genalpha_sdyn =
      dynamic_cast<const STR::TIMINT::GenAlphaDataSDyn&>(TimInt().GetDataSDyn());

  // ---------------------------------------------------------------------------
  // setup time integration parameters
  // ---------------------------------------------------------------------------
  // get a copy of the input parameters
  beta_ = genalpha_sdyn.GetBeta();
  gamma_ = genalpha_sdyn.GetGamma();
  alphaf_ = genalpha_sdyn.GetAlphaF();
  alpham_ = genalpha_sdyn.GetAlphaM();
  rhoinf_ = genalpha_sdyn.GetRhoInf();

  // ------ rho_inf specified --> calculate optimal parameters -----------------
  if (rhoinf_!=-1.)
  {
    if ( (rhoinf_ < 0.0) or (rhoinf_ > 1.0) )
      dserror("rho_inf out of range [0.0,1.0]");
    if ( (beta_!=0.25) or (gamma_!=0.5) or (alpham_!=0.5) or (alphaf_!=0.5) )
      dserror("you may only specify RHO_INF or the other four parameters");
    alpham_ = (2.0*rhoinf_-1.0)/(rhoinf_+1.0);
    alphaf_ = rhoinf_/(rhoinf_+1.0);
    beta_   = 0.25*(1.0-alpham_+alphaf_)*(1.0-alpham_+alphaf_);
    gamma_  = 0.5-alpham_+alphaf_;
  }

  // sanity checks and some screen output
  if (GlobalState().GetMyRank()==0)
  {
    if (rhoinf_>0.0)
      std::cout << "   rho = " << rhoinf_ << std::endl;
    // beta
    if ( (beta_ <= 0.0) or (beta_ > 0.5) )
      dserror("beta out of range (0.0,0.5]");
    else
      std::cout << "   beta = " << beta_ << std::endl;
    // gamma
    if ( (gamma_ <= 0.0) or (gamma_ > 1.0) )
      dserror("gamma out of range (0.0,1.0]");
    else
      std::cout << "   gamma = " << gamma_ << std::endl;
    // alpha_f
    if ( (alphaf_ < 0.0) or (alphaf_ >= 1.0) )
      dserror("alpha_f out of range [0.0,1.0)");
    else
      std::cout << "   alpha_f = " << alphaf_ << std::endl;
    // alpha_m
    if ( (alpham_ < -1.0) or (alpham_ >= 1.0) )
      dserror("alpha_m out of range [-1.0,1.0)");
    else
      std::cout << "   alpha_m = " << alpham_ << std::endl;

    /* ------ mid-averaging type -----------------------------------------------
     * In principle, there exist two mid-averaging possibilities, TR-like and
     * IMR-like, where TR-like means trapezoidal rule and IMR-like means implicit
     * mid-point rule. We used to maintain implementations of both variants, but
     * due to its significantly higher complexity, the IMR-like version has been
     * deleted (popp 02/2013). The nice thing about TR-like mid-averaging is that
     * all element (and thus also material) calls are exclusively(!) carried out
     * at the end-point t_{n+1} of each time interval, but never explicitly at
     * some generalized midpoint, such as t_{n+1-\alpha_f}. Thus, any cumbersome
     * extrapolation of history variables, etc. becomes obsolete. */
    const enum INPAR::STR::MidAverageEnum& midavg =
        genalpha_sdyn.GetMidAverageType();
    if (midavg != INPAR::STR::midavg_trlike)
      dserror("mid-averaging of internal forces only implemented TR-like");
    else
      std::cout << "   midavg = " << INPAR::STR::MidAverageString(midavg)<<std::endl;
  }

  // ---------------------------------------------------------------------------
  // setup mid-point vectors
  // ---------------------------------------------------------------------------
  const_vel_acc_update_ptr_ =
      Teuchos::rcp(new Epetra_MultiVector(*GlobalState().DofRowMapView(),2,true));

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = GlobalState().GetMutableFinertialN();
  finertianp_ptr_ = GlobalState().GetMutableFinertialNp();

  fviscon_ptr_ = GlobalState().GetMutableFviscoN();
  fvisconp_ptr_ = GlobalState().GetMutableFviscoNp();

  // Has to be set before the PostSetup() routine is called!
  issetup_ = true;

  PostSetup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::PostSetup()
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // check for applicability of classical GenAlpha scheme
  // ---------------------------------------------------------------------------
  // set the constant parameters for the element evaluation
  if (TimInt().GetDataSDyn().GetMassLinType()==INPAR::STR::ml_rotations)
  {
    dserror("MASSLIN=ml_rotations is not supported by classical GenAlpha! "
        "Choose GenAlphaLieGroup instead!");
  }

  EquilibrateInitialState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();

  if (IsPredictorState())
    return;

  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp_ptr = GlobalState().ExtractDisplEntries(x);
  GlobalState().GetMutableDisNp()->Scale(1.0,*disnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableVelNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(0),
      gamma_/(beta_*dt),*disnp_ptr,0.0);

  // ---------------------------------------------------------------------------
  // new end-point accelerations
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableAccNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(1),
      1.0/(beta_*dt*dt),*disnp_ptr,0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::UpdateConstantStateContributions()
{
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // copy the last converged state into the epetra_multivector
  // ---------------------------------------------------------------------------
  // create a copy of the last converged displacement solution and apply the
  // current DBC values --> D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp0_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().GetDisN()));
  // ---------------------------------------------------------------------------
  // (1) constant velocity update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Scale((beta_-gamma_)/beta_,
      *GlobalState().GetVelN());
  (*const_vel_acc_update_ptr_)(0)->Update((2.0*beta_-gamma_)*dt/(2.0*beta_),
      *GlobalState().GetAccN(),1.0);
  // ---------------------------------------------------------------------------
  // (2) constant acceleration update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Scale((2.0*beta_-1.0)/(2.0*beta_),
      *GlobalState().GetAccN());
  (*const_vel_acc_update_ptr_)(1)->Update(-1.0/(beta_*dt),
      *GlobalState().GetVelN(),1.0);
  // ---------------------------------------------------------------------------
  // set the current DBC values at the DBC-DoFs of the constant update vectors
  // ---------------------------------------------------------------------------
  Dbc().ApplyDirichletBC(GlobalState().GetTimeNp(),
      disnp0_ptr,                                          // displ.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(0),false),    // veloc.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(1),false),    // accel.
      false);
  /* Finally add the D_{n+1}^{0} vector to the constant velocity and
   * acceleration update vectors: */
  // ---------------------------------------------------------------------------
  // (1) velocity: V_{n+1}^{0} - gamma/(beta * dt) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Update(-gamma_/(beta_*dt),*disnp0_ptr,1.0);
  // ---------------------------------------------------------------------------
  // (2) acceleration: A_{n+1}^{0} - 1.0/(beta * dt^{2}) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Update(-1.0/(beta_*dt*dt),*disnp0_ptr,1.0);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f)
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  const bool ok = ModelEval().ApplyForce(x,f,1.0-GetIntParam());

  if (not ok)
    return ok;

  // ---------------------------------------------------------------------------
  // add the visco and mass contributions
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  const bool ok = ModelEval().ApplyStiff(x,jac,1.0-GetIntParam());

  if (not ok)
    return ok;

  // ---------------------------------------------------------------------------
  // add the visco and mass contributions
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(jac);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  const bool ok = ModelEval().ApplyForceStiff(x,f,jac,1.0-GetIntParam());

  if (not ok)
    return ok;

  // ---------------------------------------------------------------------------
  // add the visco and mass contributions
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(f);
  AddViscoMassContributions(jac);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::AssembleForce( Epetra_Vector& f,
    const std::vector<INPAR::STR::ModelType>* without_these_models ) const
{
  CheckInitSetup();

  // set the time step dependent parameters for the assembly
  const bool ok = ModelEval().AssembleForce( 1.0-GetIntParam(), f,
      without_these_models );

  if (not ok)
    return ok;

  // ---------------------------------------------------------------------------
  // add the visco and mass contributions
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::AddViscoMassContributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1-alpha_f}
  LINALG::AssembleMyVector(1.0,f,alphaf_,*fviscon_ptr_);
  LINALG::AssembleMyVector(1.0,f,1-alphaf_,*fvisconp_ptr_);
  // inertial forces at t_{n+1-alpha_m}
  LINALG::AssembleMyVector(1.0,f,1-alpham_,*finertianp_ptr_);
  LINALG::AssembleMyVector(1.0,f,alpham_,*finertian_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::AddViscoMassContributions(
    LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr =
      GlobalState().ExtractDisplBlock(jac);
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // add inertial contributions and scale the structural stiffness block
  stiff_ptr->Add(*GlobalState().GetMassMatrix(),false,
      (1.0-alpham_)/(beta_*dt*dt),1.0);
  // add Rayleigh damping contributions
  if (TimInt().GetDataSDyn().GetDampingType()==INPAR::STR::damp_rayleigh)
    stiff_ptr->Add(*GlobalState().GetDampMatrix(),false,(1.0-alphaf_)*gamma_/(beta_*dt),1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();
  // write dynamic forces
  iowriter.WriteVector("finert",finertian_ptr_);
  iowriter.WriteVector("fvisco",fviscon_ptr_);

  ModelEval().WriteRestart(iowriter,forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  ioreader.ReadVector(finertian_ptr_,"finert");
  ioreader.ReadVector(fviscon_ptr_,"fvisco");

  ModelEval().ReadRestart(ioreader);
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::GenAlpha::CalcRefNormForce(
    const enum NOX::Abstract::Vector::NormType& type)
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::GenAlpha::GetIntParam() const
{
  CheckInitSetup();
  return alphaf_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::UpdateStepState()
{
  CheckInitSetup();
  // ---------------------------------------------------------------------------
  // dynamic effects
  // ---------------------------------------------------------------------------
  // new at t_{n+1} -> t_n
  //    finertial_{n} := finertial_{n+1}
  finertian_ptr_->Scale(1.0,*GlobalState().GetFinertialNp());
  // new at t_{n+1} -> t_n
  //    fviscous_{n} := fviscous_{n+1}
  fviscon_ptr_->Scale(1.0,*fvisconp_ptr_);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState(alphaf_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::UpdateStepElement()
{
  CheckInitSetup();
  ModelEval().UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::PostUpdate()
{
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::PredictConstDisConsistVelAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // constant predictor: displacement in domain
  disnp.Scale(1.0,*disn);

  // consistent velocities following Newmark formulas
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  velnp.Update((beta_-gamma_)/beta_, *veln,
      (2.0*beta_-gamma_)*dt/(2.0*beta_), *accn,
      0.0);

  // consistent accelerations following Newmark formulas
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  accnp.Update(-1.0/(beta_*dt),*veln,
      (2.0*beta_-1.0)/(2.0*beta_),*accn,
      0.0);

  return;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::PredictConstVelConsistAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  /* In the general dynamic case there is no need to design a special start-up
   * procedure, since it is possible to prescribe an initial velocity or
   * acceleration. The corresponding accelerations are calculated in the
   * EquilibrateInitialState() routine. */

  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // extrapolated displacements based upon constant velocities
  // d_{n+1} = d_{n} + dt * v_{n}
  disnp.Update(1.0,*disn,dt,*veln,0.0);

  // consistent velocities following Newmark formulas
  velnp.Update(1.0,disnp,-1.0,*disn,0.0);
  velnp.Update((beta_-gamma_)/beta_,*veln,
      (2.*beta_-gamma_)*dt/(2.*beta_),*accn,
      gamma_/(beta_*dt));

  // consistent accelerations following Newmark formulas
  accnp.Update(1.0,disnp,-1.0,*disn,0.0);
  accnp.Update(-1.0/(beta_*dt),*veln,
      (2.0*beta_-1.0)/(2.0*beta_),*accn,
      1./(beta_*dt*dt));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::PredictConstAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  /* In the general dynamic case there is no need to design a special start-up
   * procedure, since it is possible to prescribe an initial velocity or
   * acceleration. The corresponding accelerations are calculated in the
   * EquilibrateInitialState() routine. */

  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // extrapolated displacements based upon constant accelerations
  // d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n}
  disnp.Update(1.0,*disn,dt,*veln, 0.0);
  disnp.Update(0.5*dt*dt,*accn,1.0);

  // extrapolated velocities (equal to consistent velocities)
  // v_{n+1} = v_{n} + dt * a_{n}
  velnp.Update(1.0,*veln,dt,*accn,0.0);

  // constant accelerations (equal to consistent accelerations)
  accnp.Update(1.0,*accn, 0.0);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::ResetEvalParams()
{
  // call base class
  STR::IMPLICIT::Generic::ResetEvalParams();

  // set the time step dependent parameters for the element evaluation
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  double timeintfac_dis = beta_*dt*dt;
  double timeintfac_vel = gamma_*dt;

  EvalData().SetTimIntFactorDisp(timeintfac_dis);
  EvalData().SetTimIntFactorVel(timeintfac_vel);
}
