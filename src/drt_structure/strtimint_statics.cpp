/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.cpp
\brief Statics analysis
\level 1
\maintainer Alexander Popp
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_statics.H"
#include "stru_aux.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_crack/crackUtils.H"
#include "../drt_plastic_ssn/plastic_ssn_manager.H"


/*======================================================================*/
/* constructor */
STR::TimIntStatics::TimIntStatics
(
  const Teuchos::ParameterList& timeparams,
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimIntImpl
  (
    timeparams,
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  ),
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntStatics::Init
(
    const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver
)
{
  // call Init() in base class
  STR::TimIntImpl::Init(timeparams,sdynparams,xparams,actdis,solver);

  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdynparams,"PRESTRESS");
  INPAR::STR::DynamicType dyntype = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP");

  if ( pstype != INPAR::STR::prestress_none && dyntype != INPAR::STR::dyna_statics )
  {
    dserror("Paranoia Error: PRESTRESS is only allowed in combinations with DYNAMICTYPE Statics!!");
  }

  // info to user
  if (myrank_ == 0 && bool (printscreen_))
  {
    // check if we are in prestressing mode
    if (pstype == INPAR::STR::prestress_mulf) IO::cout << "with static MULF prestress" << IO::endl;
    else if (pstype == INPAR::STR::prestress_id) IO::cout << "with static INVERSE DESIGN prestress" << IO::endl;
    else IO::cout << "with statics" << IO::endl;
  }

  // have a nice day
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntStatics::Setup()
{
  // call Setup() in base class
  STR::TimIntImpl::Setup();

  // create force vectors

  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*DofRowMapView(), true);

  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*DofRowMapView(), true);

  // internal force vector F_{int;n} at new time
  fint_ = LINALG::CreateVector(*DofRowMapView(), true);

  // external force vector F_{n} at new time
  fext_ = LINALG::CreateVector(*DofRowMapView(), true);

  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void STR::TimIntStatics::PredictConstDisConsistVelAcc()
{
  // constant predictor : displacement in domain
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  // new end-point velocities, these stay zero in static calculation
  veln_->PutScalar(0.0);

  // new end-point accelerations, these stay zero in static calculation
  accn_->PutScalar(0.0);

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* linear extrapolation of displacement field */
void STR::TimIntStatics::PredictConstVelConsistAcc()
{
  // for the first step we don't have any history to do
  // an extrapolation. Hence, we do TangDis
  if (step_==0)
  {
    PredictTangDisConsistVelAcc();
    return;
  }
  else
  {
    // Displacement increment over last time step
    Teuchos::RCP<Epetra_Vector> disp_inc = LINALG::CreateVector(*DofRowMapView(), true);
    disp_inc->Update((*dt_)[0],*(*vel_)(0),0.);
    LINALG::ApplyDirichlettoSystem(disp_inc,zeros_,*(dbcmaps_->CondMap()));
    disn_->Update(1.0, *(*dis_)(0), 0.0);
    disn_->Update(1.,*disp_inc,1.);
    veln_->Update(1.0, *(*vel_)(0), 0.0);
    accn_->Update(1.0, *(*acc_)(0), 0.0);
    return;
  }
  return;
}

/*----------------------------------------------------------------------*/
/* quadratic extrapolation of displacement field */
void STR::TimIntStatics::PredictConstAcc()
{
  // for the first step we don't have any history to do
  // an extrapolation. Hence, we do TangDis
  if (step_==0)
  {
    PredictTangDisConsistVelAcc();
    return;
  }
  else if (step_==1)
  {
    PredictConstVelConsistAcc();
    return;
  }
  else
  {
    // Displacement increment over last time step
    Teuchos::RCP<Epetra_Vector> disp_inc = LINALG::CreateVector(*DofRowMapView(), true);
    disp_inc->Update((*dt_)[0],*(*vel_)(0),0.);
    disp_inc->Update(.5*(*dt_)[0]*(*dt_)[0],*(*acc_)(0),1.);
    LINALG::ApplyDirichlettoSystem(disp_inc,zeros_,*(dbcmaps_->CondMap()));
    disn_->Update(1.0, *(*dis_)(0), 0.0);
    disn_->Update(1.,*disp_inc,1.);
    veln_->Update(1.0, *(*vel_)(0), 0.0);
    accn_->Update(1.0, *(*acc_)(0), 0.0);
    return;
  }
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its stiffness, ie derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void STR::TimIntStatics::EvaluateForceStiffResidual(Teuchos::ParameterList& params)
{
  // get info about prediction step from parameter list
  bool predict = false;
  if(params.isParameter("predict"))
    predict = params.get<bool>("predict");

  // initialize stiffness matrix to zero
  stiff_->Zero();

  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceStiffExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // ************************** (2) INTERNAL FORCES ***************************

  // initialize internal forces
  fintn_->PutScalar(0.0);

  if (HaveSemiSmoothPlasticity())
    if (plastman_->TSI())
    {
      plastman_->SetData().scale_timint_=1.;
      plastman_->SetData().dt_=(*dt_)[0];

      // pseudo-velocity for quasi-static tsi with Gough-Joule Effect
      veln_->Update(1./Dt(), *Dispnp(), -1./Dt(),*Dispn(),0.);
    }

  // ordinary internal force and stiffness
  ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_, params);

  // apply forces and stiffness due to constraints
  Teuchos::ParameterList pcon; //apply empty parameterlist, no scaling necessary
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fintn_, stiff_, pcon);

  // add forces and stiffness due to Cardiovascular0D bcs
  Teuchos::ParameterList pwindk;
  pwindk.set("time_step_size", (*dt_)[0]);
  ApplyForceStiffCardiovascular0D(timen_, disn_, fintn_, stiff_, pwindk);

  // add forces and stiffness due to spring dashpot condition
  Teuchos::ParameterList psprdash;
  ApplyForceStiffSpringDashpot(stiff_,fintn_,disn_,veln_,predict,psprdash);

  // ************************** (3) INERTIAL FORCES ***************************
  // This is statics, so there are no intertial forces.

  // ************************** (4) DAMPING FORCES ****************************
  // This is statics, so there are no viscous damping forces.

  // ******************** Finally, put everything together ********************

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update(1.0, *fintn_, 1.0);

  // build pure structural residual (only LS with EAS)
  if (fresn_str_!=Teuchos::null)
  {
    fresn_str_->Update(1.,*fintn_str_,0.);
    fresn_str_->Update(-1.,*fextn_,1.);
    LINALG::ApplyDirichlettoSystem(fresn_str_,zeros_,*(dbcmaps_->CondMap()));
  }

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = K_{T}
  // i.e. do nothing here

  // apply forces and stiffness due to beam contact
  ApplyForceStiffBeamContact(stiff_,fres_,disn_,predict);

  // apply forces and stiffness due to contact / meshtying
  ApplyForceStiffContactMeshtying(stiff_,fres_,disn_,predict);

  // close stiffness matrix
  stiff_->Complete();

  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate/define the residual force vector #fres_ for
 * relaxation solution with SolveRelaxationLinear */
void STR::TimIntStatics::EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params)
{
  // compute residual forces #fres_ and stiffness #stiff_
  EvaluateForceStiffResidual(params);

  // overwrite the residual forces #fres_ with interface load
  fres_->Update(-1.0, *fifc_, 0.0);
}

/*----------------------------------------------------------------------*/
/* Evaluate residual */
void STR::TimIntStatics::EvaluateForceResidual()
{

  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // ************************** (2) INTERNAL FORCES ***************************

  // initialize internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and stiffness
  ApplyForceInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_);

  // ************************** (3) INERTIAL FORCES ***************************
  // This is statics, so there are no inertial forces.

  // ************************** (4) DAMPING FORCES ****************************
  // This is statics, so there are no viscous damping forces.

  // ******************** Finally, put everything together ********************

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update(1.0, *fintn_, 1.0);

  // build pure structural residual (only LS with EAS)
  if (fresn_str_!=Teuchos::null)
  {
    fresn_str_->Update(1.,*fintn_str_,0.);
    fresn_str_->Update(-1.,*fextn_,1.);
    LINALG::ApplyDirichlettoSystem(fresn_str_,zeros_,*(dbcmaps_->CondMap()));
  }

  return;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for forces
 * originally by lw */
double STR::TimIntStatics::CalcRefNormForce()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = STR::AUX::CalculateVectorNorm(iternorm_, fintn_);

  // norm of the external forces
  double fextnorm = 0.0;
  fextnorm = STR::AUX::CalculateVectorNorm(iternorm_, fextn_);

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // return char norm
  return std::max(fintnorm, std::max(fextnorm, freactnorm));
}

/*----------------------------------------------------------------------*/
/* incremental iteration update of state */
void STR::TimIntStatics::UpdateIterIncrementally()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void STR::TimIntStatics::UpdateIterIteratively()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void STR::TimIntStatics::UpdateStepState()
{
  // calculate pseudo velocity and acceleration for predictor before updates
  if (pred_==INPAR::STR::pred_constvel || pred_==INPAR::STR::pred_constacc)
    veln_->Update(1./(*(*dt_)(0)),*disn_,-1./(*(*dt_)(0)),*(*dis_)(0),0.);
  if (pred_==INPAR::STR::pred_constacc)
    accn_->Update(1./(*(*dt_)(0)),*veln_,-1./(*(*dt_)(0)),*(*vel_)(0),0.);

  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->UpdateSteps(*disn_);

  //new material displacements
  if( (dismatn_!=Teuchos::null))
    dismat_->UpdateSteps(*dismatn_);

  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->UpdateSteps(*accn_);

  // update surface stress
  UpdateStepSurfstress();

  // update constraints
  UpdateStepConstraint();

  // update Cardiovascular0D
  UpdateStepCardiovascular0D();

  // update constraints
  UpdateStepSpringDashpot();

  // update contact / meshtying
  UpdateStepContactMeshtying();

  // update beam contact
  UpdateStepBeamContact();

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // look out
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntStatics::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  //p.set("alpha f", theta_);
  // action for elements
  p.set("action", "calc_struct_update_istep");
  // go to elements
  discret_->ClearState();
  discret_->SetState("displacement",(*dis_)(0));

  // Set material displacement state for ale-wear formulation
  if( (dismat_!=Teuchos::null))
    discret_->SetState("material_displacement",(*dismat_)(0));

  discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null,Teuchos::null);
  discret_->ClearState();
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntStatics::ReadRestartForce()
{
  return;
}

/*----------------------------------------------------------------------*/
/* write internal and external forces for restart */
void STR::TimIntStatics::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  output->WriteVector("fexternal", fext_);
  output->WriteVector("fint",fint_);

  //This restart output is needed in case of a static pre-simulation has to be restartet with dynamic time integration
  Teuchos::RCP<Epetra_Vector> finert= LINALG::CreateVector(*DofRowMapView(), true);  //!< inertia force at \f$t_{n}\f$
  finert->PutScalar(0.0);
  output->WriteVector("finert",finert);
  return;
}

/*-----------------------------------------------------------------------------*
 * Update all field vectors defined specific for this method,     sudhakar 12/13
 * to take into account of the new nodes introduced by crack propagation
 *----------------------------------------------------------------------------*/
void STR::TimIntStatics::updateMethodSpecificEpetraCrack( std::map<int,int>& oldnew )
{
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fintn_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fextn_, oldnew );
}

/*---------------------------------------------------------------*/
/* Apply Dirichlet boundary conditions on provided state vectors */
void STR::TimIntStatics::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc,
  bool recreatemap
)
{
  //call base ApplyDirichletBC
  STR::TimInt::ApplyDirichletBC(time,dis,vel,acc,recreatemap);

  //statics: set velocities and accelerations to zero
  if( vel!=Teuchos::null ) vel->PutScalar(0.0);
  if( acc!=Teuchos::null ) acc->PutScalar(0.0);

  return;
}
