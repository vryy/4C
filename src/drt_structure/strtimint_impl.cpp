/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.cpp
\brief Implicit time integration for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>

#include "strtimint.H"
#include "strtimint_impl.H"
#include "stru_aux.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/meshtying_abstract_strategy.H"  // needed in CmtLinearSolve (for feeding the contact solver with latest information about the contact status)
#include "../drt_contact/contact_abstract_strategy.H"  // needed in CmtLinearSolve (for feeding the contact solver with latest information about the contact status)
#include "../drt_inpar/inpar_contact.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_defines.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_potential/drt_potential_manager.H"
#include "../drt_lib/drt_locsys.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_krylov_projector.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_patspec/patspec.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"
#include "../drt_plastic_ssn/plastic_ssn_manager.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntImpl::TimIntImpl
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  ),
  pred_(DRT::INPUT::IntegralValue<INPAR::STR::PredEnum>(sdynparams,"PREDICT")),
  itertype_(DRT::INPUT::IntegralValue<INPAR::STR::NonlinSolTech>(sdynparams,"NLNSOL")),
  normtypedisi_(DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_DISP")),
  normtypefres_(DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_RESF")),
  normtypepres_(DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_PRES")),
  normtypepfres_(DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_INCO")),
  combdispre_(DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_DISPPRES")),
  combfrespfres_(DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_RESFINCO")),
  combdisifres_(DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(sdynparams,"NORMCOMBI_RESFDISP")),
  iternorm_(DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(sdynparams,"ITERNORM")),
  itermax_(sdynparams.get<int>("MAXITER")),
  itermin_(sdynparams.get<int>("MINITER")),
  divcontype_ (DRT::INPUT::IntegralValue<INPAR::STR::DivContAct>(sdynparams,"DIVERCONT")),
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
  normcon_(0.0),
  normpfres_(0.0),
  normpres_(0.0),
  normcontconstr_(0.0),  // < norm of contact constraints (saddlepoint formulation)
  normlagr_(0.0),        // < norm of lagrange multiplier increment (saddlepoint formulation)
  disi_(Teuchos::null),
  fres_(Teuchos::null),
  freact_(Teuchos::null),
  updateprojection_(false),
  stcscale_(DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdynparams,"STC_SCALING")),
  stclayer_(sdynparams.get<int>("STC_LAYER")),
  ptcdt_(sdynparams.get<double>("PTCDT")),
  dti_(1.0/ptcdt_)
{

  // verify: if system has constraints implemented with Lagrange multipliers,
  // then Uzawa-type solver is used
  if ( conman_->HaveConstraintLagr())
  {
    if ( (itertype_ != INPAR::STR::soltech_newtonuzawalin)
         and (itertype_ != INPAR::STR::soltech_newtonuzawanonlin) )
      dserror("Chosen solution technique %s does not work constrained.",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
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

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap_, true);

  //prepare matrix for scaled thickness business of thin shell structures
  stcmat_=
    Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));

  stccompl_ =false;

  // setup NOX parameter lists
  if (itertype_ == INPAR::STR::soltech_noxnewtonlinesearch)
    NoxSetup();
  else if (itertype_ == INPAR::STR::soltech_noxgeneral)
    NoxSetup(xparams.sublist("NOX"));

  // setup tolerances and binary operators for convergence check of contact/meshtying problems
  // in saddlepoint formulation
  tolcontconstr_ = tolfres_;
  tollagr_ = toldisi_;
  combfrescontconstr_ = INPAR::STR::bop_and; // default values, avoid uninitialized variables
  combdisilagr_ = INPAR::STR::bop_and;
  normtypecontconstr_ = INPAR::STR::convnorm_abs; //DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_RESF"));
  normtypeplagrincr_ = INPAR::STR::convnorm_abs; //DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdynparams,"NORM_DISP"));

  if (HaveContactMeshtying())
  {
    // extract information from parameter lists
    tolcontconstr_      = cmtman_->GetStrategy().Params().get<double>("TOLCONTCONSTR");
    tollagr_            = cmtman_->GetStrategy().Params().get<double>("TOLLAGR");
    combfrescontconstr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(cmtman_->GetStrategy().Params(),"NORMCOMBI_RESFCONTCONSTR");
    combdisilagr_       = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(cmtman_->GetStrategy().Params(),"NORMCOMBI_DISPLAGR");
  }


  // -------------------------------------------------------------------
  // setup Krylov projection if necessary
  // -------------------------------------------------------------------
  //
  // sysmat might be singular, e.g. when solid is not fully supported
  // in this case, we need a basis vector for the nullspace/kernel

  // get condition "KrylovSpaceProjection" from discretization
  std::vector<DRT::Condition*> KSPcond;
  discret_->GetCondition("KrylovSpaceProjection",KSPcond);
  int numcond = KSPcond.size();
  int numsolid = 0;

  DRT::Condition* kspcond = NULL;
  // check if for solid Krylov projection is required
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
    if (*name == "solid")
    {
      numsolid++;
      kspcond = KSPcond[icond];
    }
  }

  if (numsolid == 1)
  {
    SetupKrylovSpaceProjection(kspcond);
    if (myrank_ == 0)
      std::cout << "\nSetup of KrylovSpaceProjection in solid field\n" << std::endl;
  }
  else  if (numsolid == 0)
  {
    projector_ = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for solid field");


  // done so far
  return;
}

/*----------------------------------------------------------------------*/
/* integrate step */
int STR::TimIntImpl::IntegrateStep()
{
  int error = 0;
  Predict();
  error = Solve();
  return error;
}
/*----------------------------------------------------------------------*/
/* predict solution */
void STR::TimIntImpl::Predict()
{

  // set iteration step to 0 (predictor)
  iter_ = 0;

  // choose predictor
  if ( (pred_ == INPAR::STR::pred_constdis)
       or (pred_ == INPAR::STR::pred_constdispres) )
  {
    PredictConstDisConsistVelAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if ( (pred_ == INPAR::STR::pred_constvel))
  {
    PredictConstVelConsistAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if ( (pred_ == INPAR::STR::pred_constacc))
  {
    PredictConstAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if ( (pred_ == INPAR::STR::pred_constdisvelacc)
            or (pred_ == INPAR::STR::pred_constdisvelaccpres) )
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

  // zerofy pressure DOFs and time-derivatives
  if (pressure_ != Teuchos::null)
  {
    if ( (pred_ != INPAR::STR::pred_constdispres)
         and (pred_ != INPAR::STR::pred_constdisvelaccpres) )
    {
      pressure_->InsertCondVector(pressure_->ExtractCondVector(zeros_), disn_);
    }
    pressure_->InsertCondVector(pressure_->ExtractCondVector(zeros_), veln_);
    pressure_->InsertCondVector(pressure_->ExtractCondVector(zeros_), accn_);
  }

	// Update locals systems (which may be time dependent)
	if (locsysman_ != Teuchos::null)
  	locsysman_->Setup(timen_);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, accn_, false);

  // possibly initialise Lagrange multipliers to zero
  //  if ( (conman_->HaveConstraint())
  //       and (itertype_ == soltech_uzawalinnewton) )
  //  {
  //    conman_->ScaleLagrMult(0.0);
  //  }

  // compute residual forces fres_ and stiffness stiff_
  // (hand in boolean flag indicating that this a predictor)
  bool predict = true;
  EvaluateForceStiffResidual(predict);

  // rotate to local coordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global coordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // rotate back to global coordinate system
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
  // we set the minimum of CalcRefNormForce() and #tolfres_, because
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
/* prepare partition step */
void STR::TimIntImpl::PreparePartitionStep()
{

  // set iteration step to 0
  iter_ = 0;

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, accn_, false);

  // compute residual forces fres_ and stiffness stiff_
  // (hand in boolean flag indicating that this a predictor)
  bool predict = true;
  EvaluateForceStiffResidual(predict);

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
  // (hand in boolean flag indicating that this a predictor)
  bool predict = true;
  EvaluateForceStiffResidual(predict);

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
    locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);
  stiff_->Complete();
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

  // solve for disi_
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  if (HaveContactMeshtying())
    CmtLinearSolve();  // use contact/meshtying solver
  else
    solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, true);

  // decide which norms have to be evaluated
  bool bPressure = pressure_ != Teuchos::null;
  bool bContactSP = (HaveContactMeshtying() &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed);

  if( bPressure && bContactSP) dserror("We only support either contact/meshtying in saddlepoint formulation or structure with pressure DOFs");
  if( bPressure == false && bContactSP == false)
  {
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
  }
  if (bPressure)
  {
    Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(disi_);
    Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(disi_);
    normpres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);
  }
  if (bContactSP )
  {
    // extract subvectors
    Teuchos::RCP<Epetra_Vector> lagrincr  = cmtman_->GetStrategy().LagrMultSolveIncr();

    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    // build lagrange multiplier increment norm
    if(lagrincr!=Teuchos::null) normlagr_ = STR::AUX::CalculateVectorNorm(iternorm_, lagrincr);
    else normlagr_ = -1.0;
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

/*--------------------------------------------------------------------------*
 | setup Krylov projector including first fill                    nis Feb13 |
 *--------------------------------------------------------------------------*/
void STR::TimIntImpl::SetupKrylovSpaceProjection(DRT::Condition* kspcond)
{
  // get number of mode flags in dat-file
  const int nummodes = kspcond->GetInt("NUMMODES");

  // get rigid body mode flags - number and order as in ComputeNullspace
  // e.g. for a 3-D solid: [transx transy transz rotx roty rotz]
  const std::vector<int>* modeflags = kspcond->Get<std::vector<int> >("ONOFF");

  // get actual active mode ids given in dat-file
  std::vector<int> activemodeids;
  for(int rr=0;rr<nummodes;++rr)
  {
    if(((*modeflags)[rr])!=0)
    {
      activemodeids.push_back(rr);
    }
  }

  // get from dat-file definition how weights are to be computed
  const std::string* weighttype = kspcond->Get<std::string>("weight vector definition");

  // since we only use total Lagrange, no update necessary.
  updateprojection_ = false;

  // create the projector
  projector_ = Teuchos::rcp(new LINALG::KrylovProjector(activemodeids,weighttype,discret_->DofRowMap()));

  // update the projector
  UpdateKrylovSpaceProjection();
}

/*--------------------------------------------------------------------------*
 | update projection vectors w_ and c_ for Krylov projection      nis Feb13 |
 *--------------------------------------------------------------------------*/
void STR::TimIntImpl::UpdateKrylovSpaceProjection()
{
  const std::string* weighttype = projector_->WeightType();
  // only pointvalues are permissible for now - feel free to extend to integration!
  if(*weighttype == "integration")
  {
    dserror("option integration not implemented");
  }

  // get RCP to kernel vector of projector
  // since we are in 'pointvalue' mode, weights are changed implicitely
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  c->PutScalar(0.0);

  // We recompute the entire nullspace no matter what.
  // This is not nice yet since:
  // - translations are constant throughout the entire cmputation
  // - SAME nullspace is sometimes recomputed AGAIN for some iterative solvers
  // So here is space for optimization.

  // get number of modes and their ids
  std::vector<int> modeids = projector_->Modes();

  // RCP on vector of size 0 holding the nullspace data - resized within ComputeNullspace
  Teuchos::RCP<std::vector<double> > nullspace = Teuchos::rcp(new std::vector<double>(0));
  discret_->ComputeNullSpace(nullspace);

  //check if everything went fine
  if (nullspace->size()==0)
    dserror("nullspace not successfully computed");

  // pointer on first element of nullspace data
  double* nsdata = nullspace->data();

  // sort vector of nullspace data into kernel vector c_
  for ( size_t i=0; i < Teuchos::as<size_t>(modeids.size()); i++)
  {
    Epetra_Vector* ci = (*c)(i);
    const size_t myLength = ci->MyLength();
    for(size_t j=0; j<myLength; j++)
    {
      (*ci)[j] = nsdata[modeids[i]*myLength+j];
    }
  }

  // fillcomplete the projector to compute (w^T c)^(-1)
  projector_->FillComplete();
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ surface stresses and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::ApplyForceStiffSurfstress
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseOperator>& stiff
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
    surfstressman_->EvaluateSurfStress(p, disn, fint, stiff);
  }
  // bye bye
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate _certain_ potential forces and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::ApplyForceStiffPotential
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector>& fint,
  Teuchos::RCP<LINALG::SparseOperator>& stiff
)
{
  // potential force loads (but on internal force vector side)
  if (potman_ != Teuchos::null)
  {
    Teuchos::RCP<LINALG::SparseMatrix> mat = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff);
    ParameterList p; // create the parameters for manager
    p.set("pot_man", potman_);
    p.set("total time", time);
    potman_->EvaluatePotential(p, dis, fint, mat);
    stiff = mat;
  }

  // wooop
  return;
}


/*----------------------------------------------------------------------*/
/* TEST evaluate _certain_ potential forces and stiffness
 * evaluation happens internal-force like */
void STR::TimIntImpl::TestForceStiffPotential
(
  const double                        time,
  const Teuchos::RCP<Epetra_Vector>   dis,
  const int                           step
)
{
  // potential force loads (but on internal force vector side)
  if (potman_ != Teuchos::null)
  {
    if(potman_->ComputeAnalyticalSolution())
    {
      ParameterList p; // create the parameters for manager
      p.set("pot_man", potman_);
      p.set("total time", time);

      Teuchos::RCP<LINALG::SparseMatrix> stiff_test=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false, LINALG::SparseMatrix::FE_MATRIX));
      Teuchos::RCP<Epetra_Vector> fint_test=LINALG::CreateVector(*dofrowmap_, true);
      fint_test->PutScalar(0.0);
      stiff_test->Zero();

      potman_->TestEvaluatePotential(p, dis, fint_test, stiff_test, time, step);
    }
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
  Teuchos::RCP<LINALG::SparseOperator>& stiff,
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
/* evaluate forces and stiffness due to contact / meshtying */
void STR::TimIntImpl::ApplyForceStiffContactMeshtying
(
  Teuchos::RCP<LINALG::SparseOperator>& stiff,
  Teuchos::RCP<Epetra_Vector>& fresm,
  Teuchos::RCP<Epetra_Vector>& dis,
  bool predict
)
{
  if (HaveContactMeshtying())
  {
    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

    // contact / meshtying modifications need -fres
    fresm->Scale(-1.0);

    // make contact / meshtying modifications to lhs and rhs
    // (depending on whether this is a predictor step or not)
    cmtman_->GetStrategy().ApplyForceStiffCmt(dis,stiff,fresm,predict);

    // scaling back
    fresm->Scale(-1.0);

    // *********** time measurement ***********
    dtcmt_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // visualization of current Newton step
#ifdef MORTARGMSH2
    cmtman_->GetStrategy().VisualizeGmsh(stepn_,iter_);
#endif // #ifdef MORTARGMSH2
  }

  // wotcha
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate forces and stiffness due to beam contact */
void STR::TimIntImpl::ApplyForceStiffBeamContact
(
  Teuchos::RCP<LINALG::SparseOperator>& stiff,
  Teuchos::RCP<Epetra_Vector>& fresm,
  Teuchos::RCP<Epetra_Vector>& dis,
  bool predict
)
{
  if (HaveBeamContact())
  {
    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

    // contact / meshtying modifications need -fres
    fresm->Scale(-1.0);

    // make contact / meshtying modifications to lhs and rhs
    // (set boolean flag 'newsti' to true, which activates
    // sclaing of contact stiffness with appropriate scaling
    // factor, e.g. (1.0-alphaf), internally)
    beamcman_->Evaluate(*SystemMatrix(),*fresm,*dis,true);

    // scaling back
    fresm->Scale(-1.0);

    // *********** time measurement ***********
    dtcmt_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // visualization of current Newton step
#ifdef GMSHNEWTONSTEPS
    beamcman_->GmshOutput(*disn_,stepn_,iter_);
    beamcman_->ConsoleOutput();
#endif // #ifdef GMSHNEWTONSTEPS
  }

  // wotcha
  return;
}

void STR::TimIntImpl::ApplyForceStiffEmbedTissue
(
  Teuchos::RCP<LINALG::SparseOperator> stiff,
  Teuchos::RCP<Epetra_Vector> fint,
  Teuchos::RCP<Epetra_Vector> disn,
  bool predict
)
{
  if (pslist_ != Teuchos::null && pslist_->get<bool>("haveembedtissue") && !predict)
  {
    PATSPEC::CheckEmbeddingTissue(discret_,stiff,fint,disn);
  }

  return;
}

/*----------------------------------------------------------------------*/
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
  	convdis = normdisi_ < std::max(normchardis_*toldisi_,1e-15);
    break;
  case INPAR::STR::convnorm_mix:
  	convdis = ( (normdisi_ < toldisi_) or (normdisi_ < std::max(normchardis_*toldisi_,1e-15)) );
    break;
  default:
    dserror("Cannot check for convergence of residual displacements!");
    break;
  }

  // residual forces
  switch (normtypefres_)
  {
  case INPAR::STR::convnorm_abs:
    convfres = normfres_ < tolfres_;
    break;
  case INPAR::STR::convnorm_rel:
  	convfres = normfres_ < std::max(tolfres_*normcharforce_,1e-15);
    break;
  case INPAR::STR::convnorm_mix:
  	convfres = ( (normfres_ < tolfres_) or (normfres_ < std::max(tolfres_*normcharforce_,1e-15)) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
    break;
  }

  // check constraint
  bool cc = true;
  if (conman_->HaveConstraintLagr())
  {
    cc = normcon_ < tolcon_;
  }

  // check contact (active set)
  bool ccontact = true;
  if (HaveContactMeshtying())
  {
    // check which case (application, strategy) we are in
    INPAR::CONTACT::ApplicationType apptype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(cmtman_->GetStrategy().Params(),"APPLICATION");
    INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
    bool semismooth = DRT::INPUT::IntegralValue<int>(cmtman_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

    // only do this convergence check for semi-smooth Lagrange multiplier contact
    if (apptype == INPAR::CONTACT::app_mortarcontact && stype == INPAR::CONTACT::solution_lagmult && semismooth)
      ccontact = cmtman_->GetStrategy().ActiveSetSemiSmoothConverged();

    // add convergence check for saddlepoint formulations
    // use separate convergence checks for contact constraints and
    // LM increments
    INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM");
    if (stype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed) {
      bool convDispLagrIncr = false;
      bool convDispWIncr = false;

      switch ( normtypeplagrincr_ )
      {
      case INPAR::STR::convnorm_abs:
        convDispLagrIncr = normlagr_ < tollagr_;
        convDispWIncr = normw_ < 1e-14; // WEAR
        break;
      /*case INPAR::STR::convnorm_rel:
        convDispLagrIncr = normlagr_ < tollagr_;
        break;*/
      default:
        dserror("You should not turn up here.");
        break;
      }

      // switch between "and" and "or"
      if (combdisilagr_==INPAR::STR::bop_and)
        convdis = convdis and convDispLagrIncr and convDispWIncr;
      else if (combdisilagr_==INPAR::STR::bop_or)
        convdis = convdis or convDispLagrIncr;
      else
        dserror("Something went terribly wrong with binary operator!");

      bool convResfContConstr = false;

      switch ( normtypecontconstr_ )
      {
      case INPAR::STR::convnorm_abs:
        convResfContConstr = normcontconstr_ < tolcontconstr_;
        break;
      /*case INPAR::STR::convnorm_rel:
        //convDispLagrIncr = normcontconstr_ < std::max(tolfres_*normcharforce_,1e-15);
        convResfContConstr = normcontconstr_ < tolcontconstr_;
        break;*/
      default:
        dserror("You should not turn up here.");
        break;
      }

      // switch between "and" and "or"
      if (combfrescontconstr_==INPAR::STR::bop_and)
        convfres = convfres and convResfContConstr;
      else if (combfrescontconstr_==INPAR::STR::bop_or)
        convfres = convfres or convResfContConstr;
      else
        dserror("Something went terribly wrong with binary operator!");


    }

  }  // end HaveMeshtyingContact()

  // check convergence of plastic active set
  bool cplast=true;
  if (HaveSemiSmoothPlasticity())
  {
    cplast=plastman_->Converged();
    if (!cplast)
      if (discret_->Comm().MyPID()==0)
        std::cout << "Active plastic set has changed"<< std::endl;
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
      break;
    }

    // incompressible residual
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs:
      convfpre = normpfres_ < tolpfres_;
      break;
    default:
      dserror("Cannot check for convergence of incompressible force residuals!");
      break;
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
  return (conv and cc and ccontact and cplast);
}

/*----------------------------------------------------------------------*/
/* solve equilibrium */
int STR::TimIntImpl::Solve()
{
  int nonlin_error = 0 ;
  // special nonlinear iterations for contact / meshtying
  if (HaveContactMeshtying())
  {
    // choose solution technique in accordance with user's will
    nonlin_error = CmtNonlinearSolve();
  }

  // special nonlinear iterations for beam contact
  else if (HaveBeamContact())
  {
    // choose solution technique in accordance with user's will
    nonlin_error =  BeamContactNonlinearSolve();
  }

  // all other cases
  else
  {
    // choose solution technique in accordance with user's will
    switch (itertype_)
    {
    case INPAR::STR::soltech_newtonfull :
      nonlin_error = NewtonFull();
      break;
    case INPAR::STR::soltech_newtonuzawanonlin :
      nonlin_error = UzawaNonLinearNewtonFull();
      break;
    case INPAR::STR::soltech_newtonuzawalin :
      nonlin_error = UzawaLinearNewtonFull();
      break;
    case INPAR::STR::soltech_noxnewtonlinesearch :
    case INPAR::STR::soltech_noxgeneral :
      nonlin_error = NoxSolve();
      break;
    case INPAR::STR::soltech_ptc :
      nonlin_error = PTC();
      break;
    // catch problems
    default :
      dserror("Solution technique \"%s\" is not implemented",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
      break;
    }
  }
  return nonlin_error;
}
/*----------------------------------------------------------------------*/
/* solution with full Newton-Raphson iteration */
int STR::TimIntImpl::NewtonFull()
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
  timer_->ResetStartTime();

  int linsolve_error= 0;
  // equilibrium iteration loop
  while ( ( (not Converged() and (not linsolve_error)) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

    // STC preconditioning
    STCPreconditioning();

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

    // solve for disi_
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }

    // linear solver call (contact / meshtying case or default)
    if (HaveContactMeshtying())
      CmtLinearSolve();
    else
    {
      linsolve_error = solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, iter_==1, projector_);
      // check for problems in linear solver
      // however we only care about this if we have a fancy divcont action (meaning function will return 0 )
      linsolve_error = LinSolveErrorCheck(linsolve_error);
    }
    solver_->ResetTolerance();

    // recover standard displacements
    RecoverSTCSolution();

    // recover contact / meshtying Lagrange multipliers
    if (HaveContactMeshtying())
      cmtman_->GetStrategy().Recover(disi_);

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // update end-point displacements etc
    UpdateIter(iter_);

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

    // cancel in residual those forces that would excite rigid body modes and
    // that thus vanish in the Krylov space projection
    if (projector_!=Teuchos::null)
      projector_->ApplyPT(*fres_);

    // decide which norms have to be evaluated
    bool bPressure = pressure_ != Teuchos::null;
    bool bContactSP = (HaveContactMeshtying() &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed);

    if( bPressure && bContactSP) dserror("We only support either contact/meshtying in saddlepoint formulation or structure with pressure DOFs");
    if( bPressure == false && bContactSP == false)
    {
      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    }
    if (bPressure)
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
    if (bContactSP )
    {
      // extract subvectors
      Teuchos::RCP<Epetra_Vector> lagrincr  = cmtman_->GetStrategy().LagrMultSolveIncr();
      Teuchos::RCP<Epetra_Vector> constrrhs = cmtman_->GetStrategy().ConstrRhs();

      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
      // build residual constraint norm
      if(constrrhs!=Teuchos::null) normcontconstr_ = STR::AUX::CalculateVectorNorm(iternorm_, constrrhs);
      else normcontconstr_ = -1.0;

      // build lagrange multiplier increment norm
      if(lagrincr!=Teuchos::null) normlagr_ = STR::AUX::CalculateVectorNorm(iternorm_, lagrincr);
      else normlagr_ = -1.0;

      // for wear discretization
      INPAR::CONTACT::WearType wtype =
          DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(cmtman_->GetStrategy().Params(),"WEARTYPE");
      if(wtype==INPAR::CONTACT::wear_discr)
      {
        Teuchos::RCP<Epetra_Vector> wincr  = cmtman_->GetStrategy().WSolveIncr();
        Teuchos::RCP<Epetra_Vector> wearrhs = cmtman_->GetStrategy().WearRhs();

        if(wearrhs!=Teuchos::null) normwrhs_ = STR::AUX::CalculateVectorNorm(iternorm_, wearrhs);
        else normwrhs_ = -1.0;

        if(wincr!=Teuchos::null) normw_ = STR::AUX::CalculateVectorNorm(iternorm_, wincr);
        else normw_ = -1.0;
      }
      else
      {
        normw_=0.0;
        normwrhs_=0.0;
      }
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

  //do nonlinear solver error check
  return NewtonFullErrorCheck(linsolve_error);
}

int STR::TimIntImpl::NewtonFullErrorCheck(int linerror)
{
  // if everything is fine print to screen and return
  if (Converged())
  {
    if(myrank_ == 0)
      PrintNewtonConv();
    return 0;
  }
  // now some error checks
  // do we have a problem in the linear solver
  // only check if we want to do something fancy other wise we ignore the error in the linear solver
  if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation))
  {
    return linerror;
  }
  else
  {
    if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_stop ) )
    {
      // write restart output of last converged step before stopping
      OutputStep(true);

      dserror("Newton unconverged in %d iterations", iter_);
      return 1;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_continue ) )
    {
      if (myrank_ == 0)
        IO::cout<<"Newton unconverged in " << iter_ << " iterations, continuing" <<IO::endl;
      return 0;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation))
    {
      if (myrank_ == 0)
        IO::cout<< "Newton unconverged in " << iter_ << " iterations " << IO::endl;
      return 1;
    }
  }
  dserror("Fatal error in NonLinSolveErrorCheck, case not implemented ");
  return 0;
}


int STR::TimIntImpl::LinSolveErrorCheck(int linerror)
{
  // we only care about problems in the linear solver if we have a fancy divcont action
  if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation) )
  {
    if (myrank_ == 0)
    IO::cout<< "Linear solver is having trouble " << IO::endl;
    return linerror;
  }
  else
  {
    return 0;
  }

}

/*----------------------------------------------------------------------*/
/* do non-linear Uzawa iteration within a full NRI is called,
 * originally by tk */
int STR::TimIntImpl::UzawaNonLinearNewtonFull()
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
  int error = NewtonFull();
  if(error) return error;

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
    int error = NewtonFull();
    if(error) return error;


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
  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::UpdateStepConstraint()
{
  if (conman_ -> HaveConstraint())
    conman_->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::UpdateStepSurfstress()
{
  if (surfstressman_ -> HaveSurfStress())
    surfstressman_->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool STR::TimIntImpl::HaveConstraint()
{
  return conman_->HaveConstraintLagr();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::UpdateIterIncrConstr
(
  Teuchos::RCP<Epetra_Vector> lagrincr ///< Lagrange multiplier increment
)
{
  conman_->UpdateLagrMult(lagrincr);
}

/*----------------------------------------------------------------------*/
/* do linearised Uzawa iterations with full NRI
 * originally by tk 11/07 */
int STR::TimIntImpl::UzawaLinearNewtonFull()
{
  // allocate additional vectors and matrices
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
  timer_->ResetStartTime();

  // equilibrium iteration loop
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

//    // uncomplete stiffness matrix, so stuff can be inserted again
//    stiff_->UnComplete();

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

    // prepare residual Lagrange multiplier
    lagrincr->PutScalar(0.0);

    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

    //Use STC preconditioning on system matrix
    STCPreconditioning();

    // get constraint matrix with and without Dirichlet zeros
    Teuchos::RCP<LINALG::SparseMatrix> constr =
        (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(conman_->GetConstrMatrix()));
    Teuchos::RCP<LINALG::SparseMatrix> constrT =
        Teuchos::rcp(new LINALG::SparseMatrix (*constr));

    constr->ApplyDirichlet(*(dbcmaps_->CondMap()),false);

    // Apply STC on constraint matrices of desired
    if(stcscale_ != INPAR::STR::stc_none)
    {
      //std::cout<<"scaling constraint matrices"<<std::endl;
      constrT=LINALG::MLMultiply(*stcmat_,true,*constrT,false,false,false,true);
      if (stcscale_ == INPAR::STR::stc_currsym)
      {
        constr = LINALG::MLMultiply(*stcmat_,true,*constr,false,false,false,true);;
      }
    }
    // Call constraint solver to solve system with zeros on diagonal
    consolv_->Solve(SystemMatrix(), constr, constrT,
                    disi_, lagrincr,
                    fres_, conrhs);

    //recover unscaled solution
    RecoverSTCSolution();

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

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
	
	// no linear solver error check implemented here, always passing 0
  return UzawaLinearNewtonFullErrorCheck(0);
}

int STR::TimIntImpl::UzawaLinearNewtonFullErrorCheck(int linerror)
{
  // if everything is fine print to screen and return
  if (Converged())
  {
    	   // compute and print monitor values
    if (conman_->HaveMonitor())
    {
      conman_->ComputeMonitorValues(disn_);
    }

    // print newton message on proc 0
    if (myrank_ == 0)
      conman_->PrintMonitorValues();
      
    return 0;
  }
  // now some error checks
  // do we have a problem in the linear solver
  // only check if we want to do something fancy other wise we ignore the error in the linear solver
  if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation) )
  {
    return linerror;
  }
  else
  {
    if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_stop ) )
    {
      dserror("Newton unconverged in %d iterations", iter_);
      return 1;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_continue ) )
    {
      if (myrank_ == 0)
        IO::cout<<"Newton unconverged in " << iter_ << " iterations, continuing" <<IO::endl;
        if (conman_->HaveMonitor())
      		conman_->ComputeMonitorValues(disn_);
      return 0;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation))
    {
  		dserror("Fancy divcont actions not implemented forUzawaLinearNewtonFull ");
      return 1;
    }
  }
  dserror("Fatal error in UzawaLinearNewtonFullErrorCheck, case not implemented ");
  return 0;
}

/*----------------------------------------------------------------------*/
/* solution with nonlinear iteration for contact / meshtying */
int STR::TimIntImpl::CmtNonlinearSolve()
{
  //********************************************************************
  // get some parameters
  //********************************************************************
  // application type
  INPAR::CONTACT::ApplicationType apptype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(cmtman_->GetStrategy().Params(),"APPLICATION");

  // strategy type
  INPAR::CONTACT::SolvingStrategy soltype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");

  // semi-smooth Newton type
  bool semismooth = DRT::INPUT::IntegralValue<int>(cmtman_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

  // iteration type
  if (itertype_ != INPAR::STR::soltech_newtonfull)
    dserror("Unknown type of equilibrium iteration");

  //********************************************************************
  // Solving Strategy using Lagrangian Multipliers
  //********************************************************************
  if (soltype == INPAR::CONTACT::solution_lagmult)
  {
    //********************************************************************
    // 1) SEMI-SMOOTH NEWTON FOR CONTACT
    // The search for the correct active set (=contact nonlinearity) and
    // the large deformation linearization (=geometrical nonlinearity) are
    // merged into one semi-smooth Newton method and solved within ONE
    // iteration loop (which is then basically a standard Newton).
    //********************************************************************
    if (apptype == INPAR::CONTACT::app_mortarcontact && semismooth)
    {
      // nonlinear iteration
     int error = NewtonFull();
     if(error) return error;
    }

    //********************************************************************
    // 2) FIXED-POINT APPROACH FOR CONTACT
    // The search for the correct active set (=contact nonlinearity) is
    // represented by a fixed-point approach, whereas the large deformation
    // linearization (=geometrical nonlinearity) is treated by a standard
    // Newton scheme. This yields TWO nested iteration loops
    //********************************************************************
    else if (apptype == INPAR::CONTACT::app_mortarcontact && !semismooth)
    {
      // active set strategy
      int activeiter = 0;
      while (cmtman_->GetStrategy().ActiveSetConverged()==false)
      {
        // increase active set iteration index
        ++activeiter;

        // predictor step (except for first active set step)
        if (activeiter > 1) Predict();

        // nonlinear iteration
        int error = NewtonFull();
        if(error) return error;

        // update of active set (fixed-point)
        cmtman_->GetStrategy().UpdateActiveSet();
      }
    }

    //********************************************************************
    // 3) STANDARD NEWTON APPROACH FOR MESHTYING
    // No search for the correct active set has to be resolved for mortar
    // meshtying and mortar coupling is linear in this case. Thus, only
    // the large deformation FE problem remains to be solved as nonlinearity
    // Here, a standard Newton scheme is applied and we have ONLY ONE loop.
    //********************************************************************
    else
    {
      // nonlinear iteration
      int error = NewtonFull();
      if(error) return error;
    }
  }

  //********************************************************************
  // Solving Strategy using Regularization Techniques (Penalty Method)
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_penalty)
  {
    // nonlinear iteration
    int error = NewtonFull();
    if(error) return error;

    // update constraint norm
    cmtman_->GetStrategy().UpdateConstraintNorm();
  }

  //********************************************************************
  // Solving Strategy using Augmented Lagrange Techniques (with Uzawa)
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_auglag)
  {
    // get tolerance and maximum Uzawa steps
    double eps = cmtman_->GetStrategy().Params().get<double>("UZAWACONSTRTOL");
    int maxuzawaiter = cmtman_->GetStrategy().Params().get<int>("UZAWAMAXSTEPS");

    // Augmented Lagrangian loop (Uzawa)
    int uzawaiter=0;
    do
    {
      // increase iteration index
      ++uzawaiter;
      if (uzawaiter > maxuzawaiter) dserror("Uzawa unconverged in %d iterations",maxuzawaiter);
      if (!myrank_) std::cout << "Starting Uzawa step No. " << uzawaiter << std::endl;

      // for second, third,... Uzawa step: out-of-balance force
      if (uzawaiter>1)
      {
        fres_->Scale(-1.0);
        cmtman_->GetStrategy().InitializeUzawa(stiff_,fres_);
        fres_->Scale(-1.0);
      }

      // nonlinear iteration
      int error = NewtonFull();
      if(error) return error;

      // update constraint norm and penalty parameter
      cmtman_->GetStrategy().UpdateConstraintNorm(uzawaiter);

      // store Lagrange multipliers for next Uzawa step
      cmtman_->GetStrategy().UpdateAugmentedLagrange();
      cmtman_->GetStrategy().StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

    } while (cmtman_->GetStrategy().ConstraintNorm() >= eps);

    // reset penalty parameter
    cmtman_->GetStrategy().ResetPenalty();
  }

  return 0;
}

/*----------------------------------------------------------------------*/
/* linear solver call for contact / meshtying */
void STR::TimIntImpl::CmtLinearSolve()
{
  // adapt tolerance for contact solver
  // note: tolerance for fallback solver already adapted in NewtonFull
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normfres_;
    double wanted = tolfres_;
    contactsolver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

  // strategy and system setup types
  INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
  INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM");

  // update information about active slave dofs
  //**********************************************************************
  // feed solver/preconditioner with additional information about the contact/meshtying problem
  //**********************************************************************
  {
    // feed Aztec based solvers with contact information
    if (contactsolver_->Params().isSublist("Aztec Parameters"))
    {
      Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("Aztec Parameters");
      RCP<Epetra_Map> masterDofMap;
      RCP<Epetra_Map> slaveDofMap;
      RCP<Epetra_Map> innerDofMap;
      RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtman_->GetStrategy());
      strat->CollectMapsForPreconditioner(masterDofMap, slaveDofMap, innerDofMap, activeDofMap);
      Teuchos::ParameterList & linSystemProps = mueluParams.sublist("Linear System properties");
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact masterDofMap",masterDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact slaveDofMap",slaveDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact innerDofMap",innerDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact activeDofMap",activeDofMap);
      Teuchos::RCP<CONTACT::CoAbstractStrategy> costrat = Teuchos::rcp_dynamic_cast<CONTACT::CoAbstractStrategy>(strat);
      if (costrat != Teuchos::null) linSystemProps.set<std::string>("ProblemType", "contact");
      else                          linSystemProps.set<std::string>("ProblemType", "meshtying");
      linSystemProps.set<int>("time step",step_);
      linSystemProps.set<int>("iter",iter_);
    }
    // feed Belos based solvers with contact information
    if (contactsolver_->Params().isSublist("Belos Parameters"))
    {
      Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("Belos Parameters");
      RCP<Epetra_Map> masterDofMap;
      RCP<Epetra_Map> slaveDofMap;
      RCP<Epetra_Map> innerDofMap;
      RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtman_->GetStrategy());
      strat->CollectMapsForPreconditioner(masterDofMap, slaveDofMap, innerDofMap, activeDofMap);
      Teuchos::ParameterList & linSystemProps = mueluParams.sublist("Linear System properties");
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact masterDofMap",masterDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact slaveDofMap",slaveDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact innerDofMap",innerDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact activeDofMap",activeDofMap);
      Teuchos::RCP<CONTACT::CoAbstractStrategy> costrat = Teuchos::rcp_dynamic_cast<CONTACT::CoAbstractStrategy>(strat);
      if (costrat != Teuchos::null) linSystemProps.set<std::string>("ProblemType", "contact");
      else                          linSystemProps.set<std::string>("ProblemType", "meshtying");
      linSystemProps.set<int>("time step",step_);
      linSystemProps.set<int>("iter",iter_);
    }

  } // end: feed solver with contact/meshtying information

  // analysis of eigenvalues and condition number
#ifdef CONTACTEIG
    // global counter
    static int globindex = 0;
    ++globindex;

    // print to file in matlab format
    std::ostringstream filename;
    const std::string filebase = "sparsematrix";
    filename << "o/matlab_output/" << filebase << "_" << globindex << ".mtl";
    LINALG::PrintMatrixInMatlabFormat(filename.str().c_str(),*(SystemMatrix()->EpetraMatrix()));

    // print sparsity pattern to file
    LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
#endif // #ifdef CONTACTEIG

  //**********************************************************************
  // Solving a saddle point system
  // (1) Standard / Dual Lagrange multipliers -> SaddlePointCoupled
  // (2) Standard / Dual Lagrange multipliers -> SaddlePointSimpler
  //**********************************************************************
  if (soltype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed)
  {
    // (iter_-1 to be consistent with old time integration)
    cmtman_->GetStrategy().SaddlePointSolve(*contactsolver_,*solver_,stiff_,fres_,disi_,dirichtoggle_,iter_-1);
  }

  //**********************************************************************
  // Solving a purely displacement based system
  // (1) Dual (not Standard) Lagrange multipliers -> Condensed
  // (2) Penalty and Augmented Lagrange strategies
  //**********************************************************************
  else
  {
    // check if contact contributions are present,
    // if not we make a standard solver call to speed things up
    if (!cmtman_->GetStrategy().IsInContact() &&
        !cmtman_->GetStrategy().WasInContact() &&
        !cmtman_->GetStrategy().WasInContactLastTimeStep())
    {
      // standard solver call (fallback solver for pure structure problem)
      solver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
      return;
    }

    // solve with contact solver
    contactsolver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
    //cmtman_->GetStrategy().Solve(*contactsolver_,*solver_,stiff_,fres_,disi_,dirichtoggle_,iter_-1);
  }

  // reset tolerance for contact solver
  contactsolver_->ResetTolerance();

  return;
}

/*----------------------------------------------------------------------*/
/* solution with nonlinear iteration for beam contact */
int STR::TimIntImpl::BeamContactNonlinearSolve()
{
  //********************************************************************
  // get some parameters
  //********************************************************************
  // strategy type
  INPAR::CONTACT::SolvingStrategy soltype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(beamcman_->InputParameters(),"STRATEGY");

  // unknown types of nonlinear iteration schemes
  if (itertype_ != INPAR::STR::soltech_newtonfull)
    dserror("Unknown type of equilibrium iteration");

  //**********************************************************************
  // solving strategy using regularization with penalty method
  // (nonlinear solution approach: ordinary NEWTON)
  //**********************************************************************
  if (soltype == INPAR::CONTACT::solution_penalty)
  {
     // nonlinear iteration (Newton)
    int error = NewtonFull();
    if(error) return error;

    // update constraint norm
    beamcman_->UpdateConstrNorm();
  }
  //**********************************************************************

  //**********************************************************************
  // solving strategy using regularization with augmented Lagrange method
  // (nonlinear solution approach: nested UZAWA NEWTON)
  //**********************************************************************
  else if (soltype == INPAR::CONTACT::solution_auglag)
  {
    // get tolerance and maximum number of Uzawa steps from input file
    double eps = beamcman_->InputParameters().get<double>("UZAWACONSTRTOL");
    int maxuzawaiter = beamcman_->InputParameters().get<int>("UZAWAMAXSTEPS");

    // outer Augmented Lagrangian iteration (Uzawa)
    do
    {
      // increase iteration index by one
      beamcman_->UpdateUzawaIter();
      if (beamcman_->GetUzawaIter() > maxuzawaiter)
        dserror("Uzawa unconverged in %d iterations",maxuzawaiter);

      if (!myrank_)
        std::cout << std::endl << "Starting Uzawa step No. " << beamcman_->GetUzawaIter() << std::endl;

      // for second, third,... Uzawa step: out-of-balance force
      if (beamcman_->GetUzawaIter() > 1)
      {
        // beam contact modifications need -fres
        fres_->Scale(-1.0);

        // make contact modifications to lhs and rhs
        beamcman_->InitializeUzawa(*SystemMatrix(),*fres_,*disn_,true);

        // scaling back
        fres_->Scale(-1.0);
      }

      // inner nonlinear iteration (Newton)
      int error = NewtonFull();
      if(error) return error;
      // update constraint norm and penalty parameter
      beamcman_->UpdateConstrNorm();

      // update Uzawa Lagrange multipliers
      beamcman_->UpdateAlllmuzawa();

    } while (abs(beamcman_->GetConstrNorm()) >= eps);

    // reset penalty parameter and Uzawa index
    beamcman_->ResetCurrentpp();
    beamcman_->ResetUzawaIter();
  }

  //**********************************************************************
  // unknown solving strategy
  //**********************************************************************
  else
  {
    dserror("ERROR: Chosen strategy not yet available for beam contact");
  }

  // decide whether the original or the modified gapfunction definition is used
  bool newgapfunction =  DRT::INPUT::IntegralValue<int>(beamcman_->InputParameters(),"BEAMS_NEWGAP");

  // if the modified gap function definition is used the normal vector of the last time step has to be stored
  if (newgapfunction) {beamcman_->ShiftAllNormal();}

  return 0;
}

/*----------------------------------------------------------------------*/
/* solution with pseudo transient continuation */
int STR::TimIntImpl::PTC()
{
  // we do a PTC iteration here.
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
  timer_->ResetStartTime();

  double ptcdt     = ptcdt_;
  double nc; fres_->NormInf(&nc);
  double dti = 1/ptcdt;

  int linsolve_error= 0;
  // equilibrium iteration loop
  while ( ( ( not Converged() and (not linsolve_error) ) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

    // modify stiffness matrix with dti
    {
      RCP<Epetra_Vector> tmp = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      tmp->PutScalar(dti);
      RCP<Epetra_Vector> diag = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      SystemMatrix()->ExtractDiagonalCopy(*diag);
      diag->Update(1.0,*tmp,1.0);
      SystemMatrix()->ReplaceDiagonalValues(*diag);
    }

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

    // STC preconditioning
    STCPreconditioning();


    // solve for disi_
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    // linear solver call (contact / meshtying case or default)
    if (HaveContactMeshtying())
      CmtLinearSolve();
    else
    {
      linsolve_error = solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, iter_==1);
      // check for problems in linear solver
      // however we only care about this if we have a fancy divcont action  (meaning function will return 0 )
      linsolve_error=LinSolveErrorCheck(linsolve_error);
    }
    solver_->ResetTolerance();

    // recover standard displacements
    RecoverSTCSolution();

    // recover contact / meshtying Lagrange multipliers
    if (HaveContactMeshtying())
      cmtman_->GetStrategy().Recover(disi_);

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // update end-point displacements etc
    UpdateIter(iter_);

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

    // decide which norms have to be evaluated
    bool bPressure = pressure_ != Teuchos::null;
    bool bContactSP = (HaveContactMeshtying() &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed);

    if( bPressure && bContactSP) dserror("We only support either contact/meshtying in saddlepoint formulation or structure with pressure DOFs");
    if( bPressure == false && bContactSP == false)
    {
      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    }
    if (bPressure)
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
    if (bContactSP )
    {
      // extract subvectors
      Teuchos::RCP<Epetra_Vector> lagrincr  = cmtman_->GetStrategy().LagrMultSolveIncr();
      Teuchos::RCP<Epetra_Vector> constrrhs = cmtman_->GetStrategy().ConstrRhs();

      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
      // build residual constraint norm
      if(constrrhs!=Teuchos::null) normcontconstr_ = STR::AUX::CalculateVectorNorm(iternorm_, constrrhs);
      else normcontconstr_ = -1.0;
      // build lagrange multiplier increment norm
      if(lagrincr!=Teuchos::null) normlagr_ = STR::AUX::CalculateVectorNorm(iternorm_, lagrincr);
      else normlagr_ = -1.0;
    }

    // print stuff
    dti_ = dti;
    PrintNewtonIter();

    // update ptc
    {
      double np; fres_->NormInf(&np);
      dti *= (np/nc);
      dti = max(dti,0.0);
      nc = np;
    }
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

  //do nonlinear solver error check
  return NewtonFullErrorCheck(linsolve_error);

}


/*----------------------------------------------------------------------*/
/* Update iteration */
void STR::TimIntImpl::UpdateIter
(
  const int iter  //!< iteration counter
)
{
  // Doing UpdateIterIteratively() is not sufficient in the first Newton step
  // since the predictor might lead to velocities and accelerations that are
  // not consistently computed from the displacements based on the time
  // integration scheme.
  // Hence, in the first nonlinear iteration, we do UpdateIterIncrementally()
  // to ensure consistent velocities and accelerations across all predictors.
  //
  // From the second nonlinear iteration on, both update routines lead to
  // exactly the same results.
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

  // recover contact / meshtying Lagrange multipliers (monolithic FSI)
  // not in the case of TSI with contact
  if (DRT::Problem::Instance()->ProblemType()!=prb_tsi)
    if (HaveContactMeshtying() && disi != Teuchos::null)
      cmtman_->GetStrategy().Recover(disi_);

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
  if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0))
  {
    IO::cout << "Structural predictor "
             << INPAR::STR::PredEnumString(pred_)
             << " yields ";

    // relative check of force residual
    if ( normtypefres_ == INPAR::STR::convnorm_rel )
    {
      IO::cout << "scaled res-norm "
               << normfres_/normcharforce_
               << IO::endl;
    }
    // absolute check of force residual
    else if ( normtypefres_ == INPAR::STR::convnorm_abs )
    {
      IO::cout << "absolute res-norm "
               << normfres_
               << IO::endl;
    }
    // mixed absolute-relative check of force residual
    else if ( normtypefres_ == INPAR::STR::convnorm_mix )
    {
      IO::cout << "mixed res-norm "
               << min(normfres_, normfres_/normcharforce_)
               << IO::endl;
    }
    // default
    else
    {
      dserror("You should not turn up here.");
    }
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
  if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0) and  printiter_ )
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

/*----------------------------------------------------------------------*/
void STR::TimIntImpl::PrintNewtonIterHeader
(
  FILE* ofile
)
{
  // open outstd::stringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::STR::convnorm_rel:
    oss <<std::setw(16)<< "rel-res-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss <<std::setw(16)<< "abs-res-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss <<std::setw(16)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss <<std::setw(16)<< "abs-inco-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
    }
  }

  switch ( normtypedisi_ )
  {
  case INPAR::STR::convnorm_rel:
    oss <<std::setw(16)<< "rel-dis-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss <<std::setw(16)<< "abs-dis-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss <<std::setw(16)<< "mix-dis-norm";
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss <<std::setw(16)<< "abs-pre-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
    }
  }

  // add norms of Lagrange multiplier parts (contact/meshtying in saddlepoint formulation only)
  if(HaveContactMeshtying())
  {
    // strategy and system setup types
    INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
    INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM");
    INPAR::CONTACT::WearType        wtype   = DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(cmtman_->GetStrategy().Params(),"WEARTYPE");

    if (soltype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed)
    {
      switch ( normtypecontconstr_ )
      {
      case INPAR::STR::convnorm_rel:
        oss <<std::setw(20)<< "rel-contconstr-norm";
        break;
      case INPAR::STR::convnorm_abs :
        oss <<std::setw(20)<< "abs-contconstr-norm";
        break;
      default:
        dserror("You should not turn up here.");
        break;
      }

      switch ( normtypeplagrincr_ )
      {
      case INPAR::STR::convnorm_rel:
        oss <<std::setw(20)<< "rel-lagrincr-norm";
        break;
      case INPAR::STR::convnorm_abs :
      {
        oss <<std::setw(20)<< "abs-lagrincr-norm";
        if (wtype == INPAR::CONTACT::wear_discr)
        {
          oss <<std::setw(20)<< "abs-wearincr-norm";
          oss <<std::setw(20)<< "abs-wearcon-norm";
        }
        break;
      }
      default:
        dserror("You should not turn up here.");
        break;
      }
    }
  }

  // add constraint norm
  if (conman_->HaveConstraintLagr())
  {
    oss << std::setw(16)<< "abs-constr-norm";
  }

  if (itertype_==INPAR::STR::soltech_ptc)
  {
    oss << std::setw(16)<< "        PTC-dti";
  }

  // add solution time
  //oss << std::setw(14)<< "wct";
  oss << std::setw(13)<< "ts";
  oss << std::setw(10)<< "te";
  if (HaveContactMeshtying())
    oss << std::setw(10)<< "tc";

  // add contact set information
  if (HaveContactMeshtying())
  {
    // only print something for contact, not for meshtying
    INPAR::CONTACT::ApplicationType apptype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(cmtman_->GetStrategy().Params(),"APPLICATION");
    if (apptype == INPAR::CONTACT::app_mortarcontact)
    {
      oss << std::setw(10)<< "#active";
      if (cmtman_->GetStrategy().Friction())
        oss << std::setw(10)<< "#slip";
    }
  }

  // add plasticity information
  if (HaveSemiSmoothPlasticity())
  {
    oss << std::setw(10) << "#plastic";
  }

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
  // open outstd::stringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7)<< iter_;

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::STR::convnorm_rel:
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normfres_/normcharforce_;
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normfres_;
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << min(normfres_, normfres_/normcharforce_);
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normpfres_;
      break;
    default:
      dserror("You should not turn up here.");
      break;
    }
  }

  switch ( normtypedisi_ )
  {
  case INPAR::STR::convnorm_rel:
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_/normchardis_;
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << min(normdisi_, normdisi_/normchardis_);
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  if (pressure_ != Teuchos::null)
  {
    switch (normtypepfres_)
    {
    case INPAR::STR::convnorm_abs :
      oss << std::setw(16) << std::scientific << normpres_;
      break;
    default:
      dserror("You should not turn up here.");
      break;
    }
  }

  // add norms of Lagrange multiplier parts (contact/meshtying in saddlepoint formulation only)
  if(HaveContactMeshtying())
  {
    // strategy and system setup types
    INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
    INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtman_->GetStrategy().Params(),"SYSTEM");
    INPAR::CONTACT::WearType        wtype   = DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(cmtman_->GetStrategy().Params(),"WEARTYPE");

    if (soltype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed)
    {
      // we only support abs norms
      oss << std::setw(20) << std::setprecision(5) << std::scientific << normcontconstr_; // RHS for contact constraints
      oss << std::setw(20) << std::setprecision(5) << std::scientific << normlagr_;    // norm Lagrange multipliers

      if (wtype == INPAR::CONTACT::wear_discr)
      {
        oss << std::setw(20) << std::setprecision(5) << std::scientific << normw_;       // norm wear
        oss << std::setw(20) << std::setprecision(5) << std::scientific << normwrhs_;    // norm wear rhs
      }
    }
  }

  // add constraint norm
  if (conman_->HaveConstraintLagr())
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normcon_;
  }

  if (itertype_==INPAR::STR::soltech_ptc)
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << dti_;
  }

  // add solution time
  //oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_->ElapsedTime();
  oss << std::setw(13) << std::setprecision(2) << std::scientific << dtsolve_;
  oss << std::setw(10) << std::setprecision(2) << std::scientific << dtele_;
  if (HaveContactMeshtying())
    oss << std::setw(10) << std::setprecision(2) << std::scientific << dtcmt_;

  // add contact set information
  if (HaveContactMeshtying())
  {
    // only print something for contact, not for meshtying
    INPAR::CONTACT::ApplicationType apptype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(cmtman_->GetStrategy().Params(),"APPLICATION");
    if (apptype == INPAR::CONTACT::app_mortarcontact)
    {
      oss << std::setw(10) << cmtman_->GetStrategy().NumberOfActiveNodes();
      if (cmtman_->GetStrategy().Friction())
        oss << std::setw(10) << cmtman_->GetStrategy().NumberOfSlipNodes();
    }
  }

  // add plasticity information
  if (HaveSemiSmoothPlasticity())
  {
    oss << std::setw(10) << plastman_->NumActivePlasticGP();
  }


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
/* Export active set and characteristic calculation times into text files */
void STR::TimIntImpl::ExportContactQuantities()
{
  //add integration time contribution from every newton step
  inttime_global_+=cmtman_->GetStrategy().Inttime();

  double iteration=(double)iter_+1.0;
  double curinttime=(cmtman_->GetStrategy().Inttime())/(iteration);

  std::cout << "*** averaged inttime per newton step =  " << curinttime << std::endl;
  std::cout << "*** total inttime per time step= " << curinttime*iteration << std::endl;

  // write number of active nodes for converged newton in textfile xx x.active
  FILE* MyFile = NULL;
  std::ostringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << filebase <<".active";
  MyFile = fopen(filename.str().c_str(), "at+");

  // store active set
  if (MyFile)
  {
    fprintf(MyFile, "%d\t", cmtman_->GetStrategy().NumberOfActiveNodes());
    fprintf(MyFile, "%d\n", cmtman_->GetStrategy().NumberOfSlipNodes());
    fclose(MyFile);
  }
  else
    dserror("ERROR: File could not be opened.");


  //write required time
  FILE* MyFile2 = NULL;
  std::ostringstream filename2;
  const std::string filebase2 = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename2 << filebase2 <<".time";
  MyFile2 = fopen(filename2.str().c_str(), "at+");

  // store characteristic times
  if (MyFile2)
  {
    fprintf(MyFile2, "%g\t", dtsolve_);
    fprintf(MyFile2, "%g\t", dtele_);
    fprintf(MyFile2, "%g\t", dtcmt_);
    fprintf(MyFile2, "%g\t", curinttime);
    fprintf(MyFile2, "%g\n", curinttime*iteration);
    fclose(MyFile2);
  }
  else
    dserror("ERROR: File could not be opened.");

  return;
}

/*----------------------------------------------------------------------*/
/* print statistics of converged NRI */
void STR::TimIntImpl::PrintNewtonConv()
{
#ifdef CONTACTEXPORT
  // output integration time for contact and more...
  if (HaveContactMeshtying())
  {
    ExportContactQuantities();
  }
#endif

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
  if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0))
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
  // open outstd::stringstream
  std::ostringstream oss;

  // the text
  oss << "Finalised step " << std::setw(1) << step_;
  oss << " / " << std::setw(1) << stepmax_;
  oss << " | time " << std::setw(9) << std::setprecision(3) << std::scientific << (*time_)[0];
  oss << " | dt " << std::setw(9) << std::setprecision(3) << std::scientific << (*dt_)[0];
  oss << " | numiter " << std::setw(1) << iter_;
  oss << " | wct " << std::setw(8) << std::setprecision(2) << std::scientific << timer_->ElapsedTime();
  oss << "\n--------------------------------------------------------------------------------\n";

  // print to ofile (could be done differently...)
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // fall asleep
  return;
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
  solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, true);

  // go back
  return disi_;
}

/*----------------------------------------------------------------------*/
/* Prepare system for solving with Newton's method */
void STR::TimIntImpl::PrepareSystemForNewtonSolve()
{
  // rotate residual to local coordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global coordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);
  // blank residual at DOFs on Dirichlet BCs
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // rotate reaction forces back to global coordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(fres_);

  // make the residual negative
  fres_->Scale(-1.0);

  // transform stiff_ and fres_ to local coordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(SystemMatrix(),fres_);
  // local matrix and rhs required for correctly applying Dirichlet boundary
  // conditions: rows with inclined Dirichlet boundary condition can be blanked
  // and a '1.0' is put at the diagonal term

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));
  // final sip
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::TimIntImpl::UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
{
  // (re)allocate system matrix
  stiff_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*domainmaps,*rangemaps,81,false,true));
  mass_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*domainmaps,*rangemaps,81,false,true));
  if (damping_ != INPAR::STR::damp_none)
    damp_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*domainmaps,*rangemaps,81,false,true));

  // recalculate mass and damping matrices

  Teuchos::RCP<Epetra_Vector> fint
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  stiff_->Zero();
  mass_->Zero();

  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_nlnstiffmass");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    if (pressure_ != Teuchos::null) p.set("volume", 0.0);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement", zeros_);
    discret_->SetState("displacement", (*dis_)(0));
    if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", (*vel_)(0));
    //set state for volume coupling (e.g. for tsi and poro)
    SetCouplingState();
    discret_->Evaluate(p, stiff_, mass_, fint, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Add(*stiff_, false, dampk_, 0.0);
    damp_->Add(*mass_, false, dampm_, 1.0);
    damp_->Complete();
  }

  // in case of C0 pressure field, we need to get rid of
  // pressure equations
  if (pressure_ != Teuchos::null)
  {
    mass_->ApplyDirichlet(*(pressure_->CondMap()));
  }

  // We need to reset the stiffness matrix because its graph (topology)
  // is not finished yet in case of constraints and posssibly other side
  // effects (basically managers).
  stiff_->Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::TimIntImpl::STCPreconditioning()
{
  // print first system matrix to file in matlab format (DEBUGGING)
  #if 0
  if (iter_==1&& step_==0)
  {
    const std::string fname = "unscaled.mtl";
    if (myrank_ == 0)
      std::cout<<"Printing unscaled system matrix to file"<<std::endl;
      LINALG::PrintMatrixInMatlabFormat(fname,*((Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_))->EpetraMatrix()));
  }
  #endif

  if(stcscale_!=INPAR::STR::stc_none)
  {
    if (!stccompl_)
    {
      ComputeSTCMatrix();
      stccompl_=true;
    }

    stiff_ = MLMultiply(*(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_)),*stcmat_,false,false,true);
    if(stcscale_==INPAR::STR::stc_currsym)
    {
      stiff_ = MLMultiply(*stcmat_,true,*(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_)),false,false,false,true);
      Teuchos::RCP<Epetra_Vector> fressdc = LINALG::CreateVector(*dofrowmap_, true);
      stcmat_->Multiply(true,*fres_,*fressdc);
      fres_->Update(1.0,*fressdc,0.0);
    }

    // print first system matrix to file in matlab format (DEBUGGING)
    #if 0
    if (iter_==1&& step_==0)
    {
      const std::string fname = "scaled.mtl";
      if (myrank_ == 0)
        std::cout<<"Printing scaled system matrix to file"<<std::endl;
        LINALG::PrintMatrixInMatlabFormat(fname,*((Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_))->EpetraMatrix()));
    }
    #endif
  }

  return;
}

void STR::TimIntImpl::ComputeSTCMatrix()
{

  stcmat_->Zero();
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  discret_->  SetState("residual displacement", disi_);
  discret_->  SetState("displacement", disn_);

  const std::string action = "calc_stc_matrix";
  p.set("action", action);
  p.set<int>("stc_scaling", stcscale_);
  p.set("stc_layer",1);

  discret_-> Evaluate(p, stcmat_, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);

  stcmat_->Complete();

  #ifdef DEBUG
  if (iter_==1&& step_==0)
  {
    const std::string fname = "stcmatrix1.mtl";
    if (myrank_ == 0)
      std::cout<<"Printing stcmatrix1 to file"<<std::endl;
      LINALG::PrintMatrixInMatlabFormat(fname,*((Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stcmat_))->EpetraMatrix()));
  }
  #endif

  for (int lay = 2; lay <= stclayer_; ++lay)
  {
    Teuchos::ParameterList pe;

    pe.set("action", action);
    pe.set<int>("stc_scaling", stcscale_);
    pe.set("stc_layer", lay);

    Teuchos::RCP<LINALG::SparseMatrix> tmpstcmat=
      Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,true));
    tmpstcmat->Zero();

    discret_-> Evaluate(pe, tmpstcmat, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);
    tmpstcmat->Complete();

    #ifdef DEBUG
    if (iter_==1&& step_==0)
    {
      const std::string fname = "stcmatrix2.mtl";
      if (myrank_ == 0)
        std::cout<<"Printing stcmatrix2 to file"<<std::endl;
        LINALG::PrintMatrixInMatlabFormat(fname,*((Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(tmpstcmat))->EpetraMatrix()));
    }
    #endif

    stcmat_ = MLMultiply(*tmpstcmat,*stcmat_,false,false,true);
  }

  discret_-> ClearState();


}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::TimIntImpl::RecoverSTCSolution()
{
  if(stcscale_!=INPAR::STR::stc_none)
  {
    Teuchos::RCP<Epetra_Vector> disisdc = LINALG::CreateVector(*dofrowmap_, true);

    stcmat_->Multiply(false,*disi_,*disisdc);
    disi_->Update(1.0,*disisdc,0.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | recalculate the structural matrices for coupled problems  dano 03/10 |
 *----------------------------------------------------------------------*/
void STR::TimIntImpl::TSIMatrix()
{
  // recalculate mass and damping matrices

  Teuchos::RCP<Epetra_Vector> fint
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  stiff_->Zero();
  mass_->Zero();

  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_nlnstiffmass");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    if (pressure_ != Teuchos::null) p.set("volume", 0.0);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState(0,"residual displacement", zeros_);
    discret_->SetState(0,"displacement", (*dis_)(0));
    if (damping_ == INPAR::STR::damp_material) discret_->SetState(0,"velocity", (*vel_)(0));
    //set state for volume coupling (e.g. for tsi and poro)
    SetCouplingState();
    discret_->Evaluate(p, stiff_, mass_, fint, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Add(*stiff_, false, dampk_, 0.0);
    damp_->Add(*mass_, false, dampm_, 1.0);
    damp_->Complete();
  }

  // in case of C0 pressure field, we need to get rid of
  // pressure equations
  if (pressure_ != Teuchos::null)
  {
    mass_->ApplyDirichlet(*(pressure_->CondMap()));
  }

  // We need to reset the stiffness matrix because its graph (topology)
  // is not finished yet in case of constraints and posssibly other side
  // effects (basically managers).
  stiff_->Reset();
} // TSIMatrix


/*----------------------------------------------------------------------*/
