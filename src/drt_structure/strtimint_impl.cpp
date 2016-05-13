/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.cpp
\brief Implicit time integration for structural dynamics

\level 1

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>
#include <fenv.h>

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
#include "../drt_contact/meshtying_contact_bridge.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_wear.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_defines.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_constraint/windkessel_manager.H"
#include "../drt_constraint/springdashpot_manager.H"
#include "../drt_constraint/springdashpot.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_statmech/statmech_manager.H"
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
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex8p1j1.H"
#include "../drt_so3/so_shw6.H"
#include "../drt_so3/so_sh8p8.H"
#include "../drt_discsh3/discsh3.H"

#include "../drt_crack/crackUtils.H"

#include "../solver_nonlin/nln_operator_base.H"
#include "../solver_nonlin/nln_operator_factory.H"
#include "../solver_nonlin/nln_problem.H"
#include "../solver_nonlin/nln_problem_base.H"
#include "../solver_nonlin/nln_problem_nox.H"
#include "../solver_nonlin/nln_utils.H"
#include "strtimint_noxgroup.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntImpl::TimIntImpl
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
: TimInt
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
  toldisi_(sdynparams.get<double>("TOLDISP")),
  tolfres_(sdynparams.get<double>("TOLRES")),
  tolpfres_(sdynparams.get<double>("TOLINCO")),
  tolpres_(sdynparams.get<double>("TOLPRE")),
  uzawaparam_(sdynparams.get<double>("UZAWAPARAM")),
  uzawaitermax_(sdynparams.get<int>("UZAWAMAXITER")),
  tolcon_(sdynparams.get<double>("TOLCONSTR")),
  tolwindk_(DRT::Problem::Instance()->WindkesselStructuralParams().get<double>("TOLWINDKESSEL")),
  tolwindkdofincr_(DRT::Problem::Instance()->WindkesselStructuralParams().get<double>("TOLWINDKESSELDOFINCR")),
  iter_(-1),
  normcharforce_(0.0),
  normchardis_(0.0),
  normfres_(0.0),
  normdisi_(0.0),
  normcon_(0.0),
  normwindk_(0.0),
  normwindkdofincr_(0.0),
  normpfres_(0.0),
  normpres_(0.0),
  normcontconstr_(0.0),  // < norm of contact constraints (saddlepoint formulation)
  normlagr_(0.0),        // < norm of lagrange multiplier increment (saddlepoint formulation)
  alpha_ls_(sdynparams.get<double>("ALPHA_LS")),
  sigma_ls_(sdynparams.get<double>("SIGMA_LS")),
  ls_maxiter_(sdynparams.get<int>("LSMAXITER")),
  cond_res_(0.0),
  disi_(Teuchos::null),
  fres_(Teuchos::null),
  freact_(Teuchos::null),
  updateprojection_(false),
  stcscale_(DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdynparams, "STC_SCALING")),
  stclayer_(sdynparams.get<int>("STC_LAYER")),
  ptcdt_(sdynparams.get<double>("PTCDT")),
  dti_(1.0/ptcdt_)
{

  // general variable verifications:
  if (itermax_ < 0)
    dserror("MAXITER has to be greater than or equal to zero. Fix your input file.");

  if (itermin_ < 0)
    dserror("MINITER has to be greater than or equal to zero. Fix your input file.");

  if (toldisi_ <= 0)
    dserror("TOLDISP has to be greater than zero. Fix your input file.");

  if (tolfres_ <= 0)
    dserror("TOLRES has to be greater than zero. Fix your input file.");

  if (itermin_ > itermax_)
    dserror("ITERMIN has to be smaller than or equal to ITERMAX. Fix your input file.");

  if (tolpfres_ <= 0)
    dserror("TOLINCO has to be greater than zero. Fix your input file.");

  if (tolpres_ <= 0)
    dserror("TOLPRE has to be greater than zero. Fix your input file.");

  if (uzawaparam_ <= 0)
    dserror("UZAWAPARAM has to be greater than zero. Fix your input file.");

  if (uzawaitermax_ < 0)
    dserror("UZAWAMAXITER has to be greater than or equal to zero. Fix your input file.");

  if (tolcon_ <= 0)
    dserror("TOLCONSTR has to be greater than zero. Fix your input file.");

  if (tolwindk_ <= 0)
    dserror("TOLWINDKESSEL has to be greater than zero. Fix your input file.");

  if (tolwindkdofincr_ <= 0)
    dserror("TOLWINDKESSELDOFINCR has to be greater than zero. Fix your input file.");

  if ((alpha_ls_ <= 0) or (alpha_ls_ >= 1))
    dserror("Valid interval for ALPHA_LS is (0,1). Fix your input file.");

  if ((sigma_ls_ <= 0) or (sigma_ls_ >= 1))
    dserror("Valid interval for SIGMA_LS is (0,1). Fix your input file.");

  if (ls_maxiter_ < 0)
    dserror("LSMAXITER has to be greater than or equal to zero. Fix your input file.");

  if (ptcdt_ <= 0)
    dserror("PTCDT has to be greater than zero. Fix your input file.");


  // verify: if system has constraints implemented with Lagrange multipliers,
  // then Uzawa-type solver is used
  if (conman_->HaveConstraintLagr())
  {
    if ( (itertype_ != INPAR::STR::soltech_newtonuzawalin)
         and (itertype_ != INPAR::STR::soltech_newtonuzawanonlin) )
      dserror("Chosen solution technique %s does not work constrained.",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
  }
  else if (windkman_->HaveWindkessel())
  {
    if (itertype_ != INPAR::STR::soltech_newtonuzawalin)
      if (myrank_ == 0)
       dserror("Chosen solution technique %s does not work with Windkessel bc.",
                INPAR::STR::NonlinSolTechString(itertype_).c_str());
  }
  else if ( (itertype_ == INPAR::STR::soltech_newtonuzawalin)
            or (itertype_ == INPAR::STR::soltech_newtonuzawanonlin) )
  {
    dserror("Chosen solution technique %s does only work constrained or with Windkessel bc.",
            INPAR::STR::NonlinSolTechString(itertype_).c_str());
  }

  // Initiate Edge element for discrete shell elements
  if(DRT::Problem::Instance()->ProblemType() == prb_statmech)
  {
    const Teuchos::ParameterList&   statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
    INPAR::STATMECH::SimulationType simype  = DRT::INPUT::IntegralValue<INPAR::STATMECH::SimulationType>(statmechparams,"SIMULATION_TYPE");
    if (simype==INPAR::STATMECH::simulation_type_lipid_bilayer)
      InitializeEdgeElements();
  }

  // create empty residual force vector
  fres_ = LINALG::CreateVector(*DofRowMapView(), false);

  // create empty reaction force vector of full length
  freact_ = LINALG::CreateVector(*DofRowMapView(), false);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*DofRowMapView(), true);

  //prepare matrix for scaled thickness business of thin shell structures
  stcmat_=
    Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));

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
    tolcontconstr_      = cmtbridge_->GetStrategy().Params().get<double>("TOLCONTCONSTR");
    tollagr_            = cmtbridge_->GetStrategy().Params().get<double>("TOLLAGR");
    combfrescontconstr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(cmtbridge_->GetStrategy().Params(),"NORMCOMBI_RESFCONTCONSTR");
    combdisilagr_       = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(cmtbridge_->GetStrategy().Params(),"NORMCOMBI_DISPLAGR");
  }


  // setup binary operators for convergence check of semi-smooth plasticity problems
  combfresplconstr_ = INPAR::STR::bop_and;
  combdisiLp_       = INPAR::STR::bop_and;
  combfresEasres_   = INPAR::STR::bop_and;
  combdisiEasIncr_  = INPAR::STR::bop_and;
  if (HaveSemiSmoothPlasticity())
  {
    combfresplconstr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(*(plastman_->Params()),"NORMCOMBI_RESFPLASTCONSTR");
    combdisiLp_       = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(*(plastman_->Params()),"NORMCOMBI_DISPPLASTINCR");
    if (plastman_->EAS())
    {
      combfresEasres_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(*(plastman_->Params()),"NORMCOMBI_EASRES");
      combdisiEasIncr_= DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(*(plastman_->Params()),"NORMCOMBI_EASINCR");
    }
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

  // prepare line search
  if (itertype_==INPAR::STR::soltech_newtonls)
    PrepareLineSearch();

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

/*----------------------------------------------------------------------------*
 | Create Edges of for discrete shell elements       mukherjee (public)  04/15|
 *-----------------------------------------------------------------------------*/
void STR::TimIntImpl::InitializeEdgeElements()
{
  facediscret_ = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
  facediscret_->CreateInternalFacesExtension(true);
  facediscret_->FillCompleteFaces();
  return;
}

/*----------------------------------------------------------------------*/
/* predict solution */
void STR::TimIntImpl::Predict()
{
  // things that need to be done before Predict
  PrePredict();

  // Update locals systems (which may be time dependent)
  if (locsysman_ != Teuchos::null)
    locsysman_->Setup(timen_);

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
  else if ( pred_ == INPAR::STR::pred_constvel)
  {
    PredictConstVelConsistAcc();
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else if ( pred_ == INPAR::STR::pred_constacc)
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

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, accn_, false);

  // possibly initialise Lagrange multipliers to zero
  //  if ( (conman_->HaveConstraint())
  //       and (itertype_ == soltech_uzawalinnewton) )
  //  {
  //    conman_->ScaleLagrMult(0.0);
  //  }

  // create parameter list to hand in boolean flag indicating that this a predictor
  Teuchos::ParameterList params;
  params.set<bool>("predict",true);

  // set predictor type
  if (HaveSemiSmoothPlasticity())
  {
    plastman_->SetData().pred_type_=pred_;
    plastman_->SetData().no_recovery_=false; // recovery here also includes prediction
    plastman_->SetData().no_pl_condensation_=false;
    }
  // residual of condensed variables (e.g. EAS) for NewtonLS
  if (fresn_str_!=Teuchos::null)
  {
    params.set<double>("cond_rhs_norm",0.);
    params.set<int>("MyPID",myrank_);
  }

  // compute residual forces fres_ and stiffness stiff_
  // If we use a tangential predictor, the contact status could have been changed in contrast
  // to a constant predictor. Thus the contact status has to be reevaluated! (hiermeier 22.01.2014)
  if (pred_ == INPAR::STR::pred_tangdis) params.set<bool>("predict", false);

  // compute residual forces fres_ and stiffness stiff_
  EvaluateForceStiffResidual(params);

  // get residual of condensed variables (e.g. EAS) for NewtonLS
  if (fresn_str_!=Teuchos::null)
  {
    double loc=params.get<double>("cond_rhs_norm");
    discret_->Comm().SumAll(&loc,&cond_res_,1);
  }

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

  // things that need to be done after Predict
  PostPredict();

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

  // create parameter list to hand in boolean flag indicating that this a predictor
  Teuchos::ParameterList params;
  params.set<bool>("predict",true);

  // compute residual forces fres_ and stiffness stiff_
  EvaluateForceStiffResidual(params);

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
/* Check for LS with condensed variables and do preparations */
void STR::TimIntImpl::PrepareLineSearch()
{
  // each proc searchs through his elements
  int haveCondensationLocal=0;
  int haveCondensationGlobal=0;

  // for semi-smooth Newton plasticity we need it anyway
  if (!HaveSemiSmoothPlasticity())
  {
    // each proc searches through his elements
    for (int i=0; i<discret_->NumMyRowElements(); i++)
    {
      DRT::Element* actele = discret_->lRowElement(i);
      DRT::ELEMENTS::So_hex8* ele_hex8 = dynamic_cast<DRT::ELEMENTS::So_hex8*>(actele);
      if (
          (ele_hex8!=NULL && ele_hex8->HaveEAS()==true)
          || (actele->ElementType() == DRT::ELEMENTS::So_Hex8P1J1Type::Instance())
          || (actele->ElementType() == DRT::ELEMENTS::So_shw6Type::Instance())
      )
        haveCondensationLocal=1;
      if (actele->ElementType() == DRT::ELEMENTS::So_sh8p8Type::Instance())
        dserror("no line search for this element implemented.\n"
                "Feel free to implement similar to hex8 with EAS");
    }
    discret_->Comm().MaxAll(&haveCondensationLocal,&haveCondensationGlobal,1);
  }
  if (haveCondensationGlobal || HaveSemiSmoothPlasticity())
  {
    fresn_str_ = LINALG::CreateVector(*DofRowMapView(), true);
    fintn_str_ = LINALG::CreateVector(*DofRowMapView(), true);
  }
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
    = LINALG::CreateVector(*DofRowMapView(), true);

  // copy last converged displacements
  dbcinc->Update(1.0, *(*dis_)(0), 0.0);

  // get Dirichlet values at t_{n+1}
  ApplyDirichletBC(timen_, dbcinc, Teuchos::null, Teuchos::null, false);

  // subtract the displacements of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *(*dis_)(0), 1.0);

  // create parameter list to hand in boolean flag indicating that this a predictor
  Teuchos::ParameterList params;
  params.set<bool>("predict",true);

  // hand in flag indicating tangential predictor
  if (HaveSemiSmoothPlasticity())
  {
    plastman_->SetData().pred_type_=INPAR::STR::pred_tangdis;
    plastman_->SetData().no_pl_condensation_=true;
  }
  // compute residual forces fres_ and stiffness stiff_
  // at disn_, etc which are unchanged
  EvaluateForceStiffResidual(params);

  // add linear reaction forces to residual
  {
    // linear reactions
    Teuchos::RCP<Epetra_Vector> freact
      = LINALG::CreateVector(*DofRowMapView(), true);
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

  // recover contact / meshtying Lagrange multipliers
  if(HaveContactMeshtying())
    cmtbridge_->Recover(disi_);

  // decide which norms have to be evaluated
  bool bPressure = pressure_ != Teuchos::null;
  bool bContactSP = (HaveContactMeshtying() &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed ||
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed_lagmult));

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
    Teuchos::RCP<Epetra_Vector> lagrincr  = cmtbridge_->GetStrategy().LagrMultSolveIncr();

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
    Teuchos::ParameterList p;
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

  // get Teuchos::RCP to kernel vector of projector
  // since we are in 'pointvalue' mode, weights are changed implicitly
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  c->PutScalar(0.0);

  // We recompute the entire nullspace no matter what.
  // This is not nice yet since:
  // - translations are constant throughout the entire computation
  // - SAME nullspace is sometimes recomputed AGAIN for some iterative solvers
  // So here is space for optimization.

  // get number of modes and their ids
  std::vector<int> modeids = projector_->Modes();

  // Teuchos::RCP on vector of size 0 holding the nullspace data - resized within ComputeNullspace
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
/* evaluate external forces and its linearization at t_{n+1} */
void STR::TimIntImpl::ApplyForceStiffExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< old displacement state
  const Teuchos::RCP<Epetra_Vector> disn,  //!< new displacement state
  const Teuchos::RCP<Epetra_Vector> vel,  //!< velocity state
  Teuchos::RCP<Epetra_Vector>& fext,  //!< external force
  Teuchos::RCP<LINALG::SparseOperator>& fextlin //!<linearization of external force
)
{
  Teuchos::ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"displacement", dis);

  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);
  // get load vector
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  bool loadlin = (DRT::INPUT::IntegralValue<int>(sdyn, "LOADLIN") == 1);

  if (!loadlin)
    discret_->EvaluateNeumann(p, *fext);
  else
  {
    discret_->SetState(0,"displacement new", disn);
    discret_->EvaluateNeumann(p, fext, fextlin);
  }

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at state */
void STR::TimIntImpl::ApplyForceStiffInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,
  const Teuchos::RCP<Epetra_Vector> disi,
  const Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> fint,
  Teuchos::RCP<LINALG::SparseOperator> stiff,
  Teuchos::ParameterList& params,
  Teuchos::RCP<LINALG::SparseOperator> damp
)
{
  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********
  if(HaveStatMechBilayer())
  {
    // get reference volume
    params.set("action", "calc_struct_refvol");
    Teuchos::RCP<Epetra_SerialDenseVector> vol_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(params, vol_ref);

    // get reference CG
    params.set("action", "calc_struct_refCG");
    Teuchos::RCP<Epetra_SerialDenseVector> CG_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(params, CG_ref);

    // get reference area
    params.set("action", "calc_struct_refarea");
    Teuchos::RCP<Epetra_SerialDenseVector> area_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(params, area_ref);

    // get current volume
    discret_->SetState("displacement", disn_);
    params.set("action", "calc_struct_currvol");
    Teuchos::RCP<Epetra_SerialDenseVector> vol_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(params, vol_curr);

    // get current CG
    params.set("action", "calc_struct_currCG");
    Teuchos::RCP<Epetra_SerialDenseVector> CG_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(params, CG_curr);

    // get current area
    params.set("action", "calc_struct_currarea");
    Teuchos::RCP<Epetra_SerialDenseVector> area_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(params, area_curr);
    discret_->ClearState();

    params.set("reference volume",double((*vol_ref)(0)));
    params.set("current volume",double((*vol_curr)(0)));
    params.set<Teuchos::RCP<Epetra_SerialDenseVector> >("reference CG", CG_ref);
    params.set<Teuchos::RCP<Epetra_SerialDenseVector> >("current CG", CG_curr);
    params.set("reference area",double((*area_ref)(0)));
    params.set("current area",double((*area_curr)(0)));
   }
  // action for elements
  const std::string action = "calc_struct_nlnstiff";
  params.set("action", action);
  // other parameters that might be needed by the elements
  params.set("total time", time);
  params.set("delta time", dt);
  params.set("damping", damping_);
  params.set<int>("young_temp", young_temp_);
  if (pressure_ != Teuchos::null) params.set("volume", 0.0);

  // compute new inner radius
  discret_->ClearState();
  discret_->SetState(0,"displacement",dis);
  PATSPEC::ComputeEleInnerRadius(discret_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"residual displacement", disi);
  discret_->SetState(0,"displacement", dis);
  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector

  // Set material displacement state for ale-wear formulation
  if( (dismatn_!=Teuchos::null) )
    discret_->SetState(0,"material_displacement",dismatn_);

  // set plasticity data
  if (HaveSemiSmoothPlasticity())
  {
    plastman_->SetPlasticParams(params);
    if (DRT::Problem::Instance()->ProblemType() == prb_tsi)
      discret_->SetState(0,"velocity",vel);
    plastman_->SetData().dt_=(*dt_)[0];
  }

  /* Additionally we hand in "fint_str_"
   * This is usually Teuchos::null unless we do line search in
   * combination with elements that perform a local condensation
   * e.g. hex8 with EAS or semi-smooth Newton plasticity.
   * In such cases, fint_str_ contains the right hand side
   * without the modifications due to the local condensation procedure.
   */
  if (fintn_str_!=Teuchos::null)
    fintn_str_->PutScalar(0.);
  discret_->Evaluate(params, stiff, damp, fint, Teuchos::null, fintn_str_);
  discret_->ClearState();

  if (HaveFaceDiscret())
  {
    AssembleEdgeBasedMatandRHS(params,fint,dis,vel);
  }
  // get plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(params);

#if 0
  if (pressure_ != Teuchos::null)
    std::cout << "Total volume=" << std::scientific << p.get<double>("volume") << std::endl;
#endif

  // *********** time measurement ***********
  dtele_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate inertia force and its linearization */
void STR::TimIntImpl::ApplyForceStiffInternalAndInertial
(
  const double time,
  const double dt,
  const double timintfac_dis,
  const double timintfac_vel,
  const Teuchos::RCP<Epetra_Vector> dis,
  const Teuchos::RCP<Epetra_Vector> disi,
  const Teuchos::RCP<Epetra_Vector> vel,
  const Teuchos::RCP<Epetra_Vector> acc,
  Teuchos::RCP<Epetra_Vector> fint,
  Teuchos::RCP<Epetra_Vector> finert,
  Teuchos::RCP<LINALG::SparseOperator> stiff,
  Teuchos::RCP<LINALG::SparseOperator> mass,
  Teuchos::ParameterList& params,
  const double beta,
  const double gamma,
  const double alphaf,
  const double alpham
)
{
  // action for elements
  const std::string action = "calc_struct_nlnstiffmass";
  params.set("action", action);
  // other parameters that might be needed by the elements
  params.set("total time", time);
  params.set("delta time", dt);

  params.set("timintfac_dis", timintfac_dis);
  params.set("timintfac_vel", timintfac_vel);

  if (HaveNonlinearMass()==INPAR::STR::ml_rotations)
  {
    params.set("rot_beta", beta);
    params.set("rot_gamma", gamma);
    params.set("rot_alphaf", alphaf);
    params.set("rot_alpham", alpham);
  }

  // set plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->SetPlasticParams(params);

  // compute new inner radius
  discret_->ClearState();
  discret_->SetState(0,"displacement",dis);
  PATSPEC::ComputeEleInnerRadius(discret_);

  discret_->ClearState();
  discret_->SetState(0,"residual displacement", disi);
  discret_->SetState(0,"displacement", dis);
  discret_->SetState(0,"velocity", vel);
  discret_->SetState(0,"acceleration", acc);

  // Set material displacement state for struct-ale in cell migration
  if( (dismatn_!=Teuchos::null) )
    discret_->SetState(0,"material_displacement",dismatn_);

  /* Additionally we hand in "fint_str_"
   * This is usually Teuchos::null unless we do line search in
   * combination with elements that perform a local condensation
   * e.g. hex8 with EAS or semi-smooth Newton plasticity.
   * In such cases, fint_str_ contains the right hand side
   * without the modifications due to the local condensation procedure.
   */
  discret_->Evaluate(params, stiff, mass, fint, finert, fintn_str_);
  discret_->ClearState();

  // get plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(params);

  mass->Complete();

  return;
};

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
    Teuchos::ParameterList p;
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
    Teuchos::ParameterList p; // create the parameters for manager
    p.set("pot_man", potman_);
    p.set("total time", time);
    potman_->EvaluatePotential(p, dis, fint, mat);
    stiff = mat;
  }

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
      Teuchos::ParameterList p; // create the parameters for manager
      p.set("pot_man", potman_);
      p.set("total time", time);

      Teuchos::RCP<LINALG::SparseMatrix> stiff_test =
          Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
      Teuchos::RCP<Epetra_Vector> fint_test =
          LINALG::CreateVector(*DofRowMapView(), true);
      fint_test->PutScalar(0.0);
      stiff_test->Zero();

      potman_->TestEvaluatePotential(p, dis, fint_test, stiff_test, time, step);
    }
  }

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

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate forces due to Windkessel bcs */
void STR::TimIntImpl::ApplyForceStiffWindkessel
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis,
  const Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::ParameterList pwindk
)
{
  if (windkman_->HaveWindkessel())
  {
    windkman_->StiffnessAndInternalForces(time, dis, disn, pwindk);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate forces and stiffness due to spring dashpot BCs */
void STR::TimIntImpl::ApplyForceStiffSpringDashpot
(
  Teuchos::RCP<LINALG::SparseOperator> stiff,
  Teuchos::RCP<Epetra_Vector> fint,
  Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector> veln,
  bool predict,
  Teuchos::ParameterList psprdash
)
{
  if (springman_->HaveSpringDashpot())
  {
    springman_->StiffnessAndInternalForces(stiff,fint,disn,veln,psprdash);
  }

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

    if(cmtbridge_->HaveContact())
    {
      if (cmtbridge_->ContactManager()->GetStrategy().HasPoroNoPenetration())
      {
        // set structural velocity for poro normal no penetration
        Teuchos::RCP<Epetra_Vector> svel = Teuchos::rcp(new Epetra_Vector(*Velnp()));
        cmtbridge_->ContactManager()->GetStrategy().SetState("svelocity", svel);
      }
    }

     // make contact / meshtying modifications to lhs and rhs
    // (depending on whether this is a predictor step or not)
    if(cmtbridge_->HaveMeshtying())
      cmtbridge_->MtManager()->GetStrategy().ApplyForceStiffCmt(dis,stiff,fresm,stepn_,iter_,predict);
    if(cmtbridge_->HaveContact())
    {
      dynamic_cast<CONTACT::CoAbstractStrategy&>(cmtbridge_->ContactManager()->GetStrategy()).SetParentState("displacement",dis,discret_);
      cmtbridge_->ContactManager()->GetStrategy().ApplyForceStiffCmt(dis,stiff,fresm,stepn_,iter_,predict);
    }


    // scaling back
    fresm->Scale(-1.0);

    // *********** time measurement ***********
    dtcmt_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // visualization of current Newton step
#ifdef MORTARGMSH2
    bool gmsh = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH");
    if (gmsh) cmtbridge_->VisualizeGmsh(stepn_,iter_);
#endif // #ifdef MORTARGMSH2
  }

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

    // create empty parameter list
    Teuchos::ParameterList beamcontactparams;
    beamcontactparams.set("iter", iter_);
    beamcontactparams.set("dt", (*dt_)[0]);
    beamcontactparams.set("numstep", step_);

    // make contact / meshtying modifications to lhs and rhs
    // (set boolean flag 'newsti' to true, which activates
    // sclaing of contact stiffness with appropriate scaling
    // factor, e.g. (1.0-alphaf), internally)
    beamcman_->Evaluate(*SystemMatrix(),*fresm,*dis,beamcontactparams,true,timen_);

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

  return;
}

/*----------------------------------------------------------------------*/
/* Check residual displacement and limit it if necessary*/
void STR::TimIntImpl::LimitStepsizeBeamContact
(
  Teuchos::RCP<Epetra_Vector>& disi
)
{
  if(HaveBeamContact())
  {
    double minimal_radius = beamcman_->GetMinEleRadius();
    double maxdisiscalefac = beamcman_->BeamContactParameters().get<double>("BEAMS_MAXDISISCALEFAC");
    if(maxdisiscalefac>0)
    {
      double disi_infnorm=0.0;
      disi->NormInf(&disi_infnorm);

      while(disi_infnorm>maxdisiscalefac*minimal_radius)
      {
        if(myrank_==0)
          std::cout << "      Residual displacement scaled! (Minimal element radius: " << minimal_radius << ")" << std::endl;

        disi->Scale(0.5);
        disi->NormInf(&disi_infnorm);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Check residual displacement and limit it if necessary*/
void STR::TimIntImpl::LimitStepsizeBeam (Teuchos::RCP<Epetra_Vector>& disi)
{

  Teuchos::ParameterList StatMechParams= statmechman_->GetStatMechParams();
  double alink = StatMechParams.get<double>("ALINK", 4.751658e-06);
  if(alink==0.0)
  {
    {
      double disi_infnorm=0.0;
      disi->NormInf(&disi_infnorm);
      double minimal_radius = 0.00123;
      double maxdisiscalefac = 1000;

      while(disi_infnorm>maxdisiscalefac*minimal_radius)
      {
        std::cout << "Residual displacement scaled! (Minimal element radius: " << minimal_radius << ")" << std::endl;
        disi->Scale(0.5);
        disi->NormInf(&disi_infnorm);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for displacements
 * originally by lw */
double STR::TimIntImpl::CalcRefNormDisplacement()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  double charnormdis = 0.0;
  if (pressure_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector((*dis_)(0));
    charnormdis = STR::AUX::CalculateVectorNorm(iternorm_, disp);
  }
  else
    charnormdis = STR::AUX::CalculateVectorNorm(iternorm_, (*dis_)(0));

  // rise your hat
  return charnormdis;
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

  // check Windkessel
  bool wk = true;
  bool wkincr = true;
  if (windkman_->HaveWindkessel())
  {
    wk = normwindk_ < tolwindk_;
    wkincr = normwindkdofincr_ < tolwindkdofincr_;
  }

  // check contact (active set)
  bool ccontact = true;
  if (HaveContactMeshtying())
  {
    // check which case (application, strategy) we are in
    INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
    bool semismooth = DRT::INPUT::IntegralValue<int>(cmtbridge_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

    // only do this convergence check for semi-smooth Lagrange multiplier contact
    if (cmtbridge_->HaveContact() &&
        (stype == INPAR::CONTACT::solution_lagmult || stype == INPAR::CONTACT::solution_augmented) && semismooth)
      ccontact = cmtbridge_->GetStrategy().ActiveSetSemiSmoothConverged();

    // add convergence check for saddlepoint formulations
    // use separate convergence checks for contact constraints and
    // LM increments
    INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM");
    if ((stype==INPAR::CONTACT::solution_lagmult || stype==INPAR::CONTACT::solution_augmented) && (systype!=INPAR::CONTACT::system_condensed ||
        systype!=INPAR::CONTACT::system_condensed_lagmult)) {
      bool convDispLagrIncr = false;
      bool convDispWIncr = false;
      bool convDispWMIncr = false;

      switch ( normtypeplagrincr_ )
      {
      case INPAR::STR::convnorm_abs:
        convDispLagrIncr = normlagr_ < tollagr_;
        convDispWIncr = normw_ < 1e-12; // WEAR
        convDispWMIncr = normwm_ < 1e-12; // WEAR
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
        convdis = convdis and convDispLagrIncr and convDispWIncr and convDispWMIncr;
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

  // check convergence of plasticity
  bool cplast=true;
  if (HaveSemiSmoothPlasticity())
  {
    // check convergence of plastic active set
    cplast = cplast and plastman_->ActiveSetConverged();

    // convergence of residual
    if (combfresplconstr_==INPAR::STR::bop_and)
      convfres = convfres and plastman_->ConstraintConverged();
    else if (combfresplconstr_==INPAR::STR::bop_or)
      convfres = convfres or plastman_->ConstraintConverged();
    else
      dserror("Something went terribly wrong with binary operator!");

    // convergence of increments
    if (combdisiLp_==INPAR::STR::bop_and)
      convdis = convdis and plastman_->IncrementConverged();
    else if (combdisiLp_==INPAR::STR::bop_or)
      convdis = convdis or plastman_->IncrementConverged();
    else
      dserror("Something went terribly wrong with binary operator!");

    if (plastman_->EAS())
    {
      // convergence of residual
      if (combfresEasres_==INPAR::STR::bop_and)
        convfres = convfres and plastman_->EasResConverged();
      else if (combfresEasres_==INPAR::STR::bop_or)
        convfres = convfres or plastman_->EasResConverged();
      else
        dserror("Something went terribly wrong with binary operator!");

      // convergence of increments
      if (combdisiEasIncr_==INPAR::STR::bop_and)
        convdis = convdis and plastman_->EasIncrConverged();
      else if (combdisiEasIncr_==INPAR::STR::bop_or)
        convdis = convdis or plastman_->EasIncrConverged();
      else
        dserror("Something went terribly wrong with binary operator!");
    }
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
  return (conv and cc and wk and wkincr and ccontact and cplast);
}

/*----------------------------------------------------------------------*/
/* solve equilibrium */
INPAR::STR::ConvergenceStatus STR::TimIntImpl::Solve()
{
  // things to be done before solving
  PreSolve();

  int nonlin_error = 0 ;
  // special nonlinear iterations for contact / meshtying
  if (HaveContactMeshtying())
  {
    // check additionally if we have contact AND a Windkessel or constraint bc
    if (HaveWindkessel()) nonlin_error = CmtWindkConstrNonlinearSolve();
    else nonlin_error = CmtNonlinearSolve();
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
    case INPAR::STR::soltech_newtonls :
      nonlin_error = NewtonLS();
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
    case INPAR::STR::soltech_nlnsol :
      nonlin_error = NlnSolver();
      break;
    // catch problems
    default :
      dserror("Solution technique \"%s\" is not implemented.",
              INPAR::STR::NonlinSolTechString(itertype_).c_str());
      break;
    }
  }

  // since it is possible that the nonlinear solution fails only on some procs
  // we need to communicate the error
  int lnonlin_error = nonlin_error;
  Discretization()->Comm().MaxAll(&lnonlin_error, &nonlin_error, 1);

  INPAR::STR::ConvergenceStatus status = static_cast<INPAR::STR::ConvergenceStatus>(nonlin_error);

  // Only relevant, if the input parameter DIVERCONT is used and set to divcontype_ == adapt_step:
  // In this case, the time step size is halved as consequence of a non-converging nonlinear solver.
  // After a prescribed number of converged time steps, the time step is doubled again. The following
  // methods checks, if the time step size can be increased again.
  CheckForTimeStepIncrease(status);

  // things to be done after solving
  PostSolve();

  return status;
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

  if (outputeveryiter_)
  {
    int restart = DRT::Problem::Instance()->Restart();
    if (stepn_ == (restart + 1)) outputcounter_ = 0;
    OutputEveryIter(true);
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  int element_error= 0;
  int linsolve_error= 0;
  // equilibrium iteration loop
  while ( ( (not Converged() and (not linsolve_error) and (not element_error)) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
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

    //In beam contact applications it can be necessary to limit the Newton step size (scaled residual displacements)
    LimitStepsizeBeamContact(disi_);

    // recover standard displacements
    RecoverSTCSolution();

    // recover contact / meshtying Lagrange multipliers
    if(HaveContactMeshtying())
      cmtbridge_->Recover(disi_);

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // update end-point displacements etc
    UpdateIter(iter_);

    if (outputeveryiter_) OutputEveryIter(true);

    // create empty parameter list
    Teuchos::ParameterList params;

    // set flag for element error in form of a negative Jacobian determinant
    // in parameter list in case of potential continuation
    if (divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
    {
      params.set<bool>("tolerate_errors",true);
      params.set<bool>("eval_error",false);
    }

    // compute residual forces #fres_ and stiffness #stiff_
    // whose components are globally oriented
    EvaluateForceStiffResidual(params);

    // check for element error in form of a negative Jacobian determinant
    // in case of potential continuation
    if (divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
      element_error = ElementErrorCheck(params.get<bool>("eval_error"));

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
        ((DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
        (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed ||
         DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed_lagmult)) ||
        (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_augmented)));

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
      // extract subvectors (for mt and contact use only contact lm)
      Teuchos::RCP<Epetra_Vector> lagrincr  = cmtbridge_->GetStrategy().LagrMultSolveIncr();
      Teuchos::RCP<Epetra_Vector> constrrhs = cmtbridge_->GetStrategy().ConstrRhs();

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
      INPAR::WEAR::WearType wtype =
          DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(cmtbridge_->GetStrategy().Params(),"WEARTYPE");
      INPAR::WEAR::WearSide wside =
          DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(cmtbridge_->GetStrategy().Params(),"WEAR_SIDE");

      if(wtype==INPAR::WEAR::wear_primvar)
      {
        Teuchos::RCP<Epetra_Vector> wincr  = cmtbridge_->GetStrategy().WSolveIncr();
        Teuchos::RCP<Epetra_Vector> wearrhs = cmtbridge_->GetStrategy().WearRhs();

        if(wearrhs!=Teuchos::null) normwrhs_ = STR::AUX::CalculateVectorNorm(iternorm_, wearrhs);
        else normwrhs_ = -1.0;

        if(wincr!=Teuchos::null) normw_ = STR::AUX::CalculateVectorNorm(iternorm_, wincr);
        else normw_ = -1.0;

        if(wside==INPAR::WEAR::wear_both)
        {
          Teuchos::RCP<Epetra_Vector> wmincr  = cmtbridge_->GetStrategy().WMSolveIncr();
          Teuchos::RCP<Epetra_Vector> wearmrhs = cmtbridge_->GetStrategy().WearMRhs();

          if(wearmrhs!=Teuchos::null) normwmrhs_ = STR::AUX::CalculateVectorNorm(iternorm_, wearmrhs);
          else normwmrhs_ = -1.0;

          if(wmincr!=Teuchos::null) normwm_ = STR::AUX::CalculateVectorNorm(iternorm_, wmincr);
          else normwm_ = -1.0;
        }
        else
        {
          normwm_    = 0.0;
          normwmrhs_ = 0.0;
        }
      }
      else
      {
        normw_     = 0.0;
        normwrhs_  = 0.0;
        normwm_    = 0.0;
        normwmrhs_ = 0.0;
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
  return NewtonFullErrorCheck(linsolve_error,element_error);
}

/*----------------------------------------------------------------------*/
/* error check for full Newton problems */
int STR::TimIntImpl::NewtonFullErrorCheck(int linerror, int eleerror)
{
  // if everything is fine print to screen and return
  if (Converged())
  {
    if(myrank_ == 0)
      PrintNewtonConv();
    return 0;
  }
  // now some error checks: do we have an element problem
  // only check if we continue in this case; other wise, we ignore the error
  if (eleerror and divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
  {
    return eleerror;
  }
  // do we have a problem in the linear solver
  // only check if we want to do something fancy other wise we ignore the error in the linear solver
  else if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation or divcontype_==INPAR::STR::divcont_adapt_penaltycontact))
  {
    return linerror;
  }
  else
  {
    if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_stop ) )
    {
      // write restart output of last converged step before stopping
      Output(true);

      dserror("Newton unconverged in %d iterations", iter_);
      return 1;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_continue ) )
    {
      if (myrank_ == 0)
        IO::cout<<"Newton unconverged in " << iter_ << " iterations, continuing" <<IO::endl;
      return 0;
    }
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation or divcontype_==INPAR::STR::divcont_adapt_penaltycontact))
    {
      if (myrank_ == 0)
        IO::cout<< "Newton unconverged in " << iter_ << " iterations " << IO::endl;
      return 1;
    }
  }
  dserror("Fatal error in NonLinSolveErrorCheck, case not implemented ");
  return 0;
}

/*----------------------------------------------------------------------*/
/* error check for linear solver problems */
int STR::TimIntImpl::LinSolveErrorCheck(int linerror)
{
  // we only care about problems in the linear solver if we have a fancy divcont action
  if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation or divcontype_==INPAR::STR::divcont_adapt_penaltycontact) )
  {
    if (myrank_ == 0)
    IO::cout<< "Linear solver is having trouble " << IO::endl;
    return 2;
  }
  else
  {
    return 0;
  }

}

/*----------------------------------------------------------------------*/
/* error check for element problems in form of a negative Jacobian determinant */
int STR::TimIntImpl::ElementErrorCheck(bool evalerr)
{
  // merly care about element problems if there is a fancy divcont action
  // and element errors are considered
  if (evalerr and divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
  {
    if (myrank_ == 0)
    IO::cout<< "Element error in form of a negative Jacobian determinant " << IO::endl;
    return 3;
  }
  else
  {
    return 0;
  }

}

/*---------------------------------------------------------------------*/
/* solution with line search algorithm                  hiermeier 08/13*/
/*---------------------------------------------------------------------*/
int STR::TimIntImpl::NewtonLS()
{
  // The specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #stiff_ is the effective dynamic stiffness matrix

  int linsolve_error= 0;
  int fscontrol = 0;             // integer for a first step control (equal 1: deactivation) // fixme check if this is necessary for structural mechanics
  bool eval_error = false;       // an error occurred in the structure evaluation

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled())
    dserror("Effective stiffness matrix must be filled here");

  if (outputeveryiter_)
  {
    int restart = DRT::Problem::Instance()->Restart();
    if (stepn_ == (restart + 1)) outputcounter_ = 0;
    OutputEveryIter(true);
  }

  // initialize equilibrium loop (outer Full Newton loop)
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  // Merit function at current stage and for ls step
  std::vector<double> merit_fct (2);

  // Temporal copies of different vectors. Necessary for the sufficient decrease check.
  Teuchos::RCP<Epetra_Vector> tdisn = Teuchos::rcp(new Epetra_Vector(*disn_));
  Teuchos::RCP<Epetra_Vector> tveln = Teuchos::rcp(new Epetra_Vector(*veln_));
  Teuchos::RCP<Epetra_Vector> taccn = Teuchos::rcp(new Epetra_Vector(*accn_));

  // equilibrium iteration loop (outer full Newton loop)
  while ( ( (not Converged() and (not linsolve_error)) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // initialize the Newton line search iteration counter
    int iter_ls  = 0;
    double step_red = 1.0;

    /**************************************************************
    ***           Save successful iteration state               ***
    ***************************************************************/

    // It's necessary to save a temporal copy of the end-point displacements,
    // before any update is performed (because of the pseudo energy norm):
    tdisn->Update(1.0, *disn_, 0.0);                  // copy of the displ vector
    tveln->Update(1.0, *veln_, 0.0);                  // copy of the velocity vector
    taccn->Update(1.0, *accn_, 0.0);                  // copy of the acceleration vector

    /**************************************************************
    ***                       Solver Call                       ***
    ***************************************************************/
    linsolve_error = LsSolveNewtonStep();

    // Evaluate merit function
    if (iter_==1)
      LsEvalMeritFct(merit_fct[0]);
    else
      merit_fct[0]=merit_fct[1];

    // Check if pred_constdis is used. If yes, the first step is not controlled.
    if (pred_ == INPAR::STR::pred_constdis or pred_ == INPAR::STR::pred_constdisvelacc)
      fscontrol = 1;
    else if ((pred_==INPAR::STR::pred_tangdis || pred_==INPAR::STR::pred_constacc ||
             pred_==INPAR::STR::pred_constvel)|| (iter_ > 1))
      fscontrol=0;
    else
      dserror("The behavior of the chosen predictor is not yet tested in the line search framework.");

    /**************************************************************
    ***      Update right-hand side and stiffness matrix        ***
    ***************************************************************/
    Teuchos::ParameterList params;
    params.set<bool>("tolerate_errors",true);
    params.set<bool>("eval_error",false);
    if (fresn_str_!=Teuchos::null)
    {
      // attention: though it is called rhs_norm it actually contains sum x_i^2, i.e. the square of the L2-norm
      params.set<double>("cond_rhs_norm",0.);
      // need to know the processor id
      params.set<int>("MyPID",myrank_);
    }
    {
      int exceptcount = 0;
      fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
      EvaluateForceStiffResidual(params);
      if (fetestexcept(FE_INVALID) || fetestexcept(FE_OVERFLOW)
          || fetestexcept(FE_DIVBYZERO) || params.get<bool>("eval_error")==true)
        exceptcount  = 1;
      int tmp=0;
      discret_->Comm().SumAll(&exceptcount,&tmp,1);
      if (tmp)
        eval_error=true;
      feclearexcept(FE_ALL_EXCEPT);
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    }

    // get residual of condensed variables (e.g. EAS) for NewtonLS
    if (fresn_str_!=Teuchos::null)
    {
      double loc=params.get<double>("cond_rhs_norm");
      discret_->Comm().SumAll(&loc,&cond_res_,1);
    }

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

    /**************************************************************
    ***           merit function (current iteration)            ***
    ***************************************************************/
    int err=LsEvalMeritFct(merit_fct[1]);
    eval_error = (eval_error || err);

    if (outputeveryiter_) OutputEveryIter(true);

    /**************************************************************
    ***          1st inner LINE SEARCH loop                     ***
    ***************************************************************/

    while ((iter_-fscontrol > 0) && ((!LsConverged(& merit_fct[0], step_red) || eval_error) && (iter_ls < ls_maxiter_)))
    {
      /**************************************************************
      ***           Display line search information               ***
      ***************************************************************/
      if (iter_ls==0)
        LsPrintLineSearchIter(&merit_fct[0],iter_ls,step_red);

      // increase inner loop count
      ++iter_ls;

      /**************************************************************
      ***                   Step size control                     ***
      ***************************************************************/
      step_red *= alpha_ls_;
      // >>>> displacement, velocity, acceleration <<<<<<<<<<<<<<<
      // scale displ. increment
      disi_->Scale(alpha_ls_);
      // load old displ. vector
      disn_->Update(1.0, *tdisn, 0.0);
      // load old vel. vector
      veln_->Update(1.0, *tveln, 0.0);
      // load old acc. vector
      accn_->Update(1.0, *taccn, 0.0);

      // Update nodal displ., vel., acc., etc.
      UpdateIter(iter_);
      /**************************************************************
      ***   Update right-hand side (and part. stiffness matrix)   ***
      ***************************************************************/
      LsUpdateStructuralRHSandStiff(eval_error,merit_fct[1]);

      /**************************************************************
      ***           Display line search information               ***
      ***************************************************************/
      LsPrintLineSearchIter(&merit_fct[0],iter_ls,step_red);

      if (!(eval_error) && (outputeveryiter_)) OutputEveryIter(true, true);
    }

    if (iter_ls!=0)
    {
      if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0) and  printiter_ )
      {
        std::ostringstream oss;
        std::string dashline;
        dashline.assign(64,'-');
        oss << dashline ;
        // print to screen (could be done differently...)
        if (printerrfile_)
        {
          fprintf(errfile_, "%s\n", oss.str().c_str());
          fflush(errfile_);
        }

        fprintf(stdout, "%s\n", oss.str().c_str());
        fflush(stdout);
      }
    }

    /**************************************************************
    ***      Print Newton Step information                      ***
    ***************************************************************/

    // build residual force norm
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);

    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  } // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // call monitor
  if (conman_->HaveMonitor())
    conman_->ComputeMonitorValues(disn_);

  //do nonlinear solver error check
  return NewtonFullErrorCheck(linsolve_error,0);
}


/*----------------------------------------------------------------------*/
/*   Solver Call (line search)                          hiermeier 09/13 */
/*----------------------------------------------------------------------*/
int STR::TimIntImpl::LsSolveNewtonStep()
{
  int linsolve_error = 0;
  /**************************************************************
  ***           Prepare the solution procedure                ***
  ***************************************************************/
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

  /**************************************************************
  ***                     Solver Call                         ***
  ***************************************************************/
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

  linsolve_error = solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, iter_==1, projector_);
  // check for problems in linear solver
  // however we only care about this if we have a fancy divcont action (meaning function will return 0 )
  linsolve_error = LinSolveErrorCheck(linsolve_error);

  //In beam contact applications it can be necessary to limit the Newton step size (scaled residual displacements)
  LimitStepsizeBeamContact(disi_);

  solver_->ResetTolerance();

  // recover standard displacements
  RecoverSTCSolution();

  // *********** time measurement ***********
  dtsolve_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // update end-point displacements etc
  UpdateIter(iter_);

  return(linsolve_error);
}

/*----------------------------------------------------------------------*/
/*   Update structural RHS and stiff (line search)      hiermeier 09/13 */
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::LsUpdateStructuralRHSandStiff(bool& isexcept,double& merit_fct)
{
  // --- Checking for floating point exceptions
  fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  // compute residual forces #fres_ and stiffness #stiff_
  // whose components are globally oriented
  int exceptcount = 0;
  Teuchos::ParameterList params;
  // elements may tolerate errors usually leading to dserrors
  // in such cases the elements force the line search to reduce
  // the step size by setting "eval_error" to true
  params.set<bool>("tolerate_errors",true);
  params.set<bool>("eval_error",false);
  // condensed degrees of freedom need to know the step reduction
  params.set<double>("alpha_ls",alpha_ls_);
  // line search needs to know the residuals of additional condensed dofs
  if (fresn_str_!=Teuchos::null)
  {
    params.set<double>("cond_rhs_norm",0.);
    // need to know the processor id
    params.set<int>("MyPID",myrank_);
  }
  EvaluateForceStiffResidual(params);

  // get residual of condensed variables (e.g. EAS) for NewtonLS
  if (fresn_str_!=Teuchos::null)
  {
    double loc=params.get<double>("cond_rhs_norm");
    discret_->Comm().SumAll(&loc,&cond_res_,1);
  }

  if (fetestexcept(FE_INVALID) || fetestexcept(FE_OVERFLOW)
      || fetestexcept(FE_DIVBYZERO) || params.get<bool>("eval_error")==true)
    exceptcount  = 1;

  // synchronize the exception flag isexcept on all processors
  int exceptsum = 0;
  discret_->Comm().SumAll(& exceptcount, & exceptsum, 1);
  if (exceptsum > 0)
    isexcept = true;
  else
    isexcept=false;

  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  feclearexcept(FE_ALL_EXCEPT);
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

  /**************************************************************
  ***          merit function (current iteration)             ***
  ***************************************************************/
  int err=LsEvalMeritFct(merit_fct);
  isexcept = (isexcept || err);

  return;
}


/*----------------------------------------------------------------------*/
/*   Evaluate the merit function (line search)          hiermeier 08/13 */
/*----------------------------------------------------------------------*/
int STR::TimIntImpl::LsEvalMeritFct(double& merit_fct)
{
  fedisableexcept(FE_OVERFLOW);
  int err=0;
  // Calculate the quadratic norm of the right-hand side as merit function
  // Calculate the merit function value: (1/2) * <RHS,RHS>
  if (fresn_str_==Teuchos::null)
  {
    err=fres_->Dot(*fres_,& merit_fct);
  }
  else
  {
    merit_fct=0.;
    err=fresn_str_->Dot(*fresn_str_,&merit_fct);
    merit_fct+=cond_res_;
  }
  merit_fct *= 0.5 ;

  int exceptcount=0;
  if (fetestexcept(FE_OVERFLOW))
      exceptcount  = 1;
  int exceptsum = 0;
  discret_->Comm().SumAll(& exceptcount, & exceptsum, 1);
  if (exceptsum!=0)
    return err;
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_OVERFLOW);

  return 0;
}

/*----------------------------------------------------------------------*/
/*   Print information about the last line search step  hiermeier 09/13 */
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::LsPrintLineSearchIter(double* mf_value, int iter_ls, double step_red)
{
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
  // print to standard out
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0) and  printiter_ )
  {
    std::ostringstream oss;
    if (iter_ls== 0)
    {
      std::string dashline;
      dashline.assign(64,'-');
      oss << dashline << std::endl;
      oss << std::setw(6) << "ls_iter";
      oss << std::setw(16) << "step_scale";
      oss << std::setw(16) << "abs-dis-norm";
      oss << std::setw(16) << "merit-fct";
      oss << std::setw(10)<< "te";
      if (HaveSemiSmoothPlasticity())
        oss << std::setw(10) << "#active";
      oss << std::endl;
    }

    oss << std::setw(7) << iter_ls;
    oss << std::setw(16) << std::setprecision(5) << std::scientific << step_red;
    // build residual displacement norm
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
    if (iter_ls==0)
      oss << std::setw(16) << std::setprecision(5) << std::scientific << mf_value[0];
    else
      oss << std::setw(16) << std::setprecision(5) << std::scientific << mf_value[1];
    oss << std::setw(10) << std::setprecision(2) << std::scientific << dtele_;
    if (HaveSemiSmoothPlasticity())
      oss << std::setw(10) << std::scientific << plastman_->NumActivePlasticGP();

    // finish oss
    oss << std::ends;

      // print to screen (could be done differently...)
    if (printerrfile_)
    {
      fprintf(errfile_, "%s\n", oss.str().c_str());
      fflush(errfile_);
    }

    fprintf(stdout, "%s\n", oss.str().c_str());
    fflush(stdout);
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/*   Inner convergence check (line search)              hiermeier 08/13 */
/*----------------------------------------------------------------------*/
bool STR::TimIntImpl::LsConverged(double* mf_value,double step_red)
{
  bool check_ls_mf = false;

  /**************************************************************
  ***           Check for sufficient descent                  ***
  ***************************************************************/
  // mf_value[1]: NEW merit function value
  //            --> f(x + alpha_ls * dx)
  // mf_value[0]: OLD merit function value (initial value at the beginning of the time step
  //              or function value of the last converged iteration step. Converged means that
  //              the last step fulfilled the LsConverged test.)
  //            --> f(x)
  // The check follows to
  //            f(x + alpha_ls * dx) - f(x) <= - 2 * sigma_ls * step_red_ * f(x).
  check_ls_mf = ((mf_value[1] - mf_value[0]) <= - 2.0 * sigma_ls_ * step_red * mf_value[0]);


  return (check_ls_mf);
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
int STR::TimIntImpl::NlnSolver()
{
  // Note: fres_ has already been evaluated and comes as 'positive' residual

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled())
    dserror("Effective stiffness matrix must be filled here");

  // take care of Dirichlet boundary conditions
  ApplyDirichletBC(timen_, disn_, veln_, accn_, true);
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_, GetLocSysTrafo(), zeros_,
      *(dbcmaps_->CondMap()));

  // ---------------------------------------------------------------------------
  // Create / read parameter list for configuration of nonlinear solver
  // ---------------------------------------------------------------------------
  std::string filename =
      DRT::Problem::Instance()->NonlinearSolverParams().get<std::string>(
          "XML_FILE");

  // check for reasonable xml file
  if (filename == "none")
  {
    dserror("Seems like you forgot to set the XML file for configuration of "
        "the nonlinear solver.");
  }

  Teuchos::RCP<NLNSOL::UTILS::NlnConfig> nlnconfig =
      Teuchos::rcp(new NLNSOL::UTILS::NlnConfig());
  nlnconfig->Setup(filename);

  // ---------------------------------------------------------------------------
  // Create NOX group
  // ---------------------------------------------------------------------------
  // create initial guess vector of predictor result as CreateCopy to avoid
  // direct access to disn_
  NOX::Epetra::Vector noxSoln(disn_, NOX::Epetra::Vector::CreateCopy);

  // create NOX linear system to provide access to Jacobian
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
      NoxCreateLinearSystem(*nlnconfig->GetAllNonConstRcp(), noxSoln, noxutils_);

  /* use NOX::STR::Group to enable access to time integration
   * Note: NOX::Epetra::Group would be sufficient. */
  Teuchos::RCP<NOX::STR::Group> noxgrp =
      Teuchos::rcp(new NOX::STR::Group(*this,
          nlnconfig->GetSubListNonConst("Nonlinear Problem"), Teuchos::rcp(this, false),
          noxSoln, linSys));

  // ---------------------------------------------------------------------------
  // Create interface to nonlinear problem
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> probparams =
      Teuchos::rcp(new Teuchos::ParameterList());
  probparams->set("NOX Group",
      Teuchos::rcp_dynamic_cast<NOX::Abstract::Group>(noxgrp, true));
  probparams->set("Jacobian Operator",
      Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(stiff_));

  Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem =
      Teuchos::rcp(new NLNSOL::NlnProblemNox());
  nlnproblem->Init(Discretization()->Comm(), nlnconfig, "Nonlinear Problem",
      *probparams, DofRowMap(), Teuchos::null);
  nlnproblem->Setup();

  /* Evaluate once more to guarantee valid quantities inside of the
   * NOX::STR::Group() */
  nlnproblem->ComputeF(*disn_, *fres_);
  nlnproblem->ComputeJacobian();

  // ---------------------------------------------------------------------------
  // Create the nonlinear operator to solve the nonlinear problem
  // ---------------------------------------------------------------------------
  // use factory to create the nonlinear operator
  const std::string nlnopname = nlnconfig->GetParameter<std::string>(
      "outer nonlinear operator");

  NLNSOL::NlnOperatorFactory opfactory;
  Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnoperator = opfactory.Create(
      nlnconfig, nlnopname);
  nlnoperator->Init(Discretization()->Comm(), nlnconfig, nlnopname, nlnproblem,
      solver_, 0);
  nlnoperator->Setup();

  // ---------------------------------------------------------------------------
  // Apply the nonlinear operator
  // ---------------------------------------------------------------------------
  /* use local copy of #fres_ for ApplyInverse() to not mess with the member
   * variable #fres_ in the Evaluate call */
  Teuchos::RCP<Epetra_Vector> fres = Teuchos::rcp(new Epetra_Vector(*fres_));

  return nlnoperator->ApplyInverse(*fres, noxSoln.getEpetraVector());
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
void STR::TimIntImpl::UpdateStepWindkessel()
{
  if (windkman_ -> HaveWindkessel())
    windkman_->UpdateTimeStep();
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
bool STR::TimIntImpl::HaveWindkessel()
{
  return windkman_->HaveWindkessel();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool STR::TimIntImpl::HaveSpringDashpot()
{
  return springman_->HaveSpringDashpot();
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
/*----------------------------------------------------------------------*/
void STR::TimIntImpl::UpdateIterIncrWindkessel
(
  Teuchos::RCP<Epetra_Vector> wkdofincr ///< wk dof increment
)
{
  windkman_->UpdateWkDof(wkdofincr);
}

/*----------------------------------------------------------------------*/
/* do linearised Uzawa iterations with full NRI
 * originally by tk 11/07 */
int STR::TimIntImpl::UzawaLinearNewtonFull()
{
  int linsolve_error = 0;
  if (conman_->HaveConstraint())
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
    while ( ( (not Converged() and (not linsolve_error) ) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
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

      // create parameter list
      Teuchos::ParameterList params;

      // compute residual forces #fres_ and stiffness #stiff_
      // which contain forces and stiffness of constraints
      EvaluateForceStiffResidual(params);
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

      // why was this here? part of else statement below!!! (mhv 01/2015)
//      // build residual force norm
//      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
//      // build residual displacement norm
//      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
//      // build residual Lagrange multiplier norm
//      normcon_ = conman_->GetErrorNorm();

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
        // build residual Lagrange multiplier norm
        normcon_ = conman_->GetErrorNorm();
      }

      // print stuff
      PrintNewtonIter();

      // increment equilibrium loop index
      iter_ += 1;
    }  // end equilibrium loop

    // correct iteration counter
    iter_ -= 1;
  }
  else if (windkman_->HaveWindkessel())
  {
    // check whether we have a sanely filled stiffness matrix
    if (not stiff_->Filled())
    {
      dserror("Effective stiffness matrix must be filled here");
    }

    // initialise equilibrium loop
    iter_ = 1;
    normfres_ = CalcRefNormForce();
    // normdisi_ was already set in predictor; this is strictly >0
    normwindk_ = windkman_->GetWindkesselRHSNorm();
    normwindkdofincr_ = windkman_->GetWindkesselDofIncrNorm();
    timer_->ResetStartTime();

    // equilibrium iteration loop
    while ( ( (not Converged() and (not linsolve_error) ) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
    {
      // make negative residual
      fres_->Scale(-1.0);

      // uncomplete stiffness matrix, so stuff can be inserted again
      // stiff_->UnComplete();

      // transform to local co-ordinate systems
      if (locsysman_ != Teuchos::null)
        locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

      // apply Dirichlet BCs to system of equations
      disi_->PutScalar(0.0);  // Useful? depends on solver and more
      LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                       GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

      // *********** time measurement ***********
      double dtcpu = timer_->WallTime();
      // *********** time measurement ***********

      //Use STC preconditioning on system matrix
      STCPreconditioning();

      // linear solver call (contact / meshtying case or default)
      if (HaveContactMeshtying())
        CmtWindkConstrLinearSolve();
      else
      {
        // Call Windkessel solver to solve system
        linsolve_error = windkman_->Solve(SystemMatrix(),disi_,fres_);
      }

      // check for problems in linear solver
      // however we only care about this if we have a fancy divcont action  (meaning function will return 0)
      linsolve_error=LinSolveErrorCheck(linsolve_error);

      // recover contact / meshtying Lagrange multipliers
      if(HaveContactMeshtying())
        cmtbridge_->Recover(disi_);

      // *********** time measurement ***********
      dtsolve_ = timer_->WallTime() - dtcpu;
      // *********** time measurement ***********

      // transform back to global co-ordinate system
      if (locsysman_ != Teuchos::null)
        locsysman_->RotateLocalToGlobal(disi_);

      // update end-point displacements, velocities, accelerations
      UpdateIter(iter_);

      // create parameter list
      Teuchos::ParameterList params;

      // compute residual forces #fres_ and stiffness #stiff_
      // which contain forces and stiffness of Windkessels
      EvaluateForceStiffResidual(params);

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
        // build residual windkessel norm
        normwindk_ = windkman_->GetWindkesselRHSNorm();
        // build residual windkessel dof increment norm
        normwindkdofincr_ = windkman_->GetWindkesselDofIncrNorm();
      }

      // print stuff
      PrintNewtonIter();

      // increment equilibrium loop index
      iter_ += 1;
    }  // end equilibrium loop

    // correct iteration counter
    iter_ -= 1;
  }

  //do nonlinear solver error check
  return UzawaLinearNewtonFullErrorCheck(linsolve_error);
}

/*----------------------------------------------------------------------------*/
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

    //print Windkessel output
    if (windkman_->HaveWindkessel())
      windkman_->PrintPresFlux(false);

    return 0;
  }
  // now some error checks
  // do we have a problem in the linear solver
  // only check if we want to do something fancy other wise we ignore the error in the linear solver
  if(linerror and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation or divcontype_==INPAR::STR::divcont_adapt_penaltycontact) )
  {
    return linerror;
  }
  else
  {
    if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_stop ) )
    {
      // write restart output of last converged step before stopping
      Output(true);

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
    else if ( (iter_ >= itermax_) and (divcontype_==INPAR::STR::divcont_halve_step or divcontype_==INPAR::STR::divcont_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step or divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err or divcontype_==INPAR::STR::divcont_repeat_step or divcontype_==INPAR::STR::divcont_repeat_simulation or divcontype_==INPAR::STR::divcont_adapt_penaltycontact))
    {
      if (myrank_ == 0)
        IO::cout<< "Newton unconverged in " << iter_ << " iterations " << IO::endl;
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
  // strategy type
  INPAR::CONTACT::SolvingStrategy soltype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");

  // semi-smooth Newton type
  bool semismooth = DRT::INPUT::IntegralValue<int>(cmtbridge_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

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
    if (cmtbridge_->HaveContact() && semismooth)
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
    else if (cmtbridge_->HaveContact() && !semismooth)
    {
      // active set strategy
      int activeiter = 0;
      while (cmtbridge_->GetStrategy().ActiveSetConverged()==false)
      {
        // increase active set iteration index
        ++activeiter;

        // predictor step (except for first active set step)
        if (activeiter > 1) Predict();

        // nonlinear iteration
        int error = NewtonFull();
        if(error) return error;

        // update of active set (fixed-point)
        cmtbridge_->GetStrategy().UpdateActiveSet();
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
  // 4) AUGMENTED SEMI-SMOOTH NEWTON FOR CONTACT
  // The search for the correct active set (=contact nonlinearity) and
  // the large deformation linearization (=geometrical nonlinearity) are
  // merged into one semi-smooth Newton method and solved within ONE
  // iteration loop (which is then basically a standard Newton).
  // TODO EXTENSION TO A SMOOTH NEWTON APPROACH ARE UNDER CONSIDERATION
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_augmented)
  {
    if (cmtbridge_->HaveContact() && semismooth)
    {
      int error = 0;
      // nonlinear iteration
      if (itertype_ == INPAR::STR::soltech_newtonfull)
        error = NewtonFull();
      else if (itertype_ == INPAR::STR::soltech_newtonls)
        error = NewtonLS();

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
    cmtbridge_->GetStrategy().UpdateConstraintNorm();
  }

  //********************************************************************
  // Solving Strategy using Nitsche's method
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_nitsche)
  {
    // nonlinear iteration
    return NewtonFull();
  }

  //********************************************************************
  // Solving Strategy using Augmented Lagrange Techniques with Uzawa
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_uzawa)
  {
    // get tolerance and maximum Uzawa steps
    double eps = cmtbridge_->GetStrategy().Params().get<double>("UZAWACONSTRTOL");
    int maxuzawaiter = cmtbridge_->GetStrategy().Params().get<int>("UZAWAMAXSTEPS");

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
        cmtbridge_->GetStrategy().InitializeUzawa(stiff_,fres_);
        fres_->Scale(-1.0);
      }

      // nonlinear iteration
      int error = NewtonFull();
      if(error) return error;

      // update constraint norm and penalty parameter
      cmtbridge_->GetStrategy().UpdateConstraintNorm(uzawaiter);

      // store Lagrange multipliers for next Uzawa step
      cmtbridge_->GetStrategy().UpdateUzawaAugmentedLagrange();
      cmtbridge_->GetStrategy().StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

    } while (cmtbridge_->GetStrategy().ConstraintNorm() >= eps);

    // reset penalty parameter
    cmtbridge_->GetStrategy().ResetPenalty();
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

  INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
  INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM");

  // update information about active slave dofs
  //**********************************************************************
  // feed solver/preconditioner with additional information about the contact/meshtying problem
  //**********************************************************************
  {
    //TODO: maps for merged meshtying and contact problem !!!

    // feed Aztec based solvers with contact information
    if (contactsolver_->Params().isSublist("Aztec Parameters"))
    {
      Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("Aztec Parameters");
      Teuchos::RCP<Epetra_Map> masterDofMap;
      Teuchos::RCP<Epetra_Map> slaveDofMap;
      Teuchos::RCP<Epetra_Map> innerDofMap;
      Teuchos::RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtbridge_->GetStrategy());
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
      Teuchos::RCP<Epetra_Map> masterDofMap;
      Teuchos::RCP<Epetra_Map> slaveDofMap;
      Teuchos::RCP<Epetra_Map> innerDofMap;
      Teuchos::RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtbridge_->GetStrategy());
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
  // (1) Standard / Dual Lagrange multipliers -> SaddlePoint
  // (2) Direct Augmented Lagrange strategy
  //**********************************************************************
  if ((soltype==INPAR::CONTACT::solution_lagmult || soltype==INPAR::CONTACT::solution_augmented)
      && (systype!=INPAR::CONTACT::system_condensed && systype!=INPAR::CONTACT::system_condensed_lagmult))
  {
    // check if contact contributions are present,
    // if not we make a standard solver call to speed things up
    if (!cmtbridge_->GetStrategy().IsInContact() && !cmtbridge_->GetStrategy().WasInContact() && !cmtbridge_->GetStrategy().WasInContactLastTimeStep())
    {
      solver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
    }
    else
    {
      // otherwise, solve the saddle point linear system

      // will be removed soon
      //contactsolver_->Params().set<bool>("CONTACT",true); // for simpler precond
      //contactsolver_->Params().set<bool>("MESHTYING",true); // for simpler precond

      Teuchos::RCP<Epetra_Operator> blockMat = Teuchos::null;
      Teuchos::RCP<Epetra_Vector> blocksol = Teuchos::null;
      Teuchos::RCP<Epetra_Vector> blockrhs = Teuchos::null;

      // build the saddle point system
      cmtbridge_->GetStrategy().BuildSaddlePointSystem(stiff_, fres_, disi_, dbcmaps_, iter_-1, blockMat, blocksol, blockrhs);

      // solve the linear system
      contactsolver_->Solve(blockMat, blocksol, blockrhs, true, iter_==1);

      // split vector and update internal displacement and Lagrange multipliers
      cmtbridge_->GetStrategy().UpdateDisplacementsAndLMincrements(disi_, blocksol);
    }
  }

  //**********************************************************************
  // Solving a purely displacement based system
  // (1) Dual (not Standard) Lagrange multipliers -> Condensed
  // (2) Penalty and Uzawa Augmented Lagrange strategies
  //**********************************************************************
  else
  {
    if(cmtbridge_->HaveMeshtying())
    {
      // solve with contact solver
      contactsolver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
    }
    else if(cmtbridge_->HaveContact())
    {
      // check if contact contributions are present,
      // if not we make a standard solver call to speed things up
      if (!cmtbridge_->GetStrategy().IsInContact() &&
          !cmtbridge_->GetStrategy().WasInContact() &&
          !cmtbridge_->GetStrategy().WasInContactLastTimeStep())
      {
        // standard solver call (fallback solver for pure structure problem)
        solver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
        return;
      }

      // solve with contact solver
      contactsolver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
    }
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
  INPAR::BEAMCONTACT::Strategy strategy =
    DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcman_->BeamContactParameters(),"BEAMS_STRATEGY");

  // unknown types of nonlinear iteration schemes
  if (itertype_ != INPAR::STR::soltech_newtonfull)
    dserror("Unknown type of equilibrium iteration");

  //**********************************************************************
  // solving strategy using regularization with penalty method
  // (nonlinear solution approach: ordinary NEWTON)
  //**********************************************************************
  if (strategy == INPAR::BEAMCONTACT::bstr_penalty)
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
  else if (strategy == INPAR::BEAMCONTACT::bstr_uzawa)
  {
    // get tolerance and maximum number of Uzawa steps from input file
    double eps = beamcman_->BeamContactParameters().get<double>("BEAMS_BTBUZAWACONSTRTOL");
    int maxuzawaiter = beamcman_->BeamContactParameters().get<int>("BEAMS_BTBUZAWAMAXSTEPS");

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
        // create empty parameter list
        Teuchos::ParameterList beamcontactparams;
        beamcontactparams.set("iter", iter_);
        beamcontactparams.set("dt", (*dt_)[0]);
        beamcontactparams.set("numstep", step_);

        // make contact modifications to lhs and rhs
        beamcman_->InitializeUzawa(*SystemMatrix(),*fres_,*disn_,beamcontactparams,true);

        // scaling back
        fres_->Scale(-1.0);
      }

      // inner nonlinear iteration (Newton)
      int error = NewtonFull();
      if(error) return error;
      // update constraint norm and penalty parameter
      beamcman_->UpdateConstrNormUzawa();

      // update Uzawa Lagrange multipliers
      beamcman_->UpdateAlllmuzawa();

    } while (abs(beamcman_->GetConstrNorm()) >= eps);

    // reset penalty parameter, Uzawa index and Uzawa lagrange multiplier
    beamcman_->ResetCurrentpp();
    beamcman_->ResetUzawaIter();
    beamcman_->ResetAlllmuzawa();
    beamcman_->UpdateConstrNorm();
  }

  //**********************************************************************
  // misuse of beam contact module for GMSH output
  // (nonlinear solution approach: ordinary NEWTON)
  //**********************************************************************
  if (strategy == INPAR::BEAMCONTACT::bstr_gmshonly)
  {
     // nonlinear iteration (Newton)
    int error = NewtonFull();
    if(error) return error;
  }
  //**********************************************************************

  //**********************************************************************
  // unknown solving strategy
  //**********************************************************************
  else
  {
    dserror("ERROR: Chosen strategy not yet available for beam contact");
  }

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

  if (outputeveryiter_)
  {
    int restart = DRT::Problem::Instance()->Restart();
    if (stepn_ == (restart + 1)) outputcounter_ = 0;
    OutputEveryIter(true);
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  double ptcdt     = ptcdt_;
  double nc; fres_->NormInf(&nc);
  double dti = 1/ptcdt;

  int element_error= 0;
  int linsolve_error= 0;
  // equilibrium iteration loop
  while ( ( (not Converged() and (not linsolve_error) and (not element_error)) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // make negative residual
    fres_->Scale(-1.0);

    // transform to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

    // modify stiffness matrix with dti
    {
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      tmp->PutScalar(dti);
      Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
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
    if(HaveContactMeshtying())
      cmtbridge_->Recover(disi_);

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // update end-point displacements etc
    UpdateIter(iter_);

    if (outputeveryiter_) OutputEveryIter(true);

    // create parameter list
    Teuchos::ParameterList params;

    // set flag for element error in form of a negative Jacobian determinant
    // in parameter list in case of potential continuation
    if (divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
    {
      params.set<bool>("tolerate_errors",true);
      params.set<bool>("eval_error",false);
    }

    // compute residual forces #fres_ and stiffness #stiff_
    // whose components are globally oriented
    EvaluateForceStiffResidual(params);

    // check for element error in form of a negative Jacobian determinant
    // in case of potential continuation
    if (divcontype_==INPAR::STR::divcont_rand_adapt_step_ele_err)
      element_error = ElementErrorCheck(params.get<bool>("eval_error"));

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
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
        (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed ||
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM") != INPAR::CONTACT::system_condensed));

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
      Teuchos::RCP<Epetra_Vector> lagrincr  = cmtbridge_->GetStrategy().LagrMultSolveIncr();
      Teuchos::RCP<Epetra_Vector> constrrhs = cmtbridge_->GetStrategy().ConstrRhs();

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
      dti = std::max(dti,0.0);
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
  return NewtonFullErrorCheck(linsolve_error,element_error);

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
      cmtbridge_->Recover(disi_);

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
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0))
  {
    IO::cout << "Structural predictor for field '"<<discret_->Name()<<"' "
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
               << std::min(normfres_, normfres_/normcharforce_)
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
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0) and  printiter_ )
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
void STR::TimIntImpl::PrintNewtonIterHeader( FILE* ofile )
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
    INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
    INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM");
    INPAR::WEAR::WearType wtype   = DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(cmtbridge_->GetStrategy().Params(),"WEARTYPE");
    INPAR::WEAR::WearSide wside   = DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(cmtbridge_->GetStrategy().Params(),"WEAR_SIDE");

    if ((soltype==INPAR::CONTACT::solution_lagmult || soltype==INPAR::CONTACT::solution_augmented)
        && (systype!=INPAR::CONTACT::system_condensed && systype!=INPAR::CONTACT::system_condensed_lagmult))
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
        if (wtype == INPAR::WEAR::wear_primvar)
        {
          oss <<std::setw(20)<< "abs-wearincr-S-norm";
          oss <<std::setw(20)<< "abs-wearcon-S-norm";
          if(wside == INPAR::WEAR::wear_both)
          {
            oss <<std::setw(20)<< "abs-wearincr-M-norm";
            oss <<std::setw(20)<< "abs-wearcon-M-norm";
          }
        }
        break;
      }
      default:
        dserror("You should not turn up here.");
        break;
      }
    }
  }

  if (HaveSemiSmoothPlasticity())
  {
    oss <<std::setw(20)<< "abs-plconstr-norm";
    oss <<std::setw(20)<< "abs-dLPincr-norm";
    if (plastman_->EAS())
    {
      oss <<std::setw(20)<< "abs-easRes-norm";
      oss <<std::setw(20)<< "abs-easIncr-norm";
    }
  }

  // add constraint norm
  if (conman_->HaveConstraintLagr())
  {
    oss << std::setw(16)<< "abs-constr-norm";
  }

  // add Windkessel norm
  if (windkman_->HaveWindkessel())
  {
    oss << std::setw(16)<< "abs-wkres-norm";
    oss << std::setw(16)<< "abs-wkinc-norm";
  }

  if (itertype_==INPAR::STR::soltech_ptc)
  {
    oss << std::setw(16)<< "        PTC-dti";
  }

  // add solution time
  oss << std::setw(13)<< "ts";
  oss << std::setw(10)<< "te";
  if (HaveContactMeshtying())
    oss << std::setw(10)<< "tc";

  // add contact set information
  if (HaveContactMeshtying())
  {
    // only print something for contact, not for meshtying
    if (cmtbridge_->HaveContact())
    {
      oss << std::setw(11)<< "#active";
      if (cmtbridge_->GetStrategy().Friction())
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
void STR::TimIntImpl::PrintNewtonIterText( FILE* ofile )
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
    oss << std::setw(16) << std::setprecision(5) << std::scientific << std::min(normfres_, normfres_/normcharforce_);
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
    oss << std::setw(16) << std::setprecision(5) << std::scientific << std::min(normdisi_, normdisi_/normchardis_);
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
    INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
    INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM");
    INPAR::WEAR::WearType        wtype   = DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(cmtbridge_->GetStrategy().Params(),"WEARTYPE");
    INPAR::WEAR::WearSide        wside   = DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(cmtbridge_->GetStrategy().Params(),"WEAR_SIDE");

    if ((soltype==INPAR::CONTACT::solution_lagmult || soltype==INPAR::CONTACT::solution_augmented)
        && (systype!=INPAR::CONTACT::system_condensed && systype!=INPAR::CONTACT::system_condensed_lagmult))
    {
      // we only support abs norms
      oss << std::setw(20) << std::setprecision(5) << std::scientific << normcontconstr_; // RHS for contact constraints
      oss << std::setw(20) << std::setprecision(5) << std::scientific << normlagr_;    // norm Lagrange multipliers

      if (wtype == INPAR::WEAR::wear_primvar)
      {
        oss << std::setw(20) << std::setprecision(5) << std::scientific << normw_;       // norm wear
        oss << std::setw(20) << std::setprecision(5) << std::scientific << normwrhs_;    // norm wear rhs
        if(wside == INPAR::WEAR::wear_both)
        {
          oss << std::setw(20) << std::setprecision(5) << std::scientific << normwm_;       // norm wear
          oss << std::setw(20) << std::setprecision(5) << std::scientific << normwmrhs_;    // norm wear rhs
        }
      }
    }
  }

  if (HaveSemiSmoothPlasticity())
  {
    oss << std::setw(20) << std::setprecision(5) << std::scientific << plastman_->DeltaLp_residual_norm();
    oss << std::setw(20) << std::setprecision(5) << std::scientific << plastman_->DeltaLp_increment_norm();
    if (plastman_->EAS())
    {
      oss << std::setw(20) << std::setprecision(5) << std::scientific << plastman_->EAS_residual_norm();
      oss << std::setw(20) << std::setprecision(5) << std::scientific << plastman_->EAS_increment_norm();
    }
  }

  // add constraint norm
  if (conman_->HaveConstraintLagr())
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normcon_;
  }

  // add Windkessel norm
  if (windkman_->HaveWindkessel())
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normwindk_;
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normwindkdofincr_;
  }

  if (itertype_==INPAR::STR::soltech_ptc)
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << dti_;
  }

  // add solution time
  oss << std::setw(13) << std::setprecision(2) << std::scientific << dtsolve_;
  oss << std::setw(10) << std::setprecision(2) << std::scientific << dtele_;
  if (HaveContactMeshtying())
    oss << std::setw(10) << std::setprecision(2) << std::scientific << dtcmt_;

  // add contact set information
  if (HaveContactMeshtying())
  {
    // only print something for contact, not for meshtying
    INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
    bool semismooth =  DRT::INPUT::IntegralValue<int>(cmtbridge_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");
    if (cmtbridge_->HaveContact())
    {
      if (soltype==INPAR::CONTACT::solution_augmented && semismooth)
      {
        bool ccontact = cmtbridge_->GetStrategy().ActiveSetSemiSmoothConverged();
        // active set changed
        if (!ccontact) oss << std::setw(8) << cmtbridge_->GetStrategy().NumberOfActiveNodes() << "(c)";
        // active set didnot change
        else           oss << std::setw(8) << cmtbridge_->GetStrategy().NumberOfActiveNodes() << "(-)";
      }
      else
        oss << std::setw(11) << cmtbridge_->GetStrategy().NumberOfActiveNodes();
      if (cmtbridge_->GetStrategy().Friction())
        oss << std::setw(10) << cmtbridge_->GetStrategy().NumberOfSlipNodes();
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
  inttime_global_+=cmtbridge_->GetStrategy().Inttime();

  double iteration=(double)iter_+1.0;
  double curinttime=(cmtbridge_->GetStrategy().Inttime())/(iteration);

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
    fprintf(MyFile, "%d\t", cmtbridge_->GetStrategy().NumberOfActiveNodes());
    fprintf(MyFile, "%d\n", cmtbridge_->GetStrategy().NumberOfSlipNodes());
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
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0))
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
  // create parameter list
  Teuchos::ParameterList params;

  // Evaluate/define the residual force vector #fres_ for
  // relaxation solution with SolveRelaxationLinear
  EvaluateForceStiffResidualRelax(params);

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
void STR::TimIntImpl::PrepareSystemForNewtonSolve(const bool preparejacobian)
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

  // blank iterative increment
  disi_->PutScalar(0.0);  // Useful? depends on solver and more

  // apply Dirichlet BCs to system of equations
  if (preparejacobian)
  {
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
        GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));
  }

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
    = LINALG::CreateVector(*DofRowMapView(), true); // internal force

  stiff_->Zero();
  mass_->Zero();

  {
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // action for elements
    p.set("action", "calc_struct_nlnstiffmass");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);

    Teuchos::RCP<Epetra_Vector> finert = Teuchos::null;
    if ( HaveNonlinearMass() )
    {
      finert = LINALG::CreateVector(*DofRowMapView(), true); // intertial force
      // Note: the following parameters are just dummies, since they are only needed to calculate finert which we
      // will not use anyway
      p.set("timintfac_dis", 0.0); //dummy!
      p.set("timintfac_vel", 0.0); //dummy!
    }

    // compute new inner radius
    discret_->ClearState();
    discret_->SetState("displacement", (*dis_)(0));
    PATSPEC::ComputeEleInnerRadius(discret_);

    if (pressure_ != Teuchos::null) p.set("volume", 0.0);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement", zeros_);
    discret_->SetState("displacement", (*dis_)(0));
    discret_->SetState(0,"velocity", (*vel_)(0));
    discret_->SetState(0,"acceleration", (*acc_)(0));
    if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", (*vel_)(0));

    discret_->Evaluate(p, stiff_, mass_, fint, finert, Teuchos::null);
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
      Teuchos::RCP<Epetra_Vector> fressdc = LINALG::CreateVector(*DofRowMapView(), true);
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

/*----------------------------------------------------------------------------*/
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
    std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
    fname += ".stcmatrix1.mtl";
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
      Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(),81,true,true));
    tmpstcmat->Zero();

    discret_-> Evaluate(pe, tmpstcmat, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);
    tmpstcmat->Complete();

#ifdef DEBUG
    if (iter_==1&& step_==0)
    {
      std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
      fname += ".stcmatrix2.mtl";
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
    Teuchos::RCP<Epetra_Vector> disisdc = LINALG::CreateVector(*DofRowMapView(), true);

    stcmat_->Multiply(false,*disi_,*disisdc);
    disi_->Update(1.0,*disisdc,0.0);
  }

  return;
}


/*----------------------------------------------------------------------*/
/* solution with nonlinear iteration for contact / meshtying AND Windkessel bcs*/
int STR::TimIntImpl::CmtWindkConstrNonlinearSolve()
{
  //********************************************************************
  // get some parameters
  //********************************************************************
  // strategy type
  INPAR::CONTACT::SolvingStrategy soltype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");

  // semi-smooth Newton type
  bool semismooth = DRT::INPUT::IntegralValue<int>(cmtbridge_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

  // iteration type
  if (itertype_ != INPAR::STR::soltech_newtonuzawalin)
    dserror("Unknown type of equilibrium iteration! Choose newtonlinuzawa instead of fullnewton!");

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
    if (cmtbridge_->HaveContact() && semismooth)
    {
      // nonlinear iteration
     int error = UzawaLinearNewtonFull();
     if(error) return error;
    }

    //********************************************************************
    // 2) FIXED-POINT APPROACH FOR CONTACT
    // The search for the correct active set (=contact nonlinearity) is
    // represented by a fixed-point approach, whereas the large deformation
    // linearization (=geometrical nonlinearity) is treated by a standard
    // Newton scheme. This yields TWO nested iteration loops
    //********************************************************************
    else if (cmtbridge_->HaveContact() && !semismooth)
    {
      // active set strategy
      int activeiter = 0;
      while (cmtbridge_->GetStrategy().ActiveSetConverged()==false)
      {
        // increase active set iteration index
        ++activeiter;

        // predictor step (except for first active set step)
        if (activeiter > 1) Predict();

        // nonlinear iteration
        int error = UzawaLinearNewtonFull();
        if(error) return error;

        // update of active set (fixed-point)
        cmtbridge_->GetStrategy().UpdateActiveSet();
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
      int error = UzawaLinearNewtonFull();
      if(error) return error;
    }
  }

  //********************************************************************
  // Solving Strategy using Regularization Techniques (Penalty Method)
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_penalty)
  {
    // nonlinear iteration
    int error = UzawaLinearNewtonFull();
    if(error) return error;

    // update constraint norm
    cmtbridge_->GetStrategy().UpdateConstraintNorm();
  }

  //********************************************************************
  // Solving Strategy using Augmented Lagrange Techniques with Uzawa
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_uzawa)
  {
    // get tolerance and maximum Uzawa steps
    double eps = cmtbridge_->GetStrategy().Params().get<double>("UZAWACONSTRTOL");
    int maxuzawaiter = cmtbridge_->GetStrategy().Params().get<int>("UZAWAMAXSTEPS");

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
        cmtbridge_->GetStrategy().InitializeUzawa(stiff_,fres_);
        fres_->Scale(-1.0);
      }

      // nonlinear iteration
      int error = UzawaLinearNewtonFull();
      if(error) return error;

      // update constraint norm and penalty parameter
      cmtbridge_->GetStrategy().UpdateConstraintNorm(uzawaiter);

      // store Lagrange multipliers for next Uzawa step
      cmtbridge_->GetStrategy().UpdateUzawaAugmentedLagrange();
      cmtbridge_->GetStrategy().StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

    } while (cmtbridge_->GetStrategy().ConstraintNorm() >= eps);

    // reset penalty parameter
    cmtbridge_->GetStrategy().ResetPenalty();
  }

  return 0;
}



/*----------------------------------------------------------------------*/
/* linear solver call for contact / meshtying AND Windkessel bcs*/
void STR::TimIntImpl::CmtWindkConstrLinearSolve()
{

  // strategy and system setup types
  INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
  INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(cmtbridge_->GetStrategy().Params(),"SYSTEM");

  // update information about active slave dofs
  //**********************************************************************
  // feed solver/preconditioner with additional information about the contact/meshtying problem
  //**********************************************************************
  {
    // feed Aztec based solvers with contact information
    if (contactsolver_->Params().isSublist("Aztec Parameters"))
    {
      Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("Aztec Parameters");
      Teuchos::RCP<Epetra_Map> masterDofMap;
      Teuchos::RCP<Epetra_Map> slaveDofMap;
      Teuchos::RCP<Epetra_Map> innerDofMap;
      Teuchos::RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtbridge_->GetStrategy());
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
      Teuchos::RCP<Epetra_Map> masterDofMap;
      Teuchos::RCP<Epetra_Map> slaveDofMap;
      Teuchos::RCP<Epetra_Map> innerDofMap;
      Teuchos::RCP<Epetra_Map> activeDofMap;
      Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtbridge_->GetStrategy());
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
  // -> does not work together with constraints / Windkessel bcs
  // (1) Standard / Dual Lagrange multipliers -> SaddlePointCoupled
  // (2) Standard / Dual Lagrange multipliers -> SaddlePointSimpler
  //**********************************************************************
  if (soltype==INPAR::CONTACT::solution_lagmult &&
      (systype!=INPAR::CONTACT::system_condensed && systype!=INPAR::CONTACT::system_condensed_lagmult))
  {
    dserror("Constraints / Windkessel bcs together with saddle point contact system does not work (yet)!");
  }

  //**********************************************************************
  // Solving a purely displacement based system
  // (1) Dual (not Standard) Lagrange multipliers -> Condensed
  // (2) Penalty and Augmented Lagrange strategies
  //**********************************************************************
  else
  {
    // solve with Windkessel solver
    windkman_->Solve(SystemMatrix(),disi_,fres_);
  }

  return;
}

/*-----------------------------------------------------------------------------*
 * update the field vectors to account for new node               sudhakar 12/13
 * introduced by crack propagation
 *----------------------------------------------------------------------------*/
void STR::TimIntImpl::updateEpetraVectorsCrack( std::map<int,int>& oldnew )
{
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, disi_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fres_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, freact_, oldnew );

  updateMethodSpecificEpetraCrack( oldnew );
}

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                             meier 01/15
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void STR::TimIntImpl::CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status)
{
  const int maxnumfinestep = 4;

  if(divcontype_!=INPAR::STR::divcont_adapt_step)
    return;
  else if(status == INPAR::STR::conv_success and divconrefinementlevel_!=0)
  {
    divconnumfinestep_++;

    if(divconnumfinestep_==maxnumfinestep)
    {
      // increase the step size if the remaining number of steps is a even number
      if(((stepmax_-stepn_)%2)==0 and stepmax_!=stepn_)
      {
        IO::cout << "Nonlinear solver successful. Double timestep size!"
                 << IO::endl;

        divconrefinementlevel_--;
        divconnumfinestep_=0;

        stepmax_= stepmax_ - (stepmax_-stepn_)/2;

        // double the time step size
        (*dt_)[0]=(*dt_)[0]*2;
      }
      else  //otherwise we have to wait one more time step until the step size can be increased
      {
        divconnumfinestep_--;
      }
    }
    return;
  }
}
