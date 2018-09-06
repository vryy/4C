/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_contact.cpp

\brief Evaluation and assembly of all contact terms

\maintainer Michael Hiermeier

\date Feb 3, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "str_model_evaluator_contact.H"
#include "str_model_evaluator_data.H"
#include "str_model_evaluator.H"
#include "str_timint_base.H"
#include "str_impl_generic.H"
#include "str_dbc.H"
#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"
#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_solver_linesearchbased.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../drt_contact/contact_poro_lagrange_strategy.H"
#include "../drt_contact/contact_strategy_factory.H"
#include "../drt_contact_aug/contact_augmented_strategy.H"
#include "../drt_contact_aug/contact_aug_plot.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Contact::Contact() : strategy_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Setup()
{
  CheckInit();
  eval_contact_ptr_ = EvalData().ContactPtr();

  // ---------------------------------------------------------------------
  // create the contact factory
  // ---------------------------------------------------------------------
  CONTACT::STRATEGY::Factory factory;
  factory.Init(GStatePtr());
  factory.Setup();

  // check the problem dimension
  factory.CheckDimension();

  // create some local variables (later to be stored in strategy)
  std::vector<Teuchos::RCP<CONTACT::CoInterface>> interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  factory.ReadAndCheckInput(cparams);

  // check for FillComplete of discretization
  if (not Discret().Filled()) dserror("Discretization is not fillcomplete");

  // ---------------------------------------------------------------------
  // build the contact interfaces
  // ---------------------------------------------------------------------
  // FixMe Would be great, if we get rid of these poro parameters...
  bool poroslave = false;
  bool poromaster = false;
  factory.BuildInterfaces(cparams, interfaces, poroslave, poromaster);

  // ---------------------------------------------------------------------
  // build the solver strategy object
  // ---------------------------------------------------------------------
  strategy_ptr_ = factory.BuildStrategy(cparams, poroslave, poromaster, DofOffset(), interfaces);

  // build the search tree
  factory.BuildSearchTree(interfaces);

  // print final screen output
  factory.Print(interfaces, strategy_ptr_, cparams);

  // ---------------------------------------------------------------------
  // final touches to the contact strategy
  // ---------------------------------------------------------------------
  strategy_ptr_->StoreDirichletStatus(Int().GetDbc().GetDBCMapExtractor());
  strategy_ptr_->SetState(MORTAR::state_new_displacement, Int().GetDbc().GetZeros());
  strategy_ptr_->SaveReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->EvaluateReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->Inttime_init();
  strategy_ptr_->RedistributeContact(GState().GetDisN());
  strategy_ptr_->InitBinStrategyforTimestep(GState().GetVelN());

  CheckPseudo2D();

  PostSetup(cparams);

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostSetup(Teuchos::ParameterList& cparams)
{
  if (dynamic_cast<CONTACT::AUG::Strategy*>(strategy_ptr_.get()))
  {
    Teuchos::ParameterList& aug_params = cparams.sublist("AUGMENTED");
    Teuchos::ParameterList& plot_params = aug_params.sublist("PLOT");
    plot_params.set<const int*>("CURRENT_STEP", &GState().GetStepNp());
    plot_params.set<std::string>(
        "OUTPUT_FILE_NAME", GInOutput().GetOutputPtr()->Output()->FileName());
    plot_params.set<std::string>(
        "INPUT_FILE_NAME", GInOutput().GetOutputPtr()->Output()->InputFileName());
    plot_params.set<const DRT::DiscretizationInterface*>(
        "DISCRETIZATION", GState().GetDiscret().get());
    plot_params.set<STR::MODELEVALUATOR::Contact*>("MODELEVALUATOR", this);

    STR::IMPLICIT::Generic& impl = dynamic_cast<STR::IMPLICIT::Generic&>(Int());
    CONTACT::AUG::Plot::Create(impl.GetNoxParams(), plot_params, strategy_ptr_.get());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::CheckPseudo2D() const
{
  // print messages for multifield problems (e.g FSI)
  const PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  if ((probtype != prb_structure) and (GState().GetMyRank() == 0))
  {
    // warnings
#ifdef CONTACTPSEUDO2D
    std::cout << "WARNING: The flag CONTACTPSEUDO2D is switched on. If this "
              << "is a real 3D problem, switch it off!" << std::endl;
#else
    std::cout << "STR::MODELEVALUATOR::Contact::CheckPseudo2D -- "
              << "WARNING: \nThe flag CONTACTPSEUDO2D is switched off. If this "
              << "is a 2D problem modeled pseudo-3D, switch it on!" << std::endl;
#endif  // #ifdef CONTACTPSEUDO2D
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();
  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(2, Teuchos::null);
  eval_vec[0] = GState().GetDisNp();
  eval_vec[1] = Teuchos::rcpFromRef(x);
  EvalContact().SetActionType(MORTAR::eval_reset);

  // reset displacement state and lagrange multiplier values
  Strategy().Evaluate(EvalData().Contact(), &eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::EvaluateForce()
{
  CheckInitSetup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force);
  EvalData().SetModelEvaluator(this);
  Strategy().Evaluate(EvalData().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::EvaluateStiff()
{
  CheckInitSetup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force_stiff);
  EvalData().SetModelEvaluator(this);
  Strategy().Evaluate(EvalData().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::EvaluateForceStiff()
{
  CheckInitSetup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force_stiff);
  Strategy().Evaluate(EvalData().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PreEvaluate()
{
  EvalContact().SetActionType(MORTAR::eval_run_pre_evaluate);
  EvalData().SetModelEvaluator(this);
  Strategy().Evaluate(EvalContact());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostEvaluate()
{
  EvalContact().SetActionType(MORTAR::eval_run_post_evaluate);
  EvalData().SetModelEvaluator(this);
  Strategy().Evaluate(EvalData().Contact());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::AssembleForce(Epetra_Vector& f, const double& timefac_np) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(Strategy().Params(), "ALGORITHM") ==
      INPAR::MORTAR::algorithm_gpts)
  {
    block_vec_ptr = Strategy().GetRhsBlockPtr(DRT::UTILS::block_displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;
    LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }
  else if (Strategy().IsCondensedSystem())
  {
    block_vec_ptr = Strategy().GetCondensedRhsPtr(f, timefac_np);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;

    LINALG::AssembleMyVector(1.0, f, 1.0, *block_vec_ptr);
  }
  else if (Strategy().IsSaddlePointSystem())
  {
    // --- displ. - block ---------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(DRT::UTILS::block_displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;
    LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);

    // --- constr. - block --------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(DRT::UTILS::block_constraint);
    if (block_vec_ptr.is_null())
      dserror(
          "The constraint vector is a NULL pointer, although \n"
          "the structural part indicates, that contact contributions \n"
          "are present!");
    LINALG::AssembleMyVector(1.0, f, 1.0, *block_vec_ptr);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<LINALG::SparseMatrix> block_ptr = Teuchos::null;
  int err = 0;
  // ---------------------------------------------------------------------
  // gpts / Nitsche system: no additional/condensed dofs
  // ---------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(Strategy().Params(), "ALGORITHM") ==
      INPAR::MORTAR::algorithm_gpts)
  {
    block_ptr = Strategy().GetMatrixBlockPtr(DRT::UTILS::block_displ_displ, &EvalContact());
    Teuchos::RCP<LINALG::SparseMatrix> jac_dd = GState().ExtractDisplBlock(jac);
    jac_dd->Add(*block_ptr, false, timefac_np, 1.0);
  }
  // ---------------------------------------------------------------------
  // condensed system of equations
  // ---------------------------------------------------------------------
  else if (Strategy().IsCondensedSystem())
  {
    Teuchos::RCP<LINALG::SparseMatrix> jac_dd = GState().ExtractDisplBlock(jac);

    block_ptr = Strategy().GetCondensedMatrixBlockPtr(jac_dd, timefac_np);
    // if there are no active contact contributions, we can skip this...
    if (block_ptr.is_null()) return (err == 0);

    // here we should hand in the jac_dd matrix and modify it
    jac_dd->Add(*block_ptr, false, 1.0, 0.0);
  }
  // ---------------------------------------------------------------------
  // saddle-point system of equations or no contact contributions
  // ---------------------------------------------------------------------
  else if (Strategy().SystemType() == INPAR::CONTACT::system_saddlepoint)
  {
    // --- Kdd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(DRT::UTILS::block_displ_displ, &EvalContact());
    if (not block_ptr.is_null())
    {
      Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
      jac_dd_ptr->Add(*block_ptr, false, timefac_np, 1.0);
      // reset the block pointers, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kdz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(DRT::UTILS::block_displ_lm, &EvalContact());
    if (not block_ptr.is_null())
    {
      block_ptr->Scale(timefac_np);
      GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::block_displ_lm);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(DRT::UTILS::block_lm_displ, &EvalContact());
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::block_lm_displ);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(DRT::UTILS::block_lm_lm, &EvalContact());
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::block_lm_lm);
    }
    /* if there are no active contact contributions, we put a identity
     * matrix at the (lm,lm)-block */
    else
    {
      Teuchos::RCP<Epetra_Vector> ones =
          Teuchos::rcp(new Epetra_Vector(GState().BlockMap(Type()), false));
      err = ones->PutScalar(1.0);
      block_ptr = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
      GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::block_lm_lm);
    }
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
  }

  return (err == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // clear cache of maps due to varying vector size
  iowriter.ClearMapCache();

  // quantities to be written for restart
  std::map<std::string, Teuchos::RCP<Epetra_Vector>> restart_vectors;

  Strategy().DoWriteRestart(restart_vectors, forced_writerestart);

  // write all vectors specified by used strategy
  std::map<std::string, Teuchos::RCP<Epetra_Vector>>::const_iterator p;
  for (p = restart_vectors.begin(); p != restart_vectors.end(); ++p)
    iowriter.WriteVector(p->first, p->second);

  /* ToDo Move this stuff into the DoWriteRestart() routine of the
   * AbstractStrategy as soon as the old structural time integration
   * is gone! */
  if (Strategy().GetLagrMultN(true) != Teuchos::null)
    iowriter.WriteVector("lagrmultold", Strategy().GetLagrMultN(true));

  // since the global OutputStepState() routine is not called, if the
  // restart is written, we have to do it here manually.
  OutputStepState(iowriter);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::ReadRestart(IO::DiscretizationReader& ioreader)
{
  EvalContact().SetActionType(MORTAR::eval_force_stiff);
  // reader strategy specific stuff
  Strategy().DoReadRestart(ioreader, GState().GetDisN(), EvalData().ContactPtr());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepState(const double& timefac_n)
{
  // add the contact forces to the old structural residual state
  // vector
  Teuchos::RCP<const Epetra_Vector> strcontactrhs_ptr =
      Strategy().GetRhsBlockPtr(DRT::UTILS::block_displ);
  if (not strcontactrhs_ptr.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetMutableFstructureOld();
    fstructold_ptr->Update(timefac_n, *strcontactrhs_ptr, 1.0);
  }

  /* Note: DisN() and DisNp() have the same value at this stage, since
   * we call the structural model evaluator always in first place! */
  strategy_ptr_->Update(GState().GetDisN());

  PostUpdateStepState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostUpdateStepState()
{
  // initialize integration time for time measurement
  strategy_ptr_->Inttime_init();

  // redistribute contact
  strategy_ptr_->RedistributeContact(GState().GetDisN());

  // initialize binning strategy for new time step
  strategy_ptr_->InitBinStrategyforTimestep(GState().GetVelN());

  // setup the map extractor, since redistribute calls FillComplete
  // on the structural discretization. Though this only changes the
  // ghosted dofs while keeping row distribution fixed, the map pointers
  // in the global state are no longer valid. So we reset them.
  Int().ModelEval().SetupMultiMapExtractor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepElement()
{ /* empty */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RunPostComputeX(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  CheckInitSetup();

  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(3, Teuchos::null);
  eval_vec[0] = Teuchos::rcpFromRef(xold);
  eval_vec[1] = Teuchos::rcpFromRef(dir);
  eval_vec[2] = Teuchos::rcpFromRef(xnew);

  EvalContact().SetActionType(MORTAR::eval_run_post_compute_x);

  Strategy().Evaluate(EvalData().Contact(), &eval_vec);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineStressStrain()
{
  // evaluate contact tractions
  Strategy().OutputStresses();

  if (Strategy().WeightedWear())
  {
    /* *******************************************************************
     * We do not compute the non-weighted wear here. we just write
     * the output. the non-weighted wear will be used as dirichlet-b.
     * for the ale problem. n.w.wear will be called in
     * stru_ale_algorithm.cpp and computed in strategy.OutputWear();
     *                                                         farah 06/13
     * ******************************************************************/
    // evaluate wear (not weighted)
    Strategy().ContactWear()->PutScalar(0.0);
    Strategy().OutputWear();
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineEnergy() { return; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineOptionalQuantity() { return; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  // no output in nitsche Strategy
  if (Strategy().IsNitsche()) return;

  // *********************************************************************
  // print active set
  // *********************************************************************
  Strategy().PrintActiveSet();

  // *********************************************************************
  // active contact set and slip set
  // *********************************************************************

  // evaluate active set and slip set
  Teuchos::RCP<Epetra_Vector> activeset =
      Teuchos::rcp(new Epetra_Vector(*Strategy().ActiveRowNodes()));
  activeset->PutScalar(1.0);
  if (Strategy().Friction())
  {
    Teuchos::RCP<Epetra_Vector> slipset =
        Teuchos::rcp(new Epetra_Vector(*Strategy().SlipRowNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp =
        Teuchos::rcp(new Epetra_Vector(*Strategy().ActiveRowNodes()));
    LINALG::Export(*slipset, *slipsetexp);
    activeset->Update(1.0, *slipsetexp, 1.0);
  }

  // export to problem node row map
  Teuchos::RCP<const Epetra_Map> problemnodes = Strategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
  LINALG::Export(*activeset, *activesetexp);

  if (Strategy().WearBothDiscrete())
  {
    Teuchos::RCP<Epetra_Vector> mactiveset =
        Teuchos::rcp(new Epetra_Vector(*Strategy().MasterActiveNodes()));
    mactiveset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipset =
        Teuchos::rcp(new Epetra_Vector(*Strategy().MasterSlipNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp =
        Teuchos::rcp(new Epetra_Vector(*Strategy().MasterActiveNodes()));
    LINALG::Export(*slipset, *slipsetexp);
    mactiveset->Update(1.0, *slipsetexp, 1.0);

    Teuchos::RCP<Epetra_Vector> mactivesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
    LINALG::Export(*mactiveset, *mactivesetexp);
    activesetexp->Update(1.0, *mactivesetexp, 1.0);
  }

  iowriter.WriteVector("activeset", activesetexp, IO::nodevector);

  // *********************************************************************
  // contact tractions
  // *********************************************************************

  // export to problem dof row map
  Teuchos::RCP<const Epetra_Map> problemdofs = Strategy().ProblemDofs();

  // normal direction
  Teuchos::RCP<const Epetra_Vector> normalstresses = Strategy().ContactNorStress();
  Teuchos::RCP<Epetra_Vector> normalstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  Teuchos::RCP<const Epetra_Vector> tangentialstresses = Strategy().ContactTanStress();
  Teuchos::RCP<Epetra_Vector> tangentialstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*tangentialstresses, *tangentialstressesexp);

  // write to output
  // contact tractions in normal and tangential direction
  iowriter.WriteVector("norcontactstress", normalstressesexp);
  iowriter.WriteVector("tancontactstress", tangentialstressesexp);

#ifdef CONTACTFORCEOUTPUT
  dserror("Untested in the new structural framework!");
  // *********************************************************************
  // contact forces on slave non master side,
  // in normal and tangential direction
  // *********************************************************************
  // vectors for contact forces
  Teuchos::RCP<Epetra_Vector> fcslavenor =
      Teuchos::rcp(new Epetra_Vector(Strategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcslavetan =
      Teuchos::rcp(new Epetra_Vector(Strategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmasternor =
      Teuchos::rcp(new Epetra_Vector(Strategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertan =
      Teuchos::rcp(new Epetra_Vector(Strategy().MMatrix()->DomainMap()));

  // vectors with problem dof row map
  Teuchos::RCP<Epetra_Vector> fcslavenorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcslavetanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmasternorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmastertanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));

  // multiplication
  Strategy().DMatrix()->Multiply(true, *normalstresses, *fcslavenor);
  Strategy().DMatrix()->Multiply(true, *tangentialstresses, *fcslavetan);
  Strategy().MMatrix()->Multiply(true, *normalstresses, *fcmasternor);
  Strategy().MMatrix()->Multiply(true, *tangentialstresses, *fcmastertan);

#ifdef MASTERNODESINCONTACT
  // BEGIN: to output the global ID's of the master nodes in contact - devaal 02.2011

  int dim = DRT::Problem::Instance()->NDim();

  if (dim == 2) dserror("Only working for 3D");

  std::vector<int> lnid, gnid;

  // std::cout << "MasterNor" << fcmasternor->MyLength() << std::endl;

  for (int i = 0; i < fcmasternor->MyLength(); i = i + 3)
  {
    // check if master node in contact
    if (sqrt(((*fcmasternor)[i]) * ((*fcmasternor)[i]) +
             ((*fcmasternor)[i + 1]) * ((*fcmasternor)[i + 1]) +
             ((*fcmasternor)[i + 2]) * ((*fcmasternor)[i] + 2)) > 0.00001)
    {
      lnid.push_back((fcmasternor->Map()).GID(i) / 3);
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid, gnid, (int)allproc.size(), &allproc[0], Comm());

  // std::cout << " size of gnid:" << gnid.size() << std::endl;

  ////////////////
  ///// attempt at obtaining the nid and relative displacement u of master nodes in contact - devaal
  // define my own interface
  MORTAR::StrategyBase& myStrategy = strategy;
  CoAbstractStrategy& myContactStrategy = dynamic_cast<CoAbstractStrategy&>(myStrategy);

  std::vector<Teuchos::RCP<CONTACT::CoInterface>> myInterface = Strategy().ContactInterfaces();

  // check interface size - just doing this now for a single interface

  if (myInterface.size() != 1) dserror("Interface size should be 1");

  std::cout << "OUTPUT OF MASTER NODE IN CONTACT" << std::endl;
  // std::cout << "Master_node_in_contact x_dis y_dis z_dis" << std::endl;
  for (int i = 0; i < (int)gnid.size(); ++i)
  {
    int myGid = gnid[i];
    std::cout << gnid[i]
              << std::endl;  // << " " << myUx << " " << myUy << " " << myUz << std::endl;
  }

#endif  // MASTERNODESINCONTACT: to output the global ID's of the master nodes in contact
  // export
  LINALG::Export(*fcslavenor, *fcslavenorexp);
  LINALG::Export(*fcslavetan, *fcslavetanexp);
  LINALG::Export(*fcmasternor, *fcmasternorexp);
  LINALG::Export(*fcmastertan, *fcmastertanexp);

  // contact forces on slave and master side
  iowriter.WriteVector("norslaveforce", fcslavenorexp);
  iowriter.WriteVector("tanslaveforce", fcslavetanexp);
  iowriter.WriteVector("normasterforce", fcmasternorexp);
  iowriter.WriteVector("tanmasterforce", fcmastertanexp);

#ifdef CONTACTEXPORT
  // export averaged node forces to xxx.force
  double resultnor[fcslavenor->NumVectors()];
  double resulttan[fcslavetan->NumVectors()];
  fcslavenor->Norm2(resultnor);
  fcslavetan->Norm2(resulttan);

  if (Comm().MyPID() == 0)
  {
    std::cout << "resultnor= " << resultnor[0] << std::endl;
    std::cout << "resulttan= " << resulttan[0] << std::endl;

    FILE* MyFile = NULL;
    std::ostringstream filename;
    const std::string filebase =
        DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
    filename << filebase << ".force";
    MyFile = fopen(filename.str().c_str(), "at+");
    if (MyFile)
    {
      // fprintf(MyFile,valuename.c_str());
      fprintf(MyFile, "%g\t", resultnor[0]);
      fprintf(MyFile, "%g\n", resulttan[0]);
      fclose(MyFile);
    }
    else
      dserror("ERROR: File for Output could not be opened.");
  }
#endif  // CONTACTEXPORT
#endif  // CONTACTFORCEOUTPUT

  // *********************************************************************
  // wear with internal state variable approach
  // *********************************************************************
  if (Strategy().WeightedWear())
  {
    // write output
    Teuchos::RCP<const Epetra_Vector> wearoutput = Strategy().ContactWear();
    Teuchos::RCP<Epetra_Vector> wearoutputexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*wearoutput, *wearoutputexp);
    iowriter.WriteVector("wear", wearoutputexp);
  }

  // *********************************************************************
  // poro contact
  // *********************************************************************
  if (Strategy().HasPoroNoPenetration())
  {
    // output of poro no penetration lagrange multiplier!
    const CONTACT::PoroLagrangeStrategy& poro_strategy =
        dynamic_cast<const CONTACT::PoroLagrangeStrategy&>(Strategy());
    Teuchos::RCP<const Epetra_Vector> lambdaout = poro_strategy.LambdaNoPen();
    Teuchos::RCP<Epetra_Vector> lambdaoutexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*lambdaout, *lambdaoutexp);
    iowriter.WriteVector("poronopen_lambda", lambdaoutexp);
  }

  /// general way to write the output corresponding to the active strategy
  Strategy().WriteOutput(iowriter);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::ContactData& STR::MODELEVALUATOR::Contact::EvalContact()
{
  CheckInitSetup();
  return *eval_contact_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const STR::MODELEVALUATOR::ContactData& STR::MODELEVALUATOR::Contact::EvalContact() const
{
  CheckInitSetup();
  return *eval_contact_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<CONTACT::CoAbstractStrategy>& STR::MODELEVALUATOR::Contact::StrategyPtr()
{
  CheckInitSetup();
  return strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy()
{
  CheckInitSetup();
  return *strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const CONTACT::CoAbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy() const
{
  CheckInitSetup();
  return *strategy_ptr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Contact::GetBlockDofRowMapPtr() const
{
  DRT::Problem* problem = DRT::Problem::Instance();

  CheckInitSetup();
  if (Strategy().LMDoFRowMapPtr(false) == Teuchos::null)
    return GState().DofRowMap();
  else
  {
    enum INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(
        problem->ContactDynamicParams(), "SYSTEM");

    if (systype == INPAR::CONTACT::system_saddlepoint)
      return Strategy().LinSystemLMDoFRowMapPtr();
    else
      return GState().DofRowMap();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::GetCurrentSolutionPtr() const
{
  // TODO: this should be removed!
  DRT::Problem* problem = DRT::Problem::Instance();
  enum INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(
      problem->ContactDynamicParams(), "SYSTEM");
  if (systype == INPAR::CONTACT::system_condensed) return Teuchos::null;

  if (Strategy().GetLagrMultNp(false) != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> curr_lm_ptr =
        Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultNp(false)));
    if (not curr_lm_ptr.is_null()) curr_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));

    ExtendLagrangeMultiplierDomain(curr_lm_ptr);

    return curr_lm_ptr;
  }
  else
    return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::GetLastTimeStepSolutionPtr() const
{
  if (Strategy().GetLagrMultN(false).is_null()) return Teuchos::null;

  Teuchos::RCP<Epetra_Vector> old_lm_ptr =
      Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultN(false)));
  if (not old_lm_ptr.is_null()) old_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));

  ExtendLagrangeMultiplierDomain(old_lm_ptr);

  return old_lm_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::ExtendLagrangeMultiplierDomain(
    Teuchos::RCP<Epetra_Vector>& lm_vec) const
{
  // default case: do nothing
  if (Strategy().LMDoFRowMap(false).NumGlobalElements() ==
      GetBlockDofRowMapPtr()->NumGlobalElements())
    return;

  if (Strategy().LMDoFRowMap(false).NumGlobalElements() <
      GetBlockDofRowMapPtr()->NumGlobalElements())
  {
    Teuchos::RCP<Epetra_Vector> tmp_ptr = Teuchos::rcp(new Epetra_Vector(*GetBlockDofRowMapPtr()));
    LINALG::Export(*lm_vec, *tmp_ptr);
    lm_vec = tmp_ptr;
  }
  else
    dserror("Unconsidered case.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostOutput()
{
  CheckInitSetup();
  // empty

  return;
}  // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RunPreComputeX(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  CheckInitSetup();

  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(1, Teuchos::null);
  eval_vec[0] = Teuchos::rcpFromRef(xold);

  std::vector<Teuchos::RCP<Epetra_Vector>> eval_vec_mutable(1, Teuchos::null);
  eval_vec_mutable[0] = Teuchos::rcpFromRef(dir_mutable);

  EvalContact().SetActionType(MORTAR::eval_run_pre_compute_x);
  EvalData().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(EvalData().Contact(), &eval_vec, &eval_vec_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RunPostIterate(const NOX::Solver::Generic& solver)
{
  CheckInitSetup();

  EvalContact().SetActionType(MORTAR::eval_run_post_iterate);
  Strategy().Evaluate(EvalData().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RunPostApplyJacobianInverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  CheckInitSetup();

  EvalContact().Set(&rhs, 0);
  EvalContact().Set(&result, 1);
  EvalContact().Set(&xold, 2);
  EvalContact().Set(&grp, 3);

  EvalContact().SetActionType(MORTAR::eval_run_post_apply_jacobian_inverse);
  EvalData().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(EvalData().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> STR::MODELEVALUATOR::Contact::GetJacobianBlock(
    const DRT::UTILS::MatBlockType bt) const
{
  return GState().GetJacobianBlock(Type(), bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Contact::AssembleForceOfModels(
    const std::vector<INPAR::STR::ModelType>* without_these_models, const bool apply_dbc) const
{
  Teuchos::RCP<NOX::Epetra::Vector> force_nox = GState().CreateGlobalVector();
  Int().AssembleForce(force_nox->getEpetraVector(), without_these_models);

  // copy the vector, otherwise the storage will be freed at the end of this
  // function, resulting in a segmentation fault
  Teuchos::RCP<Epetra_Vector> force = Teuchos::rcp(new Epetra_Vector(force_nox->getEpetraVector()));

  if (apply_dbc) TimInt().GetDBC().ApplyDirichletToRhs(force);

  return force;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> STR::MODELEVALUATOR::Contact::GetAuxDisplJacobian() const
{
  std::vector<INPAR::STR::ModelType> g;
  g.push_back(INPAR::STR::ModelType::model_contact);

  Teuchos::RCP<LINALG::SparseOperator> jacaux = GState().CreateAuxJacobian();
  bool ok = Int().AssembleJac(*jacaux, &g);

  if (!ok) dserror("ERROR: CreateAuxJacobian went wrong!");

  return jacaux;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::EvaluateWeightedGapGradientError()
{
  EvalContact().SetActionType(MORTAR::eval_wgap_gradient_error);
  EvalData().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(EvalData().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::EvaluateCheapSOCRhs()
{
  CheckInitSetup();

  EvalContact().SetActionType(MORTAR::eval_static_constraint_rhs);
  EvalData().SetModelEvaluator(this);

  Strategy().Evaluate(EvalData().Contact());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::AssembleCheapSOCRhs(
    Epetra_Vector& f, const double& timefac_np) const
{
  return AssembleForce(f, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::CorrectParameters(NOX::NLN::CorrectionType type)
{
  CheckInitSetup();

  EvalContact().SetActionType(MORTAR::eval_correct_parameters);
  EvalContact().Set(&type, 0);

  Strategy().Evaluate(EvalContact());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs)
{
  CheckInitSetup();
  EvalContact().SetActionType(MORTAR::remove_condensed_contributions_from_str_rhs);

  std::vector<Teuchos::RCP<Epetra_Vector>> mutable_vec(1, Teuchos::rcpFromRef(rhs));
  Strategy().Evaluate(EvalContact(), NULL, &mutable_vec);
}
