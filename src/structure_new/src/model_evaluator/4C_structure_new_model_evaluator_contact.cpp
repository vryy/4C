/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all contact terms


\level 3
*/
/*---------------------------------------------------------------------*/

#include "4C_structure_new_model_evaluator_contact.hpp"

#include "4C_contact_aug_plot.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_lagrange_strategy_poro.hpp"
#include "4C_contact_strategy_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_solver_nonlin_nox_solver_linesearchbased.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Setup()
{
  check_init();
  eval_contact_ptr_ = eval_data().ContactPtr();

  // ---------------------------------------------------------------------
  // create the contact factory
  // ---------------------------------------------------------------------
  CONTACT::STRATEGY::Factory factory;
  factory.Init(g_state_ptr());
  factory.Setup();

  // check the problem dimension
  factory.CheckDimension();

  // create some local variables (later to be stored in strategy)
  std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  factory.read_and_check_input(cparams);

  // check for fill_complete of discretization
  if (not discret().Filled()) FOUR_C_THROW("discretization is not fillcomplete");

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
  eval_contact_ptr_->Set(&discret(), 0);
  strategy_ptr_ = factory.BuildStrategy(
      cparams, poroslave, poromaster, dof_offset(), interfaces, eval_contact_ptr_.get());
  eval_contact_ptr_->ClearEntry(CORE::GEN::AnyDataContainer::DataType::any, 0);

  // build the search tree
  factory.BuildSearchTree(interfaces);

  // print final screen output
  factory.Print(interfaces, strategy_ptr_, cparams);

  // ---------------------------------------------------------------------
  // final touches to the contact strategy
  // ---------------------------------------------------------------------
  strategy_ptr_->store_dirichlet_status(Int().GetDbc().GetDBCMapExtractor());
  strategy_ptr_->set_state(MORTAR::state_new_displacement, Int().GetDbc().GetZeros());
  strategy_ptr_->SaveReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->evaluate_reference_state();
  strategy_ptr_->Inttime_init();
  set_time_integration_info(*strategy_ptr_);
  strategy_ptr_->RedistributeContact(g_state().GetDisN(), g_state().GetVelN());

  check_pseudo2_d();

  post_setup(cparams);

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::post_setup(Teuchos::ParameterList& cparams)
{
  if (dynamic_cast<CONTACT::AUG::Strategy*>(strategy_ptr_.get()))
  {
    Teuchos::ParameterList& aug_params = cparams.sublist("AUGMENTED");
    Teuchos::ParameterList& plot_params = aug_params.sublist("PLOT");
    plot_params.set<const int*>("CURRENT_STEP", &g_state().GetStepNp());
    plot_params.set<std::string>(
        "OUTPUT_FILE_NAME", g_in_output().GetOutputPtr()->Output()->FileName());
    plot_params.set<std::string>(
        "INPUT_FILE_NAME", g_in_output().GetOutputPtr()->Output()->InputFileName());
    plot_params.set<const DRT::Discretization*>("DISCRETIZATION", g_state().GetDiscret().get());
    plot_params.set<STR::MODELEVALUATOR::Contact*>("MODELEVALUATOR", this);

    STR::IMPLICIT::Generic& impl = dynamic_cast<STR::IMPLICIT::Generic&>(integrator());
    CONTACT::AUG::Plot::Create(impl.GetNoxParams(), plot_params, strategy_ptr_.get());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::check_pseudo2_d() const
{
  // print messages for multifield problems (e.g FSI)
  const CORE::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();
  if ((probtype != CORE::ProblemType::structure) and (GState().GetMyRank() == 0))
  {
    // warnings
#ifdef CONTACTPSEUDO2D
    std::cout << "WARNING: The flag CONTACTPSEUDO2D is switched on. If this "
              << "is a real 3D problem, switch it off!" << std::endl;
#else
    std::cout << "STR::MODELEVALUATOR::Contact::check_pseudo2_d -- "
              << "WARNING: \nThe flag CONTACTPSEUDO2D is switched off. If this "
              << "is a 2D problem modeled pseudo-3D, switch it on!" << std::endl;
#endif
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::set_time_integration_info(
    CONTACT::AbstractStrategy& strategy) const
{
  const INPAR::STR::DynamicType dyntype = TimInt().GetDataSDyn().GetDynamicType();
  const double time_fac = Int().GetIntParam();

  strategy.set_time_integration_info(time_fac, dyntype);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Reset(const Epetra_Vector& x)
{
  check_init_setup();
  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(2, Teuchos::null);
  eval_vec[0] = g_state().GetDisNp();
  eval_vec[1] = Teuchos::rcpFromRef(x);
  eval_contact().SetActionType(MORTAR::eval_reset);

  // reset displacement state and lagrange multiplier values
  Strategy().Evaluate(eval_data().Contact(), &eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::evaluate_force()
{
  check_init_setup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  eval_contact().SetActionType(MORTAR::eval_force);
  eval_data().SetModelEvaluator(this);
  Strategy().Evaluate(eval_data().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::evaluate_stiff()
{
  check_init_setup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  eval_contact().SetActionType(MORTAR::eval_force_stiff);
  eval_data().SetModelEvaluator(this);
  Strategy().Evaluate(eval_data().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::evaluate_force_stiff()
{
  check_init_setup();
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  eval_contact().SetActionType(MORTAR::eval_force_stiff);
  Strategy().Evaluate(eval_data().Contact());

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::pre_evaluate()
{
  eval_contact().SetActionType(MORTAR::eval_run_pre_evaluate);
  eval_data().SetModelEvaluator(this);
  Strategy().Evaluate(eval_contact());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::post_evaluate()
{
  eval_contact().SetActionType(MORTAR::eval_run_post_evaluate);
  eval_data().SetModelEvaluator(this);
  Strategy().Evaluate(eval_data().Contact());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::assemble_force(Epetra_Vector& f, const double& timefac_np) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(Strategy().Params(), "ALGORITHM") ==
          INPAR::MORTAR::algorithm_gpts ||
      Strategy().IsPenalty() || Strategy().IsCondensedSystem())
  {
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);

    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;

    CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }
  else
  {
    // --- displ. - block ---------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;
    CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);

    // --- constr. - block --------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::constraint);
    if (block_vec_ptr.is_null()) return true;
    Epetra_Vector tmp(f.Map());
    CORE::LINALG::Export(*block_vec_ptr, tmp);
    f.Update(1., tmp, 1.);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::assemble_jacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> block_ptr = Teuchos::null;
  int err = 0;
  // ---------------------------------------------------------------------
  // Penalty / gpts / Nitsche system: no additional/condensed dofs
  // ---------------------------------------------------------------------
  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(Strategy().Params(), "ALGORITHM") ==
          INPAR::MORTAR::algorithm_gpts ||
      Strategy().IsPenalty())
  {
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ, &eval_contact());
    if (Strategy().IsPenalty() && block_ptr.is_null()) return true;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd = GState().ExtractDisplBlock(jac);
    jac_dd->Add(*block_ptr, false, timefac_np, 1.0);
  }
  // ---------------------------------------------------------------------
  // condensed system of equations
  // ---------------------------------------------------------------------
  else if (Strategy().IsCondensedSystem())
  {
    // --- Kdd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ, &eval_contact());
    if (not block_ptr.is_null())
    {
      Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
      jac_dd_ptr->Add(*block_ptr, false, timefac_np, 1.0);
      // reset the block pointers, just to be on the safe side
      block_ptr = Teuchos::null;
    }
  }
  // ---------------------------------------------------------------------
  // saddle-point system of equations or no contact contributions
  // ---------------------------------------------------------------------
  else if (Strategy().SystemType() == INPAR::CONTACT::system_saddlepoint)
  {
    // --- Kdd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ, &eval_contact());
    if (not block_ptr.is_null())
    {
      Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
      jac_dd_ptr->Add(*block_ptr, false, timefac_np, 1.0);
      // reset the block pointers, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kdz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_lm, &eval_contact());
    if (not block_ptr.is_null())
    {
      block_ptr->Scale(timefac_np);
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::displ_lm);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::lm_displ, &eval_contact());
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::lm_displ);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::lm_lm, &eval_contact());
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::lm_lm);
    }
    /* if there are no active contact contributions, we put a identity
     * matrix at the (lm,lm)-block */
    else
    {
      Teuchos::RCP<Epetra_Vector> ones =
          Teuchos::rcp(new Epetra_Vector(GState().BlockMap(Type()), false));
      err = ones->PutScalar(1.0);
      block_ptr = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*ones));
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::lm_lm);
    }
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
  }

  return (err == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::write_restart(
    CORE::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::read_restart(CORE::IO::DiscretizationReader& ioreader)
{
  eval_contact().SetActionType(MORTAR::eval_force_stiff);
  // reader strategy specific stuff
  Strategy().DoReadRestart(ioreader, g_state().GetDisN(), eval_data().ContactPtr());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepState(const double& timefac_n)
{
  // add the contact forces to the old structural residual state vector
  Teuchos::RCP<const Epetra_Vector> strcontactrhs_ptr =
      Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);
  if (not strcontactrhs_ptr.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr = g_state().GetFstructureOld();
    fstructold_ptr->Update(timefac_n, *strcontactrhs_ptr, 1.0);
  }

  /* Note: DisN() and dis_np() have the same value at this stage, since
   * we call the structural model evaluator always in first place! */
  strategy_ptr_->Update(g_state().GetDisN());

  post_update_step_state();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::post_update_step_state()
{
  // initialize integration time for time measurement
  strategy_ptr_->Inttime_init();

  // redistribute contact
  strategy_ptr_->RedistributeContact(g_state().GetDisN(), g_state().GetVelN());

  // setup the map extractor, since redistribute calls fill_complete
  // on the structural discretization. Though this only changes the
  // ghosted dofs while keeping row distribution fixed, the map pointers
  // in the global state are no longer valid. So we reset them.
  integrator().ModelEval().setup_multi_map_extractor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepElement()
{ /* empty */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  check_init_setup();

  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(3, Teuchos::null);
  eval_vec[0] = Teuchos::rcpFromRef(xold);
  eval_vec[1] = Teuchos::rcpFromRef(dir);
  eval_vec[2] = Teuchos::rcpFromRef(xnew);

  eval_contact().SetActionType(MORTAR::eval_run_post_compute_x);

  Strategy().Evaluate(eval_data().Contact(), &eval_vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::determine_stress_strain()
{
  // evaluate contact tractions
  Strategy().compute_contact_stresses();

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineEnergy() {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::determine_optional_quantity() {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::OutputStepState(CORE::IO::DiscretizationWriter& iowriter) const
{
  // no output in nitsche Strategy
  if (Strategy().IsNitsche()) return;

  // *********************************************************************
  // print summary of active set to screen
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
    CORE::LINALG::Export(*slipset, *slipsetexp);
    activeset->Update(1.0, *slipsetexp, 1.0);
  }

  // export to problem node row map
  Teuchos::RCP<const Epetra_Map> problemnodes = Strategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
  CORE::LINALG::Export(*activeset, *activesetexp);

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
    CORE::LINALG::Export(*slipset, *slipsetexp);
    mactiveset->Update(1.0, *slipsetexp, 1.0);

    Teuchos::RCP<Epetra_Vector> mactivesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
    CORE::LINALG::Export(*mactiveset, *mactivesetexp);
    activesetexp->Update(1.0, *mactivesetexp, 1.0);
  }

  iowriter.WriteVector("activeset", activesetexp, CORE::IO::nodevector);

  // *********************************************************************
  // contact tractions
  // *********************************************************************

  // export to problem dof row map
  Teuchos::RCP<const Epetra_Map> problemdofs = Strategy().ProblemDofs();

  // normal direction
  Teuchos::RCP<const Epetra_Vector> normalstresses = Strategy().ContactNorStress();
  Teuchos::RCP<Epetra_Vector> normalstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  CORE::LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  Teuchos::RCP<const Epetra_Vector> tangentialstresses = Strategy().ContactTanStress();
  Teuchos::RCP<Epetra_Vector> tangentialstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  CORE::LINALG::Export(*tangentialstresses, *tangentialstressesexp);

  // write to output
  // contact tractions in normal and tangential direction
  iowriter.WriteVector("norcontactstress", normalstressesexp);
  iowriter.WriteVector("tancontactstress", tangentialstressesexp);

#ifdef CONTACTFORCEOUTPUT
  FOUR_C_THROW("Untested in the new structural framework!");
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

  int dim = GLOBAL::Problem::Instance()->NDim();

  if (dim == 2) FOUR_C_THROW("Only working for 3D");

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
  CORE::LINALG::Gather<int>(lnid, gnid, (int)allproc.size(), allproc.data(), Comm());

  // std::cout << " size of gnid:" << gnid.size() << std::endl;

  ////////////////
  ///// attempt at obtaining the nid and relative displacement u of master nodes in contact - devaal
  // define my own interface
  MORTAR::StrategyBase& myStrategy = strategy;
  AbstractStrategy& myContactStrategy = dynamic_cast<AbstractStrategy&>(myStrategy);

  std::vector<Teuchos::RCP<CONTACT::Interface>> myInterface = Strategy().ContactInterfaces();

  // check interface size - just doing this now for a single interface

  if (myInterface.size() != 1) FOUR_C_THROW("Interface size should be 1");

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
  CORE::LINALG::Export(*fcslavenor, *fcslavenorexp);
  CORE::LINALG::Export(*fcslavetan, *fcslavetanexp);
  CORE::LINALG::Export(*fcmasternor, *fcmasternorexp);
  CORE::LINALG::Export(*fcmastertan, *fcmastertanexp);

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

    FILE* MyFile = nullptr;
    std::ostringstream filename;
    const std::string filebase =
        GLOBAL::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
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
      FOUR_C_THROW("ERROR: File for Output could not be opened.");
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
    CORE::LINALG::Export(*wearoutput, *wearoutputexp);
    iowriter.WriteVector("wear", wearoutputexp);
  }

  // *********************************************************************
  // poro contact
  // *********************************************************************
  if (Strategy().has_poro_no_penetration())
  {
    // output of poro no penetration lagrange multiplier!
    const CONTACT::LagrangeStrategyPoro& poro_strategy =
        dynamic_cast<const CONTACT::LagrangeStrategyPoro&>(Strategy());
    Teuchos::RCP<const Epetra_Vector> lambdaout = poro_strategy.LambdaNoPen();
    Teuchos::RCP<Epetra_Vector> lambdaoutexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    CORE::LINALG::Export(*lambdaout, *lambdaoutexp);
    iowriter.WriteVector("poronopen_lambda", lambdaoutexp);
  }

  /// general way to write the output corresponding to the active strategy
  Strategy().WriteOutput(iowriter);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::ResetStepState()
{
  check_init_setup();

  FOUR_C_THROW("Not yet implemented");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::ContactData& STR::MODELEVALUATOR::Contact::eval_contact()
{
  check_init_setup();
  return *eval_contact_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const STR::MODELEVALUATOR::ContactData& STR::MODELEVALUATOR::Contact::eval_contact() const
{
  check_init_setup();
  return *eval_contact_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<CONTACT::AbstractStrategy>& STR::MODELEVALUATOR::Contact::StrategyPtr()
{
  check_init_setup();
  return strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy()
{
  check_init_setup();
  return *strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const CONTACT::AbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy() const
{
  check_init_setup();
  return *strategy_ptr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Contact::get_block_dof_row_map_ptr() const
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  check_init_setup();
  if (Strategy().LMDoFRowMapPtr(false) == Teuchos::null)
    return GState().dof_row_map();
  else
  {
    enum INPAR::CONTACT::SystemType systype =
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(
            problem->contact_dynamic_params(), "SYSTEM");

    if (systype == INPAR::CONTACT::system_saddlepoint)
      return Strategy().lin_system_lm_do_f_row_map_ptr();
    else
      return GState().dof_row_map();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::get_current_solution_ptr() const
{
  // TODO: this should be removed!
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  enum INPAR::CONTACT::SystemType systype = CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(
      problem->contact_dynamic_params(), "SYSTEM");
  if (systype == INPAR::CONTACT::system_condensed) return Teuchos::null;

  if (Strategy().GetLagrMultNp(false) != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> curr_lm_ptr =
        Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultNp(false)));
    if (not curr_lm_ptr.is_null()) curr_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));

    extend_lagrange_multiplier_domain(curr_lm_ptr);

    return curr_lm_ptr;
  }
  else
    return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::get_last_time_step_solution_ptr()
    const
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  enum INPAR::CONTACT::SystemType systype = CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(
      problem->contact_dynamic_params(), "SYSTEM");
  if (systype == INPAR::CONTACT::system_condensed) return Teuchos::null;

  if (Strategy().GetLagrMultN(false).is_null()) return Teuchos::null;

  Teuchos::RCP<Epetra_Vector> old_lm_ptr =
      Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultN(false)));
  if (not old_lm_ptr.is_null()) old_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));

  extend_lagrange_multiplier_domain(old_lm_ptr);

  return old_lm_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::extend_lagrange_multiplier_domain(
    Teuchos::RCP<Epetra_Vector>& lm_vec) const
{
  // default case: do nothing
  if (Strategy().LMDoFRowMap(false).NumGlobalElements() ==
      get_block_dof_row_map_ptr()->NumGlobalElements())
    return;

  if (Strategy().LMDoFRowMap(false).NumGlobalElements() <
      get_block_dof_row_map_ptr()->NumGlobalElements())
  {
    Teuchos::RCP<Epetra_Vector> tmp_ptr =
        Teuchos::rcp(new Epetra_Vector(*get_block_dof_row_map_ptr()));
    CORE::LINALG::Export(*lm_vec, *tmp_ptr);
    lm_vec = tmp_ptr;
  }
  else
    FOUR_C_THROW("Unconsidered case.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostOutput()
{
  check_init_setup();
  // empty
}  // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::run_pre_compute_x(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  check_init_setup();

  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(1, Teuchos::null);
  eval_vec[0] = Teuchos::rcpFromRef(xold);

  std::vector<Teuchos::RCP<Epetra_Vector>> eval_vec_mutable(1, Teuchos::null);
  eval_vec_mutable[0] = Teuchos::rcpFromRef(dir_mutable);

  eval_contact().SetActionType(MORTAR::eval_run_pre_compute_x);
  eval_data().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(eval_data().Contact(), &eval_vec, &eval_vec_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::run_post_iterate(const ::NOX::Solver::Generic& solver)
{
  check_init_setup();

  const ::NOX::Epetra::Vector& nox_x =
      dynamic_cast<const ::NOX::Epetra::Vector&>(solver.getSolutionGroup().getX());

  // displacement vector after the predictor call
  Teuchos::RCP<Epetra_Vector> curr_disp = g_state().ExtractDisplEntries(nox_x.getEpetraVector());
  Teuchos::RCP<const Epetra_Vector> curr_vel = g_state().GetVelNp();

  if (Strategy().dyn_redistribute_contact(curr_disp, curr_vel, solver.getNumIterations()))
    evaluate_force();

  eval_contact().SetActionType(MORTAR::eval_run_post_iterate);
  Strategy().Evaluate(eval_data().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd = g_state().JacobianDisplBlock();
  const_cast<CONTACT::AbstractStrategy&>(Strategy())
      .run_pre_apply_jacobian_inverse(jac_dd, const_cast<Epetra_Vector&>(rhs));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RunPreSolve(const ::NOX::Solver::Generic& solver)
{
  check_init_setup();
  const ::NOX::Epetra::Vector& nox_x =
      dynamic_cast<const ::NOX::Epetra::Vector&>(solver.getSolutionGroup().getX());

  // displacement vector after the predictor call
  Teuchos::RCP<Epetra_Vector> curr_disp = g_state().ExtractDisplEntries(nox_x.getEpetraVector());

  std::vector<Teuchos::RCP<const Epetra_Vector>> eval_vec(1, Teuchos::null);
  eval_vec[0] = curr_disp;

  eval_contact().SetActionType(MORTAR::eval_run_pre_solve);
  Strategy().Evaluate(eval_data().Contact(), &eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::run_post_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  check_init_setup();

  eval_contact().Set(&rhs, 0);
  eval_contact().Set(&result, 1);
  eval_contact().Set(&xold, 2);
  eval_contact().Set(&grp, 3);

  eval_contact().SetActionType(MORTAR::eval_run_post_apply_jacobian_inverse);
  eval_data().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(eval_data().Contact());

  // clear the set any data again
  eval_contact().ClearAll(CORE::GEN::AnyDataContainer::DataType::any);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix> STR::MODELEVALUATOR::Contact::GetJacobianBlock(
    const STR::MatBlockType bt) const
{
  return GState().GetJacobianBlock(Type(), bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Contact::assemble_force_of_models(
    const std::vector<INPAR::STR::ModelType>* without_these_models, const bool apply_dbc) const
{
  Teuchos::RCP<::NOX::Epetra::Vector> force_nox = GState().CreateGlobalVector();
  Int().assemble_force(force_nox->getEpetraVector(), without_these_models);

  // copy the vector, otherwise the storage will be freed at the end of this
  // function, resulting in a segmentation fault
  Teuchos::RCP<Epetra_Vector> force = Teuchos::rcp(new Epetra_Vector(force_nox->getEpetraVector()));

  if (apply_dbc) TimInt().GetDBC().ApplyDirichletToRhs(force);

  return force;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseOperator> STR::MODELEVALUATOR::Contact::GetAuxDisplJacobian() const
{
  std::vector<INPAR::STR::ModelType> g;
  g.push_back(INPAR::STR::ModelType::model_contact);

  Teuchos::RCP<CORE::LINALG::SparseOperator> jacaux = GState().CreateAuxJacobian();
  bool ok = Int().AssembleJac(*jacaux, &g);

  if (!ok) FOUR_C_THROW("ERROR: CreateAuxJacobian went wrong!");

  return jacaux;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::evaluate_weighted_gap_gradient_error()
{
  eval_contact().SetActionType(MORTAR::eval_wgap_gradient_error);
  eval_data().SetModelEvaluator(this);

  // augment the search direction
  Strategy().Evaluate(eval_data().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::evaluate_cheap_soc_rhs()
{
  check_init_setup();

  eval_contact().SetActionType(MORTAR::eval_static_constraint_rhs);
  eval_data().SetModelEvaluator(this);

  Strategy().Evaluate(eval_data().Contact());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::assemble_cheap_soc_rhs(
    Epetra_Vector& f, const double& timefac_np) const
{
  return assemble_force(f, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::correct_parameters(NOX::NLN::CorrectionType type)
{
  check_init_setup();

  eval_contact().SetActionType(MORTAR::eval_correct_parameters);
  eval_contact().Set(&type, 0);

  Strategy().Evaluate(eval_contact());

  eval_contact().ClearEntry(CORE::GEN::AnyDataContainer::DataType::any, 0);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::remove_condensed_contributions_from_rhs(Epetra_Vector& rhs)
{
  check_init_setup();
  eval_contact().SetActionType(MORTAR::remove_condensed_contributions_from_str_rhs);

  std::vector<Teuchos::RCP<Epetra_Vector>> mutable_vec(1, Teuchos::rcpFromRef(rhs));
  Strategy().Evaluate(eval_contact(), nullptr, &mutable_vec);
}

FOUR_C_NAMESPACE_CLOSE
