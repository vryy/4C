/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all contact terms


\level 3

*/
/*---------------------------------------------------------------------*/


#include "4C_structure_new_model_evaluator_meshtying.hpp"

#include "4C_contact_aug_plot.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_lagrange_strategy_poro.hpp"
#include "4C_contact_meshtying_abstract_strategy.hpp"
#include "4C_contact_meshtying_noxinterface.hpp"
#include "4C_contact_meshtying_strategy_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_solver_nonlin_nox_solver_linesearchbased.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Meshtying::Meshtying() : strategy_ptr_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset)
{
  STR::MODELEVALUATOR::Generic::Init(
      eval_data_ptr, gstate_ptr, gio_ptr, int_ptr, timint_ptr, dof_offset);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::Setup()
{
  check_init();

  // create the meshtying factory
  MORTAR::STRATEGY::FactoryMT factory;
  factory.Init(g_state_ptr());
  factory.Setup();

  // check the problem dimension
  factory.CheckDimension();

  // create some local variables (later to be stored in strategy)
  std::vector<Teuchos::RCP<MORTAR::Interface>> interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  factory.read_and_check_input(cparams);

  // check for fill_complete of discretization
  if (not discret().Filled()) FOUR_C_THROW("discretization is not fill_complete.");

  // ---------------------------------------------------------------------
  // build the meshtying interfaces
  // ---------------------------------------------------------------------
  // FixMe Would be great, if we get rid of these poro parameters...
  bool poroslave = false;
  bool poromaster = false;
  factory.BuildInterfaces(cparams, interfaces, poroslave, poromaster);

  // ---------------------------------------------------------------------
  // build the solver strategy object
  // ---------------------------------------------------------------------
  strategy_ptr_ = factory.BuildStrategy(cparams, poroslave, poromaster, dof_offset(), interfaces);

  // build the search tree
  factory.BuildSearchTree(interfaces);

  // ---------------------------------------------------------------------
  // final touches to the meshtying strategy
  // ---------------------------------------------------------------------
  strategy_ptr_->store_dirichlet_status(Int().GetDbc().GetDBCMapExtractor());
  strategy_ptr_->set_state(MORTAR::state_new_displacement, Int().GetDbc().GetZeros());
  strategy_ptr_->SaveReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->evaluate_reference_state();
  strategy_ptr_->Inttime_init();
  set_time_integration_info(*strategy_ptr_);
  strategy_ptr_->RedistributeContact(
      Int().GetDbc().GetZerosPtr(), Int().GetDbc().GetZerosPtr());  // ToDo redistribute_meshtying??
  strategy_ptr_->MortarCoupling(Int().GetDbc().GetZerosPtr());

  strategy_ptr_->nox_interface_ptr()->Init(g_state_ptr());
  strategy_ptr_->nox_interface_ptr()->Setup();

  if (!g_state().GetRestartStep())
  {
    // perform mesh initialization if required by input parameter MESH_RELOCATION
    auto mesh_relocation_parameter = CORE::UTILS::IntegralValue<INPAR::MORTAR::MeshRelocation>(
        GLOBAL::Problem::Instance()->mortar_coupling_params(), "MESH_RELOCATION");

    if (mesh_relocation_parameter == INPAR::MORTAR::relocation_initial)
    {
      Teuchos::RCP<const Epetra_Vector> Xslavemod =
          dynamic_cast<MORTAR::StrategyBase&>(*strategy_ptr_).MeshInitialization();
      if (Xslavemod != Teuchos::null)
      {
        mesh_relocation_ = Teuchos::rcp(new Epetra_Vector(*g_state().dof_row_map()));
        for (int i = 0; i < strategy_ptr_->SlaveRowNodes()->NumMyElements(); ++i)
          for (int d = 0; d < strategy_ptr_->Dim(); ++d)
          {
            int gid = g_state().GetDiscret()->Dof(
                g_state().GetDiscret()->gNode(strategy_ptr_->SlaveRowNodes()->GID(i)), d);
            mesh_relocation_->operator[](mesh_relocation_->Map().LID(gid)) =
                g_state().GetDiscret()->gNode(strategy_ptr_->SlaveRowNodes()->GID(i))->X()[d] -
                Xslavemod->operator[](Xslavemod->Map().LID(gid));
          }
        apply_mesh_initialization(Xslavemod);
      }
    }
    else if (mesh_relocation_parameter == INPAR::MORTAR::relocation_timestep)
    {
      FOUR_C_THROW(
          "Meshtying with MESH_RELOCATION every_timestep not permitted. Change to MESH_RELOCATION "
          "initial or MESH_RELOCATION no.");
    }
  }

  issetup_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Meshtying::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;
  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(Strategy().Params(), "ALGORITHM") ==
          INPAR::MORTAR::algorithm_gpts ||
      Strategy().IsPenalty())
  {
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);
    // if there are no active contact contributions, we can skip this...
    FOUR_C_ASSERT(!block_vec_ptr.is_null(), "force not available");
    CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }
  else if (Strategy().IsCondensedSystem())
  {
    // --- displ. - block ---------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;

    CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }
  else if (Strategy().IsSaddlePointSystem())
  {
    // --- displ. - block ---------------------------------------------------
    block_vec_ptr = Strategy().GetRhsBlockPtr(CONTACT::VecBlockType::displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return true;

    CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Meshtying::assemble_jacobian(
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
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ);
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
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ);
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
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ);
    if (not block_ptr.is_null())
    {
      Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
      jac_dd_ptr->Add(*block_ptr, false, timefac_np, 1.0);
      // reset the block pointers, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kdz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::displ_lm);
    if (not block_ptr.is_null())
    {
      //      block_ptr->Scale(timefac_np);
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::displ_lm);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::lm_displ);
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::lm_displ);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }

    // --- Kzz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(CONTACT::MatBlockType::lm_lm);
    if (not block_ptr.is_null())
    {
      GState().AssignModelBlock(jac, *block_ptr, Type(), STR::MatBlockType::lm_lm);
      // reset the block pointer, just to be on the safe side
      block_ptr = Teuchos::null;
    }
  }

  return (err == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<CONTACT::MtAbstractStrategy>& STR::MODELEVALUATOR::Meshtying::StrategyPtr()
{
  check_init_setup();
  return strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::MtAbstractStrategy& STR::MODELEVALUATOR::Meshtying::Strategy()
{
  check_init_setup();
  return *strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const CONTACT::MtAbstractStrategy& STR::MODELEVALUATOR::Meshtying::Strategy() const
{
  check_init_setup();
  return *strategy_ptr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Meshtying::get_block_dof_row_map_ptr() const
{
  check_init_setup();
  if (Strategy().LMDoFRowMapPtr(true) == Teuchos::null)
    return GState().dof_row_map();
  else
  {
    enum INPAR::CONTACT::SystemType systype =
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(Strategy().Params(), "SYSTEM");

    if (systype == INPAR::CONTACT::system_saddlepoint)
      return Strategy().LMDoFRowMapPtr(true);
    else
      return GState().dof_row_map();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Meshtying::get_current_solution_ptr() const
{
  //  //TODO: this should be removed!
  //  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  //  enum INPAR::CONTACT::SystemType systype =
  //      CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(
  //          problem->contact_dynamic_params(),"SYSTEM");
  //  if (systype == INPAR::CONTACT::system_condensed)
  //    return Teuchos::null;
  //
  //  if (Strategy().GetLagrMultNp(false)!=Teuchos::null)
  //  {
  //    Teuchos::RCP<Epetra_Vector> curr_lm_ptr =
  //        Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultNp(false)));
  //    if (not curr_lm_ptr.is_null())
  //      curr_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));
  //
  //    extend_lagrange_multiplier_domain( curr_lm_ptr );
  //
  //    return curr_lm_ptr;
  //  }
  //  else
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Meshtying::get_last_time_step_solution_ptr()
    const
{
  //  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  //  enum INPAR::CONTACT::SystemType systype =
  //      CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(
  //          problem->contact_dynamic_params(),"SYSTEM");
  //  if (systype == INPAR::CONTACT::system_condensed)
  //    return Teuchos::null;
  //
  //  if (Strategy().GetLagrMultN(false).is_null())
  //    return Teuchos::null;
  //
  //  Teuchos::RCP<Epetra_Vector> old_lm_ptr =
  //      Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultN(false)));
  //  if (not old_lm_ptr.is_null())
  //    old_lm_ptr->ReplaceMap(Strategy().LMDoFRowMap(false));
  //
  //  extend_lagrange_multiplier_domain( old_lm_ptr );
  //
  //  return old_lm_ptr;
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd = g_state().JacobianDisplBlock();
  const_cast<CONTACT::MtAbstractStrategy&>(Strategy())
      .run_pre_apply_jacobian_inverse(jac_dd, const_cast<Epetra_Vector&>(rhs));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::run_post_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  const_cast<CONTACT::MtAbstractStrategy&>(Strategy()).run_post_apply_jacobian_inverse(result);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix> STR::MODELEVALUATOR::Meshtying::GetJacobianBlock(
    const STR::MatBlockType bt) const
{
  return GState().GetJacobianBlock(Type(), bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Meshtying::evaluate_force()
{
  return Strategy().evaluate_force(g_state().GetDisNp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Meshtying::evaluate_force_stiff()
{
  return Strategy().evaluate_force_stiff(g_state().GetDisNp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Meshtying::evaluate_stiff()
{
  return Strategy().evaluate_stiff(g_state().GetDisNp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::apply_mesh_initialization(
    Teuchos::RCP<const Epetra_Vector> Xslavemod)
{
  // check modified positions vector
  if (Xslavemod == Teuchos::null) return;

  // create fully overlapping slave node map
  Teuchos::RCP<Epetra_Map> slavemap = strategy_ptr_->SlaveRowNodes();
  Teuchos::RCP<Epetra_Map> allreduceslavemap = CORE::LINALG::AllreduceEMap(*slavemap);

  // export modified node positions to column map of problem discretization
  const Epetra_Map* dof_colmap = discret_ptr()->DofColMap();
  const Epetra_Map* node_colmap = discret_ptr()->NodeColMap();
  Teuchos::RCP<Epetra_Vector> Xslavemodcol = CORE::LINALG::CreateVector(*dof_colmap, false);
  CORE::LINALG::Export(*Xslavemod, *Xslavemodcol);

  const int numnode = allreduceslavemap->NumMyElements();
  const int numdim = GLOBAL::Problem::Instance()->NDim();
  const Epetra_Vector& gvector = *Xslavemodcol;

  // loop over all slave nodes (for all procs)
  for (int index = 0; index < numnode; ++index)
  {
    int gid = allreduceslavemap->GID(index);

    // only do someting for nodes in my column map
    int ilid = node_colmap->LID(gid);
    if (ilid < 0) continue;

    DRT::Node* mynode = discret_ptr()->gNode(gid);

    // get degrees of freedom associated with this fluid/structure node
    std::vector<int> nodedofs = discret_ptr()->Dof(0, mynode);
    std::vector<double> nvector(3, 0.0);

    // create new position vector
    for (int i = 0; i < numdim; ++i)
    {
      const int lid = gvector.Map().LID(nodedofs[i]);

      if (lid < 0)
        FOUR_C_THROW("ERROR: Proc %d: Cannot find gid=%d in Epetra_Vector", gvector.Comm().MyPID(),
            nodedofs[i]);

      nvector[i] += gvector[lid];
    }

    // set new reference position
    mynode->SetPos(nvector);
  }

  // re-initialize finite elements
  CORE::COMM::ParObjectFactory::Instance().initialize_elements(discret());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  check_init_setup();

  Strategy().run_post_compute_x(xold, dir, xnew);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::remove_condensed_contributions_from_rhs(Epetra_Vector& rhs)
{
  check_init_setup();

  Strategy().remove_condensed_contributions_from_rhs(rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::write_restart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  if (mesh_relocation_ != Teuchos::null)
    iowriter.WriteVector("mesh_relocation", mesh_relocation_);
  else
  {
    Teuchos::RCP<Epetra_Vector> tmp =
        Teuchos::rcp(new Epetra_Vector(*Discret().dof_row_map(), true));
    iowriter.WriteVector("mesh_relocation", tmp);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::read_restart(IO::DiscretizationReader& ioreader)
{
  mesh_relocation_ = Teuchos::rcp(new Epetra_Vector(*discret().dof_row_map(), true));
  ioreader.ReadVector(mesh_relocation_, "mesh_relocation");

  strategy_ptr_->set_state(MORTAR::state_new_displacement, *mesh_relocation_);
  strategy_ptr_->MortarCoupling(mesh_relocation_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Meshtying::set_time_integration_info(
    CONTACT::MtAbstractStrategy& strategy) const
{
  const INPAR::STR::DynamicType dyntype = TimInt().GetDataSDyn().GetDynamicType();
  const double time_fac = Int().GetIntParam();

  strategy.set_time_integration_info(time_fac, dyntype);
}

FOUR_C_NAMESPACE_CLOSE
