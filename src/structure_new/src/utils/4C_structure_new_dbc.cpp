/*-----------------------------------------------------------*/
/*! \file

\brief Wrapper for all Dirichlet boundary condition tasks.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_dbc.hpp"

#include "4C_discretization_condition_locsys.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_linearsystem_prepostoperator.hpp"
#include "4C_structure_new_timint_base.hpp"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Dbc::Dbc()
    : isinit_(false),
      issetup_(false),
      islocsys_(false),
      discret_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null),
      locsysman_ptr_(Teuchos::null),
      zeros_ptr_(Teuchos::null),
      dbcmap_ptr_(Teuchos::null),
      freact_ptr_(nullptr)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::Init(const Teuchos::RCP<DRT::Discretization>& discret_ptr,
    const Teuchos::RCP<Epetra_Vector>& freact_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  // reset the setup indicator
  issetup_ = false;

  discret_ptr_ = discret_ptr;
  freact_ptr_ = freact_ptr.get();
  timint_ptr_ = timint_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::Setup()
{
  check_init();
  // ---------------------------------------------------------------------------
  // Create Dirichlet Boundary Condition map
  // ---------------------------------------------------------------------------
  zeros_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().DofRowMapView(), true));
  Teuchos::ParameterList p;
  p.set<double>("total time", timint_ptr_->GetDataGlobalState().GetTimeNp());
  dbcmap_ptr_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
  p.set<const CORE::UTILS::FunctionManager*>(
      "function_manager", &GLOBAL::Problem::Instance()->FunctionManager());
  discret_ptr_->evaluate_dirichlet(
      p, zeros_ptr_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmap_ptr_);
  // clear the system vector of possibly inserted non-zero DBC values
  zeros_ptr_->PutScalar(0.0);

  // ---------------------------------------------------------------------------
  // Create local coordinate system manager
  // ---------------------------------------------------------------------------
  std::vector<CORE::Conditions::Condition*> locsysconditions(0);
  discret_ptr_->GetCondition("Locsys", locsysconditions);
  if (locsysconditions.size())
  {
    locsysman_ptr_ = Teuchos::rcp(
        new CORE::Conditions::LocsysManager(*discret_ptr_, GLOBAL::Problem::Instance()->NDim()));
    // in case we have no time dependent locsys conditions in our problem,
    // this is the only time where the whole setup routine is conducted.
    locsysman_ptr_->Update(-1.0, {}, GLOBAL::Problem::Instance()->FunctionManager());
    islocsys_ = true;
  }

  // ---------------------------------------------------------------------------
  // Set the new pre/post operator for the nox nln linearsystem in the parameter
  // list
  // ---------------------------------------------------------------------------
  const Teuchos::ParameterList& pnox = timint_ptr_->GetDataSDyn().GetNoxParams();
  if (pnox.sublist("Direction").isSublist("Newton"))
  {
    if (pnox.sublist("Direction").sublist("Newton").isSublist("Linear Solver"))
    {
      // get a mutable reference to the linear solver parameter list
      Teuchos::ParameterList& p_linsolver = const_cast<Teuchos::ParameterList&>(
          pnox.sublist("Direction").sublist("Newton").sublist("Linear Solver"));
      NOX::NLN::LinSystem::PrePostOperator::Map& prepostlinsystem_map =
          NOX::NLN::LinSystem::PrePostOp::GetMap(p_linsolver);
      // create the new pre/post operator for the nox nln linear system
      Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prepostdbc_ptr =
          Teuchos::rcp(new NOX::NLN::LinSystem::PrePostOp::Dbc(Teuchos::rcp(this, false)));
      // insert/replace the old pointer in the map
      prepostlinsystem_map[NOX::NLN::LinSystem::prepost_dbc] = prepostdbc_ptr;
    }
    else
      FOUR_C_THROW(
          "There is no \"[NOX]->[Direction]->[Newton]->[Linear Solver] "
          "sub-sub-sublist!");
  }
  else
    FOUR_C_THROW("There is no \"[NOX]->[Direction]->[Newton]\" sub-sublist!");

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::check_init() const { FOUR_C_ASSERT(is_init(), "Call Init() first!"); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> STR::Dbc::discret_ptr()
{
  check_init();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const DRT::Discretization> STR::Dbc::discret_ptr() const
{
  check_init();
  return Teuchos::rcp_dynamic_cast<const DRT::Discretization>(discret_ptr_, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::UpdateLocSysManager()
{
  if (!is_loc_sys()) return;

  discret_ptr()->set_state("dispnp", g_state().GetDisNp());
  locsysman_ptr_->Update(g_state().GetTimeNp(), {}, GLOBAL::Problem::Instance()->FunctionManager());
  discret_ptr()->ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::Dbc::get_dirichlet_increment()
{
  Teuchos::RCP<const Epetra_Vector> disn = timint_ptr_->GetDataGlobalState().GetDisN();
  Teuchos::RCP<Epetra_Vector> dbcincr = Teuchos::rcp(new Epetra_Vector(*disn));
  const double& timenp = g_state().GetTimeNp();

  // get the new value for the Dirichlet DOFs
  apply_dirichlet_bc(timenp, dbcincr, Teuchos::null, Teuchos::null, false);

  /* Subtract the displacements of the last converged step:
   * --> DBC-DOFs hold increments of current step
   * --> free-DOFs hold zeros. */
  dbcincr->Update(-1.0, *disn, 1.0);

  return dbcincr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::apply_dirichlet_bc(const double& time, Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> acc, bool recreatemap) const
{
  check_init_setup();
  // We have to rotate forward ...
  // ---------------------------------------------------------------------------
  if (!dis.is_null()) RotateGlobalToLocal(dis, true);
  if (!vel.is_null()) RotateGlobalToLocal(vel);
  if (!acc.is_null()) RotateGlobalToLocal(acc);

  // Apply DBCs
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList p;
  p.set("total time", time);
  p.set<const CORE::UTILS::FunctionManager*>(
      "function_manager", &GLOBAL::Problem::Instance()->FunctionManager());

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_ptr_->ClearState();
  if (recreatemap)
    discret_ptr_->evaluate_dirichlet(p, dis, vel, acc, Teuchos::null, dbcmap_ptr_);
  else
    discret_ptr_->evaluate_dirichlet(p, dis, vel, acc, Teuchos::null, Teuchos::null);

  discret_ptr_->ClearState();

  // We have to rotate back into global Cartesian frame
  // ---------------------------------------------------------------------------
  if (dis != Teuchos::null) RotateLocalToGlobal(dis, true);
  if (vel != Teuchos::null) RotateLocalToGlobal(vel);
  if (acc != Teuchos::null) RotateLocalToGlobal(acc);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::apply_dirichlet_to_local_system(
    Teuchos::RCP<CORE::LINALG::SparseOperator> A, Teuchos::RCP<Epetra_Vector>& b) const
{
  check_init_setup();
  apply_dirichlet_to_local_rhs(b);
  apply_dirichlet_to_local_jacobian(A);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::apply_dirichlet_to_vector(Teuchos::RCP<Epetra_Vector>& vec) const
{
  check_init_setup();
  // rotate the coordinate system if desired
  RotateGlobalToLocal(vec);
  // apply the dbc
  CORE::LINALG::apply_dirichlet_to_system(*vec, *zeros_ptr_, *(dbcmap_ptr_->CondMap()));
  // rotate back
  RotateLocalToGlobal(vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::apply_dirichlet_to_local_rhs(Teuchos::RCP<Epetra_Vector>& b) const
{
  check_init_setup();

  // rotate the coordinate system: global --> local
  RotateGlobalToLocal(b);

  extract_freact(b);
  CORE::LINALG::apply_dirichlet_to_system(*b, *zeros_ptr_, *(dbcmap_ptr_->CondMap()));


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToRhs(Teuchos::RCP<Epetra_Vector>& b) const
{
  check_init_setup();

  apply_dirichlet_to_local_rhs(b);

  // rotate back: local --> global
  RotateLocalToGlobal(b);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::apply_dirichlet_to_local_jacobian(Teuchos::RCP<CORE::LINALG::SparseOperator> A) const
{
  check_init_setup();
  // don't do it twice...
  /* Note: If the DBCs are applied twice this can lead to an error and strange
   * behavior during the usage of locsys. Furthermore, the consideration of
   * DBCs in an explicit way is a pretty expensive operation.
   *                                                          hiermeier 01/18 */
  if (A->IsDbcApplied(*dbcmap_ptr_->CondMap(), true, get_loc_sys_trafo().get())) return;

  if (RotateGlobalToLocal(A))
  {
    Teuchos::RCP<std::vector<CORE::LINALG::SparseMatrix*>> mats =
        g_state().extract_displ_row_of_blocks(*A);

    for (unsigned i = 0; i < mats->size(); ++i)
    {
      CORE::LINALG::SparseMatrix& mat = *(*mats)[i];

      mat.apply_dirichlet_with_trafo(
          *get_loc_sys_trafo(), *(dbcmap_ptr_->CondMap()), (i == 0), false);
    }

    if (not A->Filled()) A->Complete();
  }
  else
    A->ApplyDirichlet(*(dbcmap_ptr_->CondMap()));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(const Teuchos::RCP<Epetra_Vector>& v) const
{
  check_init_setup();
  return RotateGlobalToLocal(v, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(const Teuchos::RCP<Epetra_Vector>& v, bool offset) const
{
  check_init_setup();
  if (not is_loc_sys()) return false;

  if (g_state().MaxBlockNumber() > 1)
  {
    Epetra_Vector v_displ(*g_state().DofRowMapView());
    CORE::LINALG::ExtractMyVector(*v, v_displ);

    locsysman_ptr_->RotateGlobalToLocal(Teuchos::rcpFromRef(v_displ), offset);

    CORE::LINALG::AssembleMyVector(0.0, *v, 1.0, v_displ);
  }
  else
    locsysman_ptr_->RotateGlobalToLocal(v, offset);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(const Teuchos::RCP<CORE::LINALG::SparseOperator>& A) const
{
  check_init_setup();
  if (not is_loc_sys()) return false;

  Teuchos::RCP<std::vector<CORE::LINALG::SparseMatrix*>> mats =
      g_state().extract_displ_row_of_blocks(*A);

  for (unsigned i = 0; i < mats->size(); ++i)
    locsysman_ptr_->RotateGlobalToLocal(Teuchos::rcpFromRef(*(*mats)[i]));

  if (not A->Filled()) A->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateLocalToGlobal(const Teuchos::RCP<Epetra_Vector>& v) const
{
  check_init_setup();
  return RotateLocalToGlobal(v, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateLocalToGlobal(const Teuchos::RCP<Epetra_Vector>& v, bool offset) const
{
  check_init_setup();
  if (not is_loc_sys()) return false;

  if (g_state().MaxBlockNumber() > 1)
  {
    Epetra_Vector v_displ(*g_state().DofRowMapView());
    CORE::LINALG::ExtractMyVector(*v, v_displ);

    locsysman_ptr_->RotateLocalToGlobal(Teuchos::rcpFromRef(v_displ), offset);

    CORE::LINALG::AssembleMyVector(0.0, *v, 1.0, v_displ);
  }
  else
    locsysman_ptr_->RotateLocalToGlobal(v, offset);

  // reset flag
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix> STR::Dbc::get_loc_sys_trafo() const
{
  check_init_setup();
  if (not is_loc_sys()) return Teuchos::null;

  return locsysman_ptr_->Trafo();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::extract_freact(Teuchos::RCP<Epetra_Vector>& b) const
{
  check_init_setup();

  CORE::LINALG::ExtractMyVector(*b, freact());
  freact().Scale(-1.0);

  // put zeros on all non-DBC dofs
  insert_vector_in_non_dbc_dofs(zeros_ptr_, Teuchos::rcpFromRef(freact()));

  // turn the reaction forces back to the global coordinate system if necessary
  RotateLocalToGlobal(Teuchos::rcpFromRef(freact()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::insert_vector_in_non_dbc_dofs(
    Teuchos::RCP<const Epetra_Vector> source_ptr, Teuchos::RCP<Epetra_Vector> target_ptr) const
{
  check_init_setup();
  dbcmap_ptr_->InsertOtherVector(dbcmap_ptr_->ExtractOtherVector(source_ptr), target_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MapExtractor> STR::Dbc::GetDBCMapExtractor() const
{
  check_init_setup();
  return dbcmap_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::Conditions::LocsysManager> STR::Dbc::LocSysManagerPtr()
{
  check_init_setup();
  return locsysman_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::Dbc::GetZeros() const
{
  check_init_setup();
  return *zeros_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::Dbc::GetZerosPtr() const
{
  check_init_setup();
  return zeros_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::Dbc::freact() const
{
  FOUR_C_ASSERT(freact_ptr_, "nullptr");

  return *freact_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::Dbc::g_state() const
{
  return timint_ptr_->GetDataGlobalState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::AddDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmap_ptr_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = CORE::LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmap_ptr_ = CORE::LINALG::MapExtractor(*(discret_ptr_->dof_row_map()), condmerged);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::RemoveDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(maptoremove);
  othermaps.push_back(dbcmap_ptr_->OtherMap());
  Teuchos::RCP<Epetra_Map> othermerged = CORE::LINALG::MultiMapExtractor::MergeMaps(othermaps);
  *dbcmap_ptr_ = CORE::LINALG::MapExtractor(*(discret_ptr_->dof_row_map()), othermerged, false);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOp::Dbc::Dbc(const Teuchos::RCP<const STR::Dbc>& dbc_ptr)
    : dbc_ptr_(dbc_ptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::Dbc::run_pre_apply_jacobian_inverse(
    ::NOX::Abstract::Vector& rhs, CORE::LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  ::NOX::Epetra::Vector& rhs_epetra = dynamic_cast<::NOX::Epetra::Vector&>(rhs);
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&rhs_epetra.getEpetraVector(), false);
  Teuchos::RCP<CORE::LINALG::SparseOperator> jac_ptr = Teuchos::rcp(&jac, false);
  // apply the dirichlet condition and rotate the system if desired
  dbc_ptr_->apply_dirichlet_to_local_system(jac_ptr, rhs_ptr);
}

FOUR_C_NAMESPACE_CLOSE
