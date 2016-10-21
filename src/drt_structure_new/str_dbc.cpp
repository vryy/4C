/*-----------------------------------------------------------*/
/*!
\file str_dbc.cpp

\brief Wrapper for all Dirichlet boundary condition tasks.

\maintainer Michael Hiermeier

\date Nov 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_dbc.H"
#include "str_timint_base.H"

#include "../solver_nonlin_nox/nox_nln_linearsystem_prepostoperator.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_locsys.H"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>


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
      dbcmap_ptr_(Teuchos::null)
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
  freact_ptr_ = freact_ptr;
  timint_ptr_ = timint_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::Setup()
{
  CheckInit();
  // ---------------------------------------------------------------------------
  // Create Dirichlet Boundary Condition map
  // ---------------------------------------------------------------------------
  zeros_ptr_ = LINALG::CreateVector(*(discret_ptr_->DofRowMap()),true);
  Teuchos::ParameterList p;
  p.set("total time", timint_ptr_->GetDataGlobalState().GetTimeNp());
  dbcmap_ptr_ = Teuchos::rcp(new LINALG::MapExtractor());
  discret_ptr_->EvaluateDirichlet(p,zeros_ptr_,
      Teuchos::null,Teuchos::null,Teuchos::null,dbcmap_ptr_);
  zeros_ptr_->Scale(0.0);

  // ---------------------------------------------------------------------------
  // Create local coordinate system manager
  // ---------------------------------------------------------------------------
  std::vector<DRT::Condition*> locsysconditions(0);
  discret_ptr_->GetCondition("Locsys", locsysconditions);
  if (locsysconditions.size())
  {
    locsysman_ptr_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_ptr_));
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
          NOX::NLN::LinSystem::PrePostOp::GetMutableMap(p_linsolver);
      // create the new pre/post operator for the nox nln linear system
      Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prepostdbc_ptr =
          Teuchos::rcp(new NOX::NLN::LinSystem::PrePostOp::Dbc(
              Teuchos::rcp(this,false)));
      // insert/replace the old pointer in the map
      prepostlinsystem_map[NOX::NLN::LinSystem::prepost_dbc] =
          prepostdbc_ptr;
    }
    else
      dserror("There is no \"[NOX]->[Direction]->[Newton]->[Linear Solver] "
          "sub-sub-sublist!");
  }
  else
    dserror("There is no \"[NOX]->[Direction]->[Newton]\" sub-sublist!");

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>& STR::Dbc::DiscretPtr()
{
  CheckInit();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::UpdateLocSysManager()
{
  if (!IsLocSys())
    return;

  locsysman_ptr_->Setup(timint_ptr_->GetDataGlobalState().GetTimeNp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::Dbc::GetDirichletIncrement()
{
Teuchos::RCP<const Epetra_Vector> disn =
    timint_ptr_->GetDataGlobalState().GetDisN();
  Teuchos::RCP<Epetra_Vector> dbcincr = Teuchos::rcp(new Epetra_Vector(*disn));
  const double& timenp = timint_ptr_->GetDataGlobalState().GetTimeNp();

  // get the new value for the Dirichlet DOFs
  ApplyDirichletBC(timenp,dbcincr,Teuchos::null,Teuchos::null,false);

  /* Subtract the displacements of the last converged step:
   * --> DBC-DOFs hold increments of current step
   * --> free-DOFs hold zeros. */
  dbcincr->Update(-1.0,*disn,1.0);

  return dbcincr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletBC(
    const double& time,
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Epetra_Vector> vel,
    Teuchos::RCP<Epetra_Vector> acc,
    bool recreatemap
  )
{
  CheckInitSetup();
  // We have to rotate forward ...
  // ---------------------------------------------------------------------------
  if (!dis.is_null())
    RotateGlobalToLocal(dis,true);
  if (!vel.is_null())
    RotateGlobalToLocal(vel);
  if (!acc.is_null())
    RotateGlobalToLocal(acc);

  // Apply DBCs
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList p;
  p.set("total time",time);

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_ptr_->ClearState();
  if (recreatemap)
    discret_ptr_->EvaluateDirichlet(p, dis, vel, acc, Teuchos::null, dbcmap_ptr_);
  else
    discret_ptr_->EvaluateDirichlet(p, dis, vel, acc, Teuchos::null, Teuchos::null);

  discret_ptr_->ClearState();

  // We have to rotate back into global Cartesian frame
  // ---------------------------------------------------------------------------
  if (dis != Teuchos::null)
    RotateLocalToGlobal(dis, true);
  if (vel != Teuchos::null)
    RotateLocalToGlobal(vel);
  if (acc != Teuchos::null)
    RotateLocalToGlobal(acc);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToLocalSystem(Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<Epetra_Vector>& b) const
{
  CheckInitSetup();
  ApplyDirichletToLocalRhs(b);
  ApplyDirichletToLocalJacobian(A);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToVector(Teuchos::RCP<Epetra_Vector>& vec) const
{
  CheckInitSetup();
  // rotate the coordinate system if desired
  RotateGlobalToLocal(vec);
  // apply the dbc
  LINALG::ApplyDirichlettoSystem(vec,zeros_ptr_,*(dbcmap_ptr_->CondMap()));
  // rotate back
  RotateLocalToGlobal(vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToLocalRhs(Teuchos::RCP<Epetra_Vector>& b) const
{
  CheckInitSetup();

  // rotate the coordinate system: global --> local
  RotateGlobalToLocal(b);

  ExtractFreact(b);
  LINALG::ApplyDirichlettoSystem(b,zeros_ptr_,*(dbcmap_ptr_->CondMap()));


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToRhs(Teuchos::RCP<Epetra_Vector>& b)
{
  CheckInitSetup();

  ApplyDirichletToLocalRhs(b);

  // rotate back: local --> global
  RotateLocalToGlobal(b);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ApplyDirichletToLocalJacobian(
    Teuchos::RCP<LINALG::SparseOperator> A) const
{
  CheckInitSetup();
  RotateGlobalToLocal(A);

  Teuchos::RCP<Epetra_Vector> tmp_vec_ptr = Teuchos::null;
  LINALG::ApplyDirichlettoSystem(A,tmp_vec_ptr,tmp_vec_ptr,
      GetLocSysTrafo(),zeros_ptr_,*(dbcmap_ptr_->CondMap()));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(Teuchos::RCP<Epetra_Vector> v) const
{
  CheckInitSetup();
  return RotateGlobalToLocal(v,false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(Teuchos::RCP<Epetra_Vector> v, bool offset)
    const
{
  CheckInitSetup();
  if (not IsLocSys())
    return false;

  locsysman_ptr_->RotateGlobalToLocal(v,offset);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateGlobalToLocal(Teuchos::RCP<LINALG::SparseOperator> A)
    const
{
  CheckInitSetup();
  if (not IsLocSys())
    return false;

  Teuchos::RCP<LINALG::SparseMatrix> sysmat =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(A);
  if (sysmat.is_null())
    dserror("Dynamic cast failed!");
  locsysman_ptr_->RotateGlobalToLocal(sysmat);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> v) const
{
  CheckInitSetup();
  return RotateLocalToGlobal(v,false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Dbc::RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> v, bool offset)
    const
{
  CheckInitSetup();
  if (not IsLocSys())
    return false;

  locsysman_ptr_->RotateLocalToGlobal(v,offset);

  // reset flag
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> STR::Dbc::GetLocSysTrafo() const
{
  CheckInitSetup();
  if (not IsLocSys())
    return Teuchos::null;

  return locsysman_ptr_->Trafo();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::ExtractFreact(Teuchos::RCP<Epetra_Vector>& b) const
{
  CheckInitSetup();

  freact_ptr_->Update(-1.0,*b,0.0);

  // put zeros on all non-DBC dofs
  InsertVectorInNonDbcDofs(zeros_ptr_,freact_ptr_);

  // turn the reaction forces back to the global coordinate system if necessary
  RotateLocalToGlobal(freact_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Dbc::InsertVectorInNonDbcDofs(
    Teuchos::RCP<const Epetra_Vector> source_ptr,
    Teuchos::RCP<Epetra_Vector> target_ptr) const
{
  CheckInitSetup();
  dbcmap_ptr_->InsertOtherVector(
      dbcmap_ptr_->ExtractOtherVector(source_ptr),target_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> STR::Dbc::GetDBCMapExtractor() const
{
  CheckInitSetup();
  return dbcmap_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::LocsysManager> STR::Dbc::LocSysManagerPtr()
{
  CheckInitSetup();
  return locsysman_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::Dbc::GetZeros() const
{
  CheckInitSetup();
  return *zeros_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::Dbc::GetZerosPtr() const
{
  CheckInitSetup();
  return zeros_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOp::Dbc::Dbc(
    const Teuchos::RCP<const ::STR::Dbc>& dbc_ptr)
    : dbc_ptr_(dbc_ptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::Dbc::runPreApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  NOX::Epetra::Vector& rhs_epetra =
      dynamic_cast<NOX::Epetra::Vector&>(rhs);
  Teuchos::RCP<Epetra_Vector> rhs_ptr =
      Teuchos::rcp(&rhs_epetra.getEpetraVector(),false);
  Teuchos::RCP<LINALG::SparseOperator> jac_ptr =
      Teuchos::rcp(&jac,false);
  // apply the dirichlet condition and rotate the system if desired
  dbc_ptr_->ApplyDirichletToLocalSystem(jac_ptr,rhs_ptr);
}

