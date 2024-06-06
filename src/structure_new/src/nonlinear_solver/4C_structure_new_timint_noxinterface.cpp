/*-----------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the Jacobian, Required and
       Preconditioner %NOX::NLN interfaces.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_timint_noxinterface.hpp"

#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_constraint_group.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"

#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::NoxInterface::NoxInterface()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null),
      int_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::Init(
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<const STR::TimeInt::Base>& timint_ptr)
{
  // reset the setup flag
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;
  timint_ptr_ = timint_ptr;
  int_ptr_ = int_ptr;
  dbc_ptr_ = dbc_ptr;

  // set the initialization flag
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::Setup()
{
  check_init();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::check_init() const
{
  FOUR_C_ASSERT(is_init(), "Call Init() first!");
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator& STR::TimeInt::NoxInterface::ImplInt()
{
  check_init_setup();
  return *int_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::computeF(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  check_init_setup();

  if (not int_ptr_->ApplyForce(x, F)) return false;

  /* Apply the DBC on the right hand side, since we need the Dirichlet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&F, false);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& jac)
{
  check_init_setup();

  Core::LinAlg::SparseOperator* jac_ptr = dynamic_cast<Core::LinAlg::SparseOperator*>(&jac);
  FOUR_C_ASSERT(jac_ptr != nullptr, "Dynamic cast failed.");

  if (not int_ptr_->ApplyStiff(x, *jac_ptr)) return false;

  /* We do not consider the jacobian DBC at this point. The Dirichlet conditions
   * are applied inside the NOX::Nln::LinearSystem::applyJacobianInverse()
   * routine, instead. See the run_pre_apply_jacobian_inverse() implementation
   * for more information.                               hiermeier 01/15/2016 */

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::computeFandJacobian(
    const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac)
{
  check_init_setup();

  Core::LinAlg::SparseOperator* jac_ptr = dynamic_cast<Core::LinAlg::SparseOperator*>(&jac);
  FOUR_C_ASSERT(jac_ptr != nullptr, "Dynamic cast failed!");

  if (not int_ptr_->ApplyForceStiff(x, rhs, *jac_ptr)) return false;

  /* Apply the DBC on the right hand side, since we need the Dirichlet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&rhs, false);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

  /* We do not consider the jacobian DBC at this point. The Dirichlet conditions
   * are applied inside the NOX::Nln::LinearSystem::applyJacobianInverse()
   * routine, instead. See the run_pre_apply_jacobian_inverse() implementation
   * for more information.                               hiermeier 01/15/2016 */

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::compute_correction_system(const enum NOX::Nln::CorrectionType type,
    const ::NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
    Epetra_Operator& jac)
{
  check_init_setup();

  Core::LinAlg::SparseOperator* jac_ptr = dynamic_cast<Core::LinAlg::SparseOperator*>(&jac);
  FOUR_C_ASSERT(jac_ptr != nullptr, "Dynamic cast failed!");

  std::vector<Inpar::STR::ModelType> constraint_models;
  find_constraint_models(&grp, constraint_models);

  if (not int_ptr_->apply_correction_system(type, constraint_models, x, rhs, *jac_ptr))
    return false;

  /* Apply the DBC on the right hand side, since we need the Dirichlet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcpFromRef(rhs);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  check_init_setup();
  // currently not supported
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::GetPrimaryRHSNorms(const Epetra_Vector& F,
    const NOX::Nln::StatusTest::QuantityType& checkquantity,
    const ::NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  check_init_setup();
  double rhsnorm = -1.0;

  // convert the given quantity type to a model type
  const Inpar::STR::ModelType mt = STR::Nln::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::Nln::StatusTest::quantity_structure:
    case NOX::Nln::StatusTest::quantity_meshtying:
    case NOX::Nln::StatusTest::quantity_cardiovascular0d:
    {
      // export the model specific solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs_ptr = gstate_ptr_->ExtractModelEntries(mt, F);

      // remove entries specific to element technology
      gstate_ptr_->remove_element_technologies(rhs_ptr);

      int_ptr_->remove_condensed_contributions_from_rhs(*rhs_ptr);

      rhsnorm = calculate_norm(rhs_ptr, type, isscaled);

      break;
    }
    case NOX::Nln::StatusTest::quantity_pressure:
    {
      // export the model specific solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs_ptr = gstate_ptr_->ExtractModelEntries(mt, F);

      // extract entries specific to element technology
      gstate_ptr_->extract_element_technologies(NOX::Nln::StatusTest::quantity_pressure, rhs_ptr);

      rhsnorm = calculate_norm(rhs_ptr, type, isscaled);

      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return rhsnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::get_primary_solution_update_rms(const Epetra_Vector& xnew,
    const Epetra_Vector& xold, const double& atol, const double& rtol,
    const NOX::Nln::StatusTest::QuantityType& checkquantity,
    const bool& disable_implicit_weighting) const
{
  check_init_setup();

  double rms = -1.0;

  // convert the given quantity type to a model type
  const Inpar::STR::ModelType mt = STR::Nln::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::Nln::StatusTest::quantity_structure:
    case NOX::Nln::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->ExtractModelEntries(mt, xold)));
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr = gstate_ptr_->ExtractModelEntries(mt, xnew);

      // remove entries specific to element technology
      gstate_ptr_->remove_element_technologies(model_incr_ptr);
      gstate_ptr_->remove_element_technologies(model_xnew_ptr);

      model_incr_ptr->Update(1.0, *model_xnew_ptr, -1.0);
      rms = NOX::Nln::Aux::RootMeanSquareNorm(
          atol, rtol, model_xnew_ptr, model_incr_ptr, disable_implicit_weighting);

      break;
    }
    case NOX::Nln::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->ExtractModelEntries(mt, xold)));
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr = gstate_ptr_->ExtractModelEntries(mt, xnew);

      // extract entries specific to element technology
      gstate_ptr_->extract_element_technologies(
          NOX::Nln::StatusTest::quantity_pressure, model_incr_ptr);
      gstate_ptr_->extract_element_technologies(
          NOX::Nln::StatusTest::quantity_pressure, model_xnew_ptr);

      model_incr_ptr->Update(1.0, *model_xnew_ptr, -1.0);
      rms = NOX::Nln::Aux::RootMeanSquareNorm(
          atol, rtol, model_xnew_ptr, model_incr_ptr, disable_implicit_weighting);
      break;
    }
    case NOX::Nln::StatusTest::quantity_eas:
    case NOX::Nln::StatusTest::quantity_plasticity:
    {
      rms = int_ptr_->get_condensed_solution_update_rms(checkquantity);
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::get_primary_solution_update_norms(const Epetra_Vector& xnew,
    const Epetra_Vector& xold, const NOX::Nln::StatusTest::QuantityType& checkquantity,
    const ::NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  check_init_setup();

  double updatenorm = -1.0;

  // convert the given quantity type to a model type
  const Inpar::STR::ModelType mt = STR::Nln::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::Nln::StatusTest::quantity_structure:
    case NOX::Nln::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr = gstate_ptr_->ExtractModelEntries(mt, xold);
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr = gstate_ptr_->ExtractModelEntries(mt, xnew);

      // remove entries specific to element technology
      gstate_ptr_->remove_element_technologies(model_incr_ptr);
      gstate_ptr_->remove_element_technologies(model_xnew_ptr);

      model_incr_ptr->Update(1.0, *model_xnew_ptr, -1.0);
      updatenorm = calculate_norm(model_incr_ptr, type, isscaled);

      break;
    }
    case NOX::Nln::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr = gstate_ptr_->ExtractModelEntries(mt, xold);
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr = gstate_ptr_->ExtractModelEntries(mt, xnew);

      // extract entries specific to element technology
      gstate_ptr_->extract_element_technologies(
          NOX::Nln::StatusTest::quantity_pressure, model_incr_ptr);
      gstate_ptr_->extract_element_technologies(
          NOX::Nln::StatusTest::quantity_pressure, model_xnew_ptr);

      model_incr_ptr->Update(1.0, *model_xnew_ptr, -1.0);
      updatenorm = calculate_norm(model_incr_ptr, type, isscaled);

      break;
    }
    case NOX::Nln::StatusTest::quantity_eas:
    case NOX::Nln::StatusTest::quantity_plasticity:
    {
      // get the update norm of the condensed quantities
      updatenorm = int_ptr_->get_condensed_update_norm(checkquantity);
      // do the scaling if desired
      if (isscaled)
      {
        int gdofnumber = int_ptr_->get_condensed_dof_number(checkquantity);
        updatenorm /= static_cast<double>(gdofnumber);
      }
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::get_previous_primary_solution_norms(const Epetra_Vector& xold,
    const NOX::Nln::StatusTest::QuantityType& checkquantity,
    const ::NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  check_init_setup();

  double xoldnorm = -1.0;

  // convert the given quantity type to a model type
  const Inpar::STR::ModelType mt = STR::Nln::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::Nln::StatusTest::quantity_structure:
    case NOX::Nln::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_xold_ptr = gstate_ptr_->ExtractModelEntries(mt, xold);

      // remove entries specific to element technology
      gstate_ptr_->remove_element_technologies(model_xold_ptr);

      xoldnorm = calculate_norm(model_xold_ptr, type, isscaled);

      break;
    }
    case NOX::Nln::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_xold_ptr = gstate_ptr_->ExtractModelEntries(mt, xold);

      // extract entries specific to element technology
      gstate_ptr_->extract_element_technologies(
          NOX::Nln::StatusTest::quantity_pressure, model_xold_ptr);

      xoldnorm = calculate_norm(model_xold_ptr, type, isscaled);

      break;
    }
    case NOX::Nln::StatusTest::quantity_eas:
    case NOX::Nln::StatusTest::quantity_plasticity:
    {
      // get the update norm of the condensed quantities
      xoldnorm = int_ptr_->get_condensed_previous_sol_norm(checkquantity);
      if (isscaled)
      {
        int gdofnumber = int_ptr_->get_condensed_dof_number(checkquantity);
        xoldnorm /= static_cast<double>(gdofnumber);
      }
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return xoldnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::calculate_norm(Teuchos::RCP<Epetra_Vector> quantity,
    const ::NOX::Abstract::Vector::NormType type, const bool isscaled) const
{
  Teuchos::RCP<const ::NOX::Epetra::Vector> quantity_nox =
      Teuchos::rcp(new ::NOX::Epetra::Vector(quantity, ::NOX::Epetra::Vector::CreateView));

  double norm = quantity_nox->norm(type);
  // do the scaling if desired
  if (isscaled) norm /= static_cast<double>(quantity_nox->length());

  return norm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::GetModelValue(const Epetra_Vector& x, const Epetra_Vector& F,
    const enum NOX::Nln::MeritFunction::MeritFctName merit_func_type) const
{
  check_init_setup();

  double omval = 0.0;

  switch (merit_func_type)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    case NOX::Nln::MeritFunction::mrtfct_energy:
    {
      Core::IO::cout(Core::IO::debug) << __LINE__ << " - " << __FUNCTION__ << "\n";
      int_ptr_->get_total_mid_time_str_energy(x);
      omval = int_ptr_->GetModelValue(x);

      break;
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      // do nothing in the primary field
      break;
    }
    default:
    {
      FOUR_C_THROW("There is no objective model value for %s | %d.",
          NOX::Nln::MeritFunction::MeritFuncName2String(merit_func_type).c_str(), merit_func_type);
      exit(EXIT_FAILURE);
    }
  }

  return omval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::get_linearized_model_terms(const ::NOX::Abstract::Group* group,
    const Epetra_Vector& dir, const enum NOX::Nln::MeritFunction::MeritFctName mf_type,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  switch (mf_type)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    {
      return get_linearized_energy_model_terms(group, dir, linorder, lintype);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
      return 0.0;
    default:
    {
      FOUR_C_THROW("There is no linearization for the objective model %s | %d.",
          NOX::Nln::MeritFunction::MeritFuncName2String(mf_type).c_str(), mf_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::get_linearized_energy_model_terms(
    const ::NOX::Abstract::Group* group, const Epetra_Vector& dir,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  double lin_val = 0.0;

  switch (linorder)
  {
    case NOX::Nln::MeritFunction::linorder_first:
    case NOX::Nln::MeritFunction::linorder_all:
    {
      switch (lintype)
      {
        case NOX::Nln::MeritFunction::lin_wrt_all_dofs:
        case NOX::Nln::MeritFunction::lin_wrt_primary_dofs:
        {
          Epetra_Vector str_gradient(dir.Map(), true);

          std::vector<Inpar::STR::ModelType> constraint_models;
          find_constraint_models(group, constraint_models);

          // assemble the force and exclude all constraint models
          int_ptr_->assemble_force(str_gradient, &constraint_models);
          str_gradient.Dot(dir, &lin_val);

          Core::IO::cout(Core::IO::debug)
              << "LinEnergy   D_{d} (Energy) = " << lin_val << Core::IO::endl;

          break;
        }
        default:
        {
          /* do nothing, there are only primary dofs */
          break;
        }
      }

      break;
    }
    default:
    {
      /* do nothing, there are no high order terms */
      break;
    }
  }

  return lin_val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::find_constraint_models(
    const ::NOX::Abstract::Group* grp, std::vector<Inpar::STR::ModelType>& constraint_models) const
{
  const NOX::Nln::CONSTRAINT::Group* constr_grp =
      dynamic_cast<const NOX::Nln::CONSTRAINT::Group*>(grp);

  // direct return if this is no constraint problem
  if (not constr_grp) return;

  // find the constraint model types
  const auto& imap = constr_grp->GetConstrInterfaces();
  constraint_models.reserve(imap.size());

  for (auto cit = imap.begin(); cit != imap.end(); ++cit)
  {
    const enum NOX::Nln::SolutionType soltype = cit->first;
    const enum Inpar::STR::ModelType mtype = STR::Nln::ConvertSolType2ModelType(soltype);

    constraint_models.push_back(mtype);
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::NoxInterface::CalcRefNormForce()
{
  check_init_setup();
  const ::NOX::Epetra::Vector::NormType& nox_normtype = timint_ptr_->GetDataSDyn().GetNoxNormType();
  return int_ptr_->CalcRefNormForce(nox_normtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix>
STR::TimeInt::NoxInterface::calc_jacobian_contributions_from_element_level_for_ptc()
{
  check_init_setup();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> scalingMatrixOpPtr =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstate_ptr_->dof_row_map(), 81, true, true));
  int_ptr_->compute_jacobian_contributions_from_element_level_for_ptc(scalingMatrixOpPtr);

  return scalingMatrixOpPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::CreateBackupState(const Epetra_Vector& dir)
{
  check_init_setup();
  int_ptr_->CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::recover_from_backup_state()
{
  check_init_setup();
  int_ptr_->recover_from_backup_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::NoxInterface::compute_element_volumes(
    const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols) const
{
  check_init_setup();
  return int_ptr_->determine_element_volumes(x, ele_vols);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::NoxInterface::getDofsFromElements(
    const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const
{
  check_init_setup();

  Teuchos::RCP<const Discret::Discretization> discret_ptr = gstate_ptr_->GetDiscret();

  for (int egid : my_ele_gids)
  {
    Core::Elements::Element* ele = discret_ptr->gElement(egid);
    Core::Nodes::Node** nodes = ele->Nodes();

    for (int i = 0; i < ele->num_node(); ++i)
    {
      if (nodes[i]->Owner() != gstate_ptr_->GetComm().MyPID()) continue;

      const std::vector<int> ndofs(discret_ptr->Dof(0, nodes[i]));
      my_ele_dofs.insert(ndofs.begin(), ndofs.end());
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
