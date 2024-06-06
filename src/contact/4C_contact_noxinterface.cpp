/*---------------------------------------------------------------------*/
/*! \file
\brief Concrete mplementation of all the %NOX::Nln::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_noxinterface.hpp"

#include "4C_contact_abstract_strategy.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::NoxInterface::NoxInterface()
    : isinit_(false),
      issetup_(false),
      strategy_ptr_(Teuchos::null),
      cycling_maps_(std::vector<Teuchos::RCP<Epetra_Map>>(0))
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::NoxInterface::Init(const Teuchos::RCP<CONTACT::AbstractStrategy>& strategy_ptr)
{
  issetup_ = false;

  strategy_ptr_ = strategy_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::NoxInterface::Setup()
{
  check_init();

  // set flag at the end
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_constraint_rhs_norms(const Epetra_Vector& F,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  Teuchos::RCP<const Epetra_Vector> constrRhs =
      strategy().get_rhs_block_ptr_for_norm_check(CONTACT::VecBlockType::constraint);

  // no contact contributions present
  if (constrRhs.is_null()) return 0.0;

  // export the vector to the current redistributed map
  Teuchos::RCP<Epetra_Vector> constrRhs_red = Teuchos::null;
  // Note: PointSameAs is faster than SameAs and should do the job right here,
  // since we replace the map afterwards anyway.               hiermeier 08/17
  if (not constrRhs->Map().PointSameAs(strategy().LMDoFRowMap(true)))
  {
    constrRhs_red = Teuchos::rcp(new Epetra_Vector(strategy().LMDoFRowMap(true)));
    Core::LinAlg::Export(*constrRhs, *constrRhs_red);
  }
  else
    constrRhs_red = Teuchos::rcp(new Epetra_Vector(*constrRhs));

  // replace the map
  constrRhs_red->ReplaceMap(strategy().SlDoFRowMap(true));

  double constrNorm = -1.0;
  Teuchos::RCP<const ::NOX::Epetra::Vector> constrRhs_nox = Teuchos::null;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // create vector with redistributed slave dof row map in normal direction
      Teuchos::RCP<Epetra_Vector> nConstrRhs =
          Core::LinAlg::ExtractMyVector(*constrRhs_red, strategy().SlNormalDoFRowMap(true));


      constrRhs_nox =
          Teuchos::rcp(new ::NOX::Epetra::Vector(nConstrRhs, ::NOX::Epetra::Vector::CreateView));
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // create vector with redistributed slave dof row map in tangential directions
      Teuchos::RCP<Epetra_Vector> tConstrRhs = Core::LinAlg::ExtractMyVector(
          *constrRhs_red, strategy().sl_tangential_do_f_row_map(true));

      constrRhs_nox =
          Teuchos::rcp(new ::NOX::Epetra::Vector(tConstrRhs, ::NOX::Epetra::Vector::CreateView));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }
  constrNorm = constrRhs_nox->norm(type);
  if (isScaled) constrNorm /= static_cast<double>(constrRhs_nox->length());

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew,
    const Epetra_Vector& xOld, double aTol, double rTol,
    NOX::Nln::StatusTest::QuantityType checkQuantity, bool disable_implicit_weighting) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  double rms = -1.0;
  Teuchos::RCP<Epetra_Vector> z_ptr = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> zincr_ptr = Teuchos::null;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // extract vectors with redistributed slave dof row map in normal direction
      z_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().GetLagrMultNp(true), strategy().SlNormalDoFRowMap(true));
      zincr_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().get_lagr_mult_solve_incr(), strategy().SlNormalDoFRowMap(true));

      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // extract vectors with redistributed slave dof row map in tangential directions
      z_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().GetLagrMultNp(true), strategy().sl_tangential_do_f_row_map(true));
      zincr_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().get_lagr_mult_solve_incr(), strategy().sl_tangential_do_f_row_map(true));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }

  rms = NOX::Nln::Aux::RootMeanSquareNorm(aTol, rTol, strategy().GetLagrMultNp(true),
      strategy().get_lagr_mult_solve_incr(), disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_lagrange_multiplier_update_norms(const Epetra_Vector& xNew,
    const Epetra_Vector& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
    ::NOX::Abstract::Vector::NormType type, bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  if (strategy().GetLagrMultNp(true) == Teuchos::null) return 0.;

  double updatenorm = -1.0;
  Teuchos::RCP<Epetra_Vector> zincr_ptr = Teuchos::null;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      // extract vector with redistributed slave dof row map in normal direction
      zincr_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().get_lagr_mult_solve_incr(), strategy().SlNormalDoFRowMap(true));
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      // extract vector with redistributed slave dof row map in tangential directions
      zincr_ptr = Core::LinAlg::ExtractMyVector(
          *strategy().get_lagr_mult_solve_incr(), strategy().sl_tangential_do_f_row_map(true));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported quantity type!");
      break;
    }
  }

  Teuchos::RCP<const ::NOX::Epetra::Vector> zincr_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(zincr_ptr, ::NOX::Epetra::Vector::CreateView));

  updatenorm = zincr_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) updatenorm /= static_cast<double>(zincr_nox_ptr->length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_previous_lagrange_multiplier_norms(const Epetra_Vector& xOld,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_contact_normal and
      checkQuantity != NOX::Nln::StatusTest::quantity_contact_friction)
    return -1.0;

  double zoldnorm = -1.0;

  if (strategy().GetLagrMultNp(true) == Teuchos::null) return 0;

  /* lagrange multiplier of the previous Newton step
   * (NOT equal to zOld_, which is stored in the Strategy object!!!) */
  Teuchos::RCP<Epetra_Vector> zold_ptr =
      Teuchos::rcp(new Epetra_Vector(*strategy().GetLagrMultNp(true)));
  zold_ptr->Update(-1.0, *strategy().get_lagr_mult_solve_incr(), 1.0);
  Teuchos::RCP<::NOX::Epetra::Vector> zold_nox_ptr = Teuchos::null;
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      Teuchos::RCP<Epetra_Vector> znold_ptr =
          Core::LinAlg::ExtractMyVector(*zold_ptr, strategy().SlNormalDoFRowMap(true));

      zold_nox_ptr =
          Teuchos::rcp(new ::NOX::Epetra::Vector(znold_ptr, ::NOX::Epetra::Vector::CreateView));
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      Teuchos::RCP<Epetra_Vector> ztold_ptr =
          Core::LinAlg::ExtractMyVector(*zold_ptr, strategy().sl_tangential_do_f_row_map(true));

      zold_nox_ptr =
          Teuchos::rcp(new ::NOX::Epetra::Vector(ztold_ptr, ::NOX::Epetra::Vector::CreateView));
      break;
    }
    default:
    {
      FOUR_C_THROW("The given quantity type is unsupported!");
      break;
    }
  }

  zoldnorm = zold_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) zoldnorm /= static_cast<double>(zold_nox_ptr->length());

  // avoid very small norm values for the pure inactive case
  if (not strategy().IsInContact()) zoldnorm = std::max(zoldnorm, 1.0);

  return zoldnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum ::NOX::StatusTest::StatusType CONTACT::NoxInterface::GetActiveSetInfo(
    NOX::Nln::StatusTest::QuantityType checkQuantity, int& activesetsize) const
{
  bool semismooth = Core::UTILS::IntegralValue<int>(strategy().Params(), "SEMI_SMOOTH_NEWTON");
  if (not semismooth) FOUR_C_THROW("Currently we support only the semi-smooth Newton case!");
  // ---------------------------------------------------------------------------
  // get the number of active nodes for the given active set type
  // ---------------------------------------------------------------------------
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      activesetsize = strategy().NumberOfActiveNodes();
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      activesetsize = strategy().NumberOfSlipNodes();
      break;
    }
    default:
    {
      FOUR_C_THROW("The given quantity type is unsupported!");
      break;
    }
  }
  // ---------------------------------------------------------------------------
  // translate the active set semi-smooth Newton convergence flag
  // ---------------------------------------------------------------------------
  if (strategy().active_set_semi_smooth_converged())
    return ::NOX::StatusTest::Converged;
  else
    return ::NOX::StatusTest::Unconverged;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::NoxInterface::get_current_active_set_map(
    enum NOX::Nln::StatusTest::QuantityType checkQuantity) const
{
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      return strategy().ActiveRowNodes();
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      return strategy().SlipRowNodes();
      break;
    }
    default:
    {
      FOUR_C_THROW("The given active set type is unsupported!");
      break;
    }
  }  // switch (active_set_type)

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::NoxInterface::GetOldActiveSetMap(
    enum NOX::Nln::StatusTest::QuantityType checkQuantity) const
{
  switch (checkQuantity)
  {
    case NOX::Nln::StatusTest::quantity_contact_normal:
    {
      return strategy().get_old_active_row_nodes();
      break;
    }
    case NOX::Nln::StatusTest::quantity_contact_friction:
    {
      return strategy().GetOldSlipRowNodes();
      break;
    }
    default:
    {
      FOUR_C_THROW("The given active set type is unsupported!");
      break;
    }
  }  // switch (active_set_type)

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::GetModelValue(NOX::Nln::MeritFunction::MeritFctName name) const
{
  switch (name)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    {
      return strategy().GetPotentialValue(name);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      double val = strategy().GetPotentialValue(name);
      val = std::sqrt(val);

      return val;
    }
    case NOX::Nln::MeritFunction::mrtfct_energy:
    {
      // The energy of the primary field is considered, no contact contribution.
      return 0.0;
    }
    default:
      FOUR_C_THROW("Unsupported Merit function name! (enum = %d)", name);
      exit(EXIT_FAILURE);
  }

  FOUR_C_THROW("Impossible to reach this point.");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::NoxInterface::get_linearized_model_terms(const Epetra_Vector& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName name,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  switch (name)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    {
      return strategy().get_linearized_potential_value_terms(dir, name, linorder, lintype);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      double lin_val =
          strategy().get_linearized_potential_value_terms(dir, name, linorder, lintype);
      const double modelvalue = GetModelValue(name);
      if (modelvalue != 0.0) lin_val /= modelvalue;

      return lin_val;
    }
    default:
      FOUR_C_THROW("Unsupported Merit function name! (enum = %d)", name);
      exit(EXIT_FAILURE);
  }

  FOUR_C_THROW("Impossible to reach this point.");
  exit(EXIT_FAILURE);
}

FOUR_C_NAMESPACE_CLOSE
