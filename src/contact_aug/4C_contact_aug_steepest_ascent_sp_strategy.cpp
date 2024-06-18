/*-----------------------------------------------------------*/
/*! \file
\brief Some special solution strategy


\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_contact_aug_steepest_ascent_sp_strategy.hpp"

#include "4C_contact_aug_lagrange_multiplier_function.hpp"
#include "4C_contact_aug_penalty_update.hpp"
#include "4C_contact_aug_potential.hpp"
#include "4C_contact_aug_steepest_ascent_strategy.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_epetra_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::SteepestAscent::DataContainer::DataContainer()
{ /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof)
    : CONTACT::Aug::Strategy(
          data_ptr, dof_row_map, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  data().init_sub_data_container(Inpar::CONTACT::solution_steepest_ascent_sp);
  const Teuchos::ParameterList& sa_params =
      Params().sublist("AUGMENTED", true).sublist("STEEPESTASCENT", true);

  data().SaData().set_penalty_correction_parameter(sa_params.get<double>("CORRECTION_PARAMETER"));

  data().SaData().set_penalty_decrease_correction_parameter(
      sa_params.get<double>("DECREASE_CORRECTION_PARAMETER"));

  data().SaData().lagrange_multiplier_func_ptr() = Teuchos::rcp(new LagrangeMultiplierFunction());

  data().SaData().PenaltyUpdatePtr() = Teuchos::rcp(PenaltyUpdate::Create(sa_params));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::eval_str_contact_rhs()
{
  if (!IsInContact() and !WasInContact() and !was_in_contact_last_time_step())
  {
    data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    FOUR_C_THROW(
        "ERROR: Augmented Lagrange Formulation: Self contact is not yet "
        "considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs_exp(*ProblemDofs());
  Core::LinAlg::Export(data().SlForceLm(), augfs_exp);
  data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm_exp(*ProblemDofs());
  Core::LinAlg::Export(data().MaForceLm(), augfm_exp);
  CATCH_EPETRA_ERROR(data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::post_setup(bool redistributed, bool init)
{
  Aug::Strategy::post_setup(redistributed, init);

  if (init)
  {
    data().SaData().PenaltyUpdate().Init(this, &data());

#ifdef LAGRANGE_FUNC
    Data().SaData().lagrange_multiplier_func().Init(this, Data());
    Data().SaData().lagrange_multiplier_func().setup();
#endif
  }

  if (redistributed)
  {
#ifdef LAGRANGE_FUNC
    Data().SaData().lagrange_multiplier_func().Redistribute();
#endif
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::Nln::Group& grp)
{
  /* Note that the result vector is the result of the linear system and
   * accordingly, due to the sign convention in NOX, the negative direction
   * vector of the Newton method. Therefore, the vector is converted before and
   * after the augmentation, since the used formulas expect a direction vector.
   *                                                         hiermeier, 12/17 */
  result.Scale(-1.0);
  set_penalty_update_state(cparams, xold, result);
  result.Scale(-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::set_penalty_update_state(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  const NOX::Nln::CorrectionType corrtype = cparams.get_correction_type();

  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << __LINE__ << " -- " << CONTACT_FUNC_NAME << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << "cparams.get_correction_type() = "
                                  << NOX::Nln::CorrectionType2String(corrtype).c_str()
                                  << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << Core::IO::endl;

  /* Set the state in the penalty update object only for full second order
   * correction steps and default solution steps. Actually the only case which
   * is currently excluded is the cheap SOC step, however, the IF-condition
   * is written in this way, such that other future corrections are excluded as
   * well. If you think your correction behaves like a full step, add it here.
   *                                                         hiermeier 08/18 */
  if (corrtype != NOX::Nln::CorrectionType::soc_full and
      corrtype != NOX::Nln::CorrectionType::vague)
    return;

  data().SaData().PenaltyUpdate().set_state(cparams, xold, dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::run_post_iterate(
    const CONTACT::ParamsInterface& cparams)
{
  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << "\n";
  Core::IO::cout(Core::IO::debug) << CONTACT_FUNC_NAME << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << "is_default_step = "
                                  << (cparams.is_default_step() ? "TRUE" : "FALSE")
                                  << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << "Number of modified Newton corrections = "
                                  << cparams.get_number_of_modified_newton_corrections()
                                  << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << "\n";

  if (cparams.is_default_step() or cparams.get_number_of_modified_newton_corrections() == 0)
    update_cn(cparams);
  else
    decrease_cn(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::update_cn(
    const CONTACT::ParamsInterface& cparams)
{
  data().SaData().PenaltyUpdate().Update(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::decrease_cn(
    const CONTACT::ParamsInterface& cparams)
{
  data().SaData().PenaltyUpdate().Decrease(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::add_contributions_to_matrix_block_lm_lm(
    Core::LinAlg::SparseMatrix& kzz) const
{
  Teuchos::RCP<Epetra_Vector> active_mod_diag = get_kzz_diag_modification();
  if (Core::LinAlg::InsertMyRowDiagonalIntoUnfilledMatrix(kzz, *active_mod_diag))
  {
    Epetra_Vector kzz_diag = Epetra_Vector(kzz.RangeMap(), true);
    kzz.ExtractDiagonalCopy(kzz_diag);
    Core::LinAlg::AssembleMyVector(1.0, kzz_diag, 1.0, *active_mod_diag);

    // if the matrix is filled, we try to replace the diagonal
    if (kzz.replace_diagonal_values(kzz_diag)) FOUR_C_THROW("replace_diagonal_values failed!");
  }

  CONTACT::Aug::Strategy::add_contributions_to_matrix_block_lm_lm(kzz);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::Aug::SteepestAscentSaddlePoint::Strategy::get_kzz_diag_modification() const
{
  Teuchos::RCP<Epetra_Vector> active_mod_vec =
      Teuchos::rcp(new Epetra_Vector(*data().KappaVecPtr()));
  CATCH_EPETRA_ERROR(active_mod_vec->ReplaceMap(*data().g_active_n_dof_row_map_ptr()));

  MultiplyElementwise(*data().CnPtr(), *data().g_active_node_row_map_ptr(), *active_mod_vec, true);

  return active_mod_vec;
}

FOUR_C_NAMESPACE_CLOSE
