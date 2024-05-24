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
CONTACT::AUG::STEEPESTASCENT::DataContainer::DataContainer()
{ /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT_SP::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof)
    : CONTACT::AUG::Strategy(
          data_ptr, dof_row_map, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  Data().init_sub_data_container(INPAR::CONTACT::solution_steepest_ascent_sp);
  const Teuchos::ParameterList& sa_params =
      Params().sublist("AUGMENTED", true).sublist("STEEPESTASCENT", true);

  Data().SaData().set_penalty_correction_parameter(sa_params.get<double>("CORRECTION_PARAMETER"));

  Data().SaData().set_penalty_decrease_correction_parameter(
      sa_params.get<double>("DECREASE_CORRECTION_PARAMETER"));

  Data().SaData().lagrange_multiplier_func_ptr() = Teuchos::rcp(new LagrangeMultiplierFunction());

  Data().SaData().PenaltyUpdatePtr() = Teuchos::rcp(PenaltyUpdate::Create(sa_params));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !was_in_contact_last_time_step())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    FOUR_C_THROW(
        "ERROR: Augmented Lagrange Formulation: Self contact is not yet "
        "considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs_exp(*ProblemDofs());
  CORE::LINALG::Export(Data().SlForceLm(), augfs_exp);
  Data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm_exp(*ProblemDofs());
  CORE::LINALG::Export(Data().MaForceLm(), augfm_exp);
  CATCH_EPETRA_ERROR(Data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::post_setup(bool redistributed, bool init)
{
  AUG::Strategy::post_setup(redistributed, init);

  if (init)
  {
    Data().SaData().PenaltyUpdate().Init(this, &Data());

#ifdef LAGRANGE_FUNC
    Data().SaData().lagrange_multiplier_func().Init(this, Data());
    Data().SaData().lagrange_multiplier_func().Setup();
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
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
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
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::set_penalty_update_state(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  const NOX::NLN::CorrectionType corrtype = cparams.GetCorrectionType();

  IO::cout(IO::debug) << std::string(40, '*') << IO::endl;
  IO::cout(IO::debug) << __LINE__ << " -- " << CONTACT_FUNC_NAME << IO::endl;
  IO::cout(IO::debug) << "cparams.GetCorrectionType() = "
                      << NOX::NLN::CorrectionType2String(corrtype).c_str() << IO::endl;
  IO::cout(IO::debug) << std::string(40, '*') << IO::endl;

  /* Set the state in the penalty update object only for full second order
   * correction steps and default solution steps. Actually the only case which
   * is currently excluded is the cheap SOC step, however, the IF-condition
   * is written in this way, such that other future corrections are excluded as
   * well. If you think your correction behaves like a full step, add it here.
   *                                                         hiermeier 08/18 */
  if (corrtype != NOX::NLN::CorrectionType::soc_full and
      corrtype != NOX::NLN::CorrectionType::vague)
    return;

  Data().SaData().PenaltyUpdate().set_state(cparams, xold, dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::RunPostIterate(
    const CONTACT::ParamsInterface& cparams)
{
  IO::cout(IO::debug) << std::string(40, '*') << "\n";
  IO::cout(IO::debug) << CONTACT_FUNC_NAME << IO::endl;
  IO::cout(IO::debug) << "IsDefaultStep = " << (cparams.IsDefaultStep() ? "TRUE" : "FALSE")
                      << IO::endl;
  IO::cout(IO::debug) << "Number of modified Newton corrections = "
                      << cparams.get_number_of_modified_newton_corrections() << IO::endl;
  IO::cout(IO::debug) << std::string(40, '*') << "\n";

  if (cparams.IsDefaultStep() or cparams.get_number_of_modified_newton_corrections() == 0)
    UpdateCn(cparams);
  else
    DecreaseCn(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::UpdateCn(const CONTACT::ParamsInterface& cparams)
{
  Data().SaData().PenaltyUpdate().Update(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::DecreaseCn(const CONTACT::ParamsInterface& cparams)
{
  Data().SaData().PenaltyUpdate().Decrease(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::add_contributions_to_matrix_block_lm_lm(
    CORE::LINALG::SparseMatrix& kzz) const
{
  Teuchos::RCP<Epetra_Vector> active_mod_diag = get_kzz_diag_modification();
  if (CORE::LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(kzz, *active_mod_diag))
  {
    Epetra_Vector kzz_diag = Epetra_Vector(kzz.RangeMap(), true);
    kzz.ExtractDiagonalCopy(kzz_diag);
    CORE::LINALG::AssembleMyVector(1.0, kzz_diag, 1.0, *active_mod_diag);

    // if the matrix is filled, we try to replace the diagonal
    if (kzz.replace_diagonal_values(kzz_diag)) FOUR_C_THROW("replace_diagonal_values failed!");
  }

  CONTACT::AUG::Strategy::add_contributions_to_matrix_block_lm_lm(kzz);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::STEEPESTASCENT_SP::Strategy::get_kzz_diag_modification()
    const
{
  Teuchos::RCP<Epetra_Vector> active_mod_vec =
      Teuchos::rcp(new Epetra_Vector(*Data().KappaVecPtr()));
  CATCH_EPETRA_ERROR(active_mod_vec->ReplaceMap(*Data().g_active_n_dof_row_map_ptr()));

  MultiplyElementwise(*Data().CnPtr(), *Data().g_active_node_row_map_ptr(), *active_mod_vec, true);

  return active_mod_vec;
}

FOUR_C_NAMESPACE_CLOSE
