/*-----------------------------------------------------------*/
/*! \file
\brief Some special solution strategy


\level 3
*/
/*-----------------------------------------------------------*/

#include "contact_aug_steepest_ascent_sp_strategy.H"
#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_potential.H"
#include "contact_aug_lagrange_multiplier_function.H"
#include "contact_aug_penalty_update.H"

#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_lib/epetra_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::DataContainer::DataContainer()
{ /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT_SP::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr, const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const plain_interface_set& interfaces, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof)
    : CONTACT::AUG::Strategy(data_ptr, DofRowMap, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  Data().InitSubDataContainer(INPAR::CONTACT::solution_steepest_ascent_sp);
  const Teuchos::ParameterList& sa_params =
      Params().sublist("AUGMENTED", true).sublist("STEEPESTASCENT", true);

  Data().SaData().SetPenaltyCorrectionParameter(sa_params.get<double>("CORRECTION_PARAMETER"));

  Data().SaData().SetPenaltyDecreaseCorrectionParameter(
      sa_params.get<double>("DECREASE_CORRECTION_PARAMETER"));

  Data().SaData().LagrangeMultiplierFuncPtr() = Teuchos::rcp(new LagrangeMultiplierFunction());

  Data().SaData().PenaltyUpdatePtr() = Teuchos::rcp(PenaltyUpdate::Create(sa_params));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    dserror(
        "ERROR: Augmented Lagrange Formulation: Self contact is not yet "
        "considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs_exp(*ProblemDofs());
  LINALG::Export(Data().SlForceLm(), augfs_exp);
  Data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm_exp(*ProblemDofs());
  LINALG::Export(Data().MaForceLm(), augfm_exp);
  CATCH_EPETRA_ERROR(Data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::PostSetup(bool redistributed, bool init)
{
  AUG::Strategy::PostSetup(redistributed, init);

  if (init)
  {
    Data().SaData().PenaltyUpdate().Init(this, &Data());

#ifdef LAGRANGE_FUNC
    Data().SaData().LagrangeMultiplierFunc().Init(this, Data());
    Data().SaData().LagrangeMultiplierFunc().Setup();
#endif
  }

  if (redistributed)
  {
#ifdef LAGRANGE_FUNC
    Data().SaData().LagrangeMultiplierFunc().Redistribute();
#endif
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::RunPostApplyJacobianInverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  /* Note that the result vector is the result of the linear system and
   * accordingly, due to the sign convention in NOX, the negative direction
   * vector of the Newton method. Therefore, the vector is converted before and
   * after the augmentation, since the used formulas expect a direction vector.
   *                                                         hiermeier, 12/17 */
  result.Scale(-1.0);
  SetPenaltyUpdateState(cparams, xold, result);
  result.Scale(-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::SetPenaltyUpdateState(
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

  Data().SaData().PenaltyUpdate().SetState(cparams, xold, dir);
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
                      << cparams.GetNumberOfModifiedNewtonCorrections() << IO::endl;
  IO::cout(IO::debug) << std::string(40, '*') << "\n";

  if (cparams.IsDefaultStep() or cparams.GetNumberOfModifiedNewtonCorrections() == 0)
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
void CONTACT::AUG::STEEPESTASCENT_SP::Strategy::AddContributionsToMatrixBlockLmLm(
    LINALG::SparseMatrix& kzz) const
{
  Teuchos::RCP<Epetra_Vector> active_mod_diag = GetKzzDiagModification();
  if (LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(kzz, *active_mod_diag))
  {
    Epetra_Vector kzz_diag = Epetra_Vector(kzz.RangeMap(), true);
    kzz.ExtractDiagonalCopy(kzz_diag);
    LINALG::AssembleMyVector(1.0, kzz_diag, 1.0, *active_mod_diag);

    // if the matrix is filled, we try to replace the diagonal
    if (kzz.ReplaceDiagonalValues(kzz_diag)) dserror("ReplaceDiagonalValues failed!");
  }

  CONTACT::AUG::Strategy::AddContributionsToMatrixBlockLmLm(kzz);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::STEEPESTASCENT_SP::Strategy::GetKzzDiagModification()
    const
{
  Teuchos::RCP<Epetra_Vector> active_mod_vec =
      Teuchos::rcp(new Epetra_Vector(*Data().KappaVecPtr()));
  CATCH_EPETRA_ERROR(active_mod_vec->ReplaceMap(*Data().GActiveNDofRowMapPtr()));

  MultiplyElementwise(*Data().CnPtr(), *Data().GActiveNodeRowMapPtr(), *active_mod_vec, true);

  return active_mod_vec;
}
