/*---------------------------------------------------------------------*/
/*! \file
\brief Steepest ascent solution strategy based on the augmented contact
       formulation.

\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_aug_steepest_ascent_strategy.hpp"

#include "4C_contact_aug_lagrange_multiplier_function.hpp"
#include "4C_contact_aug_penalty_update.hpp"
#include "4C_contact_aug_potential.hpp"
#include "4C_contact_aug_steepest_ascent_interface.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_matrix_transform.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_utils_epetra_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

// #define LAGRANGE_FUNC

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr, const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const plain_interface_set& interfaces, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof)
    : CONTACT::AUG::STEEPESTASCENT_SP::Strategy(
          data_ptr, DofRowMap, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  // cast to steepest ascent interfaces
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::Interface>& interface = *cit;
    // test interfaces for the correct type
    Teuchos::rcp_dynamic_cast<CONTACT::AUG::STEEPESTASCENT::Interface>(interface, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::add_contributions_to_constr_rhs(
    Epetra_Vector& augConstrRhs) const
{
  // do nothing
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::STEEPESTASCENT::Strategy::get_rhs_block_ptr_for_norm_check(
    const enum CONTACT::VecBlockType& bt) const
{
  if (!IsInContact() and !WasInContact() and !was_in_contact_last_time_step()) return Teuchos::null;

  Teuchos::RCP<Epetra_Vector> rhs_block = Teuchos::null;

  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
    {
      FOUR_C_THROW("Unused!");

      break;
    }
    case CONTACT::VecBlockType::constraint:
    {
      rhs_block = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));

      AUG::Strategy::add_contributions_to_constr_rhs(*rhs_block);
      rhs_block->ReplaceMap(LMDoFRowMap(true));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported VecBlocktype! (enum=%d)", bt);
      exit(EXIT_FAILURE);
    }
  }

  return rhs_block;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::AUG::STEEPESTASCENT::Strategy::GetMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !was_in_contact_last_time_step()) return Teuchos::null;

  Teuchos::RCP<CORE::LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
    {
      mat_ptr = Teuchos::rcp(new CORE::LINALG::SparseMatrix(SlMaDoFRowMap(true), 100, false, true));

      // build matrix kdd
      add_contributions_to_matrix_block_displ_displ(*mat_ptr, cparams);
      mat_ptr->Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));

      // transform parallel row/column distribution
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        MORTAR::MatrixRowColTransformer& transformer = Data().matrix_row_col_transformer();
        mat_ptr = transformer.redistributed_to_unredistributed(bt, *mat_ptr);
      }

      break;
    }
    case CONTACT::MatBlockType::displ_lm:
    {
      // do nothing

      break;
    }
    case CONTACT::MatBlockType::lm_displ:
    {
      // do nothing

      break;
    }
    case CONTACT::MatBlockType::lm_lm:
    {
      // do nothing

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown STR::MatBlockType!");
      break;
    }
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::add_contributions_to_matrix_block_displ_displ(
    CORE::LINALG::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams) const
{
  //  if ( cparams and cparams->GetPredictorType() != INPAR::STR::pred_tangdis )
  kdd.Add(*Data().DGGLinMatrixPtr(), false, 1.0, 1.0);

  /* ignore the Lagrange multiplier dependent contact contributions during the
   * TangDis predictor */
  //  if ( cparams and cparams->GetPredictorType() != INPAR::STR::pred_tangdis )
  kdd.Add(*Data().DGLmLinMatrixPtr(), false, -1.0, 1.0);

  // add inactive contributions (this is not well tested)
  if (Data().add_inactiv_force_contributions())
    kdd.Add(*Data().InactiveDDMatrixPtr(), false, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  /* Note that the result vector is the result of the linear system and
   * accordingly, due to the sign convention in NOX, the negative direction
   * vector of the Newton method. Therefore, the vector is converted before and
   * after the augmentation, since the used formulas expect a direction vector.
   *                                                         hiermeier, 12/17 */
  result.Scale(-1.0);
  AugmentDirection(cparams, xold, result);
  result.Scale(-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::AugmentDirection(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  // if there are no contact contributions, do a direct return
  if (!IsInContact() && !WasInContact() && !was_in_contact_last_time_step()) return;

  // extract displ. increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_ptr =
      CORE::LINALG::ExtractMyVector(dir_mutable, *ProblemDofs());
  Teuchos::RCP<const Epetra_Vector> displ_incr_redistributed_ptr = Teuchos::null;

  if (ParRedist())
  {
    Teuchos::RCP<Epetra_Vector> tmp_exp = Teuchos::rcp(new Epetra_Vector(*gdisprowmap_));
    CORE::LINALG::Export(*displ_incr_ptr, *tmp_exp);
    displ_incr_redistributed_ptr = tmp_exp;
  }
  else
    displ_incr_redistributed_ptr = displ_incr_ptr;

  const Epetra_Vector& displ_incr = *displ_incr_redistributed_ptr;

  // --------------------------------------------------------------------------
  // extract old Lagrange multiplier from the old solution vector
  Teuchos::RCP<Epetra_Vector> zold_ptr = CORE::LINALG::ExtractMyVector(xold, LMDoFRowMap(false));
  zold_ptr->ReplaceMap(SlDoFRowMap(false));

  Teuchos::RCP<Epetra_Vector> zold_redistributed_ptr = Teuchos::null;
  if (ParRedist())
  {
    zold_redistributed_ptr = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    // export the zold vector to the zold redistributed vector
    CORE::LINALG::Export(*zold_ptr, *zold_redistributed_ptr);
  }
  else
    zold_redistributed_ptr = zold_ptr;

  const Epetra_Vector& zold = *zold_redistributed_ptr;

  // --------------------------------------------------------------------------
  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active =
      compute_active_lagrange_incr_in_normal_direction(displ_incr);

  // --------------------------------------------------------------------------
  // inactive lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive =
      compute_inactive_lagrange_incr_in_normal_direction(displ_incr, zold);

  // --------------------------------------------------------------------------
  // assemble the Lagrange multiplier contributions
  Epetra_Vector zincr_redistributed(SlDoFRowMap(true));
  CORE::LINALG::AssembleMyVector(0.0, zincr_redistributed, 1.0, *znincr_active);
  CORE::LINALG::AssembleMyVector(0.0, zincr_redistributed, 1.0, *zincr_inactive);
  zincr_redistributed.ReplaceMap(LMDoFRowMap(true));

  Epetra_Vector zincr_full(LMDoFRowMap(false));
  CORE::LINALG::Export(zincr_redistributed, zincr_full);

  CORE::LINALG::AssembleMyVector(0.0, dir_mutable, 1.0, zincr_full);

  // run at the very end...
  post_augment_direction(cparams, xold, dir_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::AUG::STEEPESTASCENT::Strategy::compute_active_lagrange_incr_in_normal_direction(
    const Epetra_Vector& displ_incr)
{
  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active_ptr =
      Teuchos::rcp(new Epetra_Vector(Data().GActiveNDofRowMap(), true));
  Epetra_Vector& znincr_active = *znincr_active_ptr;

  // nothing to do, if there are no active contributions
  if (not IsInContact()) return znincr_active_ptr;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> tempmtx12, tempmtx21, tempmtx22;

  Teuchos::RCP<CORE::LINALG::SparseMatrix> gradWGapUpdate;
  // Split dLmNWGapLinMatrix_
  CORE::LINALG::SplitMatrix2x2(Data().d_lm_nw_gap_lin_matrix_ptr(),
      Data().g_active_n_dof_row_map_ptr(), emptymap, gdisprowmap_, emptymap, gradWGapUpdate,
      tempmtx12, tempmtx21, tempmtx22);

  // calculate the Uzawa Update increment
  // *** Attention: zincr_Update has the wrong sign here! ***
  int err = gradWGapUpdate->Multiply(false, displ_incr, znincr_active);
  if (err) FOUR_C_THROW("Multiply error! (err=%d)", err);

  CATCH_EPETRA_ERROR(znincr_active.Update(1.0, Data().WGap(), 1.0));

  // Scaling of the Lagrange multiplier increment
  // --> inverse area scaling
  MultiplyElementwise(Data().KappaVec(), Data().GActiveNodeRowMap(), znincr_active, true);

  // Update the final Lagrange multiplier increment.
  // These values will also be used to update the nodal quantities during
  // the recover routine.
  MultiplyElementwise(Data().Cn(), Data().GActiveNodeRowMap(), znincr_active, false);

  // We correct the increment sign.
  znincr_active.Scale(-1.0);

  return znincr_active_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::AUG::STEEPESTASCENT::Strategy::compute_inactive_lagrange_incr_in_normal_direction(
    const Epetra_Vector& displ_incr, const Epetra_Vector& zold)
{
  Teuchos::RCP<Epetra_Map> ginactivedofs =
      CORE::LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  // inactive lagrange multiplier increment in normal and tangential direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive_ptr = Teuchos::rcp(new Epetra_Vector(*ginactivedofs));
  Epetra_Vector& zincr_inactive = *zincr_inactive_ptr;

  // extract old lagrange multipliers in normal and tangential direction
  Teuchos::RCP<const Epetra_Vector> zold_inactive_ptr =
      CORE::LINALG::ExtractMyVector(zold, *ginactivedofs);
  const Epetra_Vector& zold_inactive = *zold_inactive_ptr;

  // extract displ increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_sl_ptr =
      CORE::LINALG::ExtractMyVector(displ_incr, SlDoFRowMap(true));

  int err = Data().InactiveLinMatrix().Multiply(false, *displ_incr_sl_ptr, zincr_inactive);
  if (err) FOUR_C_THROW("Multiply error (err=%d)!", err);

  zincr_inactive.ReciprocalMultiply(-1.0, Data().InactiveDiagMatrix(), zincr_inactive, 0.0);

  CATCH_EPETRA_ERROR(zincr_inactive.Update(-1.0, zold_inactive, 1.0));

  return zincr_inactive_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::post_augment_direction(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir)
{
  set_penalty_update_state(cparams, xold, dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::remove_condensed_contributions_from_rhs(
    Epetra_Vector& str_rhs) const
{
  Epetra_Vector regforce(*Data().GSlMaDofRowMapPtr());
  CORE::LINALG::AssembleMyVector(0.0, regforce, 1.0, *Data().SlForceGPtr());
  CORE::LINALG::AssembleMyVector(1.0, regforce, 1.0, *Data().MaForceGPtr());

  Epetra_Vector regforce_exp(*ProblemDofs());
  CORE::LINALG::Export(regforce, regforce_exp);
  CATCH_EPETRA_ERROR(str_rhs.Update(-1.0, regforce_exp, 1.0));
}

FOUR_C_NAMESPACE_CLOSE
