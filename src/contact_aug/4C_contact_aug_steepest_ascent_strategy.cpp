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
CONTACT::Aug::SteepestAscent::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof)
    : CONTACT::Aug::SteepestAscentSaddlePoint::Strategy(
          data_ptr, dof_row_map, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  // cast to steepest ascent interfaces
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::Interface>& interface = *cit;
    // test interfaces for the correct type
    Teuchos::rcp_dynamic_cast<CONTACT::Aug::SteepestAscent::Interface>(interface, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscent::Strategy::add_contributions_to_constr_rhs(
    Epetra_Vector& augConstrRhs) const
{
  // do nothing
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::Aug::SteepestAscent::Strategy::get_rhs_block_ptr_for_norm_check(
    const enum CONTACT::VecBlockType& bt) const
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step())
    return Teuchos::null;

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
      rhs_block = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true), true));

      Aug::Strategy::add_contributions_to_constr_rhs(*rhs_block);
      rhs_block->ReplaceMap(lm_dof_row_map(true));

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
Teuchos::RCP<Core::LinAlg::SparseMatrix>
CONTACT::Aug::SteepestAscent::Strategy::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  // if there are no active contact contributions
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
    return Teuchos::null;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(slave_master_dof_row_map(true), 100, false, true));

      // build matrix kdd
      add_contributions_to_matrix_block_displ_displ(*mat_ptr, cparams);
      mat_ptr->Complete(slave_master_dof_row_map(true), slave_master_dof_row_map(true));

      // transform parallel row/column distribution
      // (only necessary in the parallel redistribution case)
      if (parallel_redistribution_status())
      {
        Mortar::MatrixRowColTransformer& transformer = data().matrix_row_col_transformer();
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
void CONTACT::Aug::SteepestAscent::Strategy::add_contributions_to_matrix_block_displ_displ(
    Core::LinAlg::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams) const
{
  //  if ( cparams and cparams->get_predictor_type() != Inpar::STR::pred_tangdis )
  kdd.Add(*data().DGGLinMatrixPtr(), false, 1.0, 1.0);

  /* ignore the Lagrange multiplier dependent contact contributions during the
   * TangDis predictor */
  //  if ( cparams and cparams->get_predictor_type() != Inpar::STR::pred_tangdis )
  kdd.Add(*data().DGLmLinMatrixPtr(), false, -1.0, 1.0);

  // add inactive contributions (this is not well tested)
  if (data().add_inactiv_force_contributions())
    kdd.Add(*data().InactiveDDMatrixPtr(), false, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscent::Strategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::Nln::Group& grp)
{
  /* Note that the result vector is the result of the linear system and
   * accordingly, due to the sign convention in NOX, the negative direction
   * vector of the Newton method. Therefore, the vector is converted before and
   * after the augmentation, since the used formulas expect a direction vector.
   *                                                         hiermeier, 12/17 */
  result.Scale(-1.0);
  augment_direction(cparams, xold, result);
  result.Scale(-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscent::Strategy::augment_direction(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  // if there are no contact contributions, do a direct return
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // extract displ. increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_ptr =
      Core::LinAlg::ExtractMyVector(dir_mutable, *ProblemDofs());
  Teuchos::RCP<const Epetra_Vector> displ_incr_redistributed_ptr = Teuchos::null;

  if (parallel_redistribution_status())
  {
    Teuchos::RCP<Epetra_Vector> tmp_exp = Teuchos::rcp(new Epetra_Vector(*gdisprowmap_));
    Core::LinAlg::Export(*displ_incr_ptr, *tmp_exp);
    displ_incr_redistributed_ptr = tmp_exp;
  }
  else
    displ_incr_redistributed_ptr = displ_incr_ptr;

  const Epetra_Vector& displ_incr = *displ_incr_redistributed_ptr;

  // --------------------------------------------------------------------------
  // extract old Lagrange multiplier from the old solution vector
  Teuchos::RCP<Epetra_Vector> zold_ptr = Core::LinAlg::ExtractMyVector(xold, lm_dof_row_map(false));
  zold_ptr->ReplaceMap(slave_dof_row_map(false));

  Teuchos::RCP<Epetra_Vector> zold_redistributed_ptr = Teuchos::null;
  if (parallel_redistribution_status())
  {
    zold_redistributed_ptr = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    // export the zold vector to the zold redistributed vector
    Core::LinAlg::Export(*zold_ptr, *zold_redistributed_ptr);
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
  Epetra_Vector zincr_redistributed(slave_dof_row_map(true));
  Core::LinAlg::AssembleMyVector(0.0, zincr_redistributed, 1.0, *znincr_active);
  Core::LinAlg::AssembleMyVector(0.0, zincr_redistributed, 1.0, *zincr_inactive);
  zincr_redistributed.ReplaceMap(lm_dof_row_map(true));

  Epetra_Vector zincr_full(lm_dof_row_map(false));
  Core::LinAlg::Export(zincr_redistributed, zincr_full);

  Core::LinAlg::AssembleMyVector(0.0, dir_mutable, 1.0, zincr_full);

  // run at the very end...
  post_augment_direction(cparams, xold, dir_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::Aug::SteepestAscent::Strategy::compute_active_lagrange_incr_in_normal_direction(
    const Epetra_Vector& displ_incr)
{
  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active_ptr =
      Teuchos::rcp(new Epetra_Vector(data().global_active_n_dof_row_map(), true));
  Epetra_Vector& znincr_active = *znincr_active_ptr;

  // nothing to do, if there are no active contributions
  if (not is_in_contact()) return znincr_active_ptr;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx12, tempmtx21, tempmtx22;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> gradWGapUpdate;
  // Split dLmNWGapLinMatrix_
  Core::LinAlg::SplitMatrix2x2(data().d_lm_nw_gap_lin_matrix_ptr(),
      data().global_active_n_dof_row_map_ptr(), emptymap, gdisprowmap_, emptymap, gradWGapUpdate,
      tempmtx12, tempmtx21, tempmtx22);

  // calculate the Uzawa Update increment
  // *** Attention: zincr_Update has the wrong sign here! ***
  int err = gradWGapUpdate->Multiply(false, displ_incr, znincr_active);
  if (err) FOUR_C_THROW("Multiply error! (err=%d)", err);

  CATCH_EPETRA_ERROR(znincr_active.Update(1.0, data().WGap(), 1.0));

  // Scaling of the Lagrange multiplier increment
  // --> inverse area scaling
  MultiplyElementwise(data().KappaVec(), data().global_active_node_row_map(), znincr_active, true);

  // Update the final Lagrange multiplier increment.
  // These values will also be used to update the nodal quantities during
  // the recover routine.
  MultiplyElementwise(data().Cn(), data().global_active_node_row_map(), znincr_active, false);

  // We correct the increment sign.
  znincr_active.Scale(-1.0);

  return znincr_active_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::Aug::SteepestAscent::Strategy::compute_inactive_lagrange_incr_in_normal_direction(
    const Epetra_Vector& displ_incr, const Epetra_Vector& zold)
{
  Teuchos::RCP<Epetra_Map> ginactivedofs =
      Core::LinAlg::SplitMap(slave_dof_row_map(true), data().global_active_dof_row_map());

  // inactive lagrange multiplier increment in normal and tangential direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive_ptr = Teuchos::rcp(new Epetra_Vector(*ginactivedofs));
  Epetra_Vector& zincr_inactive = *zincr_inactive_ptr;

  // extract old lagrange multipliers in normal and tangential direction
  Teuchos::RCP<const Epetra_Vector> zold_inactive_ptr =
      Core::LinAlg::ExtractMyVector(zold, *ginactivedofs);
  const Epetra_Vector& zold_inactive = *zold_inactive_ptr;

  // extract displ increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_sl_ptr =
      Core::LinAlg::ExtractMyVector(displ_incr, slave_dof_row_map(true));

  int err = data().InactiveLinMatrix().Multiply(false, *displ_incr_sl_ptr, zincr_inactive);
  if (err) FOUR_C_THROW("Multiply error (err=%d)!", err);

  zincr_inactive.ReciprocalMultiply(-1.0, data().InactiveDiagMatrix(), zincr_inactive, 0.0);

  CATCH_EPETRA_ERROR(zincr_inactive.Update(-1.0, zold_inactive, 1.0));

  return zincr_inactive_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscent::Strategy::post_augment_direction(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir)
{
  set_penalty_update_state(cparams, xold, dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::SteepestAscent::Strategy::remove_condensed_contributions_from_rhs(
    Epetra_Vector& str_rhs) const
{
  Epetra_Vector regforce(*data().global_slave_master_dof_row_map_ptr());
  Core::LinAlg::AssembleMyVector(0.0, regforce, 1.0, *data().SlForceGPtr());
  Core::LinAlg::AssembleMyVector(1.0, regforce, 1.0, *data().MaForceGPtr());

  Epetra_Vector regforce_exp(*ProblemDofs());
  Core::LinAlg::Export(regforce, regforce_exp);
  CATCH_EPETRA_ERROR(str_rhs.Update(-1.0, regforce_exp, 1.0));
}

FOUR_C_NAMESPACE_CLOSE
