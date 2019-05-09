/*---------------------------------------------------------------------*/
/*!
\brief Steepest ascent solution strategy based on the augmented contact
       formulation.

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/

#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_steepest_ascent_interface.H"
#include "contact_aug_potential.H"
#include "contact_aug_lagrange_multiplier_function.H"
#include "contact_aug_penalty_update.H"

#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_matrix_transform.H"

#include "../drt_inpar/inpar_structure.H"

#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/epetra_utils.H"

//#define LAGRANGE_FUNC

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
    const Teuchos::RCP<CONTACT::CoInterface>& interface = *cit;
    // test interfaces for the correct type
    Teuchos::rcp_dynamic_cast<CONTACT::AUG::STEEPESTASCENT::Interface>(interface, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::AddContributionsToConstrRHS(
    Epetra_Vector& augConstrRhs) const
{
  // do nothing
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::STEEPESTASCENT::Strategy::GetRhsBlockPtrForNormCheck(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep()) return Teuchos::null;

  Teuchos::RCP<Epetra_Vector> rhs_block = Teuchos::null;

  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      dserror("Unused!");

      break;
    }
    case DRT::UTILS::block_constraint:
    {
      rhs_block = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));

      AUG::Strategy::AddContributionsToConstrRHS(*rhs_block);
      rhs_block->ReplaceMap(LMDoFRowMap(true));

      break;
    }
    default:
    {
      dserror("Unsupported VecBlocktype! (enum=%d)", bt);
      exit(EXIT_FAILURE);
    }
  }

  return rhs_block;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::STEEPESTASCENT::Strategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(new LINALG::SparseMatrix(SlMaDoFRowMap(true), 100, false, true));

      // build matrix kdd
      AddContributionsToMatrixBlockDisplDispl(*mat_ptr, cparams);
      mat_ptr->Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));

      // transform parallel row/column distribution
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        MORTAR::MatrixRowColTransformer& transformer = Data().MatrixRowColTransformer();
        mat_ptr = transformer.RedistributedToUnredistributed(bt, *mat_ptr);
      }

      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      // do nothing

      break;
    }
    case DRT::UTILS::block_lm_displ:
    {
      // do nothing

      break;
    }
    case DRT::UTILS::block_lm_lm:
    {
      // do nothing

      break;
    }
    default:
    {
      dserror("Unknown STR::MatBlockType!");
      break;
    }
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::AddContributionsToMatrixBlockDisplDispl(
    LINALG::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams) const
{
  //  if ( cparams and cparams->GetPredictorType() != INPAR::STR::pred_tangdis )
  kdd.Add(*Data().DGGLinMatrixPtr(), false, 1.0, 1.0);

  /* ignore the Lagrange multiplier dependent contact contributions during the
   * TangDis predictor */
  //  if ( cparams and cparams->GetPredictorType() != INPAR::STR::pred_tangdis )
  kdd.Add(*Data().DGLmLinMatrixPtr(), false, -1.0, 1.0);

  // add inactive contributions (this is not well tested)
  if (Data().AddInactivForceContributions())
    kdd.Add(*Data().InactiveDDMatrixPtr(), false, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::RunPostApplyJacobianInverse(
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
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // extract displ. increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_ptr =
      LINALG::ExtractMyVector(dir_mutable, *ProblemDofs());
  Teuchos::RCP<const Epetra_Vector> displ_incr_redistributed_ptr = Teuchos::null;

  if (ParRedist())
  {
    Teuchos::RCP<Epetra_Vector> tmp_exp = Teuchos::rcp(new Epetra_Vector(*gdisprowmap_));
    LINALG::Export(*displ_incr_ptr, *tmp_exp);
    displ_incr_redistributed_ptr = tmp_exp;
  }
  else
    displ_incr_redistributed_ptr = displ_incr_ptr;

  const Epetra_Vector& displ_incr = *displ_incr_redistributed_ptr;

  // --------------------------------------------------------------------------
  // extract old Lagrange multiplier from the old solution vector
  Teuchos::RCP<Epetra_Vector> zold_ptr = LINALG::ExtractMyVector(xold, LMDoFRowMap(false));
  zold_ptr->ReplaceMap(SlDoFRowMap(false));

  Teuchos::RCP<Epetra_Vector> zold_redistributed_ptr = Teuchos::null;
  if (ParRedist())
  {
    zold_redistributed_ptr = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    // export the zold vector to the zold redistributed vector
    LINALG::Export(*zold_ptr, *zold_redistributed_ptr);
  }
  else
    zold_redistributed_ptr = zold_ptr;

  const Epetra_Vector& zold = *zold_redistributed_ptr;

  // --------------------------------------------------------------------------
  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active =
      ComputeActiveLagrangeIncrInNormalDirection(displ_incr);

  // --------------------------------------------------------------------------
  // inactive lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive =
      ComputeInactiveLagrangeIncrInNormalDirection(displ_incr, zold);

  // --------------------------------------------------------------------------
  // assemble the Lagrange multiplier contributions
  Epetra_Vector zincr_redistributed(SlDoFRowMap(true));
  LINALG::AssembleMyVector(0.0, zincr_redistributed, 1.0, *znincr_active);
  LINALG::AssembleMyVector(0.0, zincr_redistributed, 1.0, *zincr_inactive);
  zincr_redistributed.ReplaceMap(LMDoFRowMap(true));

  Epetra_Vector zincr_full(LMDoFRowMap(false));
  LINALG::Export(zincr_redistributed, zincr_full);

  LINALG::AssembleMyVector(0.0, dir_mutable, 1.0, zincr_full);

  // run at the very end...
  PostAugmentDirection(cparams, xold, dir_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::AUG::STEEPESTASCENT::Strategy::ComputeActiveLagrangeIncrInNormalDirection(
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
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx12, tempmtx21, tempmtx22;

  Teuchos::RCP<LINALG::SparseMatrix> gradWGapUpdate;
  // Split dLmNWGapLinMatrix_
  LINALG::SplitMatrix2x2(Data().DLmNWGapLinMatrixPtr(), Data().GActiveNDofRowMapPtr(), emptymap,
      gdisprowmap_, emptymap, gradWGapUpdate, tempmtx12, tempmtx21, tempmtx22);

  // calculate the Uzawa Update increment
  // *** Attention: zincr_Update has the wrong sign here! ***
  int err = gradWGapUpdate->Multiply(false, displ_incr, znincr_active);
  if (err) dserror("Multiply error! (err=%d)", err);

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
CONTACT::AUG::STEEPESTASCENT::Strategy::ComputeInactiveLagrangeIncrInNormalDirection(
    const Epetra_Vector& displ_incr, const Epetra_Vector& zold)
{
  Teuchos::RCP<Epetra_Map> ginactivedofs =
      LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  // inactive lagrange multiplier increment in normal and tangential direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive_ptr = Teuchos::rcp(new Epetra_Vector(*ginactivedofs));
  Epetra_Vector& zincr_inactive = *zincr_inactive_ptr;

  // extract old lagrange multipliers in normal and tangential direction
  Teuchos::RCP<const Epetra_Vector> zold_inactive_ptr =
      LINALG::ExtractMyVector(zold, *ginactivedofs);
  const Epetra_Vector& zold_inactive = *zold_inactive_ptr;

  // extract displ increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_sl_ptr =
      LINALG::ExtractMyVector(displ_incr, SlDoFRowMap(true));

  int err = Data().InactiveLinMatrix().Multiply(false, *displ_incr_sl_ptr, zincr_inactive);
  if (err) dserror("Multiply error (err=%d)!", err);

  zincr_inactive.ReciprocalMultiply(-1.0, Data().InactiveDiagMatrix(), zincr_inactive, 0.0);

  CATCH_EPETRA_ERROR(zincr_inactive.Update(-1.0, zold_inactive, 1.0));

  return zincr_inactive_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::PostAugmentDirection(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir)
{
  SetPenaltyUpdateState(cparams, xold, dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::RemoveCondensedContributionsFromRhs(
    Epetra_Vector& str_rhs) const
{
  Epetra_Vector regforce(*Data().GSlMaDofRowMapPtr());
  LINALG::AssembleMyVector(0.0, regforce, 1.0, *Data().SlForceGPtr());
  LINALG::AssembleMyVector(1.0, regforce, 1.0, *Data().MaForceGPtr());

  Epetra_Vector regforce_exp(*ProblemDofs());
  LINALG::Export(regforce, regforce_exp);
  CATCH_EPETRA_ERROR(str_rhs.Update(-1.0, regforce_exp, 1.0));
}
