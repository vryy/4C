/*--------------------------------------------------------------------------*/
/*! \file

\brief Mesh sliding for ale problems


\level 3
*/
/*--------------------------------------------------------------------------*/

#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"

#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Teuchos_TimeMonitor.hpp>

#include "ale_meshsliding.H"

ALE::Meshsliding::Meshsliding(Teuchos::RCP<DRT::Discretization> dis, LINALG::Solver& solver,
    int msht, int nsd, const UTILS::MapExtractor* surfacesplitter)
    : Meshtying(dis, solver, msht, nsd, surfacesplitter)
{
}

/*-------------------------------------------------------*/
/*  Call the constructor and the setup of the mortar     */
/*  coupling adapter                         wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::AdapterMortar(std::vector<int> coupleddof)
{
  adaptermeshsliding_ = Teuchos::rcp(new ADAPTER::CouplingNonLinMortar());

  // Setup and Output of Nonlinear meshtying adapter
  adaptermeshsliding_->Setup(discret_, discret_, coupleddof, "Mortar");
}

/*-------------------------------------------------------*/
/*  Setup mesh sliding problem               wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> ALE::Meshsliding::Setup(
    std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  Teuchos::RCP<LINALG::SparseOperator> mat = Meshtying::Setup(coupleddof, dispnp);

  lm_ = LINALG::CreateVector(*gsdofrowmap_, true);

  return mat;
}

/*-------------------------------------------------------*/
/*  Compare the size of the slave and master dof row map */
/*                                           wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::CompareNumDof()
{
  int numdofmaster = (adaptermeshsliding_->MasterDofMap())->NumGlobalElements();
  int numdofslave = (adaptermeshsliding_->SlaveDofMap())->NumGlobalElements();

  std::cout << std::endl << "number of master dof's:   " << numdofmaster << std::endl;
  std::cout << "number of slave dof's:   " << numdofslave << std::endl << std::endl;

  if (numdofmaster > numdofslave)
    std::cout << "The master side is discretized by more elements than the slave side" << std::endl;
  else
    std::cout << "The slave side is discretized by more elements than the master side" << std::endl;
}

/*-------------------------------------------------------*/
/*  Get function for the slave and master dof row map    */
/*                                           wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::DofRowMaps()
{
  // slave dof rowmap
  gsdofrowmap_ = adaptermeshsliding_->SlaveDofMap();

  // master dof rowmap
  gmdofrowmap_ = adaptermeshsliding_->MasterDofMap();
}

/*-------------------------------------------------------*/
/*  Get function for the P matrix            wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ALE::Meshsliding::GetMortarTrafo()
{
  return adaptermeshsliding_->GetMortarTrafo();
}

/*-------------------------------------------------------*/
/*  Condensation operation for a block matrix            */
/*                                          wirtz 02/16  */
/*-------------------------------------------------------*/
void ALE::Meshsliding::CondensationOperationBlockMatrix(
    Teuchos::RCP<LINALG::SparseOperator>& sysmat,  ///> sysmat established by the element routine
    Teuchos::RCP<Epetra_Vector>& residual,         ///> residual established by the element routine
    Teuchos::RCP<Epetra_Vector>& dispnp)           ///> current displacement vector
{
  // condensation operation for a block matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  SplitVector(residual, splitres);

  /**********************************************************************/
  /* Evaluate mortar matrices                                           */
  /**********************************************************************/

  // evaluate all mortar matrices and vectors necessary for mesh sliding
  adaptermeshsliding_->EvaluateSliding("displacement", dispnp, lm_);

  // get all mortar matrices necessary for mesh sliding
  Teuchos::RCP<LINALG::SparseMatrix> Aco_mm;
  Teuchos::RCP<LINALG::SparseMatrix> Aco_ms;
  Teuchos::RCP<LINALG::SparseMatrix> Aco_sm;
  Teuchos::RCP<LINALG::SparseMatrix> Aco_ss;
  Teuchos::RCP<LINALG::SparseMatrix> N_m;
  Teuchos::RCP<LINALG::SparseMatrix> N_s;
  GetMortarMatrices(Aco_mm, Aco_ms, Aco_sm, Aco_ss, N_m, N_s);
  Teuchos::RCP<LINALG::SparseMatrix> P = GetMortarTrafo();
  Teuchos::RCP<LINALG::SparseMatrix> T = adaptermeshsliding_->TMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> H = adaptermeshsliding_->HMatrix();
  Teuchos::RCP<Epetra_Vector> gap = adaptermeshsliding_->Gap();

  /**********************************************************************/
  /* Step 1: Add sliding stiffnesses to the system matrix               */
  /**********************************************************************/

  sysmatnew->Matrix(1, 1).UnComplete();
  sysmatnew->Matrix(1, 1).Add(*Aco_mm, false, 1.0, 1.0);

  sysmatnew->Matrix(1, 2).UnComplete();
  sysmatnew->Matrix(1, 2).Add(*Aco_ms, false, 1.0, 1.0);

  sysmatnew->Matrix(2, 1).UnComplete();
  sysmatnew->Matrix(2, 1).Add(*Aco_sm, false, 1.0, 1.0);

  sysmatnew->Matrix(2, 2).UnComplete();
  sysmatnew->Matrix(2, 2).Add(*Aco_ss, false, 1.0, 1.0);

  sysmatnew->Complete();

  /**********************************************************************/
  /* Store some quantities for static condensation of LM                */
  /**********************************************************************/

  A_ss_ = Teuchos::rcp(new LINALG::SparseMatrix(sysmatnew->Matrix(2, 2)));
  A_sm_ = Teuchos::rcp(new LINALG::SparseMatrix(sysmatnew->Matrix(2, 1)));
  A_sn_ = Teuchos::rcp(new LINALG::SparseMatrix(sysmatnew->Matrix(2, 0)));
  D_inv_ = adaptermeshsliding_->GetDinvMatrix();
  rs_ = Teuchos::rcp(new Epetra_Vector(*splitres[2]));

  /**********************************************************************/
  /* Step 2: condensate the system matrix                               */
  /**********************************************************************/

  // ----------------------------------------------------------------------
  // | A_nn | A_nm | A_ns |  0   |        | A_nn | A_nm | A_ns |
  // | A_mn | A_mm | A_ms | -M^T |   ->   | A_mn'| A_mm'| A_ms'|
  // | A_sn | A_sm | A_ss |  D^T |        | A_sn'| A_sm'| A_ss'|
  // |   0  |  N_m |  N_s |  0   |
  // |   0  |  0   |  H_s |  T   |

  //----------------------------------------------------------- FIRST LINE

  // do nothing

  //---------------------------------------------------------- SECOND LINE

  // compute modification for block mn       (+ P^T * A_sn)
  Teuchos::RCP<LINALG::SparseMatrix> Amn_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 0), false, false, false, true);

  // Add modification block to mn
  sysmatnew->Matrix(1, 0).UnComplete();  // sonst kann ich auf den Block nichts neues draufaddieren
  sysmatnew->Matrix(1, 0).Add(*Amn_mod, false, 1.0, 1.0);

  // compute modification for block mm       (+ P^T * A_sm)
  Teuchos::RCP<LINALG::SparseMatrix> Amm_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 1), false, false, false, true);

  // Add modification block to mm
  sysmatnew->Matrix(1, 1).UnComplete();
  sysmatnew->Matrix(1, 1).Add(*Amm_mod, false, 1.0, 1.0);

  // compute modification for block ms       (+ P^T * A_ss)
  Teuchos::RCP<LINALG::SparseMatrix> Ams_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 2), false, false, false, true);

  // Add modification block to ms
  sysmatnew->Matrix(1, 2).UnComplete();
  sysmatnew->Matrix(1, 2).Add(*Ams_mod, false, 1.0, 1.0);

  //----------------------------------------------------------- THIRD LINE

  // compute replacement for block sn      - (T * D^(-1) * A_sn)
  Teuchos::RCP<LINALG::SparseMatrix> Asn_mod_interm =
      MLMultiply(*D_inv_, false, sysmatnew->Matrix(2, 0), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> Asn_mod =
      MLMultiply(*T, false, *Asn_mod_interm, false, false, false, true);

  // Replace sn block with (negative) modification block
  sysmatnew->Matrix(2, 0).UnComplete();
  sysmatnew->Matrix(2, 0).Add(*Asn_mod, false, -1.0, 0.0);

  // compute replacement for block sm      - (T * D^(-1) * A_sm)   +  N_m
  Teuchos::RCP<LINALG::SparseMatrix> Asm_mod_interm =
      MLMultiply(*D_inv_, false, sysmatnew->Matrix(2, 1), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> Asm_mod =
      MLMultiply(*T, false, *Asm_mod_interm, false, false, false, true);

  // Replace sm block with (negative) modification block
  sysmatnew->Matrix(2, 1).UnComplete();
  sysmatnew->Matrix(2, 1).Add(*Asm_mod, false, -1.0, 0.0);
  sysmatnew->Matrix(2, 1).Add(*N_m, false, 1.0, 1.0);

  // compute replacement for block ss      (- T * D^(-1) *A_ss)   +  H  +  N_s
  Teuchos::RCP<LINALG::SparseMatrix> Ass_mod_interm =
      MLMultiply(*D_inv_, false, sysmatnew->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> Ass_mod =
      MLMultiply(*T, false, *Ass_mod_interm, false, false, false, true);

  // Replace ss block with (negative) modification block
  sysmatnew->Matrix(2, 2).UnComplete();
  sysmatnew->Matrix(2, 2).Add(*Ass_mod, false, -1.0, 0.0);
  sysmatnew->Matrix(2, 2).Add(*N_s, false, 1.0, 1.0);
  sysmatnew->Matrix(2, 2).Add(*H, false, 1.0, 1.0);

  // complete matrix
  sysmatnew->Complete();

  /**********************************************************************/
  /* Step 3: condensate the residual                                    */
  /***********************************************************************/

  //----------------------------------------------------------- FIRST LINE

  // do nothing

  //---------------------------------------------------------- SECOND LINE

  // r_m: add P^T*r_s
  Teuchos::RCP<Epetra_Vector> rm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  P->Multiply(true, *(splitres[2]), *rm_mod);

  // export and add r_m subvector to residual
  Teuchos::RCP<Epetra_Vector> rm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*rm_mod, *rm_modexp);
  residual->Update(1.0, *rm_modexp, 1.0);

  //----------------------------------------------------------- THIRD LINE

  // r_s: * 0
  Teuchos::RCP<Epetra_Vector> rs_delete = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*(splitres[2]), *rs_delete);
  residual->Update(-1.0, *rs_delete, 1.0);

  // r_s: add - T*D^(-1)*r_s
  Teuchos::RCP<Epetra_Vector> rs_mod_interm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  D_inv_->Multiply(false, *rs_, *rs_mod_interm);
  Teuchos::RCP<Epetra_Vector> rs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  T->Multiply(false, *rs_mod_interm, *rs_mod);

  // export and subtract rs_mod from residual
  Teuchos::RCP<Epetra_Vector> rs_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*rs_mod, *rs_modexp);
  residual->Update(-1.0, *rs_modexp, 1.0);

  // r_s: add gap
  Teuchos::RCP<Epetra_Vector> g_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*gap, *g_exp);
  residual->Update(1.0, *g_exp, 1.0);

  return;
}

/*-------------------------------------------------------*/
/*  Get functions for the mortar matrices    wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshsliding::GetMortarMatrices(Teuchos::RCP<LINALG::SparseMatrix>& Aco_mm,
    Teuchos::RCP<LINALG::SparseMatrix>& Aco_ms, Teuchos::RCP<LINALG::SparseMatrix>& Aco_sm,
    Teuchos::RCP<LINALG::SparseMatrix>& Aco_ss, Teuchos::RCP<LINALG::SparseMatrix>& N_m,
    Teuchos::RCP<LINALG::SparseMatrix>& N_s)
{
  Teuchos::RCP<LINALG::SparseMatrix> MLin = adaptermeshsliding_->MLinMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> DLin = adaptermeshsliding_->DLinMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> N = adaptermeshsliding_->NMatrix();

  SplitMortarMatrix(MLin, Aco_mm, Aco_ms, gmdofrowmap_);
  SplitMortarMatrix(DLin, Aco_sm, Aco_ss, gsdofrowmap_);
  SplitMortarMatrix(N, N_m, N_s, gsdofrowmap_);
}

/*-------------------------------------------------------*/
/*  Split the mortar matrix into its slave and its       */
/*  master part                              wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::SplitMortarMatrix(Teuchos::RCP<LINALG::SparseMatrix>& MortarMatrix,
    Teuchos::RCP<LINALG::SparseMatrix>& MasterMatrix,
    Teuchos::RCP<LINALG::SparseMatrix>& SlaveMatrix, Teuchos::RCP<const Epetra_Map>& dofrowmapconst)
{
  // dummy Matrices for second row and dummy map
  Teuchos::RCP<LINALG::SparseMatrix> temp21;
  Teuchos::RCP<LINALG::SparseMatrix> temp22;
  Teuchos::RCP<Epetra_Map> dummy;

  // const casts
  Teuchos::RCP<Epetra_Map> gmdofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(gmdofrowmap_);
  Teuchos::RCP<Epetra_Map> gsdofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(gsdofrowmap_);
  Teuchos::RCP<Epetra_Map> dofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(dofrowmapconst);

  // split matrix operation
  bool suc = LINALG::SplitMatrix2x2(MortarMatrix, dofrowmap, dummy, gmdofrowmap, gsdofrowmap,
      MasterMatrix, SlaveMatrix, temp21, temp22);

  if (!suc) dserror("\nCould not split Mortar Matriz!\n");

  // Complete Matrices
  MasterMatrix->Complete();
  SlaveMatrix->Complete();
}

/*-------------------------------------------------------*/
/*  Recover method for Lagrange multipliers  wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshsliding::Recover(Teuchos::RCP<Epetra_Vector>& inc)
{
  // recover lm from the condensation:
  // lm = - D^(-1) * (r_s + A_ss*d_s + A_sm*d_m + A_sn*d_n)
  // Note: A_sm is non zero because of the linearization of the D-Matrix
  // with respect to master dofs is already included

  // split displacement increment
  std::vector<Teuchos::RCP<Epetra_Vector>> splitinc(3);
  SplitVector(inc, splitinc);

  Teuchos::RCP<Epetra_Vector> lm_temp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

  // r_s
  lm_temp->Update(1.0, *rs_, 1.0);

  // + A_ss*d_s
  A_ss_->Multiply(false, *(splitinc[2]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // + A_sm*d_m
  A_sm_->Multiply(false, *(splitinc[1]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // + A_sn*d_n
  A_sn_->Multiply(false, *(splitinc[0]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // - D^(-1) *
  D_inv_->Multiply(false, *lm_temp, *lm_);
  lm_->Scale(-1.0);
}

/*-------------------------------------------------------*/
/*  Solve ALE mesh sliding problem           wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
int ALE::Meshsliding::SolveMeshtying(LINALG::Solver& solver,
    Teuchos::RCP<LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector>& disi,
    Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3)   Solve ALE mesh sliding problem");

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

  Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::null;

  mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap_, 108, false, true));

  int errorcode = 0;

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3.1)   - Preparation");

    mergedmatrix = sysmatnew->Merge();
  }

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3.2)   - Solve");

    errorcode = solver_.Solve(mergedmatrix->EpetraOperator(), disi, residual, true);

    Recover(disi);
  }
  return errorcode;
}
