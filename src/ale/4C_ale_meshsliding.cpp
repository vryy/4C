/*--------------------------------------------------------------------------*/
/*! \file

\brief Mesh sliding for ale problems


\level 3
*/
/*--------------------------------------------------------------------------*/

#include "4C_ale_meshsliding.hpp"

#include "4C_adapter_coupling_nonlin_mortar.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

ALE::Meshsliding::Meshsliding(Teuchos::RCP<Core::FE::Discretization> dis,
    Core::LinAlg::Solver& solver, int msht, int nsd, const UTILS::MapExtractor* surfacesplitter)
    : Meshtying(dis, solver, msht, nsd, surfacesplitter)
{
}


/*-------------------------------------------------------*/
/*  Call the constructor and the setup of the mortar     */
/*  coupling adapter                         wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::adapter_mortar(std::vector<int> coupleddof)
{
  adaptermeshsliding_ = Teuchos::rcp(new Adapter::CouplingNonLinMortar(
      Global::Problem::instance()->n_dim(), Global::Problem::instance()->mortar_coupling_params(),
      Global::Problem::instance()->contact_dynamic_params(),
      Global::Problem::instance()->spatial_approximation_type()));

  // Setup and Output of Nonlinear meshtying adapter
  adaptermeshsliding_->setup(discret_, discret_, coupleddof, "Mortar");
}

/*-------------------------------------------------------*/
/*  Setup mesh sliding problem               wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseOperator> ALE::Meshsliding::setup(
    std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  Teuchos::RCP<Core::LinAlg::SparseOperator> mat = Meshtying::setup(coupleddof, dispnp);

  lm_ = Core::LinAlg::CreateVector(*gsdofrowmap_, true);

  return mat;
}

/*-------------------------------------------------------*/
/*  Compare the size of the slave and master dof row map */
/*                                           wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::compare_num_dof()
{
  int numdofmaster = (adaptermeshsliding_->master_dof_map())->NumGlobalElements();
  int numdofslave = (adaptermeshsliding_->slave_dof_map())->NumGlobalElements();

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
void ALE::Meshsliding::dof_row_maps()
{
  // slave dof rowmap
  gsdofrowmap_ = adaptermeshsliding_->slave_dof_map();

  // master dof rowmap
  gmdofrowmap_ = adaptermeshsliding_->master_dof_map();
}

/*-------------------------------------------------------*/
/*  Get function for the P matrix            wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> ALE::Meshsliding::get_mortar_matrix_p()
{
  return adaptermeshsliding_->get_mortar_matrix_p();
}

/*-------------------------------------------------------*/
/*  Condensation operation for a block matrix            */
/*                                          wirtz 02/16  */
/*-------------------------------------------------------*/
void ALE::Meshsliding::condensation_operation_block_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator>&
        sysmat,                             ///> sysmat established by the element routine
    Teuchos::RCP<Epetra_Vector>& residual,  ///> residual established by the element routine
    Teuchos::RCP<Epetra_Vector>& dispnp)    ///> current displacement vector
{
  // condensation operation for a block matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);

  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  split_vector(residual, splitres);

  /**********************************************************************/
  /* Evaluate mortar matrices                                           */
  /**********************************************************************/

  // evaluate all mortar matrices and vectors necessary for mesh sliding
  adaptermeshsliding_->evaluate_sliding("displacement", dispnp, lm_);

  // get all mortar matrices necessary for mesh sliding
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Aco_mm;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Aco_ms;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Aco_sm;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Aco_ss;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> N_m;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> N_s;
  get_mortar_matrices(Aco_mm, Aco_ms, Aco_sm, Aco_ss, N_m, N_s);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = get_mortar_matrix_p();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> T = adaptermeshsliding_->t_matrix();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> H = adaptermeshsliding_->h_matrix();
  Teuchos::RCP<Epetra_Vector> gap = adaptermeshsliding_->gap();

  /**********************************************************************/
  /* Step 1: Add sliding stiffnesses to the system matrix               */
  /**********************************************************************/

  sysmatnew->matrix(1, 1).un_complete();
  sysmatnew->matrix(1, 1).add(*Aco_mm, false, 1.0, 1.0);

  sysmatnew->matrix(1, 2).un_complete();
  sysmatnew->matrix(1, 2).add(*Aco_ms, false, 1.0, 1.0);

  sysmatnew->matrix(2, 1).un_complete();
  sysmatnew->matrix(2, 1).add(*Aco_sm, false, 1.0, 1.0);

  sysmatnew->matrix(2, 2).un_complete();
  sysmatnew->matrix(2, 2).add(*Aco_ss, false, 1.0, 1.0);

  sysmatnew->complete();

  /**********************************************************************/
  /* Store some quantities for static condensation of LM                */
  /**********************************************************************/

  a_ss_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(sysmatnew->matrix(2, 2)));
  a_sm_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(sysmatnew->matrix(2, 1)));
  a_sn_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(sysmatnew->matrix(2, 0)));
  d_inv_ = adaptermeshsliding_->get_mortar_matrix_dinv();
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
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Amn_mod =
      MLMultiply(*P, true, sysmatnew->matrix(2, 0), false, false, false, true);

  // Add modification block to mn
  sysmatnew->matrix(1, 0).un_complete();  // sonst kann ich auf den Block nichts neues draufaddieren
  sysmatnew->matrix(1, 0).add(*Amn_mod, false, 1.0, 1.0);

  // compute modification for block mm       (+ P^T * A_sm)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Amm_mod =
      MLMultiply(*P, true, sysmatnew->matrix(2, 1), false, false, false, true);

  // Add modification block to mm
  sysmatnew->matrix(1, 1).un_complete();
  sysmatnew->matrix(1, 1).add(*Amm_mod, false, 1.0, 1.0);

  // compute modification for block ms       (+ P^T * A_ss)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Ams_mod =
      MLMultiply(*P, true, sysmatnew->matrix(2, 2), false, false, false, true);

  // Add modification block to ms
  sysmatnew->matrix(1, 2).un_complete();
  sysmatnew->matrix(1, 2).add(*Ams_mod, false, 1.0, 1.0);

  //----------------------------------------------------------- THIRD LINE

  // compute replacement for block sn      - (T * D^(-1) * A_sn)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Asn_mod_interm =
      MLMultiply(*d_inv_, false, sysmatnew->matrix(2, 0), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Asn_mod =
      MLMultiply(*T, false, *Asn_mod_interm, false, false, false, true);

  // Replace sn block with (negative) modification block
  sysmatnew->matrix(2, 0).un_complete();
  sysmatnew->matrix(2, 0).add(*Asn_mod, false, -1.0, 0.0);

  // compute replacement for block sm      - (T * D^(-1) * A_sm)   +  N_m
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Asm_mod_interm =
      MLMultiply(*d_inv_, false, sysmatnew->matrix(2, 1), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Asm_mod =
      MLMultiply(*T, false, *Asm_mod_interm, false, false, false, true);

  // Replace sm block with (negative) modification block
  sysmatnew->matrix(2, 1).un_complete();
  sysmatnew->matrix(2, 1).add(*Asm_mod, false, -1.0, 0.0);
  sysmatnew->matrix(2, 1).add(*N_m, false, 1.0, 1.0);

  // compute replacement for block ss      (- T * D^(-1) *A_ss)   +  H  +  N_s
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Ass_mod_interm =
      MLMultiply(*d_inv_, false, sysmatnew->matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Ass_mod =
      MLMultiply(*T, false, *Ass_mod_interm, false, false, false, true);

  // Replace ss block with (negative) modification block
  sysmatnew->matrix(2, 2).un_complete();
  sysmatnew->matrix(2, 2).add(*Ass_mod, false, -1.0, 0.0);
  sysmatnew->matrix(2, 2).add(*N_s, false, 1.0, 1.0);
  sysmatnew->matrix(2, 2).add(*H, false, 1.0, 1.0);

  // complete matrix
  sysmatnew->complete();

  /**********************************************************************/
  /* Step 3: condensate the residual                                    */
  /***********************************************************************/

  //----------------------------------------------------------- FIRST LINE

  // do nothing

  //---------------------------------------------------------- SECOND LINE

  // r_m: add P^T*r_s
  Teuchos::RCP<Epetra_Vector> rm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  P->multiply(true, *(splitres[2]), *rm_mod);

  // export and add r_m subvector to residual
  Teuchos::RCP<Epetra_Vector> rm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*rm_mod, *rm_modexp);
  residual->Update(1.0, *rm_modexp, 1.0);

  //----------------------------------------------------------- THIRD LINE

  // r_s: * 0
  Teuchos::RCP<Epetra_Vector> rs_delete = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*(splitres[2]), *rs_delete);
  residual->Update(-1.0, *rs_delete, 1.0);

  // r_s: add - T*D^(-1)*r_s
  Teuchos::RCP<Epetra_Vector> rs_mod_interm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  d_inv_->multiply(false, *rs_, *rs_mod_interm);
  Teuchos::RCP<Epetra_Vector> rs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  T->multiply(false, *rs_mod_interm, *rs_mod);

  // export and subtract rs_mod from residual
  Teuchos::RCP<Epetra_Vector> rs_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*rs_mod, *rs_modexp);
  residual->Update(-1.0, *rs_modexp, 1.0);

  // r_s: add gap
  Teuchos::RCP<Epetra_Vector> g_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*gap, *g_exp);
  residual->Update(1.0, *g_exp, 1.0);

  return;
}

/*-------------------------------------------------------*/
/*  Get functions for the mortar matrices    wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshsliding::get_mortar_matrices(Teuchos::RCP<Core::LinAlg::SparseMatrix>& Aco_mm,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& Aco_ms,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& Aco_sm,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& Aco_ss, Teuchos::RCP<Core::LinAlg::SparseMatrix>& N_m,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& N_s)
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> MLin = adaptermeshsliding_->m_lin_matrix();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> DLin = adaptermeshsliding_->d_lin_matrix();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> N = adaptermeshsliding_->n_matrix();

  split_mortar_matrix(MLin, Aco_mm, Aco_ms, gmdofrowmap_);
  split_mortar_matrix(DLin, Aco_sm, Aco_ss, gsdofrowmap_);
  split_mortar_matrix(N, N_m, N_s, gsdofrowmap_);
}

/*-------------------------------------------------------*/
/*  Split the mortar matrix into its slave and its       */
/*  master part                              wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshsliding::split_mortar_matrix(Teuchos::RCP<Core::LinAlg::SparseMatrix>& MortarMatrix,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& MasterMatrix,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& SlaveMatrix,
    Teuchos::RCP<const Epetra_Map>& dofrowmapconst)
{
  // dummy Matrices for second row and dummy map
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp21;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp22;
  Teuchos::RCP<Epetra_Map> dummy;

  // const casts
  Teuchos::RCP<Epetra_Map> gmdofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(gmdofrowmap_);
  Teuchos::RCP<Epetra_Map> gsdofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(gsdofrowmap_);
  Teuchos::RCP<Epetra_Map> dofrowmap = Teuchos::rcp_const_cast<Epetra_Map>(dofrowmapconst);

  // split matrix operation
  bool suc = Core::LinAlg::SplitMatrix2x2(MortarMatrix, dofrowmap, dummy, gmdofrowmap, gsdofrowmap,
      MasterMatrix, SlaveMatrix, temp21, temp22);

  if (!suc) FOUR_C_THROW("\nCould not split Mortar Matriz!\n");

  // Complete Matrices
  MasterMatrix->complete();
  SlaveMatrix->complete();
}

/*-------------------------------------------------------*/
/*  Recover method for Lagrange multipliers  wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshsliding::recover(Teuchos::RCP<Epetra_Vector>& inc)
{
  // recover lm from the condensation:
  // lm = - D^(-1) * (r_s + A_ss*d_s + A_sm*d_m + A_sn*d_n)
  // Note: A_sm is non zero because of the linearization of the D-Matrix
  // with respect to master dofs is already included

  // split displacement increment
  std::vector<Teuchos::RCP<Epetra_Vector>> splitinc(3);
  split_vector(inc, splitinc);

  Teuchos::RCP<Epetra_Vector> lm_temp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

  // r_s
  lm_temp->Update(1.0, *rs_, 1.0);

  // + A_ss*d_s
  a_ss_->multiply(false, *(splitinc[2]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // + A_sm*d_m
  a_sm_->multiply(false, *(splitinc[1]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // + A_sn*d_n
  a_sn_->multiply(false, *(splitinc[0]), *mod);
  lm_temp->Update(1.0, *mod, 1.0);

  // - D^(-1) *
  d_inv_->multiply(false, *lm_temp, *lm_);
  lm_->Scale(-1.0);
}

/*-------------------------------------------------------*/
/*  Solve ALE mesh sliding problem           wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
int ALE::Meshsliding::solve_meshtying(Core::LinAlg::Solver& solver,
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector>& disi,
    Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3)   Solve ALE mesh sliding problem");

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix = Teuchos::null;

  mergedmatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mergedmap_, 108, false, true));

  int errorcode = 0;

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3.1)   - Preparation");

    mergedmatrix = sysmatnew->merge();
  }

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshsliding:  3.2)   - Solve");

    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    errorcode = solver_.solve(mergedmatrix->epetra_operator(), disi, residual, solver_params);

    recover(disi);
  }
  return errorcode;
}

FOUR_C_NAMESPACE_CLOSE
