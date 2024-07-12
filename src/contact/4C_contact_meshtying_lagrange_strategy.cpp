/*----------------------------------------------------------------------------*/
/*! \file
\brief strategy handling mesh tying problems with Lagrange multipliers

\level 1


*/
/*----------------------------------------------------------------------------*/

#include "4C_contact_meshtying_lagrange_strategy.hpp"

#include "4C_contact_meshtying_defines.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"  // mesh initialization :-(
#include "4C_mortar_defines.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Epetra_SerialComm.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::MtLagrangeStrategy::MtLagrangeStrategy(const Epetra_Map* dof_row_map,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<Mortar::Interface>> interface, const int spatialDim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof)
    : MtAbstractStrategy(
          dof_row_map, NodeRowMap, params, interface, spatialDim, comm, alphaf, maxdof)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::mortar_coupling(const Teuchos::RCP<const Epetra_Vector>& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::MtLagrangeStrategy::mortar_coupling");

  // print message
  if (get_comm().MyPID() == 0)
  {
    std::cout << "Performing mortar coupling...............";
    fflush(stdout);
  }

  // time measurement
  get_comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // refer call to parent class
  MtAbstractStrategy::mortar_coupling(dis);

  //----------------------------------------------------------------------
  // Multiply Mortar matrices: m^ = inv(d) * m
  //----------------------------------------------------------------------
  invd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd_->extract_diagonal_copy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag->MyLength(); ++i)
    if ((*diag)[i] == 0.0) (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd_->replace_diagonal_values(*diag);
  if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = Core::LinAlg::MLMultiply(*invd_, false, *mmatrix_, false, false, false, true);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  // These modifications are applied once right here, thus the
  // following code (evaluate_meshtying, recover) remains unchanged.
  //----------------------------------------------------------------------
  if (dualquadslavetrafo())
  {
    // type of LM interpolation for quadratic elements
    Inpar::Mortar::LagMultQuad lagmultquad =
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

    if (lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // do nothing
    }
    else
    {
      // modify dmatrix_, invd_ and mhatmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
          Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
          Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
      dmatrix_ = temp1;
      invd_ = temp2;
      mhatmatrix_ = temp3;
    }
  }

  //----------------------------------------------------------------------
  // Build constraint matrix (containing D and M)
  //----------------------------------------------------------------------
  // case 1: saddle point system
  //    -> constraint matrix with rowmap=Problemmap, colmap=LMmap
  // case 2: condensed system
  //    -> no explicit constraint matrix needed
  //----------------------------------------------------------------------
  bool setup = true;
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");
  if (systype == Inpar::CONTACT::system_condensed ||
      systype == Inpar::CONTACT::system_condensed_lagmult)
    setup = false;

  // build constraint matrix only if necessary
  if (setup)
  {
    // first setup
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constrmt =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));
    constrmt->add(*dmatrix_, true, 1.0, 1.0);
    constrmt->add(*mmatrix_, true, -1.0, 1.0);
    constrmt->complete(*gsdofrowmap_, *gdisprowmap_);

    // transform parallel row distribution
    // (only necessary in the parallel redistribution case)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp;
    if (par_redist())
      temp = Mortar::MatrixRowTransform(constrmt, problem_dofs());
    else
      temp = constrmt;

    // always transform column GIDs of constraint matrix
    conmatrix_ = Mortar::MatrixColTransformGIDs(temp, glmdofrowmap_);
    conmatrix_->scale(1. - alphaf_);
  }

  dm_matrix_ = Teuchos::null;
  dm_matrix_t_ = Teuchos::null;
  lm_diag_matrix_ = Teuchos::null;

  // time measurement
  get_comm().Barrier();
  const double t_end = Teuchos::Time::wallTime() - t_start;
  if (get_comm().MyPID() == 0)
    std::cout << "in...." << std::scientific << std::setprecision(6) << t_end << " secs"
              << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  mesh initialization for rotational invariance             popp 12/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::MtLagrangeStrategy::mesh_initialization()
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::MtLagrangeStrategy::mesh_initialization");

  // get out of here if NTS algorithm is activated
  if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(params(), "ALGORITHM") ==
      Inpar::Mortar::algorithm_nts)
    return Teuchos::null;

  // print message
  if (get_comm().MyPID() == 0)
  {
    std::cout << "Performing mesh initialization...........";
    fflush(stdout);
  }

  // time measurement
  get_comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  //**********************************************************************
  // (1) get master positions on global level
  //**********************************************************************
  // fill Xmaster first
  Teuchos::RCP<Epetra_Vector> Xmaster = Core::LinAlg::CreateVector(*gmdofrowmap_, true);
  assemble_coords("master", true, Xmaster);

  //**********************************************************************
  // (2) solve for modified slave positions on global level
  //**********************************************************************
  // initialize modified slave positions
  Teuchos::RCP<Epetra_Vector> Xslavemod = Core::LinAlg::CreateVector(*gsdofrowmap_, true);

  // shape function type and type of LM interpolation for quadratic elements
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  // quadratic FE with dual LM
  if (dualquadslavetrafo())
  {
    if (lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // split T^-1
      Teuchos::RCP<Core::LinAlg::SparseMatrix> it_ss, it_sm, it_ms, it_mm;
      Core::LinAlg::SplitMatrix2x2(invtrafo_, gsdofrowmap_, gmdofrowmap_, gsdofrowmap_,
          gmdofrowmap_, it_ss, it_sm, it_ms, it_mm);

      // build lhs
      Teuchos::RCP<Core::LinAlg::SparseMatrix> lhs =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> direct =
          Core::LinAlg::MLMultiply(*dmatrix_, false, *it_ss, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> mixed =
          Core::LinAlg::MLMultiply(*mmatrix_, false, *it_ms, false, false, false, true);
      lhs->add(*direct, false, 1.0, 1.0);
      lhs->add(*mixed, false, -1.0, 1.0);
      lhs->complete();

      // build rhs
      Teuchos::RCP<Epetra_Vector> xm = Core::LinAlg::CreateVector(*gmdofrowmap_, true);
      assemble_coords("master", true, xm);
      Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      mmatrix_->multiply(false, *xm, *rhs);

      // solve with default solver

      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(
          solvparams, get_comm(), nullptr, Core::IO::Verbositylevel::standard);

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(lhs->epetra_operator(), Xslavemod, rhs, solver_params);
    }
    else
    {
      // this is trivial for dual Lagrange multipliers
      mhatmatrix_->multiply(false, *Xmaster, *Xslavemod);
    }
  }

  // other cases (quadratic FE with std LM, linear FE)
  else
  {
    // CASE A: DUAL LM SHAPE FUNCTIONS
    if (shapefcn == Inpar::Mortar::shape_dual)
    {
      // this is trivial for dual Lagrange multipliers
      mhatmatrix_->multiply(false, *Xmaster, *Xslavemod);
    }

    // CASE B: STANDARD LM SHAPE FUNCTIONS
    else if (shapefcn == Inpar::Mortar::shape_standard)
    {
      // create linear problem
      Teuchos::RCP<Epetra_Vector> rhs = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
      mmatrix_->multiply(false, *Xmaster, *rhs);

      // solve with default solver

      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(
          solvparams, get_comm(), nullptr, Core::IO::Verbositylevel::standard);

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(dmatrix_->epetra_operator(), Xslavemod, rhs, solver_params);
    }
  }

  //**********************************************************************
  // (3) perform mesh initialization node by node
  //**********************************************************************
  // this can be done in the AbstractStrategy now
  MtAbstractStrategy::mesh_initialization(Xslavemod);

  // time measurement
  get_comm().Barrier();
  const double t_end = Teuchos::Time::wallTime() - t_start;
  if (get_comm().MyPID() == 0)
  {
    std::cout << "in...." << std::scientific << std::setprecision(6) << t_end << " secs"
              << std::endl;
  }

  //**********************************************************************
  // blank symmetry rows in dinv
  //**********************************************************************
  invd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd_->extract_diagonal_copy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag->MyLength(); ++i)
    if ((*diag)[i] == 0.0) (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

  Teuchos::RCP<Epetra_Vector> lmDBC = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  Core::LinAlg::export_to(*pgsdirichtoggle_, *lmDBC);
  Teuchos::RCP<Epetra_Vector> tmp = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  tmp->Multiply(1., *diag, *lmDBC, 0.);
  diag->Update(-1., *tmp, 1.);

  // re-insert inverted diagonal into invd
  err = invd_->replace_diagonal_values(*diag);
  if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = Core::LinAlg::MLMultiply(*invd_, false, *mmatrix_, false, false, false, true);

  // return xslavemod for global problem
  return Xslavemod;
}

/*----------------------------------------------------------------------*
 |  evaluate meshtying (public)                               popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::evaluate_meshtying(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    Teuchos::RCP<Epetra_Vector> dis)
{
  // system type, shape function type and type of LM interpolation for quadratic elements
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_condensed ||
      systype == Inpar::CONTACT::system_condensed_lagmult)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual) FOUR_C_THROW("Condensation only for dual LM");

    // complete stiffness matrix
    // (this is a prerequisite for the Split2x2 methods to be called later)
    kteff->complete();

    /**********************************************************************/
    /* Split kteff into 3x3 block matrix                                  */
    /**********************************************************************/
    // we want to split k into 3 groups s,m,n = 9 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

    // temporarily we need the blocks ksmsm, ksmn, knsm
    // (FIXME: because a direct SplitMatrix3x3 is still missing!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ksmsm, ksmn, knsm;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

    // split into slave/master part + structure part
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);

    /**********************************************************************/
    /* Apply basis transformation to K and f                              */
    /* (currently only needed for quadratic FE with linear dual LM)       */
    /**********************************************************************/
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // basis transformation
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye = Core::LinAlg::Eye(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (par_redist())
        trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
      systrafo->add(*trafo_, false, 1.0, 1.0);
      systrafo->complete();

      // apply basis transformation to K and f
      kteffmatrix =
          Core::LinAlg::MLMultiply(*kteffmatrix, false, *systrafo, false, false, false, true);
      kteffmatrix =
          Core::LinAlg::MLMultiply(*systrafo, true, *kteffmatrix, false, false, false, true);
      systrafo->multiply(true, *feff, *feff);
    }

    if (par_redist())
    {
      // split and transform to redistributed maps
      Core::LinAlg::SplitMatrix2x2(kteffmatrix, pgsmdofrowmap_, gndofrowmap_, pgsmdofrowmap_,
          gndofrowmap_, ksmsm, ksmn, knsm, knn);
      ksmsm = Mortar::matrix_row_col_transform(ksmsm, gsmdofrowmap_, gsmdofrowmap_);
      ksmn = Mortar::MatrixRowTransform(ksmn, gsmdofrowmap_);
      knsm = Mortar::MatrixColTransform(knsm, gsmdofrowmap_);
    }
    else
    {
      // only split, no need to transform
      Core::LinAlg::SplitMatrix2x2(kteffmatrix, gsmdofrowmap_, gndofrowmap_, gsmdofrowmap_,
          gndofrowmap_, ksmsm, ksmn, knsm, knn);
    }

    // further splits into slave part + master part
    Core::LinAlg::SplitMatrix2x2(
        ksmsm, gsdofrowmap_, gmdofrowmap_, gsdofrowmap_, gmdofrowmap_, kss, ksm, kms, kmm);
    Core::LinAlg::SplitMatrix2x2(
        ksmn, gsdofrowmap_, gmdofrowmap_, gndofrowmap_, tempmap, ksn, tempmtx1, kmn, tempmtx2);
    Core::LinAlg::SplitMatrix2x2(
        knsm, gndofrowmap_, tempmap, gsdofrowmap_, gmdofrowmap_, kns, knm, tempmtx1, tempmtx2);

    /**********************************************************************/
    /* Split feff into 3 subvectors                                       */
    /**********************************************************************/
    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fs, fm, fn;

    // temporarily we need the group sm
    Teuchos::RCP<Epetra_Vector> fsm;

    // do the vector splitting smn -> sm+n
    Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);

    // we want to split fsm into 2 groups s,m
    fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

    // do the vector splitting sm -> s+m
    Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

    // store some stuff for static condensation of LM
    fs_ = fs;
    ksn_ = ksn;
    ksm_ = ksm;
    kss_ = kss;

    //***************************************************************************
    // build constraint vector
    //***************************************************************************
    // Since we enforce the meshtying constraint for the displacements u,
    // and not for the configurations x (which would also be possible in theory),
    // we avoid artificial initial stresses (+), but we might not guarantee
    // exact rotational invariance (-). However, since we always apply the
    // so-called mesh initialization procedure, we can then also guarantee
    // exact rotational invariance (+).
    //***************************************************************************

    // (nothing needs not be done here, as the constraints for meshtying are
    // LINEAR w.r.t. the displacements and in the first step, dis is zero.
    // Thus, the right hand side of the constraint rows is ALWAYS zero!)

    /**********************************************************************/
    /* Build the final K and f blocks                                     */
    /**********************************************************************/
    // knn: nothing to do

    // knm:
    Teuchos::RCP<Core::LinAlg::SparseMatrix> knmmod;
    if (systype == Inpar::CONTACT::system_condensed)
    {
      // knm: add kns*mbar
      knmmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gndofrowmap_, 100));
      knmmod->add(*knm, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> knmadd =
          Core::LinAlg::MLMultiply(*kns, false, *mhatmatrix_, false, false, false, true);
      knmmod->add(*knmadd, false, 1.0, 1.0);
      knmmod->complete(knm->domain_map(), knm->row_map());
    }

    // kmn: add T(mbar)*ksn
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmnmod->add(*kmn, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
        Core::LinAlg::MLMultiply(*mhatmatrix_, true, *ksn, false, false, false, true);
    kmnmod->add(*kmnadd, false, 1.0, 1.0);
    kmnmod->complete(kmn->domain_map(), kmn->row_map());

    // kmm: add T(mbar)*ksm
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmmmod->add(*kmm, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
        Core::LinAlg::MLMultiply(*mhatmatrix_, true, *ksm, false, false, false, true);
    kmmmod->add(*kmmadd, false, 1.0, 1.0);
    if (systype == Inpar::CONTACT::system_condensed)
    {
      // kmm: add kms*mbar + T(mbar)*kss*mbar - additionally
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd2 =
          Core::LinAlg::MLMultiply(*kms, false, *mhatmatrix_, false, false, false, true);
      kmmmod->add(*kmmadd2, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmtemp =
          Core::LinAlg::MLMultiply(*kss, false, *mhatmatrix_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd3 =
          Core::LinAlg::MLMultiply(*mhatmatrix_, true, *kmmtemp, false, false, false, true);
      kmmmod->add(*kmmadd3, false, 1.0, 1.0);
    }
    kmmmod->complete(kmm->domain_map(), kmm->row_map());

    // some modifications for kns, kms, (,ksn) ksm, kss if slave displacement increment is not
    // condensed

    // kns: nothing to do

    // kms: add T(mbar)*kss
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmsmod;
    if (systype == Inpar::CONTACT::system_condensed_lagmult)
    {
      kmsmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
      kmsmod->add(*kms, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmsadd =
          Core::LinAlg::MLMultiply(*mhatmatrix_, true, *kss, false, false, false, true);
      kmsmod->add(*kmsadd, false, 1.0, 1.0);
      kmsmod->complete(kms->domain_map(), kms->row_map());
    }

    // (ksn: do nothing as block is supposed to remain zero)

    // ksm: subtract mmatrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ksmmod;
    if (systype == Inpar::CONTACT::system_condensed_lagmult)
    {
      ksmmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100));
      ksmmod->add(*mmatrix_, false, -1.0, 1.0);  //<---- causes problems in parallel
      ksmmod->complete(ksm->domain_map(), ksm->row_map());
    }

    // kss: add dmatrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kssmod;
    if (systype == Inpar::CONTACT::system_condensed_lagmult)
    {
      kssmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100));
      kssmod->add(*dmatrix_, false, 1.0, 1.0);  //<---- causes problems in parallel
      kssmod->complete(kss->domain_map(), kss->row_map());
    }

    // fn: subtract kns*inv(D)*g
    // (nothing needs to be done, since the right hand side g is ALWAYS zero)

    // fs: subtract alphaf * old interface forces (t_n)
    Teuchos::RCP<Epetra_Vector> tempvecs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->multiply(true, *zold_, *tempvecs);
    tempvecs->Update(1.0, *fs, -alphaf_);

    // fm: add alphaf * old interface forces (t_n)
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *zold_, *tempvecm);
    fm->Update(alphaf_, *tempvecm, 1.0);

    // fm: add T(mbar)*fs
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mhatmatrix_->multiply(true, *tempvecs, *fmmod);
    fmmod->Update(1.0, *fm, 1.0);

    // fm: subtract kmsmod*inv(D)*g
    // (nothing needs to be done, since the right hand side g is ALWAYS zero)

    // RHS can remain unchanged, if slave displacement increments are not condensed
    // build identity matrix for slave dofs
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    /********************************************************************/
    /* Transform the final K blocks                                     */
    /********************************************************************/
    // The row maps of all individual matrix blocks are transformed to
    // the parallel layout of the underlying problem discretization.
    // Of course, this is only necessary in the parallel redistribution
    // case, where the meshtying interfaces have been redistributed
    // independently of the underlying problem discretization.
    if (par_redist())
    {
      kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
      kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
      onesdiag = Mortar::MatrixRowTransform(onesdiag, pgsdofrowmap_);
    }

    /**********************************************************************/
    /* Global setup of kteffnew, feffnew (including meshtying)            */
    /**********************************************************************/
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
    Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

    // add n submatrices to kteffnew
    kteffnew->add(*knn, false, 1.0, 1.0);
    if (systype == Inpar::CONTACT::system_condensed)
    {
      kteffnew->add(*knmmod, false, 1.0, 1.0);
    }
    else if (systype == Inpar::CONTACT::system_condensed_lagmult)
    {
      kteffnew->add(*knm, false, 1.0, 1.0);
      kteffnew->add(*kns, false, 1.0, 1.0);
    }

    // add m submatrices to kteffnew
    kteffnew->add(*kmnmod, false, 1.0, 1.0);
    kteffnew->add(*kmmmod, false, 1.0, 1.0);
    if (systype == Inpar::CONTACT::system_condensed_lagmult)
      kteffnew->add(*kmsmod, false, 1.0, 1.0);

    // add s submatrices to kteffnew
    if (systype == Inpar::CONTACT::system_condensed)
    {
      // add identitiy for slave increments
      kteffnew->add(*onesdiag, false, 1.0, 1.0);
    }
    else
    {
      kteffnew->add(*ksmmod, false, 1.0, 1.0);
      kteffnew->add(*kssmod, false, 1.0, 1.0);
    }

    // fill_complete kteffnew (square)
    kteffnew->complete();

    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fn, *fnexp);
    feffnew->Update(1.0, *fnexp, 1.0);

    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fmmod, *fmmodexp);
    feffnew->Update(1.0, *fmmodexp, 1.0);

    /**********************************************************************/
    /* Replace kteff and feff by kteffnew and feffnew                     */
    /**********************************************************************/
    kteff = kteffnew;
    feff = feffnew;
  }

  //**********************************************************************
  //**********************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    /**********************************************************************/
    /* Apply basis transformation to K and f                              */
    /* (currently only needed for quadratic FE with linear dual LM)       */
    /**********************************************************************/
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // basis transformation
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye = Core::LinAlg::Eye(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (par_redist())
        trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
      systrafo->add(*trafo_, false, 1.0, 1.0);
      systrafo->complete();

      // apply basis transformation to K and f
      kteff->complete();
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(
              *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
      kteffnew =
          Core::LinAlg::MLMultiply(*kteffmatrix, false, *systrafo, false, false, false, true);
      kteffnew = Core::LinAlg::MLMultiply(*systrafo, true, *kteffnew, false, false, false, true);
      kteff = kteffnew;
      systrafo->multiply(true, *feff, *feff);
    }

    // add meshtying force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->multiply(true, *z_, *fs);
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fs, *fsexp);
    feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *z_, *fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fm, *fmexp);
    feff->Update(1.0 - alphaf_, *fmexp, 1.0);

    // add old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->multiply(true, *zold_, *fsold);
    Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fsold, *fsoldexp);
    feff->Update(-alphaf_, *fsoldexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *zold_, *fmold);
    Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fmold, *fmoldexp);
    feff->Update(alphaf_, *fmoldexp, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::build_saddle_point_system(
    Teuchos::RCP<Core::LinAlg::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd,
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,
    Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
    Teuchos::RCP<Epetra_Vector>& blockrhs)
{
  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle =
      Teuchos::rcp(new Epetra_Vector(*(dbcmaps->full_map())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->cond_map())));
  temp->PutScalar(1.0);
  Core::LinAlg::export_to(*temp, *dirichtoggle);

  //**********************************************************************
  // prepare saddle point system
  //**********************************************************************
  // get system type
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");

  // the standard stiffness matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stiffmt =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kdd);

  // initialize merged system (matrix, rhs, sol)
  Teuchos::RCP<Epetra_Map> mergedmap = Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
  Teuchos::RCP<Epetra_Vector> mergedrhs = Core::LinAlg::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector> mergedsol = Core::LinAlg::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector> mergedzeros = Core::LinAlg::CreateVector(*mergedmap);

  //**********************************************************************
  // finalize matrix and vector blocks
  //**********************************************************************
  // get constraint matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> constrmt = conmatrix_;

  // build constraint rhs (=empty)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  constrrhs_ = constrrhs;  // set constraint rhs vector

  //**********************************************************************
  // build and solve saddle point system
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_saddlepoint)
  {
    // build transposed constraint matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> trconstrmt =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*glmdofrowmap_, 100, false, true));
    trconstrmt->add(*constrmt, true, 1.0, 0.0);
    trconstrmt->complete(*problem_dofs(), *glmdofrowmap_);

    // apply Dirichlet conditions to (0,1) block
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*problem_dofs(), true));
    Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
    Core::LinAlg::apply_dirichlet_to_system(*stiffmt, *sold, *rhscopy, *zeros, *dirichtoggle);
    constrmt->apply_dirichlet(*dirichtoggle, false);

    // row map (equals domain map) extractor
    Core::LinAlg::MapExtractor mapext(*mergedmap, glmdofrowmap_, problem_dofs());

    // build block matrix for SIMPLER
    blockMat =
        Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
            mapext, mapext, 81, false, false));
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> mat =
        Teuchos::rcp_dynamic_cast<
            Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(blockMat);

    mat->assign(0, 0, Core::LinAlg::View, *stiffmt);
    mat->assign(0, 1, Core::LinAlg::View, *constrmt);
    mat->assign(1, 0, Core::LinAlg::View, *trconstrmt);
    mat->complete();

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*fd, *fresmexp);
    mergedrhs->Update(1.0, *fresmexp, 1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*constrrhs, *constrexp);
    mergedrhs->Update(1.0, *constrexp, 1.0);

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*dirichtoggle, *dirichtoggleexp);
    Core::LinAlg::apply_dirichlet_to_system(*mergedsol, *mergedrhs, *mergedzeros, *dirichtoggleexp);

    // make solver SIMPLER-ready
    // solver.Params().set<bool>("MESHTYING", true); // flag makes sure that SIMPLER sets correct
    // Teuchos::null space for constraint equations

    blocksol = mergedsol;
    blockrhs = mergedrhs;
  }

  //**********************************************************************
  // invalid system types
  //**********************************************************************
  else
    FOUR_C_THROW("Invalid system type in SaddlePontSolve");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::update_displacements_and_l_mincrements(
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol)
{
  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  Teuchos::RCP<Epetra_Map> mergedmap = Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
  Core::LinAlg::MapExtractor mapext(*mergedmap, problem_dofs(), glmdofrowmap_);
  mapext.extract_cond_vector(blocksol, sold);
  mapext.extract_other_vector(blocksol, sollm);
  sollm->ReplaceMap(*gsdofrowmap_);

  zincr_->Update(1.0, *sollm, 0.0);
  z_->Update(1.0, *zincr_, 1.0);

  return;
}


/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::recover(Teuchos::RCP<Epetra_Vector> disi)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::MtLagrangeStrategy::recover");

  // system type, shape function type and type of LM interpolation for quadratic elements
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_condensed ||
      systype == Inpar::CONTACT::system_condensed_lagmult)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual) FOUR_C_THROW("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disin);

    /**********************************************************************/
    /* Update slave increment \Delta d_s                                  */
    /**********************************************************************/

    if (systype == Inpar::CONTACT::system_condensed)
    {
      mhatmatrix_->multiply(false, *disim, *disis);
      Teuchos::RCP<Epetra_Vector> disisexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*disis, *disisexp);
      disi->Update(1.0, *disisexp, 1.0);
    }

    /**********************************************************************/
    /* Undo basis transformation to solution                              */
    /* (currently only needed for quadratic FE with linear dual LM)       */
    /**********************************************************************/
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // undo basis transformation to solution
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye = Core::LinAlg::Eye(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (par_redist())
        trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
      systrafo->add(*trafo_, false, 1.0, 1.0);
      systrafo->complete();
      systrafo->multiply(false, *disi, *disi);
    }

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/

    // approximate update
    // invd_->Multiply(false,*fs_,*z_);
    // full update
    z_->Update(1.0, *fs_, 0.0);
    Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    kss_->multiply(false, *disis, *mod);
    z_->Update(-1.0, *mod, 1.0);
    ksm_->multiply(false, *disim, *mod);
    z_->Update(-1.0, *mod, 1.0);
    ksn_->multiply(false, *disin, *mod);
    z_->Update(-1.0, *mod, 1.0);
    dmatrix_->multiply(true, *zold_, *mod);
    z_->Update(-alphaf_, *mod, 1.0);
    Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
    invd_->multiply(true, *zcopy, *z_);
    z_->Scale(1 / (1 - alphaf_));
  }

  //**********************************************************************
  //**********************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // do nothing (z_ was part of solution already)

    /**********************************************************************/
    /* Undo basis transformation to solution                              */
    /* (currently only needed for quadratic FE with linear dual LM)       */
    /**********************************************************************/
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // undo basis transformation to solution
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye = Core::LinAlg::Eye(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (par_redist())
        trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
      systrafo->add(*trafo_, false, 1.0, 1.0);
      systrafo->complete();
      systrafo->multiply(false, *disi, *disi);
    }
  }

  // store updated LM into nodes
  store_nodal_quantities(Mortar::StrategyBase::lmupdate);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtLagrangeStrategy::evaluate_force(const Teuchos::RCP<const Epetra_Vector> dis)
{
  if (f_.is_null()) f_ = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  f_->PutScalar(0.);

  if (system_type() != Inpar::CONTACT::system_condensed)
  {
    // add meshtying force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (dmatrix_->multiply(true, *z_, *fs)) FOUR_C_THROW("multiply failed");
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fs, *fsexp);
    f_->Update(1.0, *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *z_, *fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fm, *fmexp);
    f_->Update(-1.0, *fmexp, 1.0);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtLagrangeStrategy::evaluate_stiff(const Teuchos::RCP<const Epetra_Vector> dis)
{
  if (!dm_matrix_.is_null() && !dm_matrix_t_.is_null() && !lm_diag_matrix_.is_null()) return true;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> constrmt =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));
  constrmt->add(*dmatrix_, true, 1.0, 1.0);
  constrmt->add(*mmatrix_, true, -1.0, 1.0);
  constrmt->complete(*gsdofrowmap_, *gdisprowmap_);

  // transform parallel row distribution
  // (only necessary in the parallel redistribution case)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp;
  if (par_redist())
    temp = Mortar::MatrixRowTransform(constrmt, problem_dofs());
  else
    temp = constrmt;


  // always transform column GIDs of constraint matrix
  dm_matrix_ = Mortar::MatrixColTransformGIDs(temp, lm_do_f_row_map_ptr(true));
  dm_matrix_t_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*lm_do_f_row_map_ptr(true), 100, false, true));
  dm_matrix_t_->add(*dm_matrix_, true, 1., 0.);
  dm_matrix_t_->complete(*problem_dofs(), *lm_do_f_row_map_ptr(true));

  dm_matrix_->scale(1. - alphaf_);
  dm_matrix_t_->scale(1. - alphaf_);

  lm_diag_matrix_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*lm_do_f_row_map_ptr(true), 1, false, true));
  lm_diag_matrix_->complete();


  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
  {
    systrafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> eye = Core::LinAlg::Eye(*gndofrowmap_);
    systrafo_->add(*eye, false, 1.0, 1.0);
    if (par_redist())
      trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
    systrafo_->add(*trafo_, false, 1.0, 1.0);
    systrafo_->complete();
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtLagrangeStrategy::evaluate_force_stiff(const Teuchos::RCP<const Epetra_Vector> dis)
{
  bool successForce = evaluate_force(dis);
  bool successStiff = evaluate_stiff(dis);

  return (successForce && successStiff);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::MtLagrangeStrategy::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  Teuchos::RCP<Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
    {
      vec_ptr = f_;
      break;
    }
    case CONTACT::VecBlockType::constraint:
    {
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown Solid::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::MtLagrangeStrategy::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt) const
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
      mat_ptr = Teuchos::null;
      break;
    case CONTACT::MatBlockType::displ_lm:
      if (dm_matrix_.is_null()) FOUR_C_THROW("matrix not available");
      mat_ptr = dm_matrix_;
      break;
    case CONTACT::MatBlockType::lm_displ:
      if (dm_matrix_t_.is_null()) FOUR_C_THROW("matrix not available");
      mat_ptr = dm_matrix_t_;
      break;
    case CONTACT::MatBlockType::lm_lm:
      if (lm_diag_matrix_.is_null()) FOUR_C_THROW("matrix not available");
      mat_ptr = lm_diag_matrix_;
      break;
    default:
      FOUR_C_THROW("Unknown Solid::MatBlockType!");
      break;
  }
  return mat_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::run_pre_apply_jacobian_inverse(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff, Epetra_Vector& rhs)
{
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");

  if (systype == Inpar::CONTACT::system_condensed)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*kteff));
    Teuchos::RCP<Epetra_Vector> r = Teuchos::rcpFromRef<Epetra_Vector>(rhs);

    Inpar::Mortar::LagMultQuad lagmultquad =
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // apply basis transformation to K and f
      k = Core::LinAlg::MLMultiply(*k, false, *systrafo_, false, false, false, true);
      k = Core::LinAlg::MLMultiply(*systrafo_, true, *k, false, false, false, true);
      systrafo_->multiply(true, *r, *r);
    }

    Mortar::UTILS::MortarMatrixCondensation(k, mhatmatrix_, mhatmatrix_);
    *kteff = *k;

    Mortar::UTILS::MortarRhsCondensation(r, mhatmatrix_);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::run_post_apply_jacobian_inverse(Epetra_Vector& result)
{
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (systype == Inpar::CONTACT::system_condensed)
  {
    Teuchos::RCP<Epetra_Vector> inc = Teuchos::rcpFromRef<Epetra_Vector>(result);
    Mortar::UTILS::MortarRecover(inc, mhatmatrix_);

    // undo basis transformation to solution
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
      systrafo_->multiply(false, *inc, *inc);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  if (system_type() != Inpar::CONTACT::system_condensed)
  {
    Teuchos::RCP<Epetra_Vector> zdir_ptr = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_, true));
    Core::LinAlg::export_to(dir, *zdir_ptr);
    zdir_ptr->ReplaceMap(*gsdofrowmap_);
    z_->Update(1., *zdir_ptr, 1.);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const
{
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (systype == Inpar::CONTACT::system_condensed)
  {
    // undo basis transformation to solution
    if (dualquadslavetrafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
      systrafo_->multiply(true, rhs, rhs);

    Teuchos::RCP<Epetra_Vector> r = Teuchos::rcpFromRef<Epetra_Vector>(rhs);
    Mortar::UTILS::MortarRhsCondensation(r, mhatmatrix_);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
