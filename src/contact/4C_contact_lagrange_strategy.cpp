/*---------------------------------------------------------------------*/
/*! \file
\brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_lagrange_strategy.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_contact_utils.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_utils_epetra_exceptions.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::LagrangeStrategy::LagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::Interface>> interface, const int spatialDim,
    Teuchos::RCP<const Epetra_Comm> comm, const double alphaf, const int maxdof)
    : AbstractStrategy(data_ptr, dof_row_map, NodeRowMap, params, spatialDim, comm, alphaf, maxdof),
      interface_(interface),
      evalForceCalled_(false),
      activesetssconv_(false),
      activesetconv_(false),
      activesetsteps_(1),
      fLTLOld_(Teuchos::null),
      fLTL_(Teuchos::null),
      fLTLn_(Teuchos::null),
      fLTLt_(Teuchos::null),
      fconservation_(Teuchos::null)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::initialize()
{
  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *gsdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *gmdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    // (re)setup global tangent matrix
    if (!friction_) tmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 3));

    // (re)setup global normal matrix
    if (regularized_) nmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 3));

    // (re)setup global matrix containing gap derivatives
    smatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 3));

    // inactive rhs for the saddle point problem
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    inactiverhs_ = Core::LinAlg::CreateVector(*gidofs, true);

    // further terms depend on friction case
    // (re)setup global matrix containing "no-friction"-derivatives
    if (!friction_)
    {
      // tangential rhs
      tangrhs_ = Core::LinAlg::CreateVector(*gactivedofs_, true);
      tderivmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 3));
      if (regularized_)
        nderivmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 3));
    }
    // (re)setup of global friction
    else
    {
      // here the calculation of gstickt is necessary
      Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);
      linstickLM_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstickdofs, 3));
      linstickDIS_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstickdofs, 3));
      linstickRHS_ = Core::LinAlg::CreateVector(*gstickdofs, true);

      linslipLM_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipdofs_, 3));
      linslipDIS_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipdofs_, 3));
      linslipRHS_ = Core::LinAlg::CreateVector(*gslipdofs_, true);
    }
  }

  else
  {
    // (re)setup global tangent matrix
    if (!friction_) tmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_, 3));

    // (re)setup global normal matrix
    if (regularized_) nmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 3));

    // (re)setup global matrix containing gap derivatives
    smatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 3));

    // inactive rhs for the saddle point problem
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    inactiverhs_ = Core::LinAlg::CreateVector(*gidofs, true);

    // further terms depend on friction case
    // (re)setup global matrix containing "no-friction"-derivatives
    if (!friction_)
    {
      // tangential rhs
      tangrhs_ = Core::LinAlg::CreateVector(*gactivet_, true);
      tderivmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_, 3));
      if (regularized_) nderivmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 3));
    }
    // (re)setup of global friction
    else
    {
      // here the calculation of gstickt is necessary
      Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
      linstickLM_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstickt, 3));
      linstickDIS_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstickt, 3));
      linstickRHS_ = Core::LinAlg::CreateVector(*gstickt, true);

      linslipLM_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipt_, 3));
      linslipDIS_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipt_, 3));
      linslipRHS_ = Core::LinAlg::CreateVector(*gslipt_, true);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate frictional contact (public)                   gitterle 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_friction(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // In case of nonsmooth contact the scenario of contacting edges (non parallel)
  // requires a penalty regularization. Here, the penalty contributions for this
  // special case are applied:
  if (nonSmoothContact_)
  {
    // add_line_to_lin_contributions(kteff,feff);
    add_line_to_lin_contributions_friction(kteff, feff);

    // FD check of weighted gap g derivatives + jump for LTL case
#ifdef CONTACTFDJUMPLTL
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckJumpDerivLTL();
    }
#endif
  }

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    gact = Core::LinAlg::CreateVector(*gactivedofs_, true);
    if (gact->GlobalLength()) Core::LinAlg::export_to(*wgap_, *gact);
  }
  else
  {
    gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  /**********************************************************************/
  /* build global matrix t with tangent vectors of active nodes         */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix linstick with derivatives of stick nodes         */
  /* and global matrix linslip with derivatives of slip nodes           */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /**********************************************************************/
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_s(*smatrix_);
    interface_[i]->assemble_lin_dm(*lindmatrix_, *linmmatrix_);
    interface_[i]->assemble_lin_stick(*linstickLM_, *linstickDIS_, *linstickRHS_);
    interface_[i]->assemble_lin_slip(*linslipLM_, *linslipDIS_, *linslipRHS_);
    if (system_type() != Inpar::CONTACT::system_condensed)
      interface_[i]->assemble_inactiverhs(*inactiverhs_);
  }
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    smatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
  }
  else
  {
    // fill_complete() global matrix S
    smatrix_->complete(*gsmdofrowmap_, *gactiven_);
  }

  // fill_complete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);

  // fill_complete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);

  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    linstickLM_->complete(*gstickdofs, *gstickdofs);
    linstickDIS_->complete(*gsmdofrowmap_, *gstickdofs);
    linslipLM_->complete(*gslipdofs_, *gslipdofs_);
    linslipDIS_->complete(*gsmdofrowmap_, *gslipdofs_);
  }
  else
  {
    linstickLM_->complete(*gstickdofs, *gstickt);
    linstickDIS_->complete(*gsmdofrowmap_, *gstickt);

    // fill_complete global Matrix linslipLM_ and linslipDIS_
    linslipLM_->complete(*gslipdofs_, *gslipt_);
    linslipDIS_->complete(*gsmdofrowmap_, *gslipt_);
  }
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify lindmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
        Core::LinAlg::MLMultiply(*invtrafo_, true, *lindmatrix_, false, false, false, true);
    lindmatrix_ = temp1;
  }

  // shape function
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype_ == Inpar::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    /********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                    */
    /********************************************************************/
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));

    // for nonsmooth contact inverting D is more complex:
    // Note: this invertation if only applicable when vertex, edge and surface nodes
    // are involved. For a falling coin (only surface and edge nodes), a special but
    // more easy implementation is needed.
    if (nonSmoothContact_)
    {
      // 1. split d matrix in vertex edge and surf part
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dss, dsev, devs, devev;
      Teuchos::RCP<Epetra_Map> gEVdofs;  // merged edge and vertex dofs

      // get dss
      Core::LinAlg::SplitMatrix2x2(
          dmatrix_, gsdofSurf_, gEVdofs, gsdofSurf_, gEVdofs, dss, dsev, devs, devev);

      // get dse and dsv
      Teuchos::RCP<Epetra_Map> temp;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dse, dsv;

      Core::LinAlg::SplitMatrix2x2(
          dsev, gsdofSurf_, temp, gsdofEdge_, gsdofVertex_, dse, dsv, tempmtx1, tempmtx2);

      // get dee dev dve dvv
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dee, dev, dve, dvv;
      Core::LinAlg::SplitMatrix2x2(
          devev, gsdofEdge_, gsdofVertex_, gsdofEdge_, gsdofVertex_, dee, dev, dve, dvv);

      // 2. invert diagonal matrices dss dee dvv
      Teuchos::RCP<Epetra_Vector> diagV = Core::LinAlg::CreateVector(*gsdofVertex_, true);
      Teuchos::RCP<Epetra_Vector> diagE = Core::LinAlg::CreateVector(*gsdofEdge_, true);
      Teuchos::RCP<Epetra_Vector> diagS = Core::LinAlg::CreateVector(*gsdofSurf_, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdV =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dvv));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdE =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dee));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdS =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dss));

      int err = 0;

      // extract diagonal of invd into diag
      invdV->extract_diagonal_copy(*diagV);
      invdE->extract_diagonal_copy(*diagE);
      invdS->extract_diagonal_copy(*diagS);

      // set zero diagonal values to dummy 1.0
      for (int i = 0; i < diagV->MyLength(); ++i)
        if (abs((*diagV)[i]) < 1e-12) (*diagV)[i] = 1.0;
      for (int i = 0; i < diagE->MyLength(); ++i)
        if (abs((*diagE)[i]) < 1e-12) (*diagE)[i] = 1.0;
      for (int i = 0; i < diagS->MyLength(); ++i)
        if (abs((*diagS)[i]) < 1e-12) (*diagS)[i] = 1.0;

      // scalar inversion of diagonal values
      err = diagV->Reciprocal(*diagV);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagE->Reciprocal(*diagE);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagS->Reciprocal(*diagS);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

      // re-insert inverted diagonal into invd
      err = invdV->replace_diagonal_values(*diagV);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdE->replace_diagonal_values(*diagE);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdS->replace_diagonal_values(*diagS);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

      // 3. multiply all sub matrices
      invd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, true, true));

      dse->scale(-1.0);
      dsv->scale(-1.0);
      dev->scale(-1.0);

      // inv_dse
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum1;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dse;
      dum1 = Core::LinAlg::MLMultiply(*dse, false, *invdE, false, false, false, true);
      dinv_dse = Core::LinAlg::MLMultiply(*invdS, false, *dum1, false, false, false, true);

      // inv_dev
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dev;
      dum2 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dinv_dev = Core::LinAlg::MLMultiply(*invdE, false, *dum2, false, false, false, true);

      // inv_dsv part1
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum3;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum4;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum5;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv1;
      dum3 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dum4 = Core::LinAlg::MLMultiply(*invdE, false, *dum3, false, false, false, true);
      dum5 = Core::LinAlg::MLMultiply(*dse, false, *dum4, false, false, false, true);
      dinv_dsv1 = Core::LinAlg::MLMultiply(*invdS, false, *dum5, false, false, false, true);

      // inv_dsv part2
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum6;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv2;
      dum6 = Core::LinAlg::MLMultiply(*dsv, false, *invdV, false, false, false, true);
      dinv_dsv2 = Core::LinAlg::MLMultiply(*invdS, false, *dum6, false, false, false, true);

      // diagonal entries
      invd->add(*invdS, false, 1.0, 1.0);
      invd->add(*invdE, false, 1.0, 1.0);
      invd->add(*invdV, false, 1.0, 1.0);

      invd->add(*dinv_dev, false, 1.0, 1.0);
      invd->add(*dinv_dse, false, 1.0, 1.0);
      invd->add(*dinv_dsv1, false, 1.0, 1.0);
      invd->add(*dinv_dsv2, false, 1.0, 1.0);

      // complete
      invd->complete();
    }
    // standard inverse diagonal matrix:
    else
    {
      Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
      int err = 0;

      // extract diagonal of invd into diag
      invd->extract_diagonal_copy(*diag);

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
      err = invd->replace_diagonal_values(*diag);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
    }

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);

    /********************************************************************/
    /* (2) Add contact stiffness terms to kteff                         */
    /********************************************************************/

    // transform if necessary
    if (parallel_redistribution_status())
    {
      lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
      linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
    }

    kteff->un_complete();
    kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->complete();

    /********************************************************************/
    /* (3) Split kteff into 3x3 matrix blocks                           */
    /********************************************************************/

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
    if (parallel_redistribution_status())
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

    /********************************************************************/
    /* (4) Split feff into 3 subvectors                                 */
    /********************************************************************/

    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fs, fm, fn;

    // temporarily we need the group sm
    Teuchos::RCP<Epetra_Vector> fsm;

    // do the vector splitting smn -> sm+n
    if (parallel_redistribution_status())
    {
      // split and transform to redistributed maps
      Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
      Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
      Core::LinAlg::export_to(*fsm, *fsmtemp);
      fsm = fsmtemp;
    }
    else
    {
      // only split, no need to transform
      Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
    }

    // abbreviations for slave and master set
    const int sset = gsdofrowmap_->NumGlobalElements();
    const int mset = gmdofrowmap_->NumGlobalElements();

    // we want to split fsm into 2 groups s,m
    fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

    // do the vector splitting sm -> s+m
    Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

    // store some stuff for static condensation of LM
    fs_ = fs;
    invd_ = invd;
    ksn_ = ksn;
    ksm_ = ksm;
    kss_ = kss;

    //--------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //--------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    // D^(-1)    ---->   T * D^(-1)
    // \hat{M}   ---->   T * \hat{M}
    //--------------------------------------------------------------------
    if (is_dual_quad_slave_trafo())
    {
      // modify dmatrix_, invd_ and mhatmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
          Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
          Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
          Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
      dmatrix_ = temp2;
      invd_ = temp3;
      mhatmatrix_ = temp4;
    }

    /********************************************************************/
    /* (5) Split slave quantities into active / inactive, stick / slip  */
    /********************************************************************/
    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;

    // do the splitting
    Core::LinAlg::SplitMatrix2x2(
        kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
    Core::LinAlg::SplitMatrix2x2(
        ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
    Core::LinAlg::SplitMatrix2x2(
        ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
    Core::LinAlg::SplitMatrix2x2(
        kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

    // we want to split kaa into 2 groups sl,st = 4 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

    // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> temp1map;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1mtx4, temp1mtx5;

    // we will get the stick rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gstdofs;

    Core::LinAlg::SplitMatrix2x2(
        kaa, gactivedofs_, gidofs, gstdofs, gslipdofs_, kast, kasl, temp1mtx4, temp1mtx5);

    // abbreviations for active and inactive set, stick and slip set
    const int aset = gactivedofs_->NumGlobalElements();
    const int iset = gidofs->NumGlobalElements();
    const int stickset = gstdofs->NumGlobalElements();
    const int slipset = gslipdofs_->NumGlobalElements();

    // we want to split fs into 2 groups a,i
    Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // do the vector splitting s -> a+i
    Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

    // we want to split fa into 2 groups sl,st
    Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
    Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

    // do the vector splitting a -> sl+st
    if (aset) Core::LinAlg::split_vector(*gactivedofs_, *fa, gslipdofs_, fsl, gstdofs, fst);

    /********************************************************************/
    /* (6) Isolate necessary parts from invd and mhatmatrix             */
    /********************************************************************/

    // active, stick and slip part of invd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda, invdsl, invdst;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);
    Core::LinAlg::SplitMatrix2x2(
        invda, gactivedofs_, gidofs, gslipdofs_, gstdofs, invdsl, tempmtx1, tempmtx2, tempmtx3);
    Core::LinAlg::SplitMatrix2x2(
        invda, gactivedofs_, gidofs, gstdofs, gslipdofs_, invdst, tempmtx1, tempmtx2, tempmtx3);

    // coupling part of dmatrix (only nonzero for 3D quadratic case!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
    Core::LinAlg::SplitMatrix2x2(
        dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

    // do the multiplication dhat = invda * dai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
    if (aset && iset)
      dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
    dhat->complete(*gidofs, *gactivedofs_);

    // active part of mmatrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
    Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
        tempmtx1, tempmtx2, tempmtx3);

    // do the multiplication mhataam = invda * mmatrixa
    // (this is only different from mhata for 3D quadratic case!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
    if (aset)
      mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
    mhataam->complete(*gmdofrowmap_, *gactivedofs_);

    // for the case without full linearization, we still need the
    // "classical" active part of mhat, which is isolated here
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
    Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
        tempmtx1, tempmtx2, tempmtx3);

    // scaling of invd and dai
    invda->scale(1 / (1 - alphaf_));
    invdsl->scale(1 / (1 - alphaf_));
    invdst->scale(1 / (1 - alphaf_));
    dai->scale(1 - alphaf_);

    /********************************************************************/
    /* (7) Build the final K blocks                                     */
    /********************************************************************/

    //--------------------------------------------------------- FIRST LINE
    // knn: nothing to do

    // knm: nothing to do

    // kns: nothing to do

    //-------------------------------------------------------- SECOND LINE
    // kmn: add T(mhataam)*kan
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmnmod->add(*kmn, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
    kmnmod->add(*kmnadd, false, 1.0, 1.0);
    kmnmod->complete(kmn->domain_map(), kmn->row_map());

    // kmm: add T(mhataam)*kam
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmmmod->add(*kmm, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
    kmmmod->add(*kmmadd, false, 1.0, 1.0);
    kmmmod->complete(kmm->domain_map(), kmm->row_map());

    // kmi: add T(mhataam)*kai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
      kmimod->add(*kmi, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
          Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
      kmimod->add(*kmiadd, false, 1.0, 1.0);
      kmimod->complete(kmi->domain_map(), kmi->row_map());
    }

    // kma: add T(mhataam)*kaa
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
      kmamod->add(*kma, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
          Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
      kmamod->add(*kmaadd, false, 1.0, 1.0);
      kmamod->complete(kma->domain_map(), kma->row_map());
    }

    //--------------------------------------------------------- THIRD LINE
    // kin: subtract T(dhat)*kan
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kinmod->add(*kin, false, 1.0, 1.0);
    if (aset && iset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
      kinmod->add(*kinadd, false, -1.0, 1.0);
    }
    kinmod->complete(kin->domain_map(), kin->row_map());

    // kim: subtract T(dhat)*kam
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kimmod->add(*kim, false, 1.0, 1.0);
    if (aset && iset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
      kimmod->add(*kimadd, false, -1.0, 1.0);
    }
    kimmod->complete(kim->domain_map(), kim->row_map());

    // kii: subtract T(dhat)*kai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
    if (iset)
    {
      kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
      kiimod->add(*kii, false, 1.0, 1.0);
      if (aset)
      {
        Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
            Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
        kiimod->add(*kiiadd, false, -1.0, 1.0);
      }
      kiimod->complete(kii->domain_map(), kii->row_map());
    }

    // kia: subtract T(dhat)*kaa
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
    if (iset && aset)
    {
      kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
      kiamod->add(*kia, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
      kiamod->add(*kiaadd, false, -1.0, 1.0);
      kiamod->complete(kia->domain_map(), kia->row_map());
    }

    //-------------------------------------------------------- FOURTH LINE

    //--------------------------------------------------------- FIFTH LINE
    // blocks for complementary conditions (stick nodes)

    // kstn: multiply with linstickLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kstnmod;
    if (stickset)
    {
      kstnmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      kstnmod = Core::LinAlg::MLMultiply(*kstnmod, false, *kan, false, false, false, true);
    }

    // kstm: multiply with linstickLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kstmmod;
    if (stickset)
    {
      kstmmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      kstmmod = Core::LinAlg::MLMultiply(*kstmmod, false, *kam, false, false, false, true);
    }

    // ksti: multiply with linstickLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kstimod;
    if (stickset && iset)
    {
      kstimod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      kstimod = Core::LinAlg::MLMultiply(*kstimod, false, *kai, false, false, false, true);
    }

    // kstsl: multiply with linstickLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kstslmod;
    if (stickset && slipset)
    {
      kstslmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      kstslmod = Core::LinAlg::MLMultiply(*kstslmod, false, *kasl, false, false, false, true);
    }

    // kststmod: multiply with linstickLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kststmod;
    if (stickset)
    {
      kststmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      kststmod = Core::LinAlg::MLMultiply(*kststmod, false, *kast, false, false, false, true);
    }

    //--------------------------------------------------------- SIXTH LINE
    // blocks for complementary conditions (slip nodes)

    // ksln: multiply with linslipLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslnmod;
    if (slipset)
    {
      kslnmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      kslnmod = Core::LinAlg::MLMultiply(*kslnmod, false, *kan, false, false, false, true);
    }

    // kslm: multiply with linslipLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslmmod;
    if (slipset)
    {
      kslmmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      kslmmod = Core::LinAlg::MLMultiply(*kslmmod, false, *kam, false, false, false, true);
    }

    // ksli: multiply with linslipLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslimod;
    if (slipset && iset)
    {
      kslimod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      kslimod = Core::LinAlg::MLMultiply(*kslimod, false, *kai, false, false, false, true);
    }

    // kslsl: multiply with linslipLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslslmod;
    if (slipset)
    {
      kslslmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      kslslmod = Core::LinAlg::MLMultiply(*kslslmod, false, *kasl, false, false, false, true);
    }

    // slstmod: multiply with linslipLM
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kslstmod;
    if (slipset && stickset)
    {
      kslstmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      kslstmod = Core::LinAlg::MLMultiply(*kslstmod, false, *kast, false, false, false, true);
    }

    /********************************************************************/
    /* (8) Build the final f blocks                                     */
    /********************************************************************/

    //--------------------------------------------------------- FIRST LINE
    // fn: nothing to do

    //---------------------------------------------------------- SECOND LINE
    // fm: add alphaf * old contact forces (t_n)
    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Mold^T * zold to fit
    if (is_self_contact())
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> tempvecm2 = Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(mold_->row_map()));
      if (mold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      mold_->multiply(true, *zoldexp, *tempvecm2);
      if (mset) Core::LinAlg::export_to(*tempvecm2, *tempvecm);
      fm->Update(alphaf_, *tempvecm, 1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->multiply(true, *zold_, *tempvecm);
      fm->Update(alphaf_, *tempvecm, 1.0);
    }

    // fs: prepare alphaf * old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Dold^T * zold to fit
    if (is_self_contact())
    {
      Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      dold_->multiply(true, *zoldexp, *tempvec);
      if (sset) Core::LinAlg::export_to(*tempvec, *fsadd);
    }
    // if there is no self contact everything is ok
    else
    {
      dold_->multiply(true, *zold_, *fsadd);
    }

    // fa: subtract alphaf * old contact forces (t_n)
    if (aset)
    {
      Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      Core::LinAlg::export_to(*fsadd, *faadd);
      fa->Update(-alphaf_, *faadd, 1.0);
    }

    // fm: add T(mhat)*fa
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (aset) mhataam->multiply(true, *fa, *fmmod);
    fmmod->Update(1.0, *fm, 1.0);

    //--------------------------------------------------------- THIRD LINE
    // fi: subtract alphaf * old contact forces (t_n)
    if (iset)
    {
      Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
      Core::LinAlg::export_to(*fsadd, *fiadd);
      fi->Update(-alphaf_, *fiadd, 1.0);
    }

    // fi: add T(dhat)*fa
    Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
    if (aset && iset) dhat->multiply(true, *fa, *fimod);
    fimod->Update(1.0, *fi, -1.0);

    //-------------------------------------------------------- FOURTH LINE

    //--------------------------------------------------------- FIFTH LINE
    Teuchos::RCP<Epetra_Map> gstickdofs =
        Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);  // get global stick dofs

    // split the lagrange multiplier vector in stick and slip part
    Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
    Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
    Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

    Core::LinAlg::split_vector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
    Core::LinAlg::split_vector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
    Teuchos::RCP<Epetra_Vector> tempvec1;

    // fst: mutliply with linstickLM
    Teuchos::RCP<Epetra_Vector> fstmod;
    if (stickset)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        fstmod = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
      else
        fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
      temp1->multiply(false, *fa, *fstmod);

      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
      else
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

      linstickLM_->multiply(false, *zst, *tempvec1);
      fstmod->Update(-1.0, *tempvec1, 1.0);
    }

    //--------------------------------------------------------- SIXTH LINE
    // fsl: mutliply with linslipLM
    Teuchos::RCP<Epetra_Vector> fslmod;
    Teuchos::RCP<Epetra_Vector> fslwmod;

    if (slipset)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        fslmod = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
      else
        fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
          Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
      temp->multiply(false, *fa, *fslmod);

      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
      else
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

      linslipLM_->multiply(false, *zsl, *tempvec1);

      fslmod->Update(-1.0, *tempvec1, 1.0);
    }

    /********************************************************************/
    /* (9) Transform the final K blocks                                 */
    /********************************************************************/
    // The row maps of all individual matrix blocks are transformed to
    // the parallel layout of the underlying problem discretization.
    // Of course, this is only necessary in the parallel redistribution
    // case, where the contact interfaces have been redistributed
    // independently of the underlying problem discretization.

    if (parallel_redistribution_status())
    {
      //----------------------------------------------------------- FIRST LINE
      // nothing to do (ndof-map independent of redistribution)

      //---------------------------------------------------------- SECOND LINE
      kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
      kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
      if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
      if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

      //----------------------------------------------------------- THIRD LINE
      if (iset)
      {
        kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
        kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
        kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
        if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
      }

      //---------------------------------------------------------- FOURTH LINE
      if (aset)
      {
        smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);
      }

      //----------------------------------------------------------- FIFTH LINE
      if (stickset)
      {
        kstnmod = Mortar::MatrixRowTransform(kstnmod, pgsdofrowmap_);
        kstmmod = Mortar::MatrixRowTransform(kstmmod, pgsdofrowmap_);
        if (iset) kstimod = Mortar::MatrixRowTransform(kstimod, pgsdofrowmap_);
        if (slipset) kstslmod = Mortar::MatrixRowTransform(kstslmod, pgsdofrowmap_);
        kststmod = Mortar::MatrixRowTransform(kststmod, pgsdofrowmap_);
        linstickDIS_ = Mortar::MatrixRowTransform(linstickDIS_, pgsdofrowmap_);
      }

      //----------------------------------------------------------- SIXTH LINE
      if (slipset)
      {
        kslnmod = Mortar::MatrixRowTransform(kslnmod, pgsdofrowmap_);
        kslmmod = Mortar::MatrixRowTransform(kslmmod, pgsdofrowmap_);
        if (iset) kslimod = Mortar::MatrixRowTransform(kslimod, pgsdofrowmap_);
        if (stickset) kslstmod = Mortar::MatrixRowTransform(kslstmod, pgsdofrowmap_);
        kslslmod = Mortar::MatrixRowTransform(kslslmod, pgsdofrowmap_);
        linslipDIS_ = Mortar::MatrixRowTransform(linslipDIS_, pgsdofrowmap_);
      }
    }

    /********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                */
    /********************************************************************/

    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
    Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

    //--------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->add(*knn, false, 1.0, 1.0);
    kteffnew->add(*knm, false, 1.0, 1.0);
    if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

    //-------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->add(*kmnmod, false, 1.0, 1.0);
    kteffnew->add(*kmmmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
    if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

    //--------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
    if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

    //-------------------------------------------------------- FOURTH LINE

    // add a submatrices to kteffnew
    if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);

    //--------------------------------------------------------- FIFTH LINE
    // add st submatrices to kteffnew
    if (stickset) kteffnew->add(*kstnmod, false, 1.0, 1.0);
    if (stickset) kteffnew->add(*kstmmod, false, 1.0, 1.0);
    if (stickset && iset) kteffnew->add(*kstimod, false, 1.0, 1.0);
    if (stickset && slipset) kteffnew->add(*kstslmod, false, 1.0, 1.0);
    if (stickset) kteffnew->add(*kststmod, false, 1.0, 1.0);

    // add terms of linearization of sick condition to kteffnew
    if (stickset) kteffnew->add(*linstickDIS_, false, -1.0, 1.0);

    //--------------------------------------------------------- SIXTH LINE
    // add sl submatrices to kteffnew
    if (slipset) kteffnew->add(*kslnmod, false, 1.0, 1.0);
    if (slipset) kteffnew->add(*kslmmod, false, 1.0, 1.0);
    if (slipset && iset) kteffnew->add(*kslimod, false, 1.0, 1.0);
    if (slipset) kteffnew->add(*kslslmod, false, 1.0, 1.0);
    if (slipset && stickset) kteffnew->add(*kslstmod, false, 1.0, 1.0);

    // add terms of linearization of slip condition to kteffnew and feffnew
    if (slipset) kteffnew->add(*linslipDIS_, false, -1.0, +1.0);

    // fill_complete kteffnew (square)
    kteffnew->complete();

    /********************************************************************/
    /* (11) Global setup of feffnew (including contact)                 */
    /********************************************************************/

    //--------------------------------------------------------- FIRST LINE
    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fn, *fnexp);
    feffnew->Update(1.0, *fnexp, 1.0);

    //-------------------------------------------------------- SECOND LINE
    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fmmod, *fmmodexp);
    feffnew->Update(1.0, *fmmodexp, 1.0);

    //--------------------------------------------------------- THIRD LINE
    // add i subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fimodexp;
    if (iset)
    {
      fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fimod, *fimodexp);
      feffnew->Update(1.0, *fimodexp, 1.0);
    }

    //-------------------------------------------------------- FOURTH LINE
    // add weighted gap vector to feffnew, if existing
    Teuchos::RCP<Epetra_Vector> gexp;
    Teuchos::RCP<Epetra_Vector> fwexp;
    Teuchos::RCP<Epetra_Vector> fgmodexp;

    if (aset)
    {
      gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*gact, *gexp);
      feffnew->Update(-1.0, *gexp, 1.0);
    }

    //--------------------------------------------------------- FIFTH LINE
    // add st subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fstmodexp;
    if (stickset)
    {
      fstmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fstmod, *fstmodexp);
      feffnew->Update(1.0, *fstmodexp, +1.0);
    }

    // add terms of linearization feffnew
    if (stickset)
    {
      Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*linstickRHS_, *linstickRHSexp);
      feffnew->Update(-1.0, *linstickRHSexp, 1.0);
    }

    //--------------------------------------------------------- SIXTH LINE

    // add a subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fslmodexp;
    Teuchos::RCP<Epetra_Vector> fwslexp;
    Teuchos::RCP<Epetra_Vector> fslwmodexp;


    if (slipset)
    {
      fslmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fslmod, *fslmodexp);
      feffnew->Update(1.0, *fslmodexp, 1.0);
    }

    if (slipset)
    {
      Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*linslipRHS_, *linslipRHSexp);
      feffnew->Update(-1.0, *linslipRHSexp, 1.0);
    }

    // finally do the replacement
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
    //----------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //----------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    //----------------------------------------------------------------------
    if (is_dual_quad_slave_trafo())
    {
      // modify dmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
          Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
      dmatrix_ = temp2;
    }

    // transform if necessary
    if (parallel_redistribution_status())
    {
      lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
      linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
    }

    // add contact stiffness
    kteff->un_complete();
    kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->complete();

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold^T * zold / D^T * z to fit
    // thus we have to export the products Mold^T * zold / M^T * z to fit
    if (is_self_contact())
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecd =
          Teuchos::rcp(new Epetra_Vector(dmatrix_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
      if (dmatrix_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*z_, *zexp);
      dmatrix_->multiply(true, *zexp, *tempvecd);
      Core::LinAlg::export_to(*tempvecd, *fsexp);
      feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecm =
          Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));
      mmatrix_->multiply(true, *zexp, *tempvecm);
      Core::LinAlg::export_to(*tempvecm, *fmexp);
      feff->Update(1.0 - alphaf_, *fmexp, 1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecdold =
          Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      dold_->multiply(true, *zoldexp, *tempvecdold);
      Core::LinAlg::export_to(*tempvecdold, *fsoldexp);
      feff->Update(-alphaf_, *fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecmold =
          Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
      mold_->multiply(true, *zoldexp, *tempvecmold);
      Core::LinAlg::export_to(*tempvecmold, *fmoldexp);
      feff->Update(alphaf_, *fmoldexp, 1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      // add contact force terms
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
      dold_->multiply(true, *zold_, *fsold);
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fsold, *fsoldexp);
      feff->Update(-alphaf_, *fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->multiply(true, *zold_, *fmold);
      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fmold, *fmoldexp);
      feff->Update(alphaf_, *fmoldexp, 1.0);

      // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
      check_conservation_laws(*fs, *fm);
#endif
    }
  }

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSLIPINCR
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->fd_check_slip_incr_deriv_txi();
    if (Dim() == 3) interface_[i]->fd_check_slip_incr_deriv_teta();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
    // FD check of stick condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      //      Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv1 = Teuchos::rcp(new
      //      Core::LinAlg::SparseMatrix(*gactivet_,81)); Teuchos::RCP<Core::LinAlg::SparseMatrix>
      //      deriv2 = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_,81));
      //
      //      deriv1->Add(*linstickLM_,false,1.0,1.0);
      //      deriv1->Complete(*gsmdofrowmap_,*gactivet_);
      //
      //      deriv2->Add(*linstickDIS_,false,1.0,1.0);
      //      deriv2->Complete(*gsmdofrowmap_,*gactivet_);
      //
      //      std::cout << "DERIV 1 *********** "<< *deriv1 << std::endl;
      //      std::cout << "DERIV 2 *********** "<< *deriv2 << std::endl;

      interface_[i]->FDCheckStickDeriv(*linstickLM_, *linstickDIS_);
    }
  }
#endif  // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
    // FD check of slip condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      //      Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv1 = Teuchos::rcp(new
      //      Core::LinAlg::SparseMatrix(*gactivet_,81)); Teuchos::RCP<Core::LinAlg::SparseMatrix>
      //      deriv2 = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_,81));
      //
      //      deriv1->Add(*linslipLM_,false,1.0,1.0);
      //      deriv1->Complete(*gsmdofrowmap_,*gslipt_);
      //
      //      deriv2->Add(*linslipDIS_,false,1.0,1.0);
      //      deriv2->Complete(*gsmdofrowmap_,*gslipt_);
      //
      //      std::cout << *deriv1 << std::endl;
      //      std::cout << *deriv2 << std::endl;

      interface_[i]->FDCheckSlipDeriv(*linslipLM_, *linslipDIS_);
    }
  }
#endif  // #ifdef CONTACTFDSLIP

  return;
}

/*----------------------------------------------------------------------*
 |  pp stresses                                              farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::compute_contact_stresses()
{
  static int step = 0;
  // call abstract function
  CONTACT::AbstractStrategy::compute_contact_stresses();

  // further scaling for nonsmooth contact
  if (nonSmoothContact_)
  {
    forcenormal_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    d_matrix()->multiply(true, *stressnormal_, *forcenormal_);
    forcetangential_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    d_matrix()->multiply(true, *stresstangential_, *forcetangential_);

    Teuchos::RCP<Epetra_Vector> forcenormal =
        Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    d_matrix()->multiply(true, *stressnormal_, *forcenormal);

    Teuchos::RCP<Epetra_Vector> forcetangential =
        Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    d_matrix()->multiply(true, *stresstangential_, *forcetangential);

    // add penalty force normal
    if (fLTLn_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> dummy = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*fLTLn_, *dummy);
      forcenormal_->Update(1.0, *dummy, 1.0);
      forcenormal->Update(1.0, *dummy, 1.0);
    }

    // add penalty force tangential
    if (fLTLt_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> dummy = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*fLTLt_, *dummy);
      forcetangential_->Update(1.0, *dummy, 1.0);
      forcetangential->Update(1.0, *dummy, 1.0);
    }

    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      // loop over all slave row nodes on the current interface
      for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_row_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cnode = dynamic_cast<Node*>(node);

        std::vector<int> locindex(n_dim());

        for (int dof = 0; dof < n_dim(); ++dof)
        {
          locindex[dof] = (forcenormal->Map()).LID(cnode->dofs()[dof]);

          if (cnode->mo_data().get_dscale() < 1e-8 and cnode->active())
          {
            std::cout << "WARNING: dscale not valid!" << std::endl;
            continue;
          }
          else if (cnode->mo_data().get_dscale() < 1e-8)
          {
            continue;
          }

          (*forcenormal)[locindex[dof]] /= cnode->mo_data().get_dscale();
          (*forcetangential)[locindex[dof]] /= cnode->mo_data().get_dscale();
        }
      }
    }
    stresstangential_->Update(1.0, *forcetangential, 0.0);
    stressnormal_->Update(1.0, *forcenormal, 0.0);

    // temporary output:
    double tangforce = 0.0;
    forcetangential_->Norm2(&tangforce);
    if (get_comm().MyPID() == 0) std::cout << "tangential force = " << tangforce << std::endl;

    double normalforce = 0.0;
    forcenormal_->Norm2(&normalforce);
    if (get_comm().MyPID() == 0) std::cout << "normal force = " << normalforce << std::endl;

    if (get_comm().MyPID() == 0)
    {
      FILE* MyFile = nullptr;
      std::ostringstream filename;
      const std::string filebase = "xxx";
      filename << filebase << ".fric";
      MyFile = fopen(filename.str().c_str(), "at+");

      // store data
      if (MyFile)
      {
        fprintf(MyFile, "%d\t", step);
        fprintf(MyFile, "%g\t", tangforce);
        fprintf(MyFile, "%g\n", normalforce);
        fclose(MyFile);
      }
      else
        FOUR_C_THROW("File could not be opened.");
    }
  }

  step++;

  return;
}

/*----------------------------------------------------------------------*
 |  add penalty terms for ltl contact                        farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::save_reference_state(Teuchos::RCP<const Epetra_Vector> dis)
{
  if (!nonSmoothContact_) return;

  // initialize the displacement field
  set_state(Mortar::state_new_displacement, *dis);

  // guarantee uniquness
  std::set<std::pair<int, int>> donebefore;

  // kappa will be the shape function integral on the slave sides
  // (1) build the nodal information
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // interface needs to be complete
    if (!interface_[i]->filled() && get_comm().MyPID() == 0)
      FOUR_C_THROW("fill_complete() not called on interface %", i);

    // reset kappa
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->master_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->master_row_nodes()->GID(j);
      Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);
      cnode->data().kappa() = 0.0;
    }

    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int j = 0; j < interface_[i]->master_col_elements()->NumMyElements(); ++j)
    {
      int gid1 = interface_[i]->master_col_elements()->GID(j);
      Core::Elements::Element* ele1 = interface_[i]->discret().g_element(gid1);
      if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
      Element* selement = dynamic_cast<Element*>(ele1);

      // loop over slave edges -> match node number for tri3/quad4
      for (int k = 0; k < selement->num_node(); ++k)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (selement->shape() == Core::FE::CellType::quad4)
        {
          if (k == 0)
          {
            nodeIds[0] = selement->node_ids()[0];
            nodeIds[1] = selement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (k == 1)
          {
            nodeIds[0] = selement->node_ids()[1];
            nodeIds[1] = selement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (k == 2)
          {
            nodeIds[0] = selement->node_ids()[2];
            nodeIds[1] = selement->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (k == 3)
          {
            nodeIds[0] = selement->node_ids()[3];
            nodeIds[1] = selement->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            FOUR_C_THROW("loop counter and edge number do not match!");
        }
        else if (selement->shape() == Core::FE::CellType::tri3)
        {
          if (k == 0)
          {
            nodeIds[0] = selement->node_ids()[0];
            nodeIds[1] = selement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (k == 1)
          {
            nodeIds[0] = selement->node_ids()[1];
            nodeIds[1] = selement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (k == 2)
          {
            nodeIds[0] = selement->node_ids()[2];
            nodeIds[1] = selement->node_ids()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }
          else
            FOUR_C_THROW("loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[0]])->is_on_edge();
        bool node1Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[1]])->is_on_edge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, selement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<Core::Nodes::Node*, 2> nodes = {
              selement->nodes()[nodeLIds[0]], selement->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes.data());

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // create integrator
          CONTACT::Integrator integrator(params(), lineEle->shape(), get_comm());

          // integrate kappe penalty
          integrator.integrate_kappa_penalty_lts(*lineEle);
        }
      }  // end edge loop
    }

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->master_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->master_row_nodes()->GID(j);
      Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // only for edge nodes!
      if (cnode->is_on_edge() and !cnode->is_on_corner())
      {
        // get nodal weighted gap
        // (this is where we stored the shape function integrals)
        double kappainv = cnode->data().kappa();

        // safety
        if (abs(kappainv) < 1e-12)
        {
          kappainv = 1.0;
          //          FOUR_C_THROW("gap is zero!");
        }
        else
        {
          // store kappa
          cnode->data().kappa() = 1.0 / kappainv;
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  add penalty terms for ltl contact                        farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::add_master_contributions(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    bool add_time_integration)
{
  // create new contact force vector for LTL contact
  Teuchos::RCP<Epetra_FEVector> fc = Teuchos::rcp(new Epetra_FEVector(feff->Map()), true);

  // create new contact stiffness matric for LTL contact
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kc = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      (dynamic_cast<Epetra_CrsMatrix*>(&(*kteff->epetra_operator())))->RowMap(), 100, true, false,
      Core::LinAlg::SparseMatrix::FE_MATRIX));

  // loop over interface and assemble force and stiffness
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // line to segment
    interface_[i]->add_lts_forces_master(fc);
    interface_[i]->add_lts_stiffness_master(kc);
    // node to segment
    interface_[i]->add_nts_forces_master(fc);
    interface_[i]->add_nts_stiffness_master(kc);
  }

  // force
  if (fc->GlobalAssemble(Add, false) != 0) FOUR_C_THROW("GlobalAssemble failed");

  // store fLTL values for time integration
  fLTL_ = Teuchos::rcp(new Epetra_Vector(fc->Map()), true);
  if (fLTL_->Update(1.0, *fc, 0.0)) FOUR_C_THROW("Update went wrong");

  if (add_time_integration)
    if (fLTLOld_ != Teuchos::null)
      if (feff->Update(alphaf_, *fLTLOld_, 1.)) FOUR_C_THROW("Update went wrong");

  double fac = 0.;
  if (add_time_integration)
    fac = 1. - alphaf_;
  else
    fac = 1.;
  if (feff->Update(fac, *fc, 1.)) FOUR_C_THROW("Update went wrong");

  // stiffness
  dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add);
  kteff->un_complete();
  kteff->add(*kc, false, fac, 1.);
  kteff->complete();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  add penalty terms for ltl contact                        farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::add_line_to_lin_contributions(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    bool add_time_integration)
{
  // create new contact force vector for LTL contact
  Teuchos::RCP<Epetra_FEVector> fc = Teuchos::rcp(new Epetra_FEVector(feff->Map()), true);

  fconservation_ = Teuchos::rcp(new Epetra_Vector(feff->Map()), true);

  // create new contact stiffness matric for LTL contact
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kc = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      (dynamic_cast<Epetra_CrsMatrix*>(&(*kteff->epetra_operator())))->RowMap(), 100, true, false,
      Core::LinAlg::SparseMatrix::FE_MATRIX));

  // loop over interface and assemble force and stiffness
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->add_ltl_forces(fc);
    interface_[i]->add_ltl_stiffness(kc);
  }

  // get info for conservation check
  fconservation_->Update(1.0, *fc, 0.0);

  // force
  if (fc->GlobalAssemble(Add, false) != 0) FOUR_C_THROW("GlobalAssemble failed");

  // store fLTL values for time integration
  fLTL_ = Teuchos::rcp(new Epetra_Vector(fc->Map()), true);
  if (fLTL_->Update(1.0, *fc, 0.0)) FOUR_C_THROW("Update went wrong");

  if (add_time_integration)
    if (fLTLOld_ != Teuchos::null)
      if (feff->Update(alphaf_, *fLTLOld_, 1.)) FOUR_C_THROW("Update went wrong");

  double fac = 0.;
  if (add_time_integration)
    fac = 1. - alphaf_;
  else
    fac = 1.;
  if (feff->Update(fac, *fc, 1.)) FOUR_C_THROW("Update went wrong");

  // stiffness
  dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add);
  kteff->un_complete();
  kteff->add(*kc, false, fac, 1.);
  kteff->complete();

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  add frictional penalty terms for ltl contact             farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::add_line_to_lin_contributions_friction(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    bool add_time_integration)
{
  // create new contact force vector for LTL contact
  Teuchos::RCP<Epetra_FEVector> fc = Teuchos::rcp(new Epetra_FEVector(feff->Map()), true);

  fconservation_ = Teuchos::rcp(new Epetra_Vector(feff->Map()), true);

  // create new contact stiffness matric for LTL contact
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kc = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      (dynamic_cast<Epetra_CrsMatrix*>(&(*kteff->epetra_operator())))->RowMap(), 100, true, false,
      Core::LinAlg::SparseMatrix::FE_MATRIX));

  // loop over interface and assemble force and stiffness
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->add_ltl_forces(fc);
    interface_[i]->add_ltl_stiffness(kc);
  }

  // store normal forces
  fLTLn_ = Teuchos::rcp(new Epetra_Vector(fc->Map()), true);
  if (fLTLn_->Update(1.0, *fc, 0.0)) FOUR_C_THROW("Update went wrong");

  // loop over interface and assemble force and stiffness
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->add_ltl_forces_friction(fc);
    interface_[i]->add_ltl_stiffness_friction(kc);
  }

  // get info for conservation check
  fconservation_->Update(1.0, *fc, 0.0);

  // store tangential forces
  fLTLt_ = Teuchos::rcp(new Epetra_Vector(fc->Map()), true);
  if (fLTLt_->Update(1.0, *fc, 0.0)) FOUR_C_THROW("Update went wrong");
  if (fLTLt_->Update(-1.0, *fLTLn_, 1.0)) FOUR_C_THROW("Update went wrong");

  // force
  if (fc->GlobalAssemble(Add, false) != 0) FOUR_C_THROW("GlobalAssemble failed");

  // store fLTL values for time integration
  fLTL_ = Teuchos::rcp(new Epetra_Vector(fc->Map()), true);
  if (fLTL_->Update(1.0, *fc, 0.0)) FOUR_C_THROW("Update went wrong");

  if (add_time_integration)
    if (fLTLOld_ != Teuchos::null)
      if (feff->Update(alphaf_, *fLTLOld_, 1.)) FOUR_C_THROW("Update went wrong");

  double fac = 0.;
  if (add_time_integration)
    fac = 1. - alphaf_;
  else
    fac = 1.;
  if (feff->Update(fac, *fc, 1.)) FOUR_C_THROW("Update went wrong");

  // stiffness
  dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add);
  kteff->un_complete();
  kteff->add(*kc, false, fac, 1.);
  kteff->complete();

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_contact(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // shape function type and type of LM interpolation for quadratic elements
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  // In case of nonsmooth contact the scenario of contacting edges (non parallel)
  // requires a penalty regularization. Here, the penalty contriutions for this
  // special case are applied:
  if (nonSmoothContact_)
  {
    // LTL contributions:
    add_line_to_lin_contributions(kteff, feff);

    // penalty support for master side quantities:
    add_master_contributions(kteff, feff);

#ifdef CONTACTFDGAPLTL
    // FD check of weighted gap g derivatives (non-penetr. condition)
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckGapDerivLTL();
    }
#endif  // #ifdef CONTACTFDGAPLTL
  }


  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    gact = Core::LinAlg::CreateVector(*gactivedofs_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
    }
  }
  else
  {
    gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
  }
  /**********************************************************************/
  /* calculate                                                          */
  /**********************************************************************/
  /* build global matrix tmatrix_ with tangent vectors of active nodes  */
  /* and global matrix nmatrix_ with normal vectors of active nodes     */
  /* and global matrix s with normal+D+M derivatives of active nodes    */
  /* and global matrix tderivmatrix_ with tangent derivatives           */
  /*     of active nodes                                                */
  /* and global matrix nderivmatrix_ with normal derivatives            */
  /*     of active nodes                                                */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /* and tangential right-hand side (incr)                              */
  /**********************************************************************/
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_tn(tmatrix_, nmatrix_);
    interface_[i]->assemble_s(*smatrix_);
    interface_[i]->assemble_t_nderiv(tderivmatrix_, nderivmatrix_);
    interface_[i]->assemble_lin_dm(*lindmatrix_, *linmmatrix_);

    if (system_type() != Inpar::CONTACT::system_condensed)
    {
      interface_[i]->assemble_inactiverhs(*inactiverhs_);
      interface_[i]->assemble_tangrhs(*tangrhs_);
    }
  }
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    // fill_complete() global matrix T
    tmatrix_->complete(*gactivedofs_, *gactivedofs_);

    // fill_complete() global matrix N
    if (nmatrix_ != Teuchos::null) nmatrix_->complete(*gactivedofs_, *gactivedofs_);
    smatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
    tderivmatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
    if (nderivmatrix_ != Teuchos::null) nderivmatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
  }
  else
  {
    // fill_complete() global matrix T
    tmatrix_->complete(*gactivedofs_, *gactivet_);

    // fill_complete() global matrix N
    if (nmatrix_ != Teuchos::null) nmatrix_->complete(*gactivedofs_, *gactiven_);

    // fill_complete() global matrix S
    smatrix_->complete(*gsmdofrowmap_, *gactiven_);

    // fill_complete() global matrix Tderiv
    // (actually gsdofrowmap_ is in general sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    tderivmatrix_->complete(*gsmdofrowmap_, *gactivet_);

    // fill_complete() global matrix Nderiv
    if (nderivmatrix_ != Teuchos::null) nderivmatrix_->complete(*gsmdofrowmap_, *gactiven_);
  }

  // fill_complete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    if (lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      if (parallel_redistribution_status())
        trafo_ = Mortar::MatrixRowTransform(trafo_, gsmdofrowmap_);
      lindmatrix_ =
          Core::LinAlg::MLMultiply(*lindmatrix_, false, *trafo_, false, false, false, true);
      linmmatrix_ =
          Core::LinAlg::MLMultiply(*linmmatrix_, false, *trafo_, false, false, false, true);
      smatrix_ = Core::LinAlg::MLMultiply(*smatrix_, false, *trafo_, false, false, false, true);
      tderivmatrix_ =
          Core::LinAlg::MLMultiply(*tderivmatrix_, false, *trafo_, false, false, false, true);
    }
    else
    {
      // modify lindmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*invtrafo_, true, *lindmatrix_, false, false, false, true);
      lindmatrix_ = temp1;
    }
  }

  // do reagularization scaling
  if (regularized_) evaluate_regularization_scaling(gact);

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype_ == Inpar::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin &&
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD") !=
            Inpar::Mortar::lagmult_const)
      FOUR_C_THROW("Condensation only for dual LM");

    /**********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                      */
    /**********************************************************************/
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));

    // for nonsmooth contact inverting D is more complex:
    // Note: this invertation if only applicable when vertex, edge and surface nodes
    // are involved. For a falling coin (only surface and edge nodes), a special but
    // more easy implementation is needed.
    if (nonSmoothContact_)
    {
      // 1. split d matrix in vertex edge and surf part
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dss, dsev, devs, devev;
      Teuchos::RCP<Epetra_Map> gEVdofs;  // merged edge and vertex dofs

      // get dss
      Core::LinAlg::SplitMatrix2x2(
          dmatrix_, gsdofSurf_, gEVdofs, gsdofSurf_, gEVdofs, dss, dsev, devs, devev);

      // get dse and dsv
      Teuchos::RCP<Epetra_Map> temp;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dse, dsv;

      Core::LinAlg::SplitMatrix2x2(
          dsev, gsdofSurf_, temp, gsdofEdge_, gsdofVertex_, dse, dsv, tempmtx1, tempmtx2);

      // get dee dev dve dvv
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dee, dev, dve, dvv;
      Core::LinAlg::SplitMatrix2x2(
          devev, gsdofEdge_, gsdofVertex_, gsdofEdge_, gsdofVertex_, dee, dev, dve, dvv);

      // 2. invert diagonal matrices dss dee dvv
      Teuchos::RCP<Epetra_Vector> diagV = Core::LinAlg::CreateVector(*gsdofVertex_, true);
      Teuchos::RCP<Epetra_Vector> diagE = Core::LinAlg::CreateVector(*gsdofEdge_, true);
      Teuchos::RCP<Epetra_Vector> diagS = Core::LinAlg::CreateVector(*gsdofSurf_, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdV =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dvv));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdE =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dee));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdS =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dss));

      int err = 0;

      // extract diagonal of invd into diag
      invdV->extract_diagonal_copy(*diagV);
      invdE->extract_diagonal_copy(*diagE);
      invdS->extract_diagonal_copy(*diagS);

      // set zero diagonal values to dummy 1.0
      for (int i = 0; i < diagV->MyLength(); ++i)
        if (abs((*diagV)[i]) < 1e-12) (*diagV)[i] = 1.0;
      for (int i = 0; i < diagE->MyLength(); ++i)
        if (abs((*diagE)[i]) < 1e-12) (*diagE)[i] = 1.0;
      for (int i = 0; i < diagS->MyLength(); ++i)
        if (abs((*diagS)[i]) < 1e-12) (*diagS)[i] = 1.0;

      // scalar inversion of diagonal values
      err = diagV->Reciprocal(*diagV);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagE->Reciprocal(*diagE);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagS->Reciprocal(*diagS);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

      // re-insert inverted diagonal into invd
      err = invdV->replace_diagonal_values(*diagV);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdE->replace_diagonal_values(*diagE);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdS->replace_diagonal_values(*diagS);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

      // 3. multiply all sub matrices
      invd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, true, true));

      dse->scale(-1.0);
      dsv->scale(-1.0);
      dev->scale(-1.0);

      // inv_dse
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum1;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dse;
      dum1 = Core::LinAlg::MLMultiply(*dse, false, *invdE, false, false, false, true);
      dinv_dse = Core::LinAlg::MLMultiply(*invdS, false, *dum1, false, false, false, true);

      // inv_dev
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dev;
      dum2 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dinv_dev = Core::LinAlg::MLMultiply(*invdE, false, *dum2, false, false, false, true);

      // inv_dsv part1
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum3;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum4;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum5;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv1;
      dum3 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dum4 = Core::LinAlg::MLMultiply(*invdE, false, *dum3, false, false, false, true);
      dum5 = Core::LinAlg::MLMultiply(*dse, false, *dum4, false, false, false, true);
      dinv_dsv1 = Core::LinAlg::MLMultiply(*invdS, false, *dum5, false, false, false, true);

      // inv_dsv part2
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum6;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv2;
      dum6 = Core::LinAlg::MLMultiply(*dsv, false, *invdV, false, false, false, true);
      dinv_dsv2 = Core::LinAlg::MLMultiply(*invdS, false, *dum6, false, false, false, true);

      // diagonal entries
      invd->add(*invdS, false, 1.0, 1.0);
      invd->add(*invdE, false, 1.0, 1.0);
      invd->add(*invdV, false, 1.0, 1.0);

      invd->add(*dinv_dev, false, 1.0, 1.0);
      invd->add(*dinv_dse, false, 1.0, 1.0);
      invd->add(*dinv_dsv1, false, 1.0, 1.0);
      invd->add(*dinv_dsv2, false, 1.0, 1.0);

      invd->complete();
    }
    // standard inverse diagonal matrix:
    else
    {
      Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
      int err = 0;

      // extract diagonal of invd into diag
      invd->extract_diagonal_copy(*diag);

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
      err = invd->replace_diagonal_values(*diag);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
    }

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);

    /**********************************************************************/
    /* (2) Add contact stiffness terms to kteff                           */
    /**********************************************************************/
    // declare sparse matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);

    if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // basis transformation
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
          Core::LinAlg::CreateIdentityMatrix(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (parallel_redistribution_status())
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

    // transform if necessary
    if (parallel_redistribution_status())
    {
      lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
      linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
    }

    kteffmatrix->un_complete();
    kteffmatrix->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
    kteffmatrix->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
    kteffmatrix->complete();

    /**********************************************************************/
    /* (3) Split kteff into 3x3 matrix blocks                             */
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
    if (parallel_redistribution_status())
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
    /* (4) Split feff into 3 subvectors                                   */
    /**********************************************************************/
    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fs, fm, fn;

    // temporarily we need the group sm
    Teuchos::RCP<Epetra_Vector> fsm;

    // do the vector splitting smn -> sm+n
    if (parallel_redistribution_status())
    {
      // split and transform to redistributed maps
      Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
      Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
      Core::LinAlg::export_to(*fsm, *fsmtemp);
      fsm = fsmtemp;
    }
    else
    {
      // only split, no need to transform
      Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
    }

    // abbreviations for slave  and master set
    const int sset = gsdofrowmap_->NumGlobalElements();
    const int mset = gmdofrowmap_->NumGlobalElements();

    // we want to split fsm into 2 groups s,m
    fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

    // do the vector splitting sm -> s+m
    Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

    // store some stuff for static condensation of LM
    fs_ = fs;
    invd_ = invd;
    ksn_ = ksn;
    ksm_ = ksm;
    kss_ = kss;

    //----------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //----------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    // D^(-1)    ---->   T * D^(-1)
    // \hat{M}   ---->   T * \hat{M}
    //----------------------------------------------------------------------
    if (is_dual_quad_slave_trafo() && lagmultquad != Inpar::Mortar::lagmult_lin)
    {
      // modify dmatrix_, invd_ and mhatmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
          Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
          Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
          Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
      dmatrix_ = temp2;
      invd_ = temp3;
      mhatmatrix_ = temp4;
    }

    /**********************************************************************/
    /* (5) Split slave quantities into active / inactive                  */
    /**********************************************************************/
    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;

    // do the splitting
    Core::LinAlg::SplitMatrix2x2(
        kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
    Core::LinAlg::SplitMatrix2x2(
        ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
    Core::LinAlg::SplitMatrix2x2(
        ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
    Core::LinAlg::SplitMatrix2x2(
        kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

    // abbreviations for active and inactive set
    const int aset = gactivedofs_->NumGlobalElements();
    const int iset = gidofs->NumGlobalElements();

    // we want to split fsmod into 2 groups a,i
    Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // do the vector splitting s -> a+i
    Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

    /**********************************************************************/
    /* (6) Isolate necessary parts from invd and mhatmatrix               */
    /**********************************************************************/
    // active part of invd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);

    // coupling part of dmatrix (only nonzero for 3D quadratic case!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
    Core::LinAlg::SplitMatrix2x2(
        dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

    // do the multiplication dhat = invda * dai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
    if (aset && iset)
      dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
    dhat->complete(*gidofs, *gactivedofs_);

    // active part of mmatrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
    Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
        tempmtx1, tempmtx2, tempmtx3);

    // do the multiplication mhataam = invda * mmatrixa
    // (this is only different from mhata for 3D quadratic case!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
    if (aset)
      mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
    mhataam->complete(*gmdofrowmap_, *gactivedofs_);

    // for the case without full linearization, we still need the
    // "classical" active part of mhat, which is isolated here
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
    Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
        tempmtx1, tempmtx2, tempmtx3);

    // scaling of invd and dai
    invda->scale(1 / (1 - alphaf_));
    dai->scale(1 - alphaf_);

    save_coupling_matrices(dhat, mhataam, invda);
    /**********************************************************************/
    /* (7) Build the final K blocks                                       */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // knn: nothing to do

    // knm: nothing to do

    // kns: nothing to do

    //---------------------------------------------------------- SECOND LINE
    // kmn: add T(mhataam)*kan
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmnmod->add(*kmn, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
    kmnmod->add(*kmnadd, false, 1.0, 1.0);
    kmnmod->complete(kmn->domain_map(), kmn->row_map());

    // kmm: add T(mhataam)*kam
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmmmod->add(*kmm, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
    kmmmod->add(*kmmadd, false, 1.0, 1.0);
    kmmmod->complete(kmm->domain_map(), kmm->row_map());

    // kmi: add T(mhataam)*kai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
      kmimod->add(*kmi, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
          Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
      kmimod->add(*kmiadd, false, 1.0, 1.0);
      kmimod->complete(kmi->domain_map(), kmi->row_map());
    }

    // kma: add T(mhataam)*kaa
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
      kmamod->add(*kma, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
          Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
      kmamod->add(*kmaadd, false, 1.0, 1.0);
      kmamod->complete(kma->domain_map(), kma->row_map());
    }

    //----------------------------------------------------------- THIRD LINE
    //------------------- FOR 3D QUADRATIC CASE ----------------------------
    // kin: subtract T(dhat)*kan --
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kinmod->add(*kin, false, 1.0, 1.0);
    if (aset && iset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
      kinmod->add(*kinadd, false, -1.0, 1.0);
    }
    kinmod->complete(kin->domain_map(), kin->row_map());

    // kim: subtract T(dhat)*kam
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kimmod->add(*kim, false, 1.0, 1.0);
    if (aset && iset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
      kimmod->add(*kimadd, false, -1.0, 1.0);
    }
    kimmod->complete(kim->domain_map(), kim->row_map());

    // kii: subtract T(dhat)*kai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
    if (iset)
    {
      kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
      kiimod->add(*kii, false, 1.0, 1.0);
      if (aset)
      {
        Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
            Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
        kiimod->add(*kiiadd, false, -1.0, 1.0);
      }
      kiimod->complete(kii->domain_map(), kii->row_map());
    }

    // kia: subtract T(dhat)*kaa
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
    if (iset && aset)
    {
      kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
      kiamod->add(*kia, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
      kiamod->add(*kiaadd, false, -1.0, 1.0);
      kiamod->complete(kia->domain_map(), kia->row_map());
    }

    //---------------------------------------------------------- FOURTH LINE
    // nothing to do - execpt its regularized contact see --> do_regularization_scaling()

    //----------------------------------------------------------- FIFTH LINE
    // kan: multiply tmatrix with invda and kan
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kanmod;
    if (aset)
    {
      kanmod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
      kanmod = Core::LinAlg::MLMultiply(*kanmod, false, *kan, false, false, false, true);
    }

    // kam: multiply tmatrix with invda and kam
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kammod;
    if (aset)
    {
      kammod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
      kammod = Core::LinAlg::MLMultiply(*kammod, false, *kam, false, false, false, true);
    }

    // kai: multiply tmatrix with invda and kai
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kaimod;
    if (aset && iset)
    {
      kaimod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
      kaimod = Core::LinAlg::MLMultiply(*kaimod, false, *kai, false, false, false, true);
    }

    // kaa: multiply tmatrix with invda and kaa
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kaamod;
    if (aset)
    {
      kaamod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
      kaamod = Core::LinAlg::MLMultiply(*kaamod, false, *kaa, false, false, false, true);
    }

    /**********************************************************************/
    /* (8) Build the final f blocks                                       */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // fn: nothing to do

    //---------------------------------------------------------- SECOND LINE
    // fm: add alphaf * old contact forces (t_n)
    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Mold^T * zold to fit
    if (is_self_contact())
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> tempvecm2 = Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(mold_->row_map()));
      if (mold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      mold_->multiply(true, *zoldexp, *tempvecm2);
      if (mset) Core::LinAlg::export_to(*tempvecm2, *tempvecm);
      fm->Update(alphaf_, *tempvecm, 1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->multiply(true, *zold_, *tempvecm);
      fm->Update(alphaf_, *tempvecm, 1.0);
    }

    // fs: prepare alphaf * old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Dold^T * zold to fit
    if (is_self_contact())
    {
      Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      dold_->multiply(true, *zoldexp, *tempvec);
      if (sset) Core::LinAlg::export_to(*tempvec, *fsadd);
    }
    // if there is no self contact everything is ok
    else
    {
      dold_->multiply(true, *zold_, *fsadd);
    }

    // fa: subtract alphaf * old contact forces (t_n)
    if (aset)
    {
      Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      Core::LinAlg::export_to(*fsadd, *faadd);
      fa->Update(-alphaf_, *faadd, 1.0);
    }

    // fm: add T(mhat)*fa
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (aset) mhataam->multiply(true, *fa, *fmmod);
    fmmod->Update(1.0, *fm, 1.0);

    //----------------------------------------------------------- THIRD LINE
    // fi: subtract alphaf * old contact forces (t_n)
    if (iset)
    {
      Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
      Core::LinAlg::export_to(*fsadd, *fiadd);
      fi->Update(-alphaf_, *fiadd, 1.0);
    }

    // fi: add T(dhat)*fa
    Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
    if (iset && aset) dhat->multiply(true, *fa, *fimod);
    fimod->Update(1.0, *fi, -1.0);

    //---------------------------------------------------------- FOURTH LINE
    // gactive: nothing to do  - execpt its regularized contact see --> do_regularization_scaling()

    //----------------------------------------------------------- FIFTH LINE
    // fa: multiply tmatrix with invda and fa
    Teuchos::RCP<Epetra_Vector> famod;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tinvda;
    if (aset)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        famod = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      else
        famod = Teuchos::rcp(new Epetra_Vector(*gactivet_));

      tinvda = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
      tinvda->multiply(false, *fa, *famod);
    }

    /********************************************************************/
    /* (9) Transform the final K blocks                                 */
    /********************************************************************/
    // The row maps of all individual matrix blocks are transformed to
    // the parallel layout of the underlying problem discretization.
    // Of course, this is only necessary in the parallel redistribution
    // case, where the contact interfaces have been redistributed
    // independently of the underlying problem discretization.

    if (parallel_redistribution_status())
    {
      //----------------------------------------------------------- FIRST LINE
      // nothing to do (ndof-map independent of redistribution)

      //---------------------------------------------------------- SECOND LINE
      kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
      kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
      if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
      if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

      //----------------------------------------------------------- THIRD LINE
      if (iset)
      {
        kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
        kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
        kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
        if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
      }

      //---------------------------------------------------------- FOURTH LINE
      if (aset) smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);

      //----------------------------------------------------------- FIFTH LINE
      if (aset)
      {
        kanmod = Mortar::MatrixRowTransform(kanmod, pgsdofrowmap_);
        kammod = Mortar::MatrixRowTransform(kammod, pgsdofrowmap_);
        kaamod = Mortar::MatrixRowTransform(kaamod, pgsdofrowmap_);
        if (iset) kaimod = Mortar::MatrixRowTransform(kaimod, pgsdofrowmap_);
        tderivmatrix_ = Mortar::MatrixRowTransform(tderivmatrix_, pgsdofrowmap_);
      }
    }

    /**********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                  */
    /**********************************************************************/

    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
    Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

    // Add all additional contributions from regularized contact to kteffnew and feffnew!
    if (regularized_)
      do_regularization_scaling(aset, iset, invda, kan, kam, kai, kaa, fa, kteffnew, feffnew);

    //----------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->add(*knn, false, 1.0, 1.0);
    kteffnew->add(*knm, false, 1.0, 1.0);
    if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

    //---------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->add(*kmnmod, false, 1.0, 1.0);
    kteffnew->add(*kmmmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
    if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

    //----------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
    if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
    if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

    //---------------------------------------------------------- FOURTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);

    //----------------------------------------------------------- FIFTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->add(*kanmod, false, 1.0, 1.0);
    if (aset) kteffnew->add(*kammod, false, 1.0, 1.0);
    if (aset && iset) kteffnew->add(*kaimod, false, 1.0, 1.0);
    if (aset) kteffnew->add(*kaamod, false, 1.0, 1.0);
    if (aset) kteffnew->add(*tderivmatrix_, false, -1.0, 1.0);

    // fill_complete kteffnew (square)
    kteffnew->complete();

    /**********************************************************************/
    /* (11) Global setup of feffnew (including contact)                   */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fn, *fnexp);
    feffnew->Update(1.0, *fnexp, 1.0);

    //---------------------------------------------------------- SECOND LINE
    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fmmod, *fmmodexp);
    feffnew->Update(1.0, *fmmodexp, 1.0);

    //----------------------------------------------------------- THIRD LINE
    // add i subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fimodexp;
    if (iset)
    {
      fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fimod, *fimodexp);
      feffnew->Update(1.0, *fimodexp, 1.0);
    }

    //---------------------------------------------------------- FOURTH LINE
    // add weighted gap vector to feffnew, if existing
    Teuchos::RCP<Epetra_Vector> gexp;
    if (aset)
    {
      gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*gact, *gexp);
      feffnew->Update(-1.0, *gexp, 1.0);
    }

    //----------------------------------------------------------- FIFTH LINE
    // add a subvector to feffnew
    Teuchos::RCP<Epetra_Vector> famodexp;
    if (aset)
    {
      famodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*famod, *famodexp);
      feffnew->Update(1.0, *famodexp, 1.0);
    }

    // finally do the replacement
    kteff = kteffnew;
    feff = feffnew;
  }

  //************************************************************************
  //************************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //************************************************************************
  //************************************************************************
  else
  {
    //----------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //----------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    //----------------------------------------------------------------------
    if (is_dual_quad_slave_trafo())
    {
      // modify dmatrix_
      // Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
      // Core::LinAlg::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true); dmatrix_    =
      // temp2;
      if (lagmultquad == Inpar::Mortar::lagmult_lin)
      {
        // basis transformation
        Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
            Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
        Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
            Core::LinAlg::CreateIdentityMatrix(*gndofrowmap_);
        systrafo->add(*eye, false, 1.0, 1.0);
        if (parallel_redistribution_status())
          trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
        systrafo->add(*trafo_, false, 1.0, 1.0);
        systrafo->complete();

        // apply basis transformation to K and f
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
    }

    // transform if necessary
    if (parallel_redistribution_status())
    {
      lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
      linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
    }

    // add contact stiffness
    kteff->un_complete();
    kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
    kteff->complete();

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold^T * zold / D^T * z to fit
    // thus we have to export the products Mold^T * zold / M^T * z to fit
    if (is_self_contact())
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecd =
          Teuchos::rcp(new Epetra_Vector(dmatrix_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
      if (dmatrix_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*z_, *zexp);
      dmatrix_->multiply(true, *zexp, *tempvecd);
      Core::LinAlg::export_to(*tempvecd, *fsexp);
      feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecm =
          Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));
      mmatrix_->multiply(true, *zexp, *tempvecm);
      Core::LinAlg::export_to(*tempvecm, *fmexp);
      feff->Update(1.0 - alphaf_, *fmexp, 1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecdold =
          Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
      Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
      dold_->multiply(true, *zoldexp, *tempvecdold);
      Core::LinAlg::export_to(*tempvecdold, *fsoldexp);
      feff->Update(-alphaf_, *fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Teuchos::RCP<Epetra_Vector> tempvecmold =
          Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
      mold_->multiply(true, *zoldexp, *tempvecmold);
      Core::LinAlg::export_to(*tempvecmold, *fmoldexp);
      feff->Update(alphaf_, *fmoldexp, 1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      // add contact force terms
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

      // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
      check_conservation_laws(*fs, *fm);
//      double normlambda = 0.0;
//      z_->Norm1(&normlambda);
//      std::cout << "normlambda =  " << normlambda << std::endl;
#endif

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dold_->multiply(true, *zold_, *fsold);
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fsold, *fsoldexp);
      feff->Update(-alphaf_, *fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->multiply(true, *zold_, *fmold);
      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fmold, *fmoldexp);
      feff->Update(alphaf_, *fmoldexp, 1.0);
    }
  }

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif  // #ifdef CONTACTFDGAP
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    std::cout << *tderivmatrix_ << std::endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif  // #ifdef CONTACTFDTANGLM

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::build_saddle_point_system(
    Teuchos::RCP<Core::LinAlg::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd,
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,
    Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
    Teuchos::RCP<Epetra_Vector>& blockrhs)
{
  // Check for saddle-point formulation
  if (system_type() != Inpar::CONTACT::system_saddlepoint)
    FOUR_C_THROW("Invalid system type! Cannot build a saddle-point system for this system type.");

  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle =
      Teuchos::rcp(new Epetra_Vector(*(dbcmaps->full_map())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->cond_map())));
  temp->PutScalar(1.0);
  Core::LinAlg::export_to(*temp, *dirichtoggle);

  // Initialize constraint matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kzz =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, true, true));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kzd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kdz =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));

  // Declare transformed constraint matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkdz = Teuchos::null;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkzd = Teuchos::null;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkzz = Teuchos::null;

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // build constraint matrix kdz
  kdz->add(*dmatrix_, true, 1.0 - alphaf_, 1.0);
  kdz->add(*mmatrix_, true, -(1.0 - alphaf_), 1.0);
  kdz->complete(*gsdofrowmap_, *gdisprowmap_);

  // *** CASE 1: FRICTIONLESS CONTACT ************************************
  if (!friction_)
  {
    // build constraint matrix kzd
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
    {
      if (gactivedofs_->NumGlobalElements())
      {
        kzd->add(*smatrix_, false, 1.0, 1.0);
        kzd->add(*tderivmatrix_, false, 1.0, 1.0);
      }
    }
    else
    {
      if (gactiven_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
      if (gactivet_->NumGlobalElements()) kzd->add(*tderivmatrix_, false, 1.0, 1.0);
    }
    kzd->complete(*gdisprowmap_, *gsdofrowmap_);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    // build constraint matrix kzz
    if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
    if (gactivet_->NumGlobalElements()) kzz->add(*tmatrix_, false, 1.0, 1.0);
    kzz->complete(*gsdofrowmap_, *gsdofrowmap_);
  }

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 2: FRICTIONAL CONTACT **************************************
  else
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
    Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);

    // build constraint matrix kzd
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
    {
      if (gactivedofs_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
      if (gstickdofs->NumGlobalElements()) kzd->add(*linstickDIS_, false, 1.0, 1.0);
      if (gslipdofs_->NumGlobalElements()) kzd->add(*linslipDIS_, false, 1.0, 1.0);
    }
    else
    {
      if (gactiven_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
      if (gstickt->NumGlobalElements()) kzd->add(*linstickDIS_, false, 1.0, 1.0);
      if (gslipt_->NumGlobalElements()) kzd->add(*linslipDIS_, false, 1.0, 1.0);
    }
    kzd->complete(*gdisprowmap_, *gsdofrowmap_);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    // build constraint matrix kzz
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
    {
      if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
      if (gstickdofs->NumGlobalElements()) kzz->add(*linstickLM_, false, 1.0, 1.0);
      if (gslipdofs_->NumGlobalElements()) kzz->add(*linslipLM_, false, 1.0, 1.0);
    }
    else
    {
      if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
      if (gstickt->NumGlobalElements()) kzz->add(*linstickLM_, false, 1.0, 1.0);
      if (gslipt_->NumGlobalElements()) kzz->add(*linslipLM_, false, 1.0, 1.0);
    }
    kzz->complete(*gsdofrowmap_, *gsdofrowmap_);
  }

  /* Step 1: Transform matrix blocks to current Lagrange multiplier dof_row_map
   *
   * This does not change the parallel layout, but only replace the slave displacement GIDs with the
   * Lagrange multiplier DOF GIDs. There is no communication involved.
   */
  {
    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = Mortar::MatrixRowTransformGIDs(kzd, glmdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (matrix_row_col_transform)
    trkzz = Mortar::MatrixRowColTransformGIDs(kzz, glmdofrowmap_, glmdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = Mortar::MatrixColTransformGIDs(kdz, glmdofrowmap_);
  }

  /* Step 2: Transform matrix blocks back to original dof_row_map as it was prior to the
   * redistribution of the interface
   *
   * Now, we keep the GID numbering, but change the parallel layout. This actually moves data
   * between MPI ranks.
   */
  if (parallel_redistribution_status())
  {
    trkzd = Mortar::matrix_row_col_transform(trkzd, pglmdofrowmap_, problem_dofs());
    trkzz = Mortar::matrix_row_col_transform(trkzz, pglmdofrowmap_, pglmdofrowmap_);
    trkdz = Mortar::matrix_row_col_transform(trkdz, problem_dofs(), pglmdofrowmap_);
  }

  // Assemble the saddle point system
  {
    // Get the standard stiffness matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiffmt =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kdd);

    // Initialize merged system (matrix, rhs, sol)
    Teuchos::RCP<Epetra_Map> mergedmap = Teuchos::null;
    if (parallel_redistribution_status())
      mergedmap = Core::LinAlg::MergeMap(problem_dofs(), pglmdofrowmap_, false);
    else
      mergedmap = Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);

    Teuchos::RCP<Epetra_Vector> mergedrhs = Core::LinAlg::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedsol = Core::LinAlg::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedzeros = Core::LinAlg::CreateVector(*mergedmap);

    /* ToDo (mayr.mt) Is this due to symmetry BCs? Basically, slave DOFs should not carry any
     * Dirichlet BCs.
     */
    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    {
      Core::LinAlg::export_to(*dirichtoggle, *dirichtoggleexp);
      Teuchos::RCP<Epetra_Vector> lmDBC = Teuchos::null;
      if (parallel_redistribution_status())
        lmDBC = Teuchos::rcp(new Epetra_Vector(*pgsdofrowmap_, true));
      else
        lmDBC = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      Core::LinAlg::export_to(*dirichtoggle, *lmDBC);

      if (parallel_redistribution_status())
        lmDBC->ReplaceMap(*pglmdofrowmap_);
      else
        lmDBC->ReplaceMap(*glmdofrowmap_);

      Teuchos::RCP<Epetra_Vector> lmDBCexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      Core::LinAlg::export_to(*lmDBC, *lmDBCexp);
      if (dirichtoggleexp->Update(1., *lmDBCexp, 1.)) FOUR_C_THROW("Update failed.");
      trkzd->apply_dirichlet(*lmDBC, false);

      trkzz->complete();
      trkzz->apply_dirichlet(*lmDBC, true);

      // apply Dirichlet conditions to (0,0) and (0,1) blocks
      Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*problem_dofs(), true));
      Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
      Core::LinAlg::apply_dirichlet_to_system(*stiffmt, *sold, *rhscopy, *zeros, *dirichtoggle);
      trkdz->apply_dirichlet(*dirichtoggle, false);
    }

    // row map (equals domain map) extractor
    Teuchos::RCP<Core::LinAlg::MapExtractor> rowmapext = Teuchos::null;
    Teuchos::RCP<Core::LinAlg::MapExtractor> dommapext = Teuchos::null;
    if (parallel_redistribution_status())
    {
      rowmapext =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*mergedmap, pglmdofrowmap_, problem_dofs()));
      dommapext =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*mergedmap, pglmdofrowmap_, problem_dofs()));
    }
    else
    {
      rowmapext =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*mergedmap, glmdofrowmap_, problem_dofs()));
      dommapext =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*mergedmap, glmdofrowmap_, problem_dofs()));
    }

    // build block matrix for SIMPLER
    blockMat =
        Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *dommapext, *rowmapext, 81, false, false));
    // blockMat is declared as an Epetra_Operator, so we need to cast it to an actual block matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> mat =
        Teuchos::rcp_dynamic_cast<
            Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(blockMat);
    mat->assign(0, 0, Core::LinAlg::View, *stiffmt);
    mat->assign(0, 1, Core::LinAlg::View, *trkdz);
    mat->assign(1, 0, Core::LinAlg::View, *trkzd);
    mat->assign(1, 1, Core::LinAlg::View, *trkzz);
    mat->complete();

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*fd, *fresmexp);
    mergedrhs->Update(1.0, *fresmexp, 1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*constrrhs_, *constrexp);
    mergedrhs->Update(1.0, *constrexp, 1.0);

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    Core::LinAlg::apply_dirichlet_to_system(*mergedsol, *mergedrhs, *mergedzeros, *dirichtoggleexp);

    blocksol = mergedsol;
    blockrhs = mergedrhs;
  }

  return;
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::update_displacements_and_l_mincrements(
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol)
{
  // Extract results for displacement and LM increments
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::null;
  if (parallel_redistribution_status())
  {
    Teuchos::RCP<Epetra_Vector> sollmOrig = Teuchos::rcp(new Epetra_Vector(*pglmdofrowmap_));
    Teuchos::RCP<Epetra_Map> mergedmapOrig =
        Core::LinAlg::MergeMap(problem_dofs(), pglmdofrowmap_, false);
    Core::LinAlg::MapExtractor mapext(*mergedmapOrig, problem_dofs(), pglmdofrowmap_);
    mapext.extract_cond_vector(blocksol, sold);
    mapext.extract_other_vector(blocksol, sollmOrig);

    sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
    Core::LinAlg::export_to(*sollmOrig, *sollm);
    sollm->ReplaceMap(*gsdofrowmap_);
  }
  else
  {
    sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
    Teuchos::RCP<Epetra_Map> mergedmap =
        Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
    Core::LinAlg::MapExtractor mapext(*mergedmap, problem_dofs(), glmdofrowmap_);
    mapext.extract_cond_vector(blocksol, sold);
    mapext.extract_other_vector(blocksol, sollm);
    sollm->ReplaceMap(*gsdofrowmap_);
  }

  /* For self contact, slave and master sets may have changed, thus we have to reinitialize the LM
   * vector map
   */
  if (is_self_contact())
  {
    zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
    Core::LinAlg::export_to(*z_, *zincr_);  // change the map of z_
    z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
    zincr_->Update(1.0, *sollm, 0.0);  // save sollm in zincr_
    z_->Update(1.0, *zincr_, 1.0);     // update z_
  }
  else
  {
    zincr_->Update(1.0, *sollm, 0.0);
    z_->Update(1.0, *zincr_, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | calculate constraint RHS entries                      hiermeier 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_constr_rhs()
{
  if (system_type() == Inpar::CONTACT::system_condensed) return;

  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
  {
    // (re)setup the vector
    constrrhs_ = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) *
  // z^(i)_(n+1)

  // export weighted gap vector
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    gact = Core::LinAlg::CreateVector(*gactivedofs_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
    }
  }
  else
  {
    gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gactiven_->NumGlobalElements())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  Teuchos::RCP<Epetra_Vector> gact_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Core::LinAlg::export_to(*gact, *gact_exp);

  constrrhs->Update(-1.0, *gact_exp, 1.0);

  // export inactive rhs
  Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Core::LinAlg::export_to(*inactiverhs_, *inactiverhsexp);

  // build constraint rhs (1)
  constrrhs->Update(1.0, *inactiverhsexp, 1.0);

  // *** CASE 1: FRICTIONLESS CONTACT *******************************************************
  if (!friction_)
  {
    // export tangential rhs
    Teuchos::RCP<Epetra_Vector> tangrhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*tangrhs_, *tangrhs_exp);

    // build constraint rhs (2)
    constrrhs->Update(1.0, *tangrhs_exp, 1.0);
  }
  // *** CASE 2: FRICTIONAL CONTACT *******************************************************
  else
  {
    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linstickRHS_, *stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linslipRHS_, *slipexp);

    // build constraint rhs
    constrrhs->Update(1.0, *stickexp, 1.0);
    constrrhs->Update(1.0, *slipexp, 1.0);
  }

  constrrhs->ReplaceMap(*glmdofrowmap_);

  // export and set constraint rhs vector
  if (parallel_redistribution_status())
  {
    constrrhs_ = Teuchos::rcp(new Epetra_Vector(lm_dof_row_map(false)));
    Core::LinAlg::export_to(*constrrhs, *constrrhs_);
  }
  else
    constrrhs_ = constrrhs;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_force(CONTACT::ParamsInterface& cparams)
{
  //---------------------------------------------------------------
  // For selfcontact the master/slave sets are updated within the -
  // contact search, see SelfBinaryTree.                          -
  // Therefore, we have to initialize the mortar matrices after   -
  // interface evaluations.                                       -
  //---------------------------------------------------------------
  if (is_self_contact())
  {
    initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
    initialize_mortar();                  // initialize mortar matrices and vectors
    assemble_mortar();                    // assemble mortar terms into global matrices
  }
  else
  {
    initialize_mortar();                  // initialize mortar matrices and vectors
    initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
    assemble_mortar();                    // assemble mortar terms into global matrices
  }

  // evaluate relative movement for friction
  if (cparams.is_predictor())
    predict_relative_movement();
  else
    evaluate_relative_movement();

  // update active set
  const bool firstTimeStepAndPredictor =
      (cparams.is_predictor() and (cparams.get_step_np() == cparams.get_restart_step() + 1));
  update_active_set_semi_smooth(firstTimeStepAndPredictor);

  // apply contact forces and stiffness
  initialize();  // init lin-matrices
  assemble_all_contact_terms();

  if (system_type() != Inpar::CONTACT::system_condensed)
  {
    eval_str_contact_rhs();  // evaluate the structure/displacement rhs
    evaluate_constr_rhs();   // evaluate the constraint rhs (saddle-point system only)

    if (constrrhs_ != Teuchos::null)
      constrrhs_->Scale(-1.0);  // scale with -1.0 --> when old structure is deleted change this!!!
  }
  else
    eval_str_contact_rhs();

  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
  {
    systrafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
        Core::LinAlg::CreateIdentityMatrix(*gndofrowmap_);
    systrafo_->add(*eye, false, 1.0, 1.0);
    if (parallel_redistribution_status())
      trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
    systrafo_->add(*trafo_, false, 1.0, 1.0);
    systrafo_->complete();

    invsystrafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
    invsystrafo_->add(*eye, false, 1.0, 1.0);
    if (parallel_redistribution_status())
      invtrafo_ = Mortar::matrix_row_col_transform(invtrafo_, pgsmdofrowmap_, pgsmdofrowmap_);
    invsystrafo_->add(*invtrafo_, false, 1.0, 1.0);
    invsystrafo_->complete();
  }

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::assemble_all_contact_terms()
{
  if (nonSmoothContact_)
  {
    nonsmooth_Penalty_force_ = Teuchos::rcp(new Epetra_Vector(*probdofs_));
    nonsmooth_Penalty_stiff_ =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*probdofs_, 100, true, true));
    Teuchos::RCP<Core::LinAlg::SparseOperator> k =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseOperator>(nonsmooth_Penalty_stiff_);
    if (!friction_)
    {
      // LTL contributions:
      add_line_to_lin_contributions(k, nonsmooth_Penalty_force_, false);

      // penalty support for master side quantities:
      add_master_contributions(k, nonsmooth_Penalty_force_, false);
    }
    else
    {
      // add_line_to_lin_contributions(kteff,feff);
      add_line_to_lin_contributions_friction(k, nonsmooth_Penalty_force_, false);
    }

    nonsmooth_Penalty_stiff_->complete();
    nonsmooth_Penalty_force_->Scale(-1.);
  }

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  if (friction_)
    assemble_all_contact_terms_friction();
  else
    assemble_all_contact_terms_frictionless();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::assemble_all_contact_terms_friction()
{
  /**********************************************************************/
  /* build global matrix t with tangent vectors of active nodes         */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix linstick with derivatives of stick nodes         */
  /* and global matrix linslip with derivatives of slip nodes           */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /**********************************************************************/
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_s(*smatrix_);
    interface_[i]->assemble_lin_dm(*lindmatrix_, *linmmatrix_);
    interface_[i]->assemble_lin_stick(*linstickLM_, *linstickDIS_, *linstickRHS_);
    interface_[i]->assemble_lin_slip(*linslipLM_, *linslipDIS_, *linslipRHS_);
    if (system_type() != Inpar::CONTACT::system_condensed)
      interface_[i]->assemble_inactiverhs(*inactiverhs_);
  }
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    smatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
  }
  else
  {
    // fill_complete() global matrix S
    smatrix_->complete(*gsmdofrowmap_, *gactiven_);
  }

  // fill_complete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);

  // fill_complete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);

  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    linstickLM_->complete(*gstickdofs, *gstickdofs);
    linstickDIS_->complete(*gsmdofrowmap_, *gstickdofs);
    linslipLM_->complete(*gslipdofs_, *gslipdofs_);
    linslipDIS_->complete(*gsmdofrowmap_, *gslipdofs_);
  }
  else
  {
    linstickLM_->complete(*gstickdofs, *gstickt);
    linstickDIS_->complete(*gsmdofrowmap_, *gstickt);

    // fill_complete global Matrix linslipLM_ and linslipDIS_
    linslipLM_->complete(*gslipdofs_, *gslipt_);
    linslipDIS_->complete(*gsmdofrowmap_, *gslipt_);
  }

  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (is_dual_quad_slave_trafo())
  {
    if (lagmultquad == Inpar::Mortar::lagmult_lin)
      FOUR_C_THROW("no linear LM interpolation for frictional contact");
    else
    {
      // modify lindmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*invtrafo_, true, *lindmatrix_, false, false, false, true);
      lindmatrix_ = temp1;
    }
  }

  if (system_type() == Inpar::CONTACT::system_condensed)
  {
    /**********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                      */
    /**********************************************************************/
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));

    // for nonsmooth contact inverting D is more complex:
    // Note: this invertation if only applicable when vertex, edge and surface nodes
    // are involved. For a falling coin (only surface and edge nodes), a special but
    // more easy implementation is needed.
    if (nonSmoothContact_)
    {
      // 1. split d matrix in vertex edge and surf part
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dss, dsev, devs, devev;
      Teuchos::RCP<Epetra_Map> gEVdofs;  // merged edge and vertex dofs

      // get dss
      Core::LinAlg::SplitMatrix2x2(
          dmatrix_, gsdofSurf_, gEVdofs, gsdofSurf_, gEVdofs, dss, dsev, devs, devev);

      // get dse and dsv
      Teuchos::RCP<Epetra_Map> temp;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dse, dsv;

      Core::LinAlg::SplitMatrix2x2(
          dsev, gsdofSurf_, temp, gsdofEdge_, gsdofVertex_, dse, dsv, tempmtx1, tempmtx2);

      // get dee dev dve dvv
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dee, dev, dve, dvv;
      Core::LinAlg::SplitMatrix2x2(
          devev, gsdofEdge_, gsdofVertex_, gsdofEdge_, gsdofVertex_, dee, dev, dve, dvv);

      // 2. invert diagonal matrices dss dee dvv
      Teuchos::RCP<Epetra_Vector> diagV = Core::LinAlg::CreateVector(*gsdofVertex_, true);
      Teuchos::RCP<Epetra_Vector> diagE = Core::LinAlg::CreateVector(*gsdofEdge_, true);
      Teuchos::RCP<Epetra_Vector> diagS = Core::LinAlg::CreateVector(*gsdofSurf_, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdV =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dvv));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdE =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dee));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdS =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dss));

      int err = 0;

      // extract diagonal of invd into diag
      invdV->extract_diagonal_copy(*diagV);
      invdE->extract_diagonal_copy(*diagE);
      invdS->extract_diagonal_copy(*diagS);

      // set zero diagonal values to dummy 1.0
      for (int i = 0; i < diagV->MyLength(); ++i)
        if (abs((*diagV)[i]) < 1e-12) (*diagV)[i] = 1.0;
      for (int i = 0; i < diagE->MyLength(); ++i)
        if (abs((*diagE)[i]) < 1e-12) (*diagE)[i] = 1.0;
      for (int i = 0; i < diagS->MyLength(); ++i)
        if (abs((*diagS)[i]) < 1e-12) (*diagS)[i] = 1.0;

      // scalar inversion of diagonal values
      err = diagV->Reciprocal(*diagV);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagE->Reciprocal(*diagE);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagS->Reciprocal(*diagS);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

      // re-insert inverted diagonal into invd
      err = invdV->replace_diagonal_values(*diagV);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdE->replace_diagonal_values(*diagE);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdS->replace_diagonal_values(*diagS);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

      // 3. multiply all sub matrices
      invd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, true, true));

      dse->scale(-1.0);
      dsv->scale(-1.0);
      dev->scale(-1.0);

      // inv_dse
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum1;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dse;
      dum1 = Core::LinAlg::MLMultiply(*dse, false, *invdE, false, false, false, true);
      dinv_dse = Core::LinAlg::MLMultiply(*invdS, false, *dum1, false, false, false, true);

      // inv_dev
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dev;
      dum2 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dinv_dev = Core::LinAlg::MLMultiply(*invdE, false, *dum2, false, false, false, true);

      // inv_dsv part1
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum3;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum4;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum5;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv1;
      dum3 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dum4 = Core::LinAlg::MLMultiply(*invdE, false, *dum3, false, false, false, true);
      dum5 = Core::LinAlg::MLMultiply(*dse, false, *dum4, false, false, false, true);
      dinv_dsv1 = Core::LinAlg::MLMultiply(*invdS, false, *dum5, false, false, false, true);

      // inv_dsv part2
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum6;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv2;
      dum6 = Core::LinAlg::MLMultiply(*dsv, false, *invdV, false, false, false, true);
      dinv_dsv2 = Core::LinAlg::MLMultiply(*invdS, false, *dum6, false, false, false, true);

      // diagonal entries
      invd->add(*invdS, false, 1.0, 1.0);
      invd->add(*invdE, false, 1.0, 1.0);
      invd->add(*invdV, false, 1.0, 1.0);

      invd->add(*dinv_dev, false, 1.0, 1.0);
      invd->add(*dinv_dse, false, 1.0, 1.0);
      invd->add(*dinv_dsv1, false, 1.0, 1.0);
      invd->add(*dinv_dsv2, false, 1.0, 1.0);

      // complete
      invd->complete();
    }
    // standard inverse diagonal matrix:
    else
    {
      Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
      int err = 0;

      // extract diagonal of invd into diag
      invd->extract_diagonal_copy(*diag);

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
      err = invd->replace_diagonal_values(*diag);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
    }

    invd_ = invd;

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);
  }

  if (is_dual_quad_slave_trafo() && lagmultquad != Inpar::Mortar::lagmult_lin)
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp2;
    if (system_type() == Inpar::CONTACT::system_condensed)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
          Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
          Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
      invd_ = temp3;
      mhatmatrix_ = temp4;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::assemble_all_contact_terms_frictionless()
{
  /**********************************************************************/
  /* calculate                                                          */
  /**********************************************************************/
  /* build global matrix tmatrix_ with tangent vectors of active nodes  */
  /* and global matrix nmatrix_ with normal vectors of active nodes     */
  /* and global matrix s with normal+D+M derivatives of active nodes    */
  /* and global matrix tderivmatrix_ with tangent derivatives           */
  /*     of active nodes                                                */
  /* and global matrix nderivmatrix_ with normal derivatives            */
  /*     of active nodes                                                */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /* and tangential right-hand side (incr)                              */
  /**********************************************************************/
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_tn(tmatrix_, nmatrix_);
    interface_[i]->assemble_s(*smatrix_);
    interface_[i]->assemble_t_nderiv(tderivmatrix_, nderivmatrix_);
    interface_[i]->assemble_lin_dm(*lindmatrix_, *linmmatrix_);

    if (system_type() != Inpar::CONTACT::system_condensed)
    {
      interface_[i]->assemble_inactiverhs(*inactiverhs_);
      interface_[i]->assemble_tangrhs(*tangrhs_);
    }
  }
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    // fill_complete() global matrix T
    tmatrix_->complete(*gactivedofs_, *gactivedofs_);

    // fill_complete() global matrix N
    if (nmatrix_ != Teuchos::null) nmatrix_->complete(*gactivedofs_, *gactivedofs_);
    smatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
    tderivmatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
    if (nderivmatrix_ != Teuchos::null) nderivmatrix_->complete(*gsmdofrowmap_, *gactivedofs_);
  }
  else
  {
    // fill_complete() global matrix T
    tmatrix_->complete(*gactivedofs_, *gactivet_);

    // fill_complete() global matrix N
    if (nmatrix_ != Teuchos::null) nmatrix_->complete(*gactivedofs_, *gactiven_);

    // fill_complete() global matrix S
    smatrix_->complete(*gsmdofrowmap_, *gactiven_);

    // fill_complete() global matrix Tderiv
    // (actually gsdofrowmap_ is in general sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    tderivmatrix_->complete(*gsmdofrowmap_, *gactivet_);

    // fill_complete() global matrix Nderiv
    if (nderivmatrix_ != Teuchos::null) nderivmatrix_->complete(*gsmdofrowmap_, *gactiven_);
  }

  // fill_complete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);

  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (is_dual_quad_slave_trafo())
  {
    if (lagmultquad == Inpar::Mortar::lagmult_lin)
    {
    }
    else
    {
      // modify lindmatrix_
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*invtrafo_, true, *lindmatrix_, false, false, false, true);
      lindmatrix_ = temp1;
    }
  }

  if (system_type() == Inpar::CONTACT::system_condensed)
  {
    /**********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                      */
    /**********************************************************************/
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));

    // for nonsmooth contact inverting D is more complex:
    // Note: this invertation if only applicable when vertex, edge and surface nodes
    // are involved. For a falling coin (only surface and edge nodes), a special but
    // more easy implementation is needed.
    if (nonSmoothContact_)
    {
      // 1. split d matrix in vertex edge and surf part
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dss, dsev, devs, devev;
      Teuchos::RCP<Epetra_Map> gEVdofs;  // merged edge and vertex dofs

      // get dss
      Core::LinAlg::SplitMatrix2x2(
          dmatrix_, gsdofSurf_, gEVdofs, gsdofSurf_, gEVdofs, dss, dsev, devs, devev);

      // get dse and dsv
      Teuchos::RCP<Epetra_Map> temp;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dse, dsv;

      Core::LinAlg::SplitMatrix2x2(
          dsev, gsdofSurf_, temp, gsdofEdge_, gsdofVertex_, dse, dsv, tempmtx1, tempmtx2);

      // get dee dev dve dvv
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dee, dev, dve, dvv;
      Core::LinAlg::SplitMatrix2x2(
          devev, gsdofEdge_, gsdofVertex_, gsdofEdge_, gsdofVertex_, dee, dev, dve, dvv);

      // 2. invert diagonal matrices dss dee dvv
      Teuchos::RCP<Epetra_Vector> diagV = Core::LinAlg::CreateVector(*gsdofVertex_, true);
      Teuchos::RCP<Epetra_Vector> diagE = Core::LinAlg::CreateVector(*gsdofEdge_, true);
      Teuchos::RCP<Epetra_Vector> diagS = Core::LinAlg::CreateVector(*gsdofSurf_, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdV =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dvv));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdE =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dee));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> invdS =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dss));

      int err = 0;

      // extract diagonal of invd into diag
      invdV->extract_diagonal_copy(*diagV);
      invdE->extract_diagonal_copy(*diagE);
      invdS->extract_diagonal_copy(*diagS);

      // set zero diagonal values to dummy 1.0
      for (int i = 0; i < diagV->MyLength(); ++i)
        if (abs((*diagV)[i]) < 1e-12) (*diagV)[i] = 1.0;
      for (int i = 0; i < diagE->MyLength(); ++i)
        if (abs((*diagE)[i]) < 1e-12) (*diagE)[i] = 1.0;
      for (int i = 0; i < diagS->MyLength(); ++i)
        if (abs((*diagS)[i]) < 1e-12) (*diagS)[i] = 1.0;

      // scalar inversion of diagonal values
      err = diagV->Reciprocal(*diagV);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagE->Reciprocal(*diagE);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");
      err = diagS->Reciprocal(*diagS);
      if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

      // re-insert inverted diagonal into invd
      err = invdV->replace_diagonal_values(*diagV);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdE->replace_diagonal_values(*diagE);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
      err = invdS->replace_diagonal_values(*diagS);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

      // 3. multiply all sub matrices
      invd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, true, true));

      dse->scale(-1.0);
      dsv->scale(-1.0);
      dev->scale(-1.0);

      // inv_dse
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum1;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dse;
      dum1 = Core::LinAlg::MLMultiply(*dse, false, *invdE, false, false, false, true);
      dinv_dse = Core::LinAlg::MLMultiply(*invdS, false, *dum1, false, false, false, true);

      // inv_dev
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum2;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dev;
      dum2 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dinv_dev = Core::LinAlg::MLMultiply(*invdE, false, *dum2, false, false, false, true);

      // inv_dsv part1
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum3;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum4;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum5;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv1;
      dum3 = Core::LinAlg::MLMultiply(*dev, false, *invdV, false, false, false, true);
      dum4 = Core::LinAlg::MLMultiply(*invdE, false, *dum3, false, false, false, true);
      dum5 = Core::LinAlg::MLMultiply(*dse, false, *dum4, false, false, false, true);
      dinv_dsv1 = Core::LinAlg::MLMultiply(*invdS, false, *dum5, false, false, false, true);

      // inv_dsv part2
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dum6;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> dinv_dsv2;
      dum6 = Core::LinAlg::MLMultiply(*dsv, false, *invdV, false, false, false, true);
      dinv_dsv2 = Core::LinAlg::MLMultiply(*invdS, false, *dum6, false, false, false, true);

      // diagonal entries
      invd->add(*invdS, false, 1.0, 1.0);
      invd->add(*invdE, false, 1.0, 1.0);
      invd->add(*invdV, false, 1.0, 1.0);

      invd->add(*dinv_dev, false, 1.0, 1.0);
      invd->add(*dinv_dse, false, 1.0, 1.0);
      invd->add(*dinv_dsv1, false, 1.0, 1.0);
      invd->add(*dinv_dsv2, false, 1.0, 1.0);

      invd->complete();
    }
    // standard inverse diagonal matrix:
    else
    {
      Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
      int err = 0;

      // extract diagonal of invd into diag
      invd->extract_diagonal_copy(*diag);

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
      err = invd->replace_diagonal_values(*diag);
      if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);
    }

    invd_ = invd;

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);
  }

  if (is_dual_quad_slave_trafo() && lagmultquad != Inpar::Mortar::lagmult_lin)
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp2;
    if (system_type() == Inpar::CONTACT::system_condensed)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
          Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
          Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
      invd_ = temp3;
      mhatmatrix_ = temp4;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::assemble_contact_rhs()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    if (system_type() != Inpar::CONTACT::system_condensed)
    {
      interface_[i]->assemble_inactiverhs(*inactiverhs_);
      if (!is_friction()) interface_[i]->assemble_tangrhs(*tangrhs_);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::eval_str_contact_rhs()
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step())
  {
    strcontactrhs_ = Teuchos::null;
    return;
  }

  strcontactrhs_ = Teuchos::rcp(new Epetra_Vector(*problem_dofs(), true));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (is_self_contact())
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecd = Teuchos::rcp(new Epetra_Vector(dmatrix_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
    if (dmatrix_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*z_, *zexp);
    dmatrix_->multiply(true, *zexp, *tempvecd);
    Core::LinAlg::export_to(*tempvecd, *fsexp);
    strcontactrhs_->Update(1., *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));
    mmatrix_->multiply(true, *zexp, *tempvecm);
    Core::LinAlg::export_to(*tempvecm, *fmexp);
    strcontactrhs_->Update(-1., *fmexp, 1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->multiply(true, *z_, *fs);
    //    Core::LinAlg::MLMultiply(*lindmatrix_,false,*trafo_,false,false,false,true);
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fs, *fsexp);
    strcontactrhs_->Update(+1., *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *z_, *fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fm, *fmexp);
    strcontactrhs_->Update(-1., *fmexp, 1.0);

    // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
    check_conservation_laws(*fs, *fm);
#endif
  }


  return;
}

/*----------------------------------------------------------------------*
 | set force evaluation flag before evaluation step          farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::pre_evaluate(CONTACT::ParamsInterface& cparams)
{
  const enum Mortar::ActionType& act = cparams.get_action_type();

  switch (act)
  {
      // -------------------------------------------------------------------
      // reset force evaluation flag for predictor step
      // -------------------------------------------------------------------
    case Mortar::eval_force_stiff:
    {
      if (cparams.is_predictor()) evalForceCalled_ = false;
      break;
    }
    // -------------------------------------------------------------------
    // default
    // -------------------------------------------------------------------
    default:
    {
      // do nothing
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set force evaluation flag after evaluation                farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::post_evaluate(CONTACT::ParamsInterface& cparams)
{
  const enum Mortar::ActionType& act = cparams.get_action_type();

  switch (act)
  {
    // -------------------------------------------------------------------
    // set flag to false after force stiff evaluation
    // -------------------------------------------------------------------
    case Mortar::eval_force_stiff:
    {
      evalForceCalled_ = false;
      break;
    }
    // -------------------------------------------------------------------
    // set flag for force evaluation to true
    // -------------------------------------------------------------------
    case Mortar::eval_force:
    {
      evalForceCalled_ = true;
      break;
    }
    // -------------------------------------------------------------------
    // default
    // -------------------------------------------------------------------
    default:
    {
      // do nothing
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_force_stiff(CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine if not done before
  if (!evalForceCalled_) evaluate_force(cparams);

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::LagrangeStrategy::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, const ParamsInterface* cparams) const
{
  // if there are no active LM contact contributions
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
  {
    if (nonSmoothContact_ && bt == CONTACT::MatBlockType::displ_displ)
      return nonsmooth_Penalty_stiff_;
    else
      return Teuchos::null;
  }
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(slave_master_dof_row_map(true), 100, false, true));

      // build matrix kdd
      mat_ptr->add(*lindmatrix_, false, 1.0, 1.0);
      mat_ptr->add(*linmmatrix_, false, 1.0, 1.0);
      if (nonSmoothContact_ && !nonsmooth_Penalty_stiff_.is_null())
        mat_ptr->add(*nonsmooth_Penalty_stiff_, false, 1.0, 1.0);
      mat_ptr->complete();

      // transform parallel row/column distribution of matrix kdd
      // (only necessary in the parallel redistribution case)
      if (parallel_redistribution_status())
        mat_ptr = Mortar::matrix_row_col_transform(
            mat_ptr, slave_master_dof_row_map_ptr(false), slave_master_dof_row_map_ptr(false));

      Teuchos::RCP<Core::LinAlg::SparseMatrix> full_mat_ptr =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      full_mat_ptr->add(*mat_ptr, false, 1., 1.);
      full_mat_ptr->complete();
      if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
        full_mat_ptr =
            Core::LinAlg::MLMultiply(*invsystrafo_, true, *full_mat_ptr, false, false, false, true);

      mat_ptr = full_mat_ptr;

      break;
    }
    case CONTACT::MatBlockType::displ_lm:
    {
      // build constraint matrix kdz
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kdz_ptr =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));

      kdz_ptr->add(*dmatrix_, true, 1.0, 1.0);
      kdz_ptr->add(*mmatrix_, true, -1.0, 1.0);
      kdz_ptr->complete(*gsdofrowmap_, *gdisprowmap_);

      // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
      mat_ptr = Mortar::MatrixColTransformGIDs(kdz_ptr, lm_dof_row_map_ptr(true));

      // transform parallel row/column distribution of matrix kdz
      // (only necessary in the parallel redistribution case)
      if (parallel_redistribution_status() or is_self_contact())
        mat_ptr = Mortar::matrix_row_col_transform(
            mat_ptr, problem_dofs(), lin_system_lm_dof_row_map_ptr());

      break;
    }
    case CONTACT::MatBlockType::lm_displ:
    {
      // build constraint matrix kzd
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kzd_ptr =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100, false, true));

      // build constraint matrix kzd
      if (gactiven_->NumGlobalElements())
      {
        kzd_ptr->add(*smatrix_, false, 1.0, 1.0);

        // frictionless contact
        if (!is_friction()) kzd_ptr->add(*tderivmatrix_, false, 1.0, 1.0);

        // frictional contact
        else
        {
          //          if (gslipnodes_->NumGlobalElements())
          kzd_ptr->add(*linslipDIS_, false, 1.0, 1.0);
          //          if (gslipnodes_->NumGlobalElements()!=gactivenodes_->NumGlobalElements())
          kzd_ptr->add(*linstickDIS_, false, 1.0, 1.0);
        }
      }
      kzd_ptr->complete(*gdisprowmap_, *gsdofrowmap_);

      // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
      mat_ptr = Mortar::MatrixRowTransformGIDs(kzd_ptr, lm_dof_row_map_ptr(true));

      // transform parallel row/column distribution of matrix kzd
      // (only necessary in the parallel redistribution case)
      if (parallel_redistribution_status() or is_self_contact())
        mat_ptr = Mortar::matrix_row_col_transform(
            mat_ptr, lin_system_lm_dof_row_map_ptr(), problem_dofs());
      break;
    }
    case CONTACT::MatBlockType::lm_lm:
    {
      // build constraint matrix kzz
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kzz_ptr = Teuchos::null;
      if (is_self_contact())
      {
        kzz_ptr = Teuchos::rcp(
            new Core::LinAlg::SparseMatrix(global_self_contact_ref_map(), 100, false, true));

        Teuchos::RCP<Epetra_Map> unused_lmdofs =
            Core::LinAlg::SplitMap(global_self_contact_ref_map(), *gsdofrowmap_);
        Epetra_Vector ones = Epetra_Vector(*unused_lmdofs, false);
        ones.PutScalar(1.0);
        if (Core::LinAlg::InsertMyRowDiagonalIntoUnfilledMatrix(*kzz_ptr, ones))
          FOUR_C_THROW("Unexpected error!");
      }
      else
        kzz_ptr =
            Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100, false, true));

      // build unity matrix for inactive dofs
      Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
      Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
      ones->PutScalar(1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
      onesdiag->complete();

      // build constraint matrix kzz
      if (gidofs->NumGlobalElements()) kzz_ptr->add(*onesdiag, false, 1.0, 1.0);

      if (!is_friction())
      {
        if (gactivet_->NumGlobalElements()) kzz_ptr->add(*tmatrix_, false, 1.0, 1.0);
      }
      else
      {
        kzz_ptr->add(*linslipLM_, false, 1., 1.);
        kzz_ptr->add(*linstickLM_, false, 1., 1.);
      }

      // transform constraint matrix kzz to lmdofmap
      if (is_self_contact())
      {
        kzz_ptr->complete(*gsmdofrowmap_, *gsmdofrowmap_);
        mat_ptr = Mortar::MatrixRowColTransformGIDs(
            kzz_ptr, lin_system_lm_dof_row_map_ptr(), lin_system_lm_dof_row_map_ptr());
      }
      else
      {
        kzz_ptr->complete(*gsdofrowmap_, *gsdofrowmap_);
        mat_ptr = Mortar::MatrixRowColTransformGIDs(
            kzz_ptr, lm_dof_row_map_ptr(true), lm_dof_row_map_ptr(true));
      }

      // transform parallel row/column distribution of matrix kzz
      // (only necessary in the parallel redistribution case)
      if (parallel_redistribution_status())
        mat_ptr = Mortar::matrix_row_col_transform(
            mat_ptr, lin_system_lm_dof_row_map_ptr(), lin_system_lm_dof_row_map_ptr());

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown Solid::MatBlockType!");
      break;
    }
  }

  return mat_ptr;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::run_post_compute_x(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  if (system_type() != Inpar::CONTACT::system_condensed)
  {
    if (lm_dof_row_map(true).NumGlobalElements() > 0)
    {
      Teuchos::RCP<Epetra_Vector> zdir_ptr =
          Teuchos::rcp(new Epetra_Vector(lm_dof_row_map(true), true));
      Core::LinAlg::export_to(dir, *zdir_ptr);
      // get the current step length
      const double stepLength = cparams.get_step_length();
      // ---------------------------------------------------------------------
      // store the SCALED Lagrange multiplier increment in the contact
      // strategy
      // ---------------------------------------------------------------------
      zdir_ptr->ReplaceMap(zincr_->Map());
      zincr_->Scale(stepLength, *zdir_ptr);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::LagrangeStrategy::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  // if there are no active LM contact contributions
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
  {
    if (nonSmoothContact_ && bt == CONTACT::VecBlockType::displ)
    {
      return nonsmooth_Penalty_force_;
    }
    else
      return Teuchos::null;
  }

  Teuchos::RCP<Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
    {
      if (nonSmoothContact_ && nonsmooth_Penalty_force_ != Teuchos::null)
      {
        vec_ptr = Teuchos::rcp(new Epetra_Vector(*nonsmooth_Penalty_force_));
      }

      if (vec_ptr == Teuchos::null)
        vec_ptr = strcontactrhs_;
      else if (strcontactrhs_ != Teuchos::null)
        vec_ptr->Update(1., *strcontactrhs_, 1.);

      Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      if (is_dual_quad_slave_trafo() && Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(
                                            params(), "LM_QUAD") == Inpar::Mortar::lagmult_lin)
      {
        invsystrafo_->multiply(true, *vec_ptr, *tmp);
        vec_ptr = tmp;
      }

      break;
    }
    case CONTACT::VecBlockType::constraint:
    {
      vec_ptr = constrrhs_;
      if (is_self_contact() && !is_condensed_system())
      {
        static Teuchos::RCP<Epetra_Vector> tmp_ptr =
            Teuchos::rcp(new Epetra_Vector(lin_system_lm_dof_row_map(), false));
        tmp_ptr->PutScalar(0.0);
        Core::LinAlg::export_to(*vec_ptr, *tmp_ptr);
        vec_ptr = tmp_ptr;
      }

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
void CONTACT::LagrangeStrategy::reset_lagrange_multipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  if (system_type() != Inpar::CONTACT::system_condensed)
  {
    if (lm_dof_row_map(true).NumGlobalElements() == 0) return;

    Teuchos::RCP<Epetra_Vector> znew_ptr =
        Teuchos::rcp(new Epetra_Vector(lm_dof_row_map(true), true));
    Core::LinAlg::export_to(xnew, *znew_ptr);
    // ---------------------------------------------------------------------
    // Update the current lagrange multiplier
    // ---------------------------------------------------------------------
    znew_ptr->ReplaceMap(z_->Map());

    z_->Scale(1.0, *znew_ptr);

    // ---------------------------------------------------------------------
    // store the new Lagrange multiplier in the nodes
    // ---------------------------------------------------------------------
    store_nodal_quantities(Mortar::StrategyBase::lmupdate);
  }
  return;
}


/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::recover(Teuchos::RCP<Epetra_Vector> disi)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::LagrangeStrategy::recover");

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // shape function type and type of LM interpolation for quadratic elements
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (system_type() == Inpar::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if ((shapefcn != Inpar::Mortar::shape_dual &&
            shapefcn != Inpar::Mortar::shape_petrovgalerkin) &&
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD") !=
            Inpar::Mortar::lagmult_const)
      FOUR_C_THROW("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, tempmap, gactivedofs_, tempmap, invda, tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invdmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    invdmod->add(*invda, false, 1.0, 1.0);
    invdmod->complete();

    /**********************************************************************/
    /* Undo basis transformation to solution                              */
    /* (currently only needed for quadratic FE with linear dual LM)       */
    /**********************************************************************/
    if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // undo basis transformation to solution
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
          Core::LinAlg::CreateIdentityMatrix(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (parallel_redistribution_status())
        trafo_ = Mortar::matrix_row_col_transform(trafo_, pgsmdofrowmap_, pgsmdofrowmap_);
      systrafo->add(*trafo_, false, 1.0, 1.0);
      systrafo->complete();
      systrafo->multiply(false, *disi, *disi);
    }

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/
    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold * zold and Mold^T * zold to fit
    if (is_self_contact())
    {
      // approximate update
      // z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      // invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new Epetra_Vector((dold_->row_map())));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *mod2);
      Teuchos::RCP<Epetra_Vector> mod3 = Teuchos::rcp(new Epetra_Vector((dold_->row_map())));
      dold_->multiply(true, *mod2, *mod3);
      Teuchos::RCP<Epetra_Vector> mod4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*mod3, *mod4);
      z_->Update(-alphaf_, *mod4, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }
    else
    {
      // approximate update
      // invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      dold_->multiply(true, *zold_, *mod);
      z_->Update(-alphaf_, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }
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
    if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      // undo basis transformation to solution
      Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*problem_dofs(), 100, false, true));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
          Core::LinAlg::CreateIdentityMatrix(*gndofrowmap_);
      systrafo->add(*eye, false, 1.0, 1.0);
      if (parallel_redistribution_status())
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
 |  Update active set and check for convergence               popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::update_active_set()
{
  // get input parameter ftype
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(params(), "FRICTION");

  // assume that active set has converged and check for opposite
  activesetconv_ = true;

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // compute weighted gap
      double wgap = (*wgap_)[wgap_->Map().LID(gid)];

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        nz += cnode->mo_data().n()[k] * cnode->mo_data().lm()[k];
      }

      // friction
      double tz = 0.0;
      double tjump = 0.0;

      if (friction_)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute tangential part of Lagrange multiplier
        tz = frinode->data().txi()[0] * frinode->mo_data().lm()[0] +
             frinode->data().txi()[1] * frinode->mo_data().lm()[1];

        // compute tangential part of jump FIXME -- also the teta component should be considered
        if (Core::UTILS::IntegralValue<int>(params(), "GP_SLIP_INCR") == true)
          tjump = frinode->fri_data().jump_var()[0];
        else
          tjump = frinode->data().txi()[0] * frinode->fri_data().jump()[0] +
                  frinode->data().txi()[1] * frinode->fri_data().jump()[1];
      }

      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->active() == false)
      {
        // check for penetration
        if (wgap < 0)
        {
          cnode->active() = true;
          activesetconv_ = false;
        }
      }

      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        // check for tensile contact forces
        if (nz <= 0)  // no averaging of Lagrange multipliers
        // if (0.5*nz+0.5*nzold <= 0) // averaging of Lagrange multipliers
        {
          cnode->active() = false;
          activesetconv_ = false;

          // friction
          if (friction_) dynamic_cast<FriNode*>(cnode)->fri_data().slip() = false;
        }

        // only do something for friction
        else
        {
          // friction tresca
          if (ftype == Inpar::CONTACT::friction_tresca)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);
            double ct =
                interface_[i]->get_ct_ref()[interface_[i]->get_ct_ref().Map().LID(frinode->id())];

            // CAREFUL: friction bound is now interface-local (popp 08/2012)
            double frbound = interface_[i]->interface_params().get<double>("FRBOUND");

            if (frinode->fri_data().slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if (abs(tz + ct * tjump) - frbound <= 0)
              {
              }
              // do nothing (stick was correct)
              else
              {
                frinode->fri_data().slip() = true;
                activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if (abs(tz + ct * tjump) - frbound > 0)
              {
              }
              // do nothing (slip was correct)
              else
              {
                frinode->fri_data().slip() = false;
                activesetconv_ = false;
              }
            }
          }  // if (ftype == Inpar::CONTACT::friction_tresca)

          // friction coulomb
          if (ftype == Inpar::CONTACT::friction_coulomb)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);
            double ct =
                interface_[i]->get_ct_ref()[interface_[i]->get_ct_ref().Map().LID(frinode->id())];

            // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
            double frcoeff = interface_[i]->interface_params().get<double>("FRCOEFF");

            if (frinode->fri_data().slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if (abs(tz + ct * tjump) - frcoeff * nz <= 0)
              {
              }
              // do nothing (stick was correct)
              else
              {
                frinode->fri_data().slip() = true;
                activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if (abs(tz + ct * tjump) - frcoeff * nz > 0)
              {
              }
              // do nothing (slip was correct)
              else
              {
                frinode->fri_data().slip() = false;
                activesetconv_ = false;
              }
            }
          }  // if (ftype == Inpar::CONTACT::friction_coulomb)
        }    // if (nz <= 0)
      }      // if (cnode->Active()==false)
    }        // loop over all slave nodes
  }          // loop over all interfaces

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  get_comm().SumAll(&localcheck, &convcheck, 1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck != get_comm().NumProc())
  {
    activesetconv_ = false;
    activesetsteps_ += 1;
  }

  // update zig-zagging history (shift by one)
  if (zigzagtwo_ != Teuchos::null) zigzagthree_ = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_ != Teuchos::null) zigzagtwo_ = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_ != Teuchos::null) zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  ginactivenodes_ = Teuchos::null;
  ginactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->build_active_set();
    gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interface_[i]->active_nodes(), false);
    gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interface_[i]->active_dofs(), false);

    ginactivenodes_ =
        Core::LinAlg::MergeMap(ginactivenodes_, interface_[i]->in_active_nodes(), false);
    ginactivedofs_ = Core::LinAlg::MergeMap(ginactivedofs_, interface_[i]->in_active_dofs(), false);

    gactiven_ = Core::LinAlg::MergeMap(gactiven_, interface_[i]->active_n_dofs(), false);
    gactivet_ = Core::LinAlg::MergeMap(gactivet_, interface_[i]->active_t_dofs(), false);

    if (friction_)
    {
      gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interface_[i]->slip_nodes(), false);
      gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interface_[i]->slip_dofs(), false);
      gslipt_ = Core::LinAlg::MergeMap(gslipt_, interface_[i]->slip_t_dofs(), false);
    }
  }

  // CHECK FOR ZIG-ZAGGING / JAMMING OF THE ACTIVE SET
  // *********************************************************************
  // A problem of the active set strategy which sometimes arises is known
  // from optimization literature as jamming or zig-zagging. This means
  // that within a load/time-step the algorithm can have more than one
  // solution due to the fact that the active set is not unique. Hence the
  // algorithm jumps between the solutions of the active set. The non-
  // uniquenesss results either from highly curved contact surfaces or
  // from the FE discretization, Thus the uniqueness of the closest-point-
  // projection cannot be guaranteed.
  // *********************************************************************
  // To overcome this problem we monitor the development of the active
  // set scheme in our contact algorithms. We can identify zig-zagging by
  // comparing the current active set with the active set of the second-
  // and third-last iteration. If an identity occurs, we consider the
  // active set strategy as converged instantly, accepting the current
  // version of the active set and proceeding with the next time/load step.
  // This very simple approach helps stabilizing the contact algorithm!
  // *********************************************************************
  bool zigzagging = false;
  // FIXGIT: For tresca friction zig-zagging is not eliminated
  if (ftype != Inpar::CONTACT::friction_tresca && ftype != Inpar::CONTACT::friction_coulomb)
  {
    // frictionless contact
    if (active_set_steps() > 2)
    {
      if (zigzagtwo_ != Teuchos::null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (get_comm().MyPID() == 0)
            std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << std::endl;
        }
      }

      if (zigzagthree_ != Teuchos::null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (get_comm().MyPID() == 0)
            std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << std::endl;
        }
      }
    }
  }  // if (ftype != Inpar::CONTACT::friction_tresca && ftype != Inpar::CONTACT::friction_coulomb)


  // reset zig-zagging history
  if (activesetconv_ == true)
  {
    zigzagone_ = Teuchos::null;
    zigzagtwo_ = Teuchos::null;
    zigzagthree_ = Teuchos::null;
  }

  // output of active set status to screen
  if (get_comm().MyPID() == 0 && activesetconv_ == false)
    std::cout << "ACTIVE SET ITERATION " << active_set_steps() - 1
              << " NOT CONVERGED - REPEAT TIME STEP................." << std::endl;
  else if (get_comm().MyPID() == 0 && activesetconv_ == true)
    std::cout << "ACTIVE SET CONVERGED IN " << active_set_steps() - zigzagging
              << " STEP(S)................." << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
  }
  else
    isincontact_ = false;

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::update_active_set_semi_smooth(const bool firstStepPredictor)
{
  // FIXME: Here we do not consider zig-zagging yet!
  //  print_active_set();

  // get out gof here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = Core::UTILS::IntegralValue<int>(params(), "SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      // loop over all slave nodes on the current interface
      for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_row_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cnode = dynamic_cast<Node*>(node);

        // The nested active set strategy cannot deal with the case of
        // active nodes that have no integration segments/cells attached,
        // as this leads to zero rows in D and M and thus to singular systems.
        // However, this case might possibly happen when slave nodes slide
        // over the edge of a master body within one fixed active set step.
        // (Remark: Semi-smooth Newton has no problems in this case, as it
        // updates the active set after EACH Newton step, see below, and thus
        // would always set the corresponding nodes to INACTIVE.)
        if (cnode->active() && !cnode->has_segment() && !cnode->is_on_boundor_ce())
          FOUR_C_THROW("Active node %i without any segment/cell attached", cnode->id());
      }
    }
    return;
  }

  // get input parameter ftype
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(params(), "FRICTION");

  // assume that active set has converged and check for opposite
  activesetconv_ = true;

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    bool my_check = interface_[i]->update_active_set_semi_smooth();
    activesetconv_ = activesetconv_ and my_check;
  }  // loop over all interfaces

  // Overwrite active set with input file information in predictor step of first time step
  if (firstStepPredictor)
  {
    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->update_active_set_initial_status();
    }  // loop over all interfaces
  }

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  get_comm().SumAll(&localcheck, &convcheck, 1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck != get_comm().NumProc())
  {
    activesetconv_ = false;
    activesetsteps_ += 1;
  }

  // only if it's a full Newton step...
  // store the previous active set
  if (gactivenodes_ != Teuchos::null)
  {
    gOldActiveSlaveNodes_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));
    if (friction_) gOldslipnodes_ = Teuchos::rcp(new Epetra_Map(*gslipnodes_));
  }
  else
  {
    gOldActiveSlaveNodes_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    if (friction_) gOldslipnodes_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
  }

  // also update special flag for semi-smooth Newton convergence
  activesetssconv_ = activesetconv_;

  // update zig-zagging history (shift by one)
  if (zigzagtwo_ != Teuchos::null) zigzagthree_ = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_ != Teuchos::null) zigzagtwo_ = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_ != Teuchos::null) zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  ginactivenodes_ = Teuchos::null;
  ginactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->build_active_set();
    gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interface_[i]->active_nodes(), false);
    gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interface_[i]->active_dofs(), false);
    ginactivenodes_ =
        Core::LinAlg::MergeMap(ginactivenodes_, interface_[i]->in_active_nodes(), false);
    ginactivedofs_ = Core::LinAlg::MergeMap(ginactivedofs_, interface_[i]->in_active_dofs(), false);
    gactiven_ = Core::LinAlg::MergeMap(gactiven_, interface_[i]->active_n_dofs(), false);
    gactivet_ = Core::LinAlg::MergeMap(gactivet_, interface_[i]->active_t_dofs(), false);

    if (friction_)
    {
      gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interface_[i]->slip_nodes(), false);
      gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interface_[i]->slip_dofs(), false);
      gslipt_ = Core::LinAlg::MergeMap(gslipt_, interface_[i]->slip_t_dofs(), false);
    }
  }

  // CHECK FOR ZIG-ZAGGING / JAMMING OF THE ACTIVE SET
  // *********************************************************************
  // A problem of the active set strategy which sometimes arises is known
  // from optimization literature as jamming or zig-zagging. This means
  // that within a load/time-step the semi-smooth Newton algorithm can get
  // stuck between more than one intermediate solution due to the fact that
  // the active set decision is a discrete decision. Hence the semi-smooth
  // Newton algorithm fails to converge. The non-uniquenesss results either
  // from highly curved contact surfaces or from the FE discretization.
  // *********************************************************************
  // To overcome this problem we monitor the development of the active
  // set scheme in our contact algorithms. We can identify zig-zagging by
  // comparing the current active set with the active set of the second-
  // and third-last iteration. If an identity occurs, we interfere and
  // let the semi-smooth Newton algorithm restart from another active set
  // (e.g. intermediate set between the two problematic candidates), thus
  // leading to some kind of damped / modified semi-smooth Newton method.
  // This very simple approach helps stabilizing the contact algorithm!
  // *********************************************************************
  int zigzagging = 0;
  // FIXGIT: For friction zig-zagging is not eliminated
  if (ftype != Inpar::CONTACT::friction_tresca && ftype != Inpar::CONTACT::friction_coulomb)
  {
    // frictionless contact
    if (active_set_steps() > 2)
    {
      if (zigzagtwo_ != Teuchos::null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 1;
        }
      }

      if (zigzagthree_ != Teuchos::null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 2;
        }
      }
    }
  }  // if (ftype != Inpar::CONTACT::friction_tresca && ftype != Inpar::CONTACT::friction_coulomb)

  // output to screen
  if (get_comm().MyPID() == 0)
  {
    if (zigzagging == 1)
    {
      std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << std::endl;
    }
    else if (zigzagging == 2)
    {
      std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << std::endl;
    }
    else
    {
      // do nothing, no zig-zagging
    }
  }

  // reset zig-zagging history
  if (activesetconv_ == true)
  {
    zigzagone_ = Teuchos::null;
    zigzagtwo_ = Teuchos::null;
    zigzagthree_ = Teuchos::null;
  }

  // output of active set status to screen
  if (get_comm().MyPID() == 0 && activesetconv_ == false)
    std::cout << "ACTIVE CONTACT SET HAS CHANGED... CHANGE No. " << active_set_steps() - 1
              << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
  }
  else
    isincontact_ = false;

  return;
}

/*----------------------------------------------------------------------*
 |  update routine for ltl forces                            farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::update(Teuchos::RCP<const Epetra_Vector> dis)
{
  if (fLTL_ != Teuchos::null)
  {
    // store fLTL values for time integration
    fLTLOld_ = Teuchos::rcp(new Epetra_Vector(fLTL_->Map()), true);
    if (fLTLOld_->Update(1.0, *fLTL_, 0.0)) FOUR_C_THROW("Update went wrong");
  }

  // abstract routine
  CONTACT::AbstractStrategy::update(dis);

  if (fconservation_ == Teuchos::null) return;

  // *****************************************************************
  // This is output functionality for conservation properties of LTL
  // penalty contact
  // *****************************************************************

  //  Teuchos::RCP<Epetra_Vector> fconservationS =
  //      Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)),true);
  //  Teuchos::RCP<Epetra_Vector> fconservationM =
  //      Teuchos::rcp(new Epetra_Vector(master_dof_row_map(true)),true);
  //
  //  Core::LinAlg::export_to(*fconservation_,*fconservationS);
  //  Core::LinAlg::export_to(*fconservation_,*fconservationM);
  ////  fconservationM->Scale(-1.0);
  //  // check conservation properties:
  //  interface_[0]->EvalResultantMoment(*fconservationS, *fconservationM);
  //
  //
  //  double lssum = 0.0;   // local slave sum
  //  double gssum = 0.0;   // global slave sum
  //  double lmsum = 0.0;   // local master sum
  //  double gmsum = 0.0;   // global master sum
  //  double gcsum = 0.0;   // global complete sum
  //  // slave
  ////  for (int i=0;i<fconservationS->MyLength();++i)
  ////  {
  ////    lssum+=(*fconservationS)[i];
  ////  }
  ////  Comm().SumAll(&lssum,&gssum,1);
  ////  // master
  ////  for (int i=0;i<fconservationM->MyLength();++i)
  ////  {
  ////    lmsum+=(*fconservationM)[i];
  ////  }
  ////  Comm().SumAll(&lmsum,&gmsum,1);
  //
  //
  //
  //  // complete balance check
  //  gcsum = gssum+gmsum;
  //  if (abs(gcsum)>1.0e-11)
  //    FOUR_C_THROW("Conservation of linear momentum is not fulfilled!");
  //  if (Comm().MyPID()==0)
  //  {
  //    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  //    std::cout << ">>      Linear Momentum Conservation      <<" << std::endl;
  //    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  //    std::cout << ">>      Standard terms (lm)               <<" << std::endl;
  //    std::cout << "SLAVE:   " << std::setw(14) << gssum<< std::endl;
  //    std::cout << "MASTER:  " << std::setw(14) << gmsum << std::endl;
  //    std::cout << "Balance: " << std::setw(14) << gcsum << std::endl;
  //    std::cout << "--------------------------------------------" << std::endl;
  //  }
  //
  //
  //  FILE* MyFile = nullptr;
  //  std::ostringstream filename;
  //  const std::string filebase = "xxx";
  //  filename << filebase <<".lmom";
  //  MyFile = fopen(filename.str().c_str(), "at+");
  //
  //  // store data
  //  if (MyFile)
  //  {
  //    fprintf(MyFile, "%d\t", step);
  //    fprintf(MyFile, "%g\t", gssum);
  //    fprintf(MyFile, "%g\t", gmsum);
  //    fprintf(MyFile, "%g\n", gcsum);
  //    fclose(MyFile);
  //  }
  //  else
  //    FOUR_C_THROW("File could not be opened.");
  //
  //  ++step;

  return;
}
/*----------------------------------------------------------------------*
 |  Check conservation laws                              hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::check_conservation_laws(
    const Epetra_Vector& fs, const Epetra_Vector& fm)
{
  // change sign (adapt sign convention from the augmented Lagrange formulation)
  Teuchos::RCP<Epetra_Vector> tmpFm = Teuchos::rcp(new Epetra_Vector(fm));
  tmpFm->Scale(-1.0);

  double lssum = 0.0;  // local slave sum
  double gssum = 0.0;  // global slave sum
  double lmsum = 0.0;  // local master sum
  double gmsum = 0.0;  // global master sum
  double gcsum = 0.0;  // global complete sum
  // slave
  for (int i = 0; i < fs.MyLength(); ++i) lssum += (fs)[i];
  get_comm().SumAll(&lssum, &gssum, 1);
  // master
  for (int i = 0; i < tmpFm->MyLength(); ++i) lmsum += (*tmpFm)[i];
  get_comm().SumAll(&lmsum, &gmsum, 1);
  // complete balance check
  gcsum = gssum + gmsum;
  if (abs(gcsum) > 1.0e-11) FOUR_C_THROW("Conservation of linear momentum is not fulfilled!");
  if (get_comm().MyPID() == 0)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    std::cout << ">>      Linear Momentum Conservation      <<" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    std::cout << ">>      Standard terms (lm)               <<" << std::endl;
    std::cout << "SLAVE:   " << std::setw(14) << gssum << std::endl;
    std::cout << "MASTER:  " << std::setw(14) << gmsum << std::endl;
    std::cout << "Balance: " << std::setw(14) << gcsum << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
  }



  if (get_comm().MyPID() == 0)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    std::cout << ">>      Angular Momentum Conservation     <<" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  for (std::size_t i = 0; i < interface_.size(); ++i)
  {
    if (get_comm().MyPID() == 0)
    {
      std::cout << ">>----- Interface " << std::setw(2) << i;
      std::cout << " ---------------------<<" << std::endl;
      std::cout << ">>      Standard terms (lm)               <<" << std::endl;
    }
    interface_[i]->eval_resultant_moment(fs, *tmpFm);
  }
  if (get_comm().MyPID() == 0)
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;

  return;
}

/*-----------------------------------------------------------------------------*
 | do additional matrix manipulations for regularization scaling     ager 07/15|
 *----------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::do_regularization_scaling(bool aset, bool iset,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& invda, Teuchos::RCP<Core::LinAlg::SparseMatrix>& kan,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& kam, Teuchos::RCP<Core::LinAlg::SparseMatrix>& kai,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& kaa, Teuchos::RCP<Epetra_Vector>& fa,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteffnew, Teuchos::RCP<Epetra_Vector>& feffnew)
{
  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kanmod_n;
  if (aset)
  {
    kanmod_n = Core::LinAlg::MLMultiply(*nmatrix_, false, *invda, true, false, false, true);
    kanmod_n = Core::LinAlg::MLMultiply(*kanmod_n, false, *kan, false, false, false, true);
  }

  // kam: multiply tmatrix with invda and kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kammod_n;
  if (aset)
  {
    kammod_n = Core::LinAlg::MLMultiply(*nmatrix_, false, *invda, true, false, false, true);
    kammod_n = Core::LinAlg::MLMultiply(*kammod_n, false, *kam, false, false, false, true);
  }

  // kai: multiply tmatrix with invda and kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaimod_n;
  if (aset && iset)
  {
    kaimod_n = Core::LinAlg::MLMultiply(*nmatrix_, false, *invda, true, false, false, true);
    kaimod_n = Core::LinAlg::MLMultiply(*kaimod_n, false, *kai, false, false, false, true);
  }

  // kaa: multiply tmatrix with invda and kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaamod_n;
  if (aset)
  {
    kaamod_n = Core::LinAlg::MLMultiply(*nmatrix_, false, *invda, true, false, false, true);
    kaamod_n = Core::LinAlg::MLMultiply(*kaamod_n, false, *kaa, false, false, false, true);
  }

  /**********************************************************************/
  /* (8) Build the final f blocks                                       */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE

  Teuchos::RCP<Epetra_Vector> famod_n;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ninvda;
  if (aset)
  {
    famod_n = Teuchos::rcp(new Epetra_Vector(*gactiven_));
    ninvda = Core::LinAlg::MLMultiply(*nmatrix_, false, *invda, true, false, false, true);
    ninvda->multiply(false, *fa, *famod_n);
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  if (parallel_redistribution_status())
  {
    if (aset) nderivmatrix_ = Mortar::MatrixRowTransform(nderivmatrix_, pgsdofrowmap_);
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  if (aset) kteffnew->add(*nderivmatrix_, false, 1.0, 1.0);

  if (aset) kteffnew->add(*kanmod_n, false, -1.0, 1.0);
  if (aset) kteffnew->add(*kammod_n, false, -1.0, 1.0);
  if (aset && iset) kteffnew->add(*kaimod_n, false, -1.0, 1.0);
  if (aset) kteffnew->add(*kaamod_n, false, -1.0, 1.0);

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  Teuchos::RCP<Epetra_Vector> famodexp_n;
  if (aset)
  {
    famodexp_n = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*famod_n, *famodexp_n);
    feffnew->Update(-1.0, *famodexp_n, 1.0);
  }
}

/*----------------------------------------------------------------------*
 | calculate regularization scaling and apply it to matrixes   ager 06/15|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::evaluate_regularization_scaling(Teuchos::RCP<Epetra_Vector> gact)
{
  /********************************************************************/
  /* (0) Get Nodal Gap Scaling                                        */
  /********************************************************************/
  // vector with all values alpha
  Teuchos::RCP<Epetra_Vector> scalevec = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));

  // vecor with (alpha + dalpha/dgap*scaling + d(1-alpha)/dgap)
  Teuchos::RCP<Epetra_Vector> derivscalevec = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));

  // vector with all values (1-alpha)
  Teuchos::RCP<Epetra_Vector> scalevec2 = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));
  scalevec->PutScalar(1.0);

  {
    // setting of the parameters is still under investigation Ager Chr.
    // Global::Problem* problem = Global::Problem::instance();
    // const Teuchos::ParameterList& structdyn   = problem->structural_dynamic_params(); //just for
    // now!

    double scaling = 1e6;  // structdyn.get<double>("TOLRES"); //1e6
    //
    // sign of the term lambda*n (std is +1.0)?
    static double sign = 1.0;
    // scaling also with (1-alpha) of lambda*n term?
    static double scal1ma = 1.0;

    double epsilon = 1e-6;  // structdyn.get<double>("TOLINCO");//1e-5;+

    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      // loop over all slave nodes on the current interface
      for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_row_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

        Node* cnode = dynamic_cast<Node*>(node);

        if (cnode->active())
        {
          // normal lagrange multiplier! ... replace by multiplication
          double lmval = cnode->mo_data().n()[0] * cnode->mo_data().lm()[0] +
                         cnode->mo_data().n()[1] * cnode->mo_data().lm()[1] +
                         cnode->mo_data().n()[2] * cnode->mo_data().lm()[2];
          double alpha;
          double alphaderiv;

          std::cout << "cnode->Data().Getg(): " << cnode->data().getg() << "( eps = " << epsilon
                    << ", scal = " << scaling << " ) / lambdan: " << lmval << std::endl;


          if (cnode->data().getg() >= 0)
          {
            alpha = 0.0;
            alphaderiv = 0.0;
          }
          else
          {
            double scaledgsum = (5 * cnode->data().getg() / epsilon + 0.5);
            alpha = 0.5 * (1 - tanh(scaledgsum));
            alphaderiv = -0.5 * 5 / (cosh(scaledgsum) * cosh(scaledgsum) * epsilon);

            alphaderiv = alphaderiv * (scaling * cnode->data().getg() - sign * lmval);
          }

          int lid = gactiven_->LID(cnode->dofs()[0]);
          if (lid != -1)
          {
            (*scalevec)[lid] = alpha;
            (*derivscalevec)[lid] = alphaderiv;
          }
          else
          {
            FOUR_C_THROW("lid = -1!");
          }

          std::cout << std::setprecision(32) << "cnode->PoroData().GetnScale() for node "
                    << cnode->id() << " is: " << alpha << std::endl;
        }
      }  // loop over all slave nodes
    }    // loop over all interfaces

    // Evaluate final scaling vetors!
    {
      scalevec2->PutScalar(sign);
      scalevec2->Update(-sign * scal1ma, *scalevec, 1.0);  // 1-alpha

      scalevec->Scale(scaling);
      derivscalevec->Update(1.0, *scalevec, 1.0);
    }
  }

  // scale all matrixes!!!
  int err = smatrix_->left_scale(
      *derivscalevec);  // scale with (alpha + dalpha/dgap*scaling + d(1-alpha)/dgap*lambda*n
  if (err) FOUR_C_THROW("LeftScale failed!");

  err = nderivmatrix_->left_scale(*scalevec2);  // scale with 1.0-alpha
  if (err) FOUR_C_THROW("LeftScale failed!");

  err = nmatrix_->left_scale(*scalevec2);  // scale with 1.0-alpha
  if (err) FOUR_C_THROW("LeftScale failed!");

  Teuchos::RCP<Epetra_Vector> tmpgact =
      Teuchos::rcp(new Epetra_Vector(*gact));  // check later if this is required!!!
  gact->Multiply(1.0, *scalevec, *tmpgact, 0.0);

  return;
}

void CONTACT::LagrangeStrategy::condense_friction(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff, Epetra_Vector& rhs)
{
  Teuchos::RCP<Epetra_Vector> feff = Teuchos::rcpFromRef<Epetra_Vector>(rhs);

  feff->Update(-1. + alphaf_, *strcontactrhs_, 1.);
  feff->Scale(-1.);

  Teuchos::RCP<Core::LinAlg::SparseOperator> kteff_op =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseOperator>(kteff);

  // In case of nonsmooth contact the scenario of contacting edges (non parallel)
  // requires a penalty regularization. Here, the penalty contriutions for this
  // special case are applied:
  if (nonSmoothContact_)
  {
    // add_line_to_lin_contributions(kteff,feff);
    add_line_to_lin_contributions_friction(kteff_op, feff);

    // FD check of weighted gap g derivatives + jump for LTL case
#ifdef CONTACTFDJUMPLTL
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckJumpDerivLTL();
    }
#endif  // #ifdef CONTACTFDJUMPLTL
  }

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    gact = Core::LinAlg::CreateVector(*gactivedofs_, true);
    if (gact->GlobalLength()) Core::LinAlg::export_to(*wgap_, *gact);
  }
  else
  {
    gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  // fill_complete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);

  // shape function
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Condensation only for dual LM");

  /********************************************************************/
  /* (1) Multiply Mortar matrices: m^ = inv(d) * m                    */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invd = invd_;

  /********************************************************************/
  /* (2) Add contact stiffness terms to kteff                         */
  /********************************************************************/

  /********************************************************************/
  /* (3) Split kteff into 3x3 matrix blocks                           */
  /********************************************************************/

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
  if (parallel_redistribution_status())
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

  /********************************************************************/
  /* (4) Split feff into 3 subvectors                                 */
  /********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    Core::LinAlg::export_to(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
  }

  // abbreviations for slave set
  const int sset = gsdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

  // store some stuff for static condensation of LM
  fs_ = fs;
  invd_ = invd;
  ksn_ = ksn;
  ksm_ = ksm;
  kss_ = kss;

  /********************************************************************/
  /* (5) Split slave quantities into active / inactive, stick / slip  */
  /********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  Core::LinAlg::SplitMatrix2x2(kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
  Core::LinAlg::SplitMatrix2x2(
      ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

  // we want to split kaa into 2 groups sl,st = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> temp1map;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1mtx4, temp1mtx5;

  // we will get the stick rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gstdofs;

  Core::LinAlg::SplitMatrix2x2(
      kaa, gactivedofs_, gidofs, gstdofs, gslipdofs_, kast, kasl, temp1mtx4, temp1mtx5);

  // abbreviations for active and inactive set, stick and slip set
  const int aset = gactivedofs_->NumGlobalElements();
  const int iset = gidofs->NumGlobalElements();
  const int stickset = gstdofs->NumGlobalElements();
  const int slipset = gslipdofs_->NumGlobalElements();

  // we want to split fs into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

  // we want to split fa into 2 groups sl,st
  Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
  Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

  // do the vector splitting a -> sl+st
  if (aset) Core::LinAlg::split_vector(*gactivedofs_, *fa, gslipdofs_, fsl, gstdofs, fst);

  /********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix             */
  /********************************************************************/

  // active, stick and slip part of invd
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invda, invdsl, invdst;
  Core::LinAlg::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gslipdofs_, gstdofs, invdsl, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gstdofs, gslipdofs_, invdst, tempmtx1, tempmtx2, tempmtx3);

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
  Core::LinAlg::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  // do the multiplication dhat = invda * dai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset && iset) dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
  dhat->complete(*gidofs, *gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
  Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
      tempmtx1, tempmtx2, tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset) mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
  mhataam->complete(*gmdofrowmap_, *gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
  Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
      tempmtx1, tempmtx2, tempmtx3);

  // scaling of invd and dai
  invda->scale(1 / (1 - alphaf_));
  invdsl->scale(1 / (1 - alphaf_));
  invdst->scale(1 / (1 - alphaf_));
  dai->scale(1 - alphaf_);

  /********************************************************************/
  /* (7) Build the final K blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //-------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmnmod->add(*kmn, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
  kmnmod->add(*kmnadd, false, 1.0, 1.0);
  kmnmod->complete(kmn->domain_map(), kmn->row_map());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmmmod->add(*kmm, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
  kmmmod->add(*kmmadd, false, 1.0, 1.0);
  kmmmod->complete(kmm->domain_map(), kmm->row_map());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmimod->add(*kmi, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
    kmimod->add(*kmiadd, false, 1.0, 1.0);
    kmimod->complete(kmi->domain_map(), kmi->row_map());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmamod->add(*kma, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
    kmamod->add(*kmaadd, false, 1.0, 1.0);
    kmamod->complete(kma->domain_map(), kma->row_map());
  }

  //--------------------------------------------------------- THIRD LINE
  // kin: subtract T(dhat)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kinmod->add(*kin, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
    kinmod->add(*kinadd, false, -1.0, 1.0);
  }
  kinmod->complete(kin->domain_map(), kin->row_map());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kimmod->add(*kim, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
    kimmod->add(*kimadd, false, -1.0, 1.0);
  }
  kimmod->complete(kim->domain_map(), kim->row_map());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiimod->add(*kii, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
      kiimod->add(*kiiadd, false, -1.0, 1.0);
    }
    kiimod->complete(kii->domain_map(), kii->row_map());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiamod->add(*kia, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
    kiamod->add(*kiaadd, false, -1.0, 1.0);
    kiamod->complete(kia->domain_map(), kia->row_map());
  }

  //-------------------------------------------------------- FOURTH LINE

  //--------------------------------------------------------- FIFTH LINE
  // blocks for complementary conditions (stick nodes)

  // kstn: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstnmod;
  if (stickset)
  {
    kstnmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstnmod = Core::LinAlg::MLMultiply(*kstnmod, false, *kan, false, false, false, true);
  }

  // kstm: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstmmod;
  if (stickset)
  {
    kstmmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstmmod = Core::LinAlg::MLMultiply(*kstmmod, false, *kam, false, false, false, true);
  }

  // ksti: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstimod;
  if (stickset && iset)
  {
    kstimod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstimod = Core::LinAlg::MLMultiply(*kstimod, false, *kai, false, false, false, true);
  }

  // kstsl: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstslmod;
  if (stickset && slipset)
  {
    kstslmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstslmod = Core::LinAlg::MLMultiply(*kstslmod, false, *kasl, false, false, false, true);
  }

  // kststmod: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kststmod;
  if (stickset)
  {
    kststmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kststmod = Core::LinAlg::MLMultiply(*kststmod, false, *kast, false, false, false, true);
  }

  //--------------------------------------------------------- SIXTH LINE
  // blocks for complementary conditions (slip nodes)

  // ksln: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslnmod;
  if (slipset)
  {
    kslnmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslnmod = Core::LinAlg::MLMultiply(*kslnmod, false, *kan, false, false, false, true);
  }

  // kslm: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslmmod;
  if (slipset)
  {
    kslmmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslmmod = Core::LinAlg::MLMultiply(*kslmmod, false, *kam, false, false, false, true);
  }

  // ksli: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslimod;
  if (slipset && iset)
  {
    kslimod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslimod = Core::LinAlg::MLMultiply(*kslimod, false, *kai, false, false, false, true);
  }

  // kslsl: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslslmod;
  if (slipset)
  {
    kslslmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslslmod = Core::LinAlg::MLMultiply(*kslslmod, false, *kasl, false, false, false, true);
  }

  // slstmod: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslstmod;
  if (slipset && stickset)
  {
    kslstmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslstmod = Core::LinAlg::MLMultiply(*kslstmod, false, *kast, false, false, false, true);
  }

  /********************************************************************/
  /* (8) Build the final f blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->multiply(true, *fa, *fmmod);
  fmmod->Update(1.0, *fm, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (aset && iset) dhat->multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //-------------------------------------------------------- FOURTH LINE

  //--------------------------------------------------------- FIFTH LINE
  // split the lagrange multiplier vector in stick and slip part
  Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
  Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
  Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

  Core::LinAlg::split_vector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
  Core::LinAlg::split_vector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
  Teuchos::RCP<Epetra_Vector> tempvec1;

  // fst: mutliply with linstickLM
  Teuchos::RCP<Epetra_Vector> fstmod;
  if (stickset)
  {
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
      fstmod = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
    else
      fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
        Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    temp1->multiply(false, *fa, *fstmod);

    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
      tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
    else
      tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

    linstickLM_->multiply(false, *zst, *tempvec1);
    fstmod->Update(-1.0, *tempvec1, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE
  // fsl: mutliply with linslipLM
  Teuchos::RCP<Epetra_Vector> fslmod;
  Teuchos::RCP<Epetra_Vector> fslwmod;

  if (slipset)
  {
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
      fslmod = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
    else
      fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    temp->multiply(false, *fa, *fslmod);

    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
      tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
    else
      tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

    linslipLM_->multiply(false, *zsl, *tempvec1);

    fslmod->Update(-1.0, *tempvec1, 1.0);
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (parallel_redistribution_status())
  {
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
    kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
    if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
    if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
      kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
      kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
      if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset)
    {
      smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);
    }

    //----------------------------------------------------------- FIFTH LINE
    if (stickset)
    {
      kstnmod = Mortar::MatrixRowTransform(kstnmod, pgsdofrowmap_);
      kstmmod = Mortar::MatrixRowTransform(kstmmod, pgsdofrowmap_);
      if (iset) kstimod = Mortar::MatrixRowTransform(kstimod, pgsdofrowmap_);
      if (slipset) kstslmod = Mortar::MatrixRowTransform(kstslmod, pgsdofrowmap_);
      kststmod = Mortar::MatrixRowTransform(kststmod, pgsdofrowmap_);
      linstickDIS_ = Mortar::MatrixRowTransform(linstickDIS_, pgsdofrowmap_);
    }

    //----------------------------------------------------------- SIXTH LINE
    if (slipset)
    {
      kslnmod = Mortar::MatrixRowTransform(kslnmod, pgsdofrowmap_);
      kslmmod = Mortar::MatrixRowTransform(kslmmod, pgsdofrowmap_);
      if (iset) kslimod = Mortar::MatrixRowTransform(kslimod, pgsdofrowmap_);
      if (stickset) kslstmod = Mortar::MatrixRowTransform(kslstmod, pgsdofrowmap_);
      kslslmod = Mortar::MatrixRowTransform(kslslmod, pgsdofrowmap_);
      linslipDIS_ = Mortar::MatrixRowTransform(linslipDIS_, pgsdofrowmap_);
    }
  }

  /********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                */
  /********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

  //--------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->add(*knn, false, 1.0, 1.0);
  kteffnew->add(*knm, false, 1.0, 1.0);
  if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->add(*kmnmod, false, 1.0, 1.0);
  kteffnew->add(*kmmmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
  if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

  //-------------------------------------------------------- FOURTH LINE

  // add a submatrices to kteffnew
  if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);

  //--------------------------------------------------------- FIFTH LINE
  // add st submatrices to kteffnew
  if (stickset) kteffnew->add(*kstnmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kstmmod, false, 1.0, 1.0);
  if (stickset && iset) kteffnew->add(*kstimod, false, 1.0, 1.0);
  if (stickset && slipset) kteffnew->add(*kstslmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kststmod, false, 1.0, 1.0);

  // add terms of linearization of sick condition to kteffnew
  if (stickset) kteffnew->add(*linstickDIS_, false, -1.0, 1.0);

  //--------------------------------------------------------- SIXTH LINE
  // add sl submatrices to kteffnew
  if (slipset) kteffnew->add(*kslnmod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslmmod, false, 1.0, 1.0);
  if (slipset && iset) kteffnew->add(*kslimod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslslmod, false, 1.0, 1.0);
  if (slipset && stickset) kteffnew->add(*kslstmod, false, 1.0, 1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (slipset) kteffnew->add(*linslipDIS_, false, -1.0, +1.0);

  // add diagonal entries to sparsity pattern for dbc
  for (int i = 0; i < kteffnew->row_map().NumMyElements(); ++i)
  {
    int gid = kteffnew->row_map().GID(i);
    kteffnew->assemble(0., gid, gid);
  }

  // fill_complete kteffnew (square)
  kteffnew->complete();

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fn, *fnexp);
  feffnew->Update(1.0, *fnexp, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fmmod, *fmmodexp);
  feffnew->Update(1.0, *fmmodexp, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //-------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  Teuchos::RCP<Epetra_Vector> fwexp;
  Teuchos::RCP<Epetra_Vector> fgmodexp;

  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*gact, *gexp);
    feffnew->Update(-1.0, *gexp, 1.0);
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fstmodexp;
  if (stickset)
  {
    fstmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fstmod, *fstmodexp);
    feffnew->Update(1.0, *fstmodexp, +1.0);
  }

  // add terms of linearization feffnew
  if (stickset)
  {
    Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linstickRHS_, *linstickRHSexp);
    feffnew->Update(-1.0, *linstickRHSexp, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE

  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fslmodexp;
  Teuchos::RCP<Epetra_Vector> fwslexp;
  Teuchos::RCP<Epetra_Vector> fslwmodexp;


  if (slipset)
  {
    fslmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fslmod, *fslmodexp);
    feffnew->Update(1.0, *fslmodexp, 1.0);
  }

  if (slipset)
  {
    Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linslipRHS_, *linslipRHSexp);
    feffnew->Update(-1.0, *linslipRHSexp, 1.0);
  }

  // finally do the replacement
  *kteff = *kteffnew;
  *feff = *feffnew;

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSLIPINCR
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->fd_check_slip_incr_deriv_txi();
    if (Dim() == 3) interface_[i]->fd_check_slip_incr_deriv_teta();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
    // FD check of stick condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      //      Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv1 = Teuchos::rcp(new
      //      Core::LinAlg::SparseMatrix(*gactivet_,81)); Teuchos::RCP<Core::LinAlg::SparseMatrix>
      //      deriv2 = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_,81));
      //
      //      deriv1->Add(*linstickLM_,false,1.0,1.0);
      //      deriv1->Complete(*gsmdofrowmap_,*gactivet_);
      //
      //      deriv2->Add(*linstickDIS_,false,1.0,1.0);
      //      deriv2->Complete(*gsmdofrowmap_,*gactivet_);
      //
      //      std::cout << "DERIV 1 *********** "<< *deriv1 << std::endl;
      //      std::cout << "DERIV 2 *********** "<< *deriv2 << std::endl;

      interface_[i]->FDCheckStickDeriv(*linstickLM_, *linstickDIS_);
    }
  }
#endif  // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
    // FD check of slip condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      //      Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv1 = Teuchos::rcp(new
      //      Core::LinAlg::SparseMatrix(*gactivet_,81)); Teuchos::RCP<Core::LinAlg::SparseMatrix>
      //      deriv2 = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_,81));
      //
      //      deriv1->Add(*linslipLM_,false,1.0,1.0);
      //      deriv1->Complete(*gsmdofrowmap_,*gslipt_);
      //
      //      deriv2->Add(*linslipDIS_,false,1.0,1.0);
      //      deriv2->Complete(*gsmdofrowmap_,*gslipt_);
      //
      //      std::cout << *deriv1 << std::endl;
      //      std::cout << *deriv2 << std::endl;

      interface_[i]->FDCheckSlipDeriv(*linslipLM_, *linslipDIS_);
    }
  }
#endif  // #ifdef CONTACTFDSLIP

  feff->Scale(-1.);
  return;
}

void CONTACT::LagrangeStrategy::condense_frictionless(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff, Epetra_Vector& rhs)
{
  Teuchos::RCP<Epetra_Vector> feff = Teuchos::rcpFromRef<Epetra_Vector>(rhs);

  feff->Update(-1. + alphaf_, *strcontactrhs_, 1.);
  feff->Scale(-1.);
  Teuchos::RCP<Core::LinAlg::SparseOperator> kteff_op =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseOperator>(kteff);

  // shape function type and type of LM interpolation for quadratic elements
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  // In case of nonsmooth contact the scenario of contacting edges (non parallel)
  // requires a penalty regularization. Here, the penalty contriutions for this
  // special case are applied:
  if (nonSmoothContact_)
  {
    // LTL contributions:
    add_line_to_lin_contributions(kteff_op, feff);

    // penalty support for master side quantities:
    add_master_contributions(kteff_op, feff);

#ifdef CONTACTFDGAPLTL
    // FD check of weighted gap g derivatives (non-penetr. condition)
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckGapDerivLTL();
    }
#endif  // #ifdef CONTACTFDGAPLTL
  }

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
  {
    gact = Core::LinAlg::CreateVector(*gactivedofs_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
    }
  }
  else
  {
    gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gact->GlobalLength())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  // do reagularization scaling
  if (regularized_) evaluate_regularization_scaling(gact);

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin &&
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD") !=
          Inpar::Mortar::lagmult_const)
    FOUR_C_THROW("Condensation only for dual LM");

  /**********************************************************************/
  /* (1) Multiply Mortar matrices: m^ = inv(d) * m                      */
  /**********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invd = invd_;

  /**********************************************************************/
  /* (2) Add contact stiffness terms to kteff                           */
  /**********************************************************************/
  // declare sparse matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);

  /**********************************************************************/
  /* (3) Split kteff into 3x3 matrix blocks                             */
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
  if (parallel_redistribution_status())
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
  /* (4) Split feff into 3 subvectors                                   */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    Core::LinAlg::export_to(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
  }

  // abbreviations for slave set
  const int sset = gsdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

  // store some stuff for static condensation of LM
  fs_ = fs;
  invd_ = invd;
  ksn_ = ksn;
  ksm_ = ksm;
  kss_ = kss;

  /**********************************************************************/
  /* (5) Split slave quantities into active / inactive                  */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  Core::LinAlg::SplitMatrix2x2(kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
  Core::LinAlg::SplitMatrix2x2(
      ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

  // abbreviations for active and inactive set
  const int aset = gactivedofs_->NumGlobalElements();
  const int iset = gidofs->NumGlobalElements();

  // we want to split fsmod into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

  /**********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix               */
  /**********************************************************************/
  // active part of invd
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
  Core::LinAlg::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
  Core::LinAlg::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  // do the multiplication dhat = invda * dai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset && iset) dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
  dhat->complete(*gidofs, *gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
  Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
      tempmtx1, tempmtx2, tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset) mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
  mhataam->complete(*gmdofrowmap_, *gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
  Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
      tempmtx1, tempmtx2, tempmtx3);

  // scaling of invd and dai
  invda->scale(1 / (1 - alphaf_));
  dai->scale(1 - alphaf_);

  save_coupling_matrices(dhat, mhataam, invda);
  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmnmod->add(*kmn, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
  kmnmod->add(*kmnadd, false, 1.0, 1.0);
  kmnmod->complete(kmn->domain_map(), kmn->row_map());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmmmod->add(*kmm, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
  kmmmod->add(*kmmadd, false, 1.0, 1.0);
  kmmmod->complete(kmm->domain_map(), kmm->row_map());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmimod->add(*kmi, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
    kmimod->add(*kmiadd, false, 1.0, 1.0);
    kmimod->complete(kmi->domain_map(), kmi->row_map());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmamod->add(*kma, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
    kmamod->add(*kmaadd, false, 1.0, 1.0);
    kmamod->complete(kma->domain_map(), kma->row_map());
  }

  //----------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------
  // kin: subtract T(dhat)*kan --
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kinmod->add(*kin, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
    kinmod->add(*kinadd, false, -1.0, 1.0);
  }
  kinmod->complete(kin->domain_map(), kin->row_map());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kimmod->add(*kim, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
    kimmod->add(*kimadd, false, -1.0, 1.0);
  }
  kimmod->complete(kim->domain_map(), kim->row_map());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiimod->add(*kii, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
      kiimod->add(*kiiadd, false, -1.0, 1.0);
    }
    kiimod->complete(kii->domain_map(), kii->row_map());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiamod->add(*kia, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
    kiamod->add(*kiaadd, false, -1.0, 1.0);
    kiamod->complete(kia->domain_map(), kia->row_map());
  }

  //---------------------------------------------------------- FOURTH LINE
  // nothing to do - execpt its regularized contact see --> do_regularization_scaling()

  //----------------------------------------------------------- FIFTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kanmod;
  if (aset)
  {
    kanmod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
    kanmod = Core::LinAlg::MLMultiply(*kanmod, false, *kan, false, false, false, true);
  }

  // kam: multiply tmatrix with invda and kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kammod;
  if (aset)
  {
    kammod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
    kammod = Core::LinAlg::MLMultiply(*kammod, false, *kam, false, false, false, true);
  }

  // kai: multiply tmatrix with invda and kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaimod;
  if (aset && iset)
  {
    kaimod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
    kaimod = Core::LinAlg::MLMultiply(*kaimod, false, *kai, false, false, false, true);
  }

  // kaa: multiply tmatrix with invda and kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaamod;
  if (aset)
  {
    kaamod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
    kaamod = Core::LinAlg::MLMultiply(*kaamod, false, *kaa, false, false, false, true);
  }

  /**********************************************************************/
  /* (8) Build the final f blocks                                       */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->multiply(true, *fa, *fmmod);
  fmmod->Update(1.0, *fm, 1.0);

  //----------------------------------------------------------- THIRD LINE
  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (iset && aset) dhat->multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //---------------------------------------------------------- FOURTH LINE
  // gactive: nothing to do  - execpt its regularized contact see --> do_regularization_scaling()

  //----------------------------------------------------------- FIFTH LINE
  // fa: multiply tmatrix with invda and fa
  Teuchos::RCP<Epetra_Vector> famod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tinvda;
  if (aset)
  {
    if (constr_direction_ == Inpar::CONTACT::constr_xyz)
      famod = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    else
      famod = Teuchos::rcp(new Epetra_Vector(*gactivet_));

    tinvda = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda, true, false, false, true);
    tinvda->multiply(false, *fa, *famod);
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (parallel_redistribution_status())
  {
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
    kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
    if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
    if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
      kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
      kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
      if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset) smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);

    //----------------------------------------------------------- FIFTH LINE
    if (aset)
    {
      kanmod = Mortar::MatrixRowTransform(kanmod, pgsdofrowmap_);
      kammod = Mortar::MatrixRowTransform(kammod, pgsdofrowmap_);
      kaamod = Mortar::MatrixRowTransform(kaamod, pgsdofrowmap_);
      if (iset) kaimod = Mortar::MatrixRowTransform(kaimod, pgsdofrowmap_);
      tderivmatrix_ = Mortar::MatrixRowTransform(tderivmatrix_, pgsdofrowmap_);
    }
  }

  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

  // Add all additional contributions from regularized contact to kteffnew and feffnew!
  if (regularized_)
    do_regularization_scaling(aset, iset, invda, kan, kam, kai, kaa, fa, kteffnew, feffnew);

  //----------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->add(*knn, false, 1.0, 1.0);
  kteffnew->add(*knm, false, 1.0, 1.0);
  if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->add(*kmnmod, false, 1.0, 1.0);
  kteffnew->add(*kmmmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

  //----------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
  if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

  //---------------------------------------------------------- FOURTH LINE
  // add a submatrices to kteffnew
  if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);

  //----------------------------------------------------------- FIFTH LINE
  // add a submatrices to kteffnew
  if (aset) kteffnew->add(*kanmod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kammod, false, 1.0, 1.0);
  if (aset && iset) kteffnew->add(*kaimod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kaamod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*tderivmatrix_, false, -1.0, 1.0);

  // fill_complete kteffnew (square)
  kteffnew->complete();

  /**********************************************************************/
  /* (11) Global setup of feffnew (including contact)                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fn, *fnexp);
  feffnew->Update(1.0, *fnexp, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fmmod, *fmmodexp);
  feffnew->Update(1.0, *fmmodexp, 1.0);

  //----------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //---------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*gact, *gexp);
    feffnew->Update(-1.0, *gexp, 1.0);
  }

  //----------------------------------------------------------- FIFTH LINE
  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*famod, *famodexp);
    feffnew->Update(1.0, *famodexp, 1.0);
  }

  // finally do the replacement
  *kteff = *kteffnew;
  *feff = *feffnew;


#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif  // #ifdef CONTACTFDGAP
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    std::cout << *tderivmatrix_ << std::endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif  // #ifdef CONTACTFDTANGLM


  feff->Scale(-1.);
}

void CONTACT::LagrangeStrategy::run_pre_apply_jacobian_inverse(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff, Epetra_Vector& rhs)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if (is_dual_quad_slave_trafo() && lagmultquad == Inpar::Mortar::lagmult_lin)
  {
    if (not(systrafo_->domain_map().SameAs(kteff->domain_map()))) FOUR_C_THROW("stop");

    *kteff = *Core::LinAlg::Multiply(*systrafo_, true, *kteff, false, false, false, true);
    Epetra_Vector rhs_str(*problem_dofs());
    Epetra_Vector rhs_str2(*problem_dofs());
    Core::LinAlg::export_to(rhs, rhs_str);
    if (systrafo_->multiply(true, rhs_str, rhs_str2)) FOUR_C_THROW("multiply failed");
    for (int i = 0; i < rhs_str2.Map().NumMyElements(); ++i)
      rhs[rhs.Map().LID(rhs_str2.Map().GID(i))] = rhs_str2[i];
  }


  if (systype_ != Inpar::CONTACT::system_condensed) return;

  if (friction_)
    condense_friction(kteff, rhs);
  else
    condense_frictionless(kteff, rhs);

  return;
}

void CONTACT::LagrangeStrategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::Nln::Group& grp)
{
  if (system_type() != Inpar::CONTACT::system_condensed) return;

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcpFromRef<Epetra_Vector>(result);
  disi->Scale(-1.);
  {
    // shape function type and type of LM interpolation for quadratic elements
    Inpar::Mortar::ShapeFcn shapefcn =
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

    // double-check if this is a dual LM system
    if ((shapefcn != Inpar::Mortar::shape_dual &&
            shapefcn != Inpar::Mortar::shape_petrovgalerkin) &&
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD") !=
            Inpar::Mortar::lagmult_const)
      FOUR_C_THROW("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, tempmap, gactivedofs_, tempmap, invda, tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invdmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    invdmod->add(*invda, false, 1.0, 1.0);
    invdmod->complete();

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/
    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold * zold and Mold^T * zold to fit
    if (is_self_contact())
    {
      // full update
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }
    else
    {
      // full update
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }

    // store updated LM into nodes
    store_nodal_quantities(Mortar::StrategyBase::lmupdate);
  }
  disi->Scale(-1.);
}

FOUR_C_NAMESPACE_CLOSE
