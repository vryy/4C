/*---------------------------------------------------------------------*/
/*! \file
\brief Contact Strategy handling the porous no penetraction condition on the active contact
interface

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_lagrange_strategy_poro.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_utils.hpp"

#include <Epetra_SerialComm.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                              ager 08/14|
 *----------------------------------------------------------------------*/
CONTACT::LagrangeStrategyPoro::LagrangeStrategyPoro(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
    Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof, bool poroslave, bool poromaster)
    : MonoCoupledLagrangeStrategy(
          data_ptr, dof_row_map, NodeRowMap, params, interface, dim, comm, alphaf, maxdof),
      no_penetration_(params.get<bool>("CONTACTNOPEN")),
      nopenalpha_(0.0),
      poroslave_(poroslave),
      poromaster_(poromaster)
{
  if (!poroslave_ and !poromaster_)
    FOUR_C_THROW(
        "you called a poroelastic meshtying method without participating poroelastic domains on "
        "your interface");
  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      ager 12/16|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::do_read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  Teuchos::RCP<Core::FE::Discretization> discret =
      Global::Problem::instance()->get_dis("structure");
  if (discret == Teuchos::null) FOUR_C_THROW("didn't get my discretization");

  Teuchos::RCP<Epetra_Vector> global =
      Teuchos::rcp(new Epetra_Vector(*discret->dof_col_map(), true));
  // it's clear that we get some zeros here ... but poroelast monolithic fixes this a little bit
  // later by doing the same thing with correct displacements again :-)
  Core::LinAlg::Export(*dis, *global);
  set_parent_state("displacement", global, discret);

  // Call (nearly absolute)Base Class
  CONTACT::AbstractStrategy::do_read_restart(reader, dis, cparams_ptr);
  return;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                ager 12/16 |
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::setup(bool redistributed, bool init)
{
  // Call Base Class
  CONTACT::AbstractStrategy::setup(redistributed, init);

  if (no_penetration_) setup_no_penetration_condition();
}

/*----------------------------------------------------------------------*
 | Activate No-Penetration for the contact surface (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::setup_no_penetration_condition()
{
  lambda_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  if (!poroslave_ and poromaster_)
    FOUR_C_THROW("poroelastic meshtying/contact method needs the slave side to be poroelastic");
  /*
   *  The first error may also occur when there is no single element on the interface (it is empty)
   * but POROELASTICITY DYNAMIC coupling algorithmus (COUPALGO) is chosen as
   * poro_monolithicmeshtying that means that the method creates an interface but has nothing to
   * fill in
   *
   *  the second error shows that there are no elements with PoroCoupling on the slave side of the
   * interface and the master side is fully poroelastic
   *
   *  having the condition with structure slave and poro master the problem is that there is no
   * diagonal D Matrix, from the porofluid lagrange multiplier equality condition, that can be
   * inverted easily to condense out the respective lagrange multipliers
   *
   *  there are two alternatives:
   *  1. do not condense out the lagrange multiplier for the porofluid meshtying condition in
   *  this constellation (master/slave) and solve the saddlepoint system
   *
   *  2. invert the M Matrix to condense out the Lagrange multiplier
   *  - which is only applicable with adequate costs for small problems as the M Matrix is not at
   * all diagonal, except for matching meshes - for matching meshes it doesn't hurt to chose the
   * poro side as the slave side. Being applicable only for small problems this alternative is too
   * restricted in its application overall to be implemented.
   */
}

/*----------------------------------------------------------------------*
 | initialize global poro contact variables                             |
 |                            for next Newton step (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::poro_initialize(
    Core::Adapter::Coupling& coupfs, Teuchos::RCP<const Epetra_Map> fluiddofs, bool fullinit)
{
  if (fullinit)  // fullinit is true by default, but needed when this method is called for
                 // meshtying, as the maps and matrix mapping stay the same for meshtying. would
                 // work without this but would do things repeatedly unnecessarily
  {
    if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
    {
      //  (1)                                                          //
      //      Get required fluid maps from structural maps             //
      //                                                               //
      fgsdofrowmap_ = coupfs.master_to_slave_map(gsdofrowmap_);
      fgmdofrowmap_ = coupfs.master_to_slave_map(gmdofrowmap_);
      fgsmdofrowmap_ = coupfs.master_to_slave_map(gsmdofrowmap_);
      fgndofrowmap_ = Core::LinAlg::SplitMap(*fluiddofs,
          *fgsmdofrowmap_);  // Not equal to transforming gndofrowmap_ (pressure dofs missing!)
      fgactivedofs_ = coupfs.master_to_slave_map(gactivedofs_);
      falldofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(*fluiddofs));
      fgactiven_ = coupfs.master_to_slave_map(gactiven_);
      fgactivet_ = coupfs.master_to_slave_map(gactivet_);
    }
  }
  //  (2)                                                          //
  //      Initialize Matrices                                      //
  //                                                               //
  if (fullinit)
  {
    if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
    {
      // (re)setup global nCoup Vector
      NCoup_ = Teuchos::rcp(new Epetra_Vector(*gactiven_));

      // (re)setup global linearisation matrices of nCoup
      NCoup_lindisp_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 10));
      NCoup_linvel_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 10));

      // (re)setup global tangential and lin(tangential)*lambda matrices
      Tangential_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_, 10));
      linTangentiallambda_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_, 10));

      // (re)setup global lin of D & M * lambda - Matrix
      porolindmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gsdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      porolinmmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    }
  }
  else
  {
    // In the case of meshtying resetting the matrices is sufficient as they retain their size
    NCoup_->PutScalar(0.0);

    // (re)setup global linearisation matrices of nCoup
    NCoup_lindisp_->zero();
    NCoup_linvel_->zero();

    // (re)setup global tangential and lin(tangential)*lambda matrices
    Tangential_->zero();
    linTangentiallambda_->zero();

    // (re)setup global lin of D & M * lambda - Matrix
    porolindmatrix_->zero();
    porolinmmatrix_->zero();
  }
  //  (3)                                                          //
  //      Assemble Matrices                                        //
  //                                                               //
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
    {
      interface_[i]->assemble_n_coup(*NCoup_);

      interface_[i]->assemble_n_coup_lin(*NCoup_lindisp_, coupfs);
      interface_[i]->assemble_n_coup_lin(*NCoup_linvel_, coupfs, true);

      interface_[i]->assemble_tn(Tangential_, Teuchos::null);
      interface_[i]->assemble_t_nderiv(linTangentiallambda_, Teuchos::null,
          true);  // use lambda(n +1) for tangential condition!!!

      interface_[i]->assemble_lin_dm(*porolindmatrix_, *porolinmmatrix_, true);
    }
  }

  //  (4)                                                          //
  //      Complete Matrices                                        //
  //                                                               //
  if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
  {
    NCoup_lindisp_->complete(*problem_dofs(), *gactiven_);
    NCoup_linvel_->complete(*fluiddofs, *gactiven_);

    Tangential_->complete(*gactivedofs_, *gactivet_);
    linTangentiallambda_->complete(*gsmdofrowmap_, *gactivet_);

    porolindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
    porolinmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);
  }

  //  (5)                                                          //
  //      Reset Matrix Transform Objects                           //
  //                                                               //
  if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
  {
    linncoupveltransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    linncoupdisptransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    tanginvtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
    lintangentlambdatransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    porolindmatrixtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    porolinmmatrixtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    mhataamtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);  // h.Willmann
    dhattransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    doldtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    moldtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
    invDatransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  }

  //  (6)                                                          //
  //      Transform Matrices from structural dofs to fluid dofs    //
  //                                                               //
  if (no_penetration_ && (is_in_contact() || was_in_contact() || was_in_contact_last_time_step()))
  // eventually a coupling object just on the mortar interface would make sense!!! ChrAg
  {
    // transform matrices coming from contact to fluid maps, as they are all in structure maps!
    //
    // A generell problem here is that we would need to update coupling objects in everey newton
    // step if the active set changes. To avoid this, a 'bigger' coupling object is used - but
    // therefore now the Row-Maps of the created Sparse Matrixes are to big! --- Leads to problems
    // for Matrix - Multiplications where the Row - Map is used!
    //
    // At the moment this is solved by using SplitMatrix2x2 to cut out just the relevant part of the
    // matrix, but this shouldn't be the final solution!
    //
    //************************************************************************************************
    //
    Teuchos::RCP<Epetra_Vector> tmpfullncoup =
        Teuchos::rcp(new Epetra_Vector(*coupfs.master_dof_map()));
    Core::LinAlg::Export(*NCoup_, *tmpfullncoup);
    tmpfullncoup = coupfs.master_to_slave(tmpfullncoup);
    fNCoup_ = Teuchos::rcp(new Epetra_Vector(*fgactiven_));
    Core::LinAlg::Export(*tmpfullncoup, *fNCoup_);
    //
    //************************************************************************************************
    //
    fdoldtransp_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
        new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*doldtransform_)(*dold_->transpose(), 1.0, Core::Adapter::CouplingMasterConverter(coupfs),
        *fdoldtransp_, false);
    fdoldtransp_->complete(dold_->row_map(), *fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    if (poromaster_)
    {
      fmoldtransp_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
          new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*moldtransform_)(*mold_->transpose(), 1.0, Core::Adapter::CouplingMasterConverter(coupfs),
          *fmoldtransp_, false);
      fmoldtransp_->complete(mold_->row_map(), *fgmdofrowmap_);
    }
    //
    //************************************************************************************************
    //
    fporolindmatrix_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
        new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*porolindmatrixtransform_)(*porolindmatrix_, 1.0,
        Core::Adapter::CouplingMasterConverter(coupfs), *fporolindmatrix_, false);
    fporolindmatrix_->complete(porolindmatrix_->domain_map(), *fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    if (poromaster_)
    {
      fporolinmmatrix_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
          new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*porolinmmatrixtransform_)(*porolinmmatrix_, 1.0,
          Core::Adapter::CouplingMasterConverter(coupfs), *fporolinmmatrix_, false);
      fporolinmmatrix_->complete(porolinmmatrix_->domain_map(), *fgmdofrowmap_);
    }
    // porolinmmatrixtransform_ is no longer missing as twosided poro meshtying is considered!
    //
    //************************************************************************************************
    if (poromaster_)
    //
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpfmhataam =
          Teuchos::rcp<Core::LinAlg::SparseMatrix>(
              new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));

      fmhataam_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
          new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 1, true, false));
      (*mhataamtransform_)(*mhataam_, 1.0, Core::Adapter::CouplingMasterConverter(coupfs),
          Core::Adapter::CouplingMasterConverter(coupfs), *tmpfmhataam, false, false);

      tmpfmhataam->complete(*fgmdofrowmap_, *falldofrowmap_);

      // better solution to get maps as wanted? -- for this matrix map as important as there will be
      // a matrix-matrix multiplication

      Teuchos::RCP<Epetra_Map> restfgmdofrowmap, restfgactivedofs;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpm1, tmpm2, tmpm3;

      // This should just be a temporary solution to change the row map of the matrix ...
      Core::LinAlg::SplitMatrix2x2(tmpfmhataam, fgactivedofs_, restfgactivedofs, fgmdofrowmap_,
          restfgmdofrowmap, fmhataam_, tmpm1, tmpm2, tmpm3);
      fmhataam_->complete(*fgmdofrowmap_, *fgactivedofs_);
    }
    //
    //************************************************************************************************
    //
    fdhat_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
        new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*dhattransform_)(*dhat_, 1.0, Core::Adapter::CouplingMasterConverter(coupfs), *fdhat_, false);
    fdhat_->complete(dhat_->domain_map(), *fgactivedofs_);
    // fdhat is expected to be zero in this method but still used for condensation
    //
    //************************************************************************************************
    //
    if (gactivedofs_->NumGlobalElements())
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tanginvD =
          Core::LinAlg::MLMultiply(*Tangential_, false, *invda_, true, false, false, true);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpftanginvD =
          Teuchos::rcp<Core::LinAlg::SparseMatrix>(
              new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*tanginvtransform_)(*tanginvD, 1.0, Core::Adapter::CouplingMasterConverter(coupfs),
          Core::Adapter::CouplingMasterConverter(coupfs), *tmpftanginvD, false, false);
      tmpftanginvD->complete(*fgactivedofs_, *falldofrowmap_);
      // better solution to get maps as wanted? -- for this matrix map as important as there will be
      // a matrix-matrix multiplication

      Teuchos::RCP<Epetra_Map> restfgactivet, restfgactivedofs;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpm1, tmpm2, tmpm3;

      // This should just be a temporary solution to change the row map of the matrix ...
      Core::LinAlg::SplitMatrix2x2(tmpftanginvD, fgactivet_, restfgactivet, fgactivedofs_,
          restfgactivedofs, ftanginvD_, tmpm1, tmpm2, tmpm3);

      //
      //************************************************************************************************
      //
      fNCoup_linvel_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*falldofrowmap_, 108, false));

      (*linncoupveltransform_)(*NCoup_linvel_, 1.0, Core::Adapter::CouplingMasterConverter(coupfs),
          *fNCoup_linvel_, false);
      fNCoup_linvel_->complete(*falldofrowmap_, *fgactiven_);
      //
      //************************************************************************************************
      //
      fNCoup_lindisp_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*falldofrowmap_, 81, false));
      (*linncoupdisptransform_)(*NCoup_lindisp_, 1.0,
          Core::Adapter::CouplingMasterConverter(coupfs), *fNCoup_lindisp_, false);
      fNCoup_lindisp_->complete(*problem_dofs(), *fgactiven_);
      //
      //************************************************************************************************
      //
      flinTangentiallambda_ =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*falldofrowmap_, 81, false));
      (*lintangentlambdatransform_)(*linTangentiallambda_, 1.0,
          Core::Adapter::CouplingMasterConverter(coupfs), *flinTangentiallambda_, false);
      flinTangentiallambda_->complete(*gsdofrowmap_, *fgactivet_);

      //
      //************************************************************************************************
      //
      finvda_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(
          new Core::LinAlg::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*invDatransform_)(
          *invda_, 1.0, Core::Adapter::CouplingMasterConverter(coupfs), *finvda_, false);
      finvda_->complete(invda_->domain_map(), *fgactivedofs_);
    }
  }
  return;
}  // CONTACT::LagrangeStrategyPoro::PoroInitialize

/*-------------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration             |
 |  condition on contact surface (pure porous problem)(public)   ager 07/15|
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::evaluate_poro_no_pen_contact(
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& k_fseff,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& Feff, Teuchos::RCP<Epetra_Vector>& feff)
{
  evaluate_mat_poro_no_pen(k_fseff, feff);

  evaluate_other_mat_poro_no_pen(Feff, 0);
}

/*-------------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration             |
 |  condition on contact surface(public)                         ager 07/15|
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::evaluate_poro_no_pen_contact(
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& k_fseff,
    std::map<int, Teuchos::RCP<Core::LinAlg::SparseMatrix>*>& Feff,
    Teuchos::RCP<Epetra_Vector>& feff)
{
  evaluate_mat_poro_no_pen(k_fseff, feff);

  // Take care of the alternative condensation of the off-diagonal blocks!!!
  std::map<int, Teuchos::RCP<Core::LinAlg::SparseMatrix>*>::iterator matiter;
  for (matiter = Feff.begin(); matiter != Feff.end(); ++matiter)
  {
    evaluate_other_mat_poro_no_pen(*(matiter->second), matiter->first);
  }
}

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration          |
 |  condition on contact surface (public)                     ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::evaluate_mat_poro_no_pen(
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& k_fseff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ ||
      (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()))
    return;
  // h.Willmann this method should be renamed as it handles twosided poro meshtying now aswell and
  // is able to handle fluid coupling for twosided contact

  nopenalpha_ = alphaf_;  // to use different alpha for nopen condition (not used at the moment)

  // shape function
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Condensation only for dual LM");

  /**********************************************************************/
  /* (2) Add contact stiffness terms to kteff                           */
  /**********************************************************************/

  // transform if necessary
  if (parallel_redistribution_status())
  {
    FOUR_C_THROW("CHECK ME!!!");
    lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
    linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
  }

  k_fseff->un_complete();
  k_fseff->add(*fporolindmatrix_, false, (1.0 - nopenalpha_) * 1.0, 1.0);

  if (poromaster_)
    k_fseff->add(*fporolinmmatrix_, false, (1.0 - nopenalpha_) * 1.0,
        1.0);  // is needed only for twosided poro contact or meshtying

  k_fseff->complete(*problem_dofs(), *falldofrowmap_);  // gets bigger because of linearisation
                                                        // w.r.t. to pure structural displacements!
  /**********************************************************************/
  /* (3) Split k_fseff and Feff into 3x3 matrix blocks                             */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_ss, k_fs_sm, k_fs_sn, k_fs_ms, k_fs_mm, k_fs_mn,
      k_fs_ns, k_fs_nm, k_fs_nn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_smsm, k_fs_smn, k_fs_nsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<Epetra_Map> ftempmap1, ftempmap2, ftempmap3;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;

  // split into slave/master part + structure part
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fseffmatrix =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(k_fseff);
  if (parallel_redistribution_status())
  {
    FOUR_C_THROW("CHECK ME!");
    // split and transform to redistributed maps
    //        Core::LinAlg::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    //        ksmsm = Mortar::matrix_row_col_transform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    //        ksmn  = Mortar::MatrixRowTransform(ksmn,gsmdofrowmap_);
    //        knsm  = Mortar::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::SplitMatrix2x2(k_fseffmatrix, fgsmdofrowmap_, fgndofrowmap_, gsmdofrowmap_,
        gndofrowmap_, k_fs_smsm, k_fs_smn, k_fs_nsm, k_fs_nn);
  }

  // further splits into slave part + master part
  Core::LinAlg::SplitMatrix2x2(k_fs_smsm, fgsdofrowmap_, fgmdofrowmap_, gsdofrowmap_, gmdofrowmap_,
      k_fs_ss, k_fs_sm, k_fs_ms, k_fs_mm);
  Core::LinAlg::SplitMatrix2x2(k_fs_smn, fgsdofrowmap_, fgmdofrowmap_, gndofrowmap_, tempmap,
      k_fs_sn, tempmtx1, k_fs_mn, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(k_fs_nsm, fgndofrowmap_, ftempmap1, gsdofrowmap_, gmdofrowmap_,
      k_fs_ns, k_fs_nm, tempmtx1, tempmtx2);

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
    FOUR_C_THROW("CHECK ME!");
    // split and transform to redistributed maps
    Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    Core::LinAlg::Export(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::split_vector(*falldofrowmap_, *feff, fgsmdofrowmap_, fsm, fgndofrowmap_, fn);
  }

  // abbreviations for slave  and master set
  int sset = fgsdofrowmap_->NumGlobalElements();
  int mset = fgmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*fgsmdofrowmap_, *fsm, fgsdofrowmap_, fs, fgmdofrowmap_, fm);

  // store some stuff for static condensation of poro no pen. LM

  ffs_ = fs;

  cfssn_ = k_fs_sn;
  cfssm_ = k_fs_sm;
  cfsss_ = k_fs_ss;

  /**********************************************************************/
  /* (5) Split slave quantities into active / inactive                  */
  /**********************************************************************/

  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_aa, k_fs_ai, k_fs_ia, k_fs_ii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_an, k_fs_in, k_fs_am, k_fs_im, k_fs_ma, k_fs_mi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;
  Teuchos::RCP<Epetra_Map> fgidofs;

  // some more temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap1, tempmap2;
  Teuchos::RCP<Epetra_Map> ftempmap4, ftempmap5, ftempmap6, ftempmap7;

  // do the splitting
  Core::LinAlg::SplitMatrix2x2(
      k_fs_ss, fgactivedofs_, fgidofs, gactivedofs_, gidofs, k_fs_aa, k_fs_ai, k_fs_ia, k_fs_ii);
  Core::LinAlg::SplitMatrix2x2(k_fs_sn, fgactivedofs_, fgidofs, gndofrowmap_, tempmap1, k_fs_an,
      tempmtx1, k_fs_in, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(k_fs_sm, fgactivedofs_, fgidofs, gmdofrowmap_, tempmap2, k_fs_am,
      tempmtx1, k_fs_im, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(k_fs_ms, fgmdofrowmap_, ftempmap4, gactivedofs_, gidofs, k_fs_ma,
      k_fs_mi, tempmtx1, tempmtx2);

  // abbreviations for active and inactive set
  int aset = fgactivedofs_->NumGlobalElements();
  int iset = fgidofs->NumGlobalElements();

  // we want to split fsmod into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*fgidofs));

  // do the vector splitting s -> a+i
  Core::LinAlg::split_vector(*fgsdofrowmap_, *fs, fgactivedofs_, fa, fgidofs, fi);

  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /* where K stands for k_fs and F!!!                                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do
  //---------------------------------------------------------- SECOND LINE --- Will just exist when
  // starting with two sided poro contact!!!
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mnmod;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mmmod;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mimod;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mamod;
  if (mset)
  {
    // kmn: add T(mhataam)*kan
    k_fs_mnmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 100));
    k_fs_mnmod->add(*k_fs_mn, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mnadd =
        Core::LinAlg::MLMultiply(*fmhataam_, true, *k_fs_an, false, false, false, true);
    k_fs_mnmod->add(*k_fs_mnadd, false, 1.0, 1.0);
    k_fs_mnmod->complete(k_fs_mn->domain_map(), k_fs_mn->row_map());

    // kmm: add T(mhataam)*kam
    k_fs_mmmod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 100));
    k_fs_mmmod->add(*k_fs_mm, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_mmadd =
        Core::LinAlg::MLMultiply(*fmhataam_, true, *k_fs_am, false, false, false, true);
    k_fs_mmmod->add(*k_fs_mmadd, false, 1.0, 1.0);
    k_fs_mmmod->complete(k_fs_mm->domain_map(), k_fs_mm->row_map());

    // kmi: add T(mhataam)*kai
    if (iset)
    {
      k_fs_mimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 100));
      k_fs_mimod->add(*k_fs_mi, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_miadd =
          Core::LinAlg::MLMultiply(*fmhataam_, true, *k_fs_ai, false, false, false, true);
      k_fs_mimod->add(*k_fs_miadd, false, 1.0, 1.0);
      k_fs_mimod->complete(k_fs_mi->domain_map(), k_fs_mi->row_map());
    }

    // kma: add T(mhataam)*kaa
    if (aset)
    {
      k_fs_mamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 100));
      k_fs_mamod->add(*k_fs_ma, false, 1.0, 1.0);
      Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_maadd =
          Core::LinAlg::MLMultiply(*fmhataam_, true, *k_fs_aa, false, false, false, true);
      k_fs_mamod->add(*k_fs_maadd, false, 1.0, 1.0);
      k_fs_mamod->complete(k_fs_ma->domain_map(), k_fs_ma->row_map());
    }
  }

  //----------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------
  // fdhat is expected to be zero here but still used for condensation
  // kin: subtract T(dhat)*kan --
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_inmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgidofs, 100));
  k_fs_inmod->add(*k_fs_in, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_inadd =
      Core::LinAlg::MLMultiply(*fdhat_, true, *k_fs_an, false, false, false, true);
  k_fs_inmod->add(*k_fs_inadd, false, -1.0, 1.0);
  k_fs_inmod->complete(k_fs_in->domain_map(), k_fs_in->row_map());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_immod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgidofs, 100));
  k_fs_immod->add(*k_fs_im, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_imadd =
      Core::LinAlg::MLMultiply(*fdhat_, true, *k_fs_am, false, false, false, true);
  k_fs_immod->add(*k_fs_imadd, false, -1.0, 1.0);
  k_fs_immod->complete(k_fs_im->domain_map(), k_fs_im->row_map());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_iimod;
  if (iset)
  {
    k_fs_iimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgidofs, 100));
    k_fs_iimod->add(*k_fs_ii, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_iiadd =
        Core::LinAlg::MLMultiply(*fdhat_, true, *k_fs_ai, false, false, false, true);
    k_fs_iimod->add(*k_fs_iiadd, false, -1.0, 1.0);
    k_fs_iimod->complete(k_fs_ii->domain_map(), k_fs_ii->row_map());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_iamod;
  if (iset && aset)
  {
    k_fs_iamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgidofs, 100));
    k_fs_iamod->add(*k_fs_ia, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_iaadd =
        Core::LinAlg::MLMultiply(*fdhat_, true, *k_fs_aa, false, false, false, true);
    k_fs_iamod->add(*k_fs_iaadd, false, -1.0, 1.0);
    k_fs_iamod->complete(k_fs_ia->domain_map(), k_fs_ia->row_map());
  }

  //---------------------------------------------------------- FOURTH LINE
  // nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_anmod;
  if (aset)
    k_fs_anmod = Core::LinAlg::MLMultiply(*ftanginvD_, false, *k_fs_an, false, false, false, true);

  // kam: multiply tmatrix with invda and kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_ammod;
  if (aset)
    k_fs_ammod = Core::LinAlg::MLMultiply(*ftanginvD_, false, *k_fs_am, false, false, false, true);

  // kai: multiply tmatrix with invda and kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_aimod;
  if (aset && iset)
    k_fs_aimod = Core::LinAlg::MLMultiply(*ftanginvD_, false, *k_fs_ai, false, false, false, true);

  // kaa: multiply tmatrix with invda and kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_aamod;
  if (aset)
    k_fs_aamod = Core::LinAlg::MLMultiply(*ftanginvD_, false, *k_fs_aa, false, false, false, true);

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
    FOUR_C_THROW("CHECK ME!");
    //        Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    //        Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new
    //        Epetra_Vector(mold_->DomainMap())); Teuchos::RCP<Epetra_Vector> zoldexp  =
    //        Teuchos::rcp(new Epetra_Vector(mold_->RowMap())); if
    //        (mold_->RowMap().NumGlobalElements()) Core::LinAlg::Export(*zold_,*zoldexp);
    //        mold_->Multiply(true,*zoldexp,*tempvecm2);
    //        if (mset) Core::LinAlg::Export(*tempvecm2,*tempvecm);
    //        fm->Update(alphaf_,*tempvecm,1.0);
  }
  // if there is no self contact everything is ok
  else if (poromaster_)
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    fmoldtransp_->multiply(false, *lambdaold_, *tempvecm);
    fm->Update(nopenalpha_, *tempvecm, 1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (is_self_contact())
  {
    FOUR_C_THROW("CHECK ME!");
    //        Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new
    //        Epetra_Vector(dold_->DomainMap())); Teuchos::RCP<Epetra_Vector> zoldexp  =
    //        Teuchos::rcp(new Epetra_Vector(dold_->RowMap())); if
    //        (dold_->RowMap().NumGlobalElements()) Core::LinAlg::Export(*zold_,*zoldexp);
    //        dold_->Multiply(true,*zoldexp,*tempvec);
    //        if (sset) Core::LinAlg::Export(*tempvec,*fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    fdoldtransp_->multiply(false, *lambdaold_, *fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
    Core::LinAlg::Export(*fsadd, *faadd);

    fa->Update(-nopenalpha_, *faadd, 1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod;
  if (mset)
  {
    fmmod = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    if (aset) fmhataam_->multiply(true, *fa, *fmmod);
    fmmod->Update(1.0, *fm, 1.0);
  }

  //----------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*fgidofs));
    Core::LinAlg::Export(*fsadd, *fiadd);
    fi->Update(-nopenalpha_, *fiadd, 1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*fgidofs));
  if (aset) fdhat_->multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //---------------------------------------------------------- FOURTH LINE
  // gactive: nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // fa: multiply tmatrix with invda and fa
  Teuchos::RCP<Epetra_Vector> famod;
  if (aset)
  {
    famod = Teuchos::rcp(new Epetra_Vector(*fgactivet_));
    ftanginvD_->multiply(false, *fa, *famod);
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
    FOUR_C_THROW("CHECK ME!");
    //        //----------------------------------------------------------- FIRST LINE
    //        // nothing to do (ndof-map independent of redistribution)
    //
    //        //---------------------------------------------------------- SECOND LINE
    //        kmnmod = Mortar::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    //        kmmmod = Mortar::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    //        if (iset) kmimod = Mortar::MatrixRowTransform(kmimod,pgmdofrowmap_);
    //        if (aset) kmamod = Mortar::MatrixRowTransform(kmamod,pgmdofrowmap_);
    //
    //        //----------------------------------------------------------- THIRD LINE
    //        if (iset)
    //        {
    //          kinmod = Mortar::MatrixRowTransform(kinmod,pgsdofrowmap_);
    //          kimmod = Mortar::MatrixRowTransform(kimmod,pgsdofrowmap_);
    //          kiimod = Mortar::MatrixRowTransform(kiimod,pgsdofrowmap_);
    //          if (aset) kiamod = Mortar::MatrixRowTransform(kiamod,pgsdofrowmap_);
    //        }
    //
    //        //---------------------------------------------------------- FOURTH LINE
    //        if (aset) smatrix_ = Mortar::MatrixRowTransform(smatrix_,pgsdofrowmap_);
    //
    //        //----------------------------------------------------------- FIFTH LINE
    //        if (aset)
    //        {
    //          kanmod = Mortar::MatrixRowTransform(kanmod,pgsdofrowmap_);
    //          kammod = Mortar::MatrixRowTransform(kammod,pgsdofrowmap_);
    //          kaamod = Mortar::MatrixRowTransform(kaamod,pgsdofrowmap_);
    //          if (iset) kaimod = Mortar::MatrixRowTransform(kaimod,pgsdofrowmap_);
    //          tderivmatrix_ = Mortar::MatrixRowTransform(tderivmatrix_,pgsdofrowmap_);
    //        }
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs_effnew = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*falldofrowmap_, 81, true, false, k_fseff->get_matrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*falldofrowmap_);

  //----------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  k_fs_effnew->add(*k_fs_nn, false, 1.0, 1.0);
  k_fs_effnew->add(*k_fs_nm, false, 1.0, 1.0);
  if (sset) k_fs_effnew->add(*k_fs_ns, false, 1.0, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  if (mset)
  {
    k_fs_effnew->add(*k_fs_mnmod, false, 1.0, 1.0);
    k_fs_effnew->add(*k_fs_mmmod, false, 1.0, 1.0);
    if (iset) k_fs_effnew->add(*k_fs_mimod, false, 1.0, 1.0);
    if (aset) k_fs_effnew->add(*k_fs_mamod, false, 1.0, 1.0);
  }

  //----------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) k_fs_effnew->add(*k_fs_inmod, false, 1.0, 1.0);
  if (iset) k_fs_effnew->add(*k_fs_immod, false, 1.0, 1.0);
  if (iset) k_fs_effnew->add(*k_fs_iimod, false, 1.0, 1.0);
  if (iset && aset) k_fs_effnew->add(*k_fs_iamod, false, 1.0, 1.0);

  //---------------------------------------------------------- FOURTH LINE
  // add a submatrices to kteffnew
  if (aset) k_fs_effnew->add(*fNCoup_lindisp_, false, 1.0, 1.0);

  //----------------------------------------------------------- FIFTH LINE
  // add a submatrices to kteffnew
  if (aset) k_fs_effnew->add(*k_fs_anmod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->add(*k_fs_ammod, false, -1.0, 1.0);
  if (aset && iset) k_fs_effnew->add(*k_fs_aimod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->add(*k_fs_aamod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->add(*flinTangentiallambda_, false, 1.0, 1.0);

  // fill_complete kteffnew (square)
  k_fs_effnew->complete(*problem_dofs(), *falldofrowmap_);
  /**********************************************************************/
  /* (11) Global setup of feffnew (including contact)                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));

  Core::LinAlg::Export(*fn, *fnexp);

  feffnew->Update(1.0, *fnexp, 1.0);
  //---------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  if (mset)
  {
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    Core::LinAlg::Export(*fmmod, *fmmodexp);
    feffnew->Update(1.0, *fmmodexp, 1.0);
  }
  //----------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    Core::LinAlg::Export(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //---------------------------------------------------------- FOURTH LINE
  // add weighted nCoup vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> nCoupexp;
  if (aset)
  {
    nCoupexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    Core::LinAlg::Export(*fNCoup_, *nCoupexp);
    feffnew->Update(-1.0, *nCoupexp, 1.0);
  }

  //----------------------------------------------------------- FIFTH LINE
  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    Core::LinAlg::Export(*famod, *famodexp);
    feffnew->Update(-1.0, *famodexp, 1.0);
  }

  // finally do the replacement
  k_fseff = k_fs_effnew;
  feff = feffnew;
  return;
}  // CONTACT::LagrangeStrategyPoro::EvaluatePoroNoPen

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration          |
 |  condition on contact surface (public)                     ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::evaluate_other_mat_poro_no_pen(
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& Feff, int Column_Block_Id)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ ||
      (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()))
    return;
  // this method should be renamed as it handles twosided poro meshtying now aswell and is
  // able to handle fluid coupling for twosided contact

  nopenalpha_ = alphaf_;  // to use different alpha for nopen condition (not used at the moment)

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  Feff->complete();

  Teuchos::RCP<Epetra_Map> domainmap = Teuchos::rcp(new Epetra_Map(Feff->domain_map()));

  // shape function
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Condensation only for dual LM");

  /**********************************************************************/
  /* (3) Split k_fseff and Feff into 3x3 matrix blocks                             */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_s, F_m, F_n;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_sm, F_sm0, F_n0, F_m0, F_s0;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap0;
  Teuchos::RCP<Epetra_Map> tempmap1;
  Teuchos::RCP<Epetra_Map> ftempmap;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<Core::LinAlg::SparseMatrix> Feffmatrix =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(Feff);
  if (parallel_redistribution_status())
  {
    FOUR_C_THROW("CHECK ME!");
    // split and transform to redistributed maps
    //        Core::LinAlg::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    //        ksmsm = Mortar::matrix_row_col_transform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    //        ksmn  = Mortar::MatrixRowTransform(ksmn,gsmdofrowmap_);
    //        knsm  = Mortar::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::SplitMatrix2x2(
        Feffmatrix, fgsmdofrowmap_, fgndofrowmap_, domainmap, tempmap0, F_sm, F_sm0, F_n, F_n0);
  }

  // further splits into slave part + master part
  Core::LinAlg::SplitMatrix2x2(
      F_sm, fgsdofrowmap_, fgmdofrowmap_, domainmap, tempmap0, F_s, F_s0, F_m, F_m0);

  // store some stuff for static condensation of LM
  cfx_s_.insert(std::pair<int, Teuchos::RCP<Core::LinAlg::SparseMatrix>>(Column_Block_Id, F_s));

  /**********************************************************************/
  /* (5) Split slave quantities into active / inactive                  */
  /**********************************************************************/

  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_a, F_a0, F_i, F_i0;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_an, F_in, F_am, F_im, F_ma, F_mi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> fgidofs;

  Core::LinAlg::SplitMatrix2x2(
      F_s, fgactivedofs_, fgidofs, domainmap, tempmap1, F_a, F_a0, F_i, F_i0);

  // abbreviations for active and inactive set
  int aset = fgactivedofs_->NumGlobalElements();
  int iset = fgidofs->NumGlobalElements();

  // abbreviations for slave  and master set
  // int sset = fgsdofrowmap_->NumGlobalElements(); //usally slave should anyway exist!
  int mset = fgmdofrowmap_->NumGlobalElements();

  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /* where K stands for k_fs and F!!!                                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // kn: nothing to do
  //---------------------------------------------------------- SECOND LINE --- Will just exist when
  // starting with two sided poro contact!!!
  // km: add T(mhataam)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_mmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgmdofrowmap_, 100));
  F_mmod->add(*F_m, false, 1.0, 1.0);
  if (aset && mset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> F_madd =
        Core::LinAlg::MLMultiply(*fmhataam_, true, *F_a, false, false, false, true);
    F_mmod->add(*F_madd, false, 1.0, 1.0);
  }
  F_mmod->complete(F_m->domain_map(), F_m->row_map());

  //----------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------

  //--- For using non diagonal D-Matrix, it should be checked if this assumtion isn't anywhere
  // else!!!

  // kin: subtract T(dhat)*kan --
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_imod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgidofs, 100));
  F_imod->add(*F_i, false, 1.0, 1.0);
  if (aset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> F_iadd =
        Core::LinAlg::MLMultiply(*fdhat_, true, *F_a, false, false, false, true);
    F_imod->add(*F_iadd, false, -1.0, 1.0);
  }
  F_imod->complete(F_i->domain_map(), F_i->row_map());

  //---------------------------------------------------------- FOURTH LINE
  // nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_amod;
  if (aset)
  {
    F_amod = Core::LinAlg::MLMultiply(*ftanginvD_, false, *F_a, false, false, false, true);
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
    FOUR_C_THROW("CHECK ME!");
    //        //----------------------------------------------------------- FIRST LINE
    //        // nothing to do (ndof-map independent of redistribution)
    //
    //        //---------------------------------------------------------- SECOND LINE
    //        kmnmod = Mortar::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    //        kmmmod = Mortar::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    //        if (iset) kmimod = Mortar::MatrixRowTransform(kmimod,pgmdofrowmap_);
    //        if (aset) kmamod = Mortar::MatrixRowTransform(kmamod,pgmdofrowmap_);
    //
    //        //----------------------------------------------------------- THIRD LINE
    //        if (iset)
    //        {
    //          kinmod = Mortar::MatrixRowTransform(kinmod,pgsdofrowmap_);
    //          kimmod = Mortar::MatrixRowTransform(kimmod,pgsdofrowmap_);
    //          kiimod = Mortar::MatrixRowTransform(kiimod,pgsdofrowmap_);
    //          if (aset) kiamod = Mortar::MatrixRowTransform(kiamod,pgsdofrowmap_);
    //        }
    //
    //        //---------------------------------------------------------- FOURTH LINE
    //        if (aset) smatrix_ = Mortar::MatrixRowTransform(smatrix_,pgsdofrowmap_);
    //
    //        //----------------------------------------------------------- FIFTH LINE
    //        if (aset)
    //        {
    //          kanmod = Mortar::MatrixRowTransform(kanmod,pgsdofrowmap_);
    //          kammod = Mortar::MatrixRowTransform(kammod,pgsdofrowmap_);
    //          kaamod = Mortar::MatrixRowTransform(kaamod,pgsdofrowmap_);
    //          if (iset) kaimod = Mortar::MatrixRowTransform(kaimod,pgsdofrowmap_);
    //          pmatrix_ = Mortar::MatrixRowTransform(tderivmatrix_,pgsdofrowmap_);
    //        }
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> F_effnew = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*falldofrowmap_, 108, true, false, Feff->get_matrixtype()));

  //----------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  F_effnew->add(*F_n, false, 1.0, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  if (mset) F_effnew->add(*F_mmod, false, 1.0, 1.0);
  //----------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) F_effnew->add(*F_imod, false, 1.0, 1.0);
  //---------------------------------------------------------- FOURTH LINE
  // add a submatrices to kteffnew //assume that Column_Block_Id==0 is the porofluid coupling
  // block!!!
  if (Column_Block_Id == 0 && aset) F_effnew->add(*fNCoup_linvel_, false, 1.0, 1.0);
  //----------------------------------------------------------- FIFTH LINE
  // add a submatrices to kteffnew
  if (aset) F_effnew->add(*F_amod, false, -1.0, 1.0);

  // fill_complete kteffnew (square)
  F_effnew->complete(*domainmap, *falldofrowmap_);

  // finally do the replacement
  Feff = F_effnew;
  return;
}  // CONTACT::LagrangeStrategyPoro::evaluate_other_mat_poro_no_pen

/*----------------------------------------------------------------------------*
 | Poro Recovery method for no penetration LM (pure porous problem) ager 08/14|
 *---------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::recover_poro_no_pen(
    Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> inc)
{
  std::map<int, Teuchos::RCP<Epetra_Vector>> incm;
  incm.insert(std::pair<int, Teuchos::RCP<Epetra_Vector>>(0, inc));

  recover_poro_no_pen(disi, incm);
  return;
}

/*----------------------------------------------------------------------*
 | Poro Recovery method for no penetration LM                 ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::recover_poro_no_pen(
    Teuchos::RCP<Epetra_Vector> disi, std::map<int, Teuchos::RCP<Epetra_Vector>> inc)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ ||
      (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()))
    return;

  // shape function and system types
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::Export(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::Export(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::Export(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> finvda;
    Teuchos::RCP<Epetra_Map> tempmap1, tempmap2;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    Core::LinAlg::SplitMatrix2x2(finvda_, fgactivedofs_, tempmap1, gactivedofs_, tempmap2, finvda,
        tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> finvdmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*fgsdofrowmap_, 10));
    finvdmod->add(*finvda, false, 1.0, 1.0);
    finvdmod->complete(*gsdofrowmap_, *fgsdofrowmap_);

    /**********************************************************************/
    /* Update Lagrange multipliers lambda_n+1                                  */
    /**********************************************************************/
    {
      Teuchos::RCP<Epetra_Vector> flambda = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_, true));

      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

      cfssn_->multiply(false, *disin, *mod);
      flambda->Update(-1.0, *mod, 0.0);
      cfssm_->multiply(false, *disim, *mod);
      flambda->Update(-1.0, *mod, 1.0);
      cfsss_->multiply(false, *disis, *mod);
      flambda->Update(-1.0, *mod, 1.0);

      // loop over all offdiag blocks!!!
      std::map<int, Teuchos::RCP<Core::LinAlg::SparseOperator>>::iterator matiter;
      std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator inciter;
      for (matiter = cfx_s_.begin(); matiter != cfx_s_.end(); ++matiter)
      {
        inciter = inc.find(matiter->first);
        if (inciter == inc.end())
          FOUR_C_THROW(
              "CONTACT::LagrangeStrategyPoro::RecoverPoroNoPen: Couldn't find increment block %d "
              "for recovery of the lagrange multiplier!",
              matiter->first);

        matiter->second->multiply(false, *inciter->second, *mod);
        flambda->Update(-1.0, *mod, 1.0);
      }

      flambda->Update(1.0, *ffs_, 1.0);



      fdoldtransp_->multiply(false, *lambdaold_, *mod);
      flambda->Update(-nopenalpha_, *mod, 1.0);

      Teuchos::RCP<Epetra_Vector> lambdacopy = Teuchos::rcp(new Epetra_Vector(*flambda));

      finvdmod->multiply(true, *lambdacopy, *lambda_);  // should be lambda_ at the end!!!

      lambda_->Scale(
          (1 - alphaf_) / (1 - nopenalpha_));  //-- is already scaled by this factor by scaling
                                               // invda_!!! --- scale it back to with nopenalpha_...
    }
  }
  // store updated LM into nodes
  set_state(Mortar::state_lagrange_multiplier, *lambda_);

  return;
}

/*------------------------------------------------------------------------------*
 |  Additional update and output poro contact at end of time step     ager 08/14|
 *-----------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::update_poro_contact()
{
  if (no_penetration_)
  {
    // std::cout << "print lambda: " << *lambda_ << std::endl;
    lambdaold_->Update(1.0, *lambda_, 0.0);
  }
  return;
}

/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                          ager 08/14|
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::set_state(
    const enum Mortar::StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case Mortar::state_fvelocity:
    case Mortar::state_svelocity:
    case Mortar::state_lagrange_multiplier:
    case Mortar::state_fpressure:
    {
      // set state on interfaces
      for (int i = 0; i < (int)interface_.size(); ++i)
      {
        // interface_[i]->set_state(statename, vec);
        Core::FE::Discretization& idiscret_ = interface_[i]->discret();

        switch (statetype)
        {
          case Mortar::state_fvelocity:
          case Mortar::state_svelocity:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.dof_col_map(), true));

            Core::LinAlg::Export(vec, *global);

            // loop over all nodes to set current velocity
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.num_my_col_nodes(); ++i)
            {
              CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_.l_col_node(i));
              const int numdof = node->num_dof();
              std::vector<double> myvel(numdof);
              std::vector<int> lm(numdof);

              for (int j = 0; j < numdof; ++j) lm[j] = node->dofs()[j];
              Core::FE::ExtractMyValues(*global, myvel, lm);

              // add myvel[2]=0 for 2D problems
              if (myvel.size() < 3) myvel.resize(3);
              // set current configuration
              for (int j = 0; j < 3; ++j)
              {
                if (statetype == Mortar::state_fvelocity)
                  node->poro_data().fvel()[j] = myvel[j];
                else
                  node->poro_data().svel()[j] = myvel[j];
              }
            }
            break;
          }
          case Mortar::state_lagrange_multiplier:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.dof_col_map(), true));
            Core::LinAlg::Export(vec, *global);

            // loop over all nodes to set current velocity
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.num_my_col_nodes(); ++i)
            {
              CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_.l_col_node(i));

              const int numdof = node->num_dof();
              std::vector<double> mylm(numdof);
              std::vector<int> lm(numdof);

              for (int j = 0; j < numdof; ++j) lm[j] = node->dofs()[j];

              Core::FE::ExtractMyValues(*global, mylm, lm);

              // add myvel[2]=0 for 2D problems
              if (mylm.size() < 3) mylm.resize(3);
              // set current configuration
              if (node->is_slave())
              {
                for (int j = 0; j < 3; ++j)
                {
                  node->poro_data().poro_lm()[j] = mylm[j];
                }
              }
            }
            break;
          }
          case Mortar::state_fpressure:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.dof_col_map(), true));
            Core::LinAlg::Export(vec, *global);

            // loop over all nodes to set current pressure
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.num_my_col_nodes(); ++i)
            {
              CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_.l_col_node(i));

              double myfpres;
              int fpres;

              fpres = node->dofs()[0];  // here get ids of first component of node

              myfpres = global->Values()[global->Map().LID(fpres)];

              *node->poro_data().fpres() = myfpres;
            }
            break;
          }
          default:
          {
            FOUR_C_THROW("Shouldn't happen!");
            break;
          }
        }  // end inner switch statement
      }    // end loop over all interfaces
      break;
    }
    default:
    {
      CONTACT::AbstractStrategy::set_state(statetype, vec);
      break;
    }
  }  // end outer switch statement
  return;
}


// this should add the displacement of the parent element to the contact element into the
// datacontainer...
//...and it should additionally add the dof IDs to its Data container
/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                          ager 10/14|
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::set_parent_state(const std::string& statename,
    const Teuchos::RCP<Epetra_Vector> vec, const Teuchos::RCP<Core::FE::Discretization> dis)
{
  if (statename == "displacement")
  {
    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->dof_col_map(), true));
    Core::LinAlg::Export(*vec, *global);

    // set state on interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      Core::FE::Discretization& idiscret_ = interface_[i]->discret();

      if (poroslave_)
      {
        for (int j = 0; j < interface_[i]->slave_col_elements()->NumMyElements();
             ++j)  // will just work for onesided poro contact as the porosity is just on slave
                   // side!!!
        {
          int gid = interface_[i]->slave_col_elements()->GID(j);

          Mortar::Element* ele = dynamic_cast<Mortar::Element*>(idiscret_.g_element(gid));

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          // this gets values in local order
          ele->parent_element()->location_vector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          Core::FE::ExtractMyValues(*global, myval, lm);

          ele->mo_data().parent_disp() = myval;
          ele->mo_data().parent_dof() = lm;
        }
      }
      if (poromaster_)  // add master parent element displacements
      {
        for (int j = 0; j < interface_[i]->master_col_elements()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->master_col_elements()->GID(j);

          Mortar::Element* mele = dynamic_cast<Mortar::Element*>(idiscret_.g_element(gid));

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          // this gets values in local order
          mele->parent_element()->location_vector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          Core::FE::ExtractMyValues(*global, myval, lm);

          mele->mo_data().parent_disp() = myval;
          mele->mo_data().parent_dof() = lm;
        }
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | Initialize poro meshtying matrices with corresponding maps  h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::poro_mt_initialize()
{
  // (re)setup global lin of D & M * lambda - Matrix
  porolindmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *gsdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
  porolinmmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *gmdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

  // (re)setup global Mortar Core::LinAlg::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
  mmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100));

  // this is needed because on the mesh tying way to build this strategy it is not done elsewhere
  setup_no_penetration_condition();

  isincontact_ = true;      // simply set true for meshtying
  wasincontact_ = true;     // as meshtying interfaces stay the same and are fully in contact
  wasincontactlts_ = true;  // this is necessary for other methods in this strategy
  return;
}  // CONTACT::LagrangeStrategyPoro::PoroLinkDM

/*------------------------------------------------------------------------*
 | assemble porofluid meshtying matrices                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::poro_mt_prepare_fluid_coupling()
{
  // reset D and M matrices for new Newton iteration

  dmatrix_->zero();
  mmatrix_->zero();

  // Assemble D and M matrices
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_dm(*dmatrix_, *mmatrix_);
  }

  // complete D and M matrices
  dmatrix_->complete();
  mmatrix_->complete(*gmdofrowmap_, *gsdofrowmap_);

  // as mhataam-, dhat_ and invda_ are not computed in poro - meshtying before this point it is
  // necessary here
  poro_mt_set_coupling_matrices();

  return;
}  // CONTACT::LagrangeStrategyPoro::PoroAssembleFluidCoupling()

/*------------------------------------------------------------------------*
 | assemble porofluid meshtying matrices                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::poro_mt_set_coupling_matrices()
{
  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

  Teuchos::RCP<Epetra_Map> gidofs;

  int aset = gactivedofs_->NumGlobalElements();

  // invert dmatrix to invd
  invd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd_->extract_diagonal_copy(*diag);

  // set zero diagonal values to dummy 1.0
  //  for (int i = 0; i < diag->MyLength(); ++i)
  //    if ((*diag)[i] == 0.0)
  //      (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err > 0) FOUR_C_THROW("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd_->replace_diagonal_values(*diag);
  // inversion end

  // active part of invd
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
  Core::LinAlg::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)(not considered for poro
  // mt/contact yet-> ==0)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
  Core::LinAlg::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  int iset = gidofs->NumGlobalElements();
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

  // scaling of invd and dai
  invda->scale(1 / (1 - alphaf_));

  save_coupling_matrices(dhat, mhataam, invda);

  return;
}  // CONTACT::LagrangeStrategyPoro::poro_mt_prepare_fluid_coupling

/*------------------------------------------------------------------------*
 | update old meshtying matrices and LMP                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::LagrangeStrategyPoro::poro_mt_update()
{
  // set dold_ mold_ and lambdaold_ after every time step
  dold_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(new Core::LinAlg::SparseMatrix(*dmatrix_));
  mold_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(new Core::LinAlg::SparseMatrix(*mmatrix_));
  dold_->complete(dmatrix_->domain_map(), dmatrix_->range_map());
  mold_->complete(mmatrix_->domain_map(), mmatrix_->range_map());
  update_poro_contact();
  return;
}  // CONTACT::LagrangeStrategyPoro::PoroMtUpdate()

FOUR_C_NAMESPACE_CLOSE
