/*---------------------------------------------------------------------*/
/*! \file
\brief This class provides the functionality to use contact with Lagrangian
multipliers for monolithical coupled multifield problems!
Therefore apply_force_stiff_cmt() & recover() are overloaded by this class and
do nothing, as they are called directly in the structure. To use the contact
the additional methods apply_force_stiff_cmt_coupled() & RecoverCoupled() have
to be called!


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_monocoupled_lagrange_strategy.hpp"

#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_SerialComm.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                              ager 02/15|
 *----------------------------------------------------------------------*/
CONTACT::MonoCoupledLagrangeStrategy::MonoCoupledLagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
    Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
    : LagrangeStrategy(
          data_ptr, dof_row_map, NodeRowMap, params, interface, dim, comm, alphaf, maxdof),
      has_to_evaluate_(false),
      has_to_recover_(false)
{
  // do some security checks ...
  return;
}

/*----------------------------------------------------------------------*
 | structural contact global evaluation method called from time         |
 | integrator + condensation of offdiagonal Matrixes (public) ager 02/15|
 *----------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::apply_force_stiff_cmt_coupled(
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Core::LinAlg::SparseOperator>& k_ss,
    std::map<int, Teuchos::RCP<Core::LinAlg::SparseOperator>*> k_sx,
    Teuchos::RCP<Epetra_Vector>& rhs_s, const int step, const int iter, bool predictor)
{
  // call the main routine for contact!!!
  CONTACT::AbstractStrategy::apply_force_stiff_cmt(dis, k_ss, rhs_s, step, iter, predictor);

  // Take care of the alternative condensation of the off-diagonal blocks!!!
  std::map<int, Teuchos::RCP<Core::LinAlg::SparseOperator>*>::iterator matiter;
  for (matiter = k_sx.begin(); matiter != k_sx.end(); ++matiter)
  {
    evaluate_off_diag_contact(*(matiter->second), matiter->first);
  }
  has_to_evaluate_ = false;
  return;
}

/*-----------------------------------------------------------------------------*
 | structural contact global evaluation method called from time                |
 | integrator + condensation of one!!! offdiagonal Matrixes (public) ager 02/15|
 *----------------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::apply_force_stiff_cmt_coupled(
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Core::LinAlg::SparseOperator>& k_ss,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& k_sx, Teuchos::RCP<Epetra_Vector>& rhs_s,
    const int step, const int iter, bool predictor)
{
  // call the main routine for contact!!!
  CONTACT::AbstractStrategy::apply_force_stiff_cmt(dis, k_ss, rhs_s, step, iter, predictor);

  // Take care of the alternative condensation of the off-diagonal blocks!!!
  evaluate_off_diag_contact(k_sx, 0);

  has_to_evaluate_ = false;
  return;
}

/*------------------------------------------------------------------------*
 |  condense off-diagonal blocks                      (public)  ager 02/15|
 *-----------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::evaluate_off_diag_contact(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, int Column_Block_Id)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();

  Teuchos::RCP<Epetra_Map> domainmap = Teuchos::rcp(new Epetra_Map(kteff->DomainMap()));

  // system type
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(Params(), "SYSTEM");

  // shape function
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(Params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    /**********************************************************************/
    /* (3) Split kteff into 3x3 matrix blocks                             */  // just split the rows
                                                                              // !!!
    /**********************************************************************/

    // we want to split k into 3 groups s,m,n = 9 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ks, km, kn;

    // temporarily we need the blocks ksmsm, ksmn, knsm
    // (FIXME: because a direct SplitMatrix3x3 is still missing!)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ksm, ksm0, kn0, km0, ks0;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap0;
    Teuchos::RCP<Epetra_Map> tempmap1;
    Teuchos::RCP<Epetra_Map> ftempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

    // split into slave/master part + structure part
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);

    if (parallel_redistribution_status())  // TODO Check if how to modify
    {
      FOUR_C_THROW("parallel_redistribution_status(): CHECK ME!");
      // split and transform to redistributed maps
      //      Core::LinAlg::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
      //      ksmsm = Mortar::matrix_row_col_transform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
      //      ksmn  = Mortar::MatrixRowTransform(ksmn,gsmdofrowmap_);
      //      knsm  = Mortar::MatrixColTransform(knsm,gsmdofrowmap_);
    }
    else
    {
      // only split, no need to transform
      Core::LinAlg::SplitMatrix2x2(
          kteffmatrix, gsmdofrowmap_, gndofrowmap_, domainmap, tempmap0, ksm, ksm0, kn, kn0);
    }

    // further splits into slave part + master part
    Core::LinAlg::SplitMatrix2x2(
        ksm, gsdofrowmap_, gmdofrowmap_, domainmap, tempmap0, ks, ks0, km, km0);

    // store some stuff for static condensation of LM
    csx_s_.insert(std::pair<int, Teuchos::RCP<Core::LinAlg::SparseMatrix>>(Column_Block_Id, ks));


    /**********************************************************************/
    /* (5) Split slave quantities into active / inactive                  */
    /**********************************************************************/

    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ka, ka0, ki, ki0;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;
    Teuchos::RCP<Epetra_Map> fgidofs;

    // do the splitting
    Core::LinAlg::SplitMatrix2x2(ks, gactivedofs_, gidofs, domainmap, tempmap1, ka, ka0, ki, ki0);

    // abbreviations for master, active and inactive set
    int aset = gactivedofs_->NumGlobalElements();
    int iset = gidofs->NumGlobalElements();

    /**********************************************************************/
    /* (7) Build the final K blocks                                       */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // kn: nothing to do

    //---------------------------------------------------------- SECOND LINE
    // km: add T(mhataam)*kan
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmmod->Add(*km, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kmadd =
          Core::LinAlg::MLMultiply(*mhataam_, true, *ka, false, false, false, true);
      kmmod->Add(*kmadd, false, 1.0, 1.0);
    }
    kmmod->Complete(kteff->DomainMap(), km->RowMap());

    //----------------------------------------------------------- THIRD LINE
    //------------------- FOR 3D QUADRATIC CASE ----------------------------

    //--- For using non diagonal D-Matrix, it should be checked if this assumtion isn't anywhere
    // else!!!

    // kin: subtract T(dhat)*kan --
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kimod->Add(*ki, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiadd =
          Core::LinAlg::MLMultiply(*dhat_, true, *ka, false, false, false, true);
      kimod->Add(*kiadd, false, -1.0, 1.0);
    }
    kimod->Complete(kteff->DomainMap(), ki->RowMap());

    //---------------------------------------------------------- FOURTH LINE
    // nothing to do

    //----------------------------------------------------------- FIFTH LINE
    // ka: multiply tmatrix with invda and ka
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kamod;
    if (aset)
    {
      kamod = Core::LinAlg::MLMultiply(*tmatrix_, false, *invda_, true, false, false, true);
      kamod = Core::LinAlg::MLMultiply(*kamod, false, *ka, false, false, false, true);
    }

    /********************************************************************/
    /* (9) Transform the final K blocks                                 */
    /********************************************************************/
    // The row maps of all individual matrix blocks are transformed to
    // the parallel layout of the underlying problem discretization.
    // Of course, this is only necessary in the parallel redistribution
    // case, where the contact interfaces have been redistributed
    // independently of the underlying problem discretization.

    if (parallel_redistribution_status())  // check what to do
    {
      FOUR_C_THROW("not checked so far!!!");
      //----------------------------------------------------------- FIRST LINE
      // nothing to do (ndof-map independent of redistribution)

      //      //---------------------------------------------------------- SECOND LINE
      //      kmnmod = Mortar::MatrixRowTransform(kmnmod,pgmdofrowmap_);
      //      kmmmod = Mortar::MatrixRowTransform(kmmmod,pgmdofrowmap_);
      //      if (iset) kmimod = Mortar::MatrixRowTransform(kmimod,pgmdofrowmap_);
      //      if (aset) kmamod = Mortar::MatrixRowTransform(kmamod,pgmdofrowmap_);
      //
      //      //----------------------------------------------------------- THIRD LINE
      //      if (iset)
      //      {
      //        kinmod = Mortar::MatrixRowTransform(kinmod,pgsdofrowmap_);
      //        kimmod = Mortar::MatrixRowTransform(kimmod,pgsdofrowmap_);
      //        kiimod = Mortar::MatrixRowTransform(kiimod,pgsdofrowmap_);
      //        if (aset) kiamod = Mortar::MatrixRowTransform(kiamod,pgsdofrowmap_);
      //      }
      //
      //      //---------------------------------------------------------- FOURTH LINE
      //      if (aset) smatrix_ = Mortar::MatrixRowTransform(smatrix_,pgsdofrowmap_);
      //
      //      //----------------------------------------------------------- FIFTH LINE
      //      if (aset)
      //      {
      //        kanmod = Mortar::MatrixRowTransform(kanmod,pgsdofrowmap_);
      //        kammod = Mortar::MatrixRowTransform(kammod,pgsdofrowmap_);
      //        kaamod = Mortar::MatrixRowTransform(kaamod,pgsdofrowmap_);
      //        if (iset) kaimod = Mortar::MatrixRowTransform(kaimod,pgsdofrowmap_);
      //        pmatrix_ = Mortar::MatrixRowTransform(pmatrix_,pgsdofrowmap_);
      //      }
    }

    /**********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                  */
    /**********************************************************************/

    Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *gdisprowmap_, 81, true, false, kteffmatrix->GetMatrixtype()));

    //----------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->Add(*kn, false, 1.0, 1.0);
    //---------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->Add(*kmmod, false, 1.0, 1.0);
    //----------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->Add(*kimod, false, 1.0, 1.0);

    //---------------------------------------------------------- FOURTH LINE
    // for off diag blocks this line is empty (weighted normal = f(disp))

    //----------------------------------------------------------- FIFTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->Add(*kamod, false, 1.0, 1.0);

    // fill_complete kteffnew (square)
    kteffnew->Complete(*domainmap, *gdisprowmap_);

    // finally do the replacement
    kteff = kteffnew;
  }
  else
  {
    FOUR_C_THROW("Trying to use not condensed form --- Feel Free to implement!");
  }
  return;
}

/*------------------------------------------------------------------------*
 | Coupled Recovery method for contact LM                       ager 02/15|
 *-----------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::RecoverCoupled(
    Teuchos::RCP<Epetra_Vector> disi, std::map<int, Teuchos::RCP<Epetra_Vector>> inc)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  LagrangeStrategy::recover(
      disi);  // Update Structural Part! --> Here just Part from Coupling Matrix will be added!

  if (inc.size() == 0 && csx_s_.size() == 0)
    return;  // done already here if there are no off-diag blocks

  // shape function and system types
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(Params(), "LM_SHAPEFCN");
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(Params(), "SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    if (inc.size() != csx_s_.size())
      FOUR_C_THROW(
          "CONTACT::MonoCoupledLagrangeStrategy::RecoverCoupled: For Recovery the same number of "
          "off-diagonal increment blocks is required! %d != %d !",
          inc.size(), csx_s_.size());

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
    invdmod->Add(*invda, false, 1.0, 1.0);
    invdmod->Complete();

    std::map<int, Teuchos::RCP<Core::LinAlg::SparseOperator>>::iterator matiter;
    std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator inciter;

    // loop over all offdiag blocks!!!
    for (matiter = csx_s_.begin(); matiter != csx_s_.end(); ++matiter)
    {
      inciter = inc.find(matiter->first);
      if (inciter == inc.end())
        FOUR_C_THROW(
            "CONTACT::MonoCoupledLagrangeStrategy::RecoverCoupled: Couldn't find increment block "
            "%d for recovery of the lagrange multiplier!",
            matiter->first);

      /**********************************************************************/
      /* Update Lagrange multipliers z_n+1                                  */
      /**********************************************************************/
      // for self contact, slave and master sets may have changed,
      // thus we have to export the products Dold * zold and Mold^T * zold to fit
      if (is_self_contact())  // is not considered yet!
      {
        FOUR_C_THROW(
            "Trying to make coupled selfcontact condensation... Check if this makes any sense!!!");
        // approximate update
        // z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
        // invdmod->Multiply(false,*fs_,*z_);

        // full update
        //      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
        //      z_->Update(1.0,*fs_,0.0);
        //      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
        //      kss_->Multiply(false,*disis,*mod);
        //      z_->Update(-1.0,*mod,1.0);
        //      ksm_->Multiply(false,*disim,*mod);
        //      z_->Update(-1.0,*mod,1.0);
        //      ksn_->Multiply(false,*disin,*mod);
        //      z_->Update(-1.0,*mod,1.0);
        //      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new
        //      Epetra_Vector((dold_->RowMap()))); if (dold_->RowMap().NumGlobalElements())
        //      Core::LinAlg::Export(*zold_,*mod2); Teuchos::RCP<Epetra_Vector> mod3 =
        //      Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
        //      dold_->Multiply(true,*mod2,*mod3); Teuchos::RCP<Epetra_Vector> mod4 =
        //      Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_)); if
        //      (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::Export(*mod3,*mod4);
        //      z_->Update(-alphaf_,*mod4,1.0);
        //      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
        //      invdmod->Multiply(true,*zcopy,*z_);
        //      z_->Scale(1/(1-alphaf_));
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> zfluid = Teuchos::rcp(new Epetra_Vector(z_->Map(), true));

        Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
        matiter->second->Multiply(false, *inciter->second, *mod);
        zfluid->Update(-1.0, *mod, 0.0);
        Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*zfluid));
        invdmod->Multiply(true, *zcopy, *zfluid);
        zfluid->Scale(1 / (1 - alphaf_));

        z_->Update(1.0, *zfluid, 1.0);  // Add Offdiag  -  Coupling Contribution to LM!!!
      }
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
  }

  // store updated LM into nodes
  store_nodal_quantities(Mortar::StrategyBase::lmupdate);  // Here done twice: already in structural
                                                           // contact --> not wanted

  has_to_recover_ = false;
  return;
}


/*-------------------------------------------------------------------------*
 | Coupled Recovery method for contact LM with one offdiag block ager 02/15|
 *------------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::RecoverCoupled(
    Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> inc)
{
  std::map<int, Teuchos::RCP<Epetra_Vector>> incm;
  incm.insert(std::pair<int, Teuchos::RCP<Epetra_Vector>>(0, inc));

  RecoverCoupled(disi, incm);
  return;
}

/*---------------------------------------------------------------------------*
 | Save mortar coupling matrices for evaluation of off diag terms! ager 08/14|
 *--------------------------------------------------------------------------*/
void CONTACT::MonoCoupledLagrangeStrategy::save_coupling_matrices(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat, Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda)
{
  dhat_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(new Core::LinAlg::SparseMatrix(*dhat));
  mhataam_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(new Core::LinAlg::SparseMatrix(*mhataam));
  invda_ = Teuchos::rcp<Core::LinAlg::SparseMatrix>(new Core::LinAlg::SparseMatrix(*invda));
}

FOUR_C_NAMESPACE_CLOSE
