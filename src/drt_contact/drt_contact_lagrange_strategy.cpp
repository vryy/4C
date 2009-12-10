/*!----------------------------------------------------------------------
\file drt_contact_lagrange_strategy.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_lagrange_strategy.H"
#include "drt_contact_abstract_strategy.H"

#include "Epetra_SerialComm.h"
#include "drt_contact_manager_base.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::LagrangeStrategy::LagrangeStrategy(RCP<Epetra_Map> problemrowmap,
                                            Teuchos::ParameterList params,
                                            vector<RCP<CONTACT::Interface> > interface, int dim,
                                            RCP<Epetra_Comm> comm, double alphaf) :
AbstractStrategy(problemrowmap, params, interface, dim, comm, alphaf)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::Initialize()
{
  // (re)setup global normal and tangent matrices
  nmatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  tmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));

  // (re)setup global matrix containing gap derivatives
  smatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));

  // further terms depend on friction case
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

  // (re)setup global matrix containing "no-friction"-derivatives
  if (ftype == INPAR::CONTACT::friction_none)
  {
    pmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  }

  // (re)setup global Tresca friction / perfect stick / MeshTying
  if (ftype == INPAR::CONTACT::friction_tresca ||
      ftype == INPAR::CONTACT::friction_coulomb ||
      ftype == INPAR::CONTACT::friction_stick)
  {
    lmatrix_ = rcp(new LINALG::SparseMatrix(*gslipt_,10));
    r_       = LINALG::CreateVector(*gslipt_,true);

    // here the calculation of gstickt is necessary
    RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
    linstickLM_ = rcp(new LINALG::SparseMatrix(*gstickt,3));
    linstickDIS_ = rcp(new LINALG::SparseMatrix(*gstickt,3));
    linstickRHS_ = LINALG::CreateVector(*gstickt,true);

    linslipLM_ = rcp(new LINALG::SparseMatrix(*gslipt_,3));
    linslipDIS_ = rcp(new LINALG::SparseMatrix(*gslipt_,3));
    linslipRHS_ = LINALG::CreateVector(*gslipt_,true);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate frictional contact (public)                    gitterle 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::EvaluateFriction(RCP<LINALG::SparseMatrix> kteff,
                                                  RCP<Epetra_Vector> feff)
{
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...

  // input parameters
  bool fulllin = Teuchos::getIntegralValue<int>(Params(),"FULL_LINEARIZATION");

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }

  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix linstick with derivatives of stick nodes         */
  /* and global matrix linslip with derivatives of slip nodes           */
  /**********************************************************************/
  // here and for the splitting later, we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
    interface_[i]->AssembleLinStick(*linstickLM_,*linstickDIS_,*linstickRHS_);
    interface_[i]->AssembleLinSlip(*linslipLM_,*linslipDIS_,*linslipRHS_);
  }

  // FillComplete() global matrices N and T and L
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);

  // FillComplete() global matrix  L
  lmatrix_->Complete(*gslipt_,*gslipt_);

  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofs,*gactiven_);

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofs,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofs,*gmdofrowmap_);

  // FillComplete global Matrix LinStick
  RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
  RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);
  linstickLM_->Complete(*gstickdofs,*gstickt);
  linstickDIS_->Complete(*gsmdofs,*gstickt);

  // FillComplete global Matrix linslipLM and linslipDIS
  linslipLM_->Complete(*gslipdofs_,*gslipt_);
  linslipDIS_->Complete(*gsmdofs,*gslipt_);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);

  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);

  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;

  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split LinD and LinM into blocks                                    */
  /**********************************************************************/
  // we want to split lindmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> lindai, lindaa, lindam, lindas;

  // we want to split linmmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> linmmi, linmma, linmmm, linmms;

  if (fulllin)
  {
    // do the splitting
    LINALG::SplitMatrix2x2(lindmatrix_,gactivedofs_,gidofs,gmdofrowmap_,gsdofrowmap_,lindam,lindas,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(lindas,gactivedofs_,tempmap,gactivedofs_,gidofs,lindaa,lindai,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmmatrix_,gmdofrowmap_,tempmap,gmdofrowmap_,gsdofrowmap_,linmmm,linmms,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,linmma,linmmi,tempmtx1,tempmtx2);

    // modification of kai, kaa, kam
    // (this has to be done first as they are needed below)
    // (note, that kai, kaa, kam have to be UNcompleted again first!!!)
    kai->UnComplete();
    kaa->UnComplete();
    kam->UnComplete();
    kai->Add(*lindai,false,1.0-alphaf_,1.0);
    kaa->Add(*lindaa,false,1.0-alphaf_,1.0);
    kam->Add(*lindam,false,1.0-alphaf_,1.0);
    kai->Complete(*gidofs,*gactivedofs_);
    kaa->Complete();
    kam->Complete(*gmdofrowmap_,*gactivedofs_);
  }

  /**********************************************************************/
  /* Split active quantities into slip / stick                          */
  /**********************************************************************/

  // we want to split kaa into 2 groups sl,st = 4 blocks
  RCP<LINALG::SparseMatrix> kslsl, kslst, kstsl, kstst;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  RCP<LINALG::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary RCPs
  RCP<Epetra_Map> temp1map;
  RCP<LINALG::SparseMatrix> temp1mtx1;
  RCP<LINALG::SparseMatrix> temp1mtx2;
  RCP<LINALG::SparseMatrix> temp1mtx3;
  RCP<LINALG::SparseMatrix> temp1mtx4;

  // we will get the stick rowmap as a by-product
  RCP<Epetra_Map> gstdofs;

  // do the splitting
  LINALG::SplitMatrix2x2(kaa,gslipdofs_,gstdofs,gslipdofs_,gstdofs,kslsl,kslst,kstsl,kstst);
  LINALG::SplitMatrix2x2(kan,gslipdofs_,gstdofs,gndofrowmap_,temp1map,ksln,temp1mtx1,kstn,temp1mtx2);
  LINALG::SplitMatrix2x2(kam,gslipdofs_,gstdofs,gmdofrowmap_,temp1map,kslm,temp1mtx1,kstm,temp1mtx2);
  LINALG::SplitMatrix2x2(kai,gslipdofs_,gstdofs,gidofs,temp1map,ksli,temp1mtx1,ksti,temp1mtx2);

  // we want to split fs into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*gactivedofs_));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    *fa = *fs;
  else if (!gactivedofs_->NumGlobalElements())
    *fi = *fs;
  else
  {
    LINALG::SplitVector(*fs,*gactivedofs_,fa,*gidofs,fi);
  }

  /**********************************************************************/
  /* Isolate active and slip part from mhat, invd and dold              */
  /* Also isolate slip part form dmatrix_, mmatrix_, dold_ and mold_    */
  /* Isolate slip part from T                                           */
  /**********************************************************************/

  RCP<LINALG::SparseMatrix> mhata, mhatst;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);
  mhata_=mhata;

  RCP<LINALG::SparseMatrix> invda, invdsl, invdst;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invd_,gslipdofs_,gstdofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,invdst);
  invda->Scale(1/(1-alphaf_));
  invdsl->Scale(1/(1-alphaf_));
  invdst->Scale(1/(1-alphaf_));

  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);

  RCP<LINALG::SparseMatrix> dmatrixsl, doldsl, dmatrixst, doldst, mmatrixsl, mmatrixst, moldsl, moldst;
  LINALG::SplitMatrix2x2(dmatrix_,gslipdofs_,gstdofs,gslipdofs_,gstdofs,dmatrixsl,tempmtx1,tempmtx2,dmatrixst);
  LINALG::SplitMatrix2x2(dold_,gslipdofs_,gstdofs,gslipdofs_,gstdofs,doldsl,tempmtx1,tempmtx2,doldst);
  LINALG::SplitMatrix2x2(mmatrix_,gslipdofs_,gstdofs,gmdofrowmap_,tempmap,mmatrixsl,tempmtx2,mmatrixst,tempmtx3);
  LINALG::SplitMatrix2x2(mold_,gslipdofs_,gstdofs,gmdofrowmap_,tempmap,moldsl,tempmtx2,moldst,tempmtx3);

  // FIXGIT: Is this scaling really necessary
  dmatrixsl->Scale(1/(1-alphaf_));
  doldsl->Scale(1/(1-alphaf_));
  mmatrixsl->Scale(1/(1-alphaf_));
  moldsl->Scale(1/(1-alphaf_));

  // temporary RCPs
  RCP<Epetra_Map> tmap;

  // temporary RCPs
  RCP<LINALG::SparseMatrix> tm1, tm2;

  // we want to split the tmatrix_ into 2 groups
  RCP<LINALG::SparseMatrix> tslmatrix, tstmatrix;
  LINALG::SplitMatrix2x2(tmatrix_,gslipt_,gstickt,gslipdofs_,tmap,tslmatrix,tm1,tm2,tstmatrix);

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhata,true,*kan,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhata,true,*kam,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  if (fulllin) kmmmod->Add(*linmmm,false,1.0-alphaf_,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod = LINALG::Multiply(*mhata,true,*kai,false,false);
  kmimod->Add(*kmi,false,1.0,1.0);
  if (fulllin) kmimod->Add(*linmmi,false,1.0-alphaf_,1.0);
  kmimod->Complete(kmi->DomainMap(),kmi->RowMap());

  // kma: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod = LINALG::Multiply(*mhata,true,*kaa,false,false);
  kmamod->Add(*kma,false,1.0,1.0);
  if (fulllin) kmamod->Add(*linmma,false,1.0-alphaf_,1.0);
  kmamod->Complete(kma->DomainMap(),kma->RowMap());

  // kin: nothing to do

  // kim: nothing to do

  // kii: nothing to do

  // kisl: nothing to do

  // kist: nothing to do

  // n*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata,false,true);

 // nmatrix: nothing to do

 // blocks for complementary conditions (stick nodes) - from LM

  // kstn: multiply with linstickLM
  RCP<LINALG::SparseMatrix> kstnmod = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);
  kstnmod = LINALG::Multiply(*kstnmod,false,*kstn,false,true);
  kstnmod->Complete(kstn->DomainMap(),kstn->RowMap());

  // kstm: multiply with linstickLM
  RCP<LINALG::SparseMatrix> kstmmod = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);
  kstmmod = LINALG::Multiply(*kstmmod,false,*kstm,false,false);
  kstmmod->Complete(kstm->DomainMap(),kstm->RowMap());

  // ksti: multiply with linstickLM
  RCP<LINALG::SparseMatrix> kstimod = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);
  kstimod = LINALG::Multiply(*kstimod,false,*ksti,false,true);
  kstimod->Complete(ksti->DomainMap(),ksti->RowMap());

  // kstsl: multiply with linstickLM
  RCP<LINALG::SparseMatrix> kstslmod = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);
  kstslmod = LINALG::Multiply(*kstslmod,false,*kstsl,false,true);
  kstslmod->Complete(kstsl->DomainMap(),kstsl->RowMap());

  // kststmod: multiply with linstickLM
  RCP<LINALG::SparseMatrix> kststmod = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);
  kststmod = LINALG::Multiply(*kststmod,false,*kstst,false,true);
  kststmod->Complete(kstst->DomainMap(),kstst->RowMap());

  // blocks for complementary conditions (slip nodes) - from LM

  // ksln: multiply with linslipLM
  RCP<LINALG::SparseMatrix> kslnmod = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);
  kslnmod = LINALG::Multiply(*kslnmod,false,*ksln,false,true);
  kslnmod->Complete(ksln->DomainMap(),ksln->RowMap());

  // kslm: multiply with linslipLM
  RCP<LINALG::SparseMatrix> kslmmod = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);
  kslmmod = LINALG::Multiply(*kslmmod,false,*kslm,false,false);
  kslmmod->Complete(kslm->DomainMap(),kslm->RowMap());

  // ksli: multiply with linslipLM
  RCP<LINALG::SparseMatrix> kslimod = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);
  kslimod = LINALG::Multiply(*kslimod,false,*ksli,false,true);
  kslimod->Complete(ksli->DomainMap(),ksli->RowMap());

  // kslsl: multiply with linslipLM
  RCP<LINALG::SparseMatrix> kslslmod = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);
  kslslmod = LINALG::Multiply(*kslslmod,false,*kslsl,false,true);
  kslslmod->Complete(kslsl->DomainMap(),kslsl->RowMap());

  // slstmod: multiply with linslipLM
  RCP<LINALG::SparseMatrix> kslstmod = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);
  kslstmod = LINALG::Multiply(*kslstmod,false,*kslst,false,true);
  kslstmod->Complete(kslst->DomainMap(),kslst->RowMap());

  // fn: nothing to do

  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fi->Update(-alphaf_,*tempveci,1.0);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    fa->Update(-alphaf_,*tempvec,1.0);
  }

  // we want to split famod into 2 groups sl,st
  RCP<Epetra_Vector> fsl, fst;

  // do the vector splitting a -> sl+st
  if(gactivedofs_->NumGlobalElements())
  {
    if (!gstdofs->NumGlobalElements())
      fsl = rcp(new Epetra_Vector(*fa));
    else if (!gslipdofs_->NumGlobalElements())
      fst = rcp(new Epetra_Vector(*fa));
    else
    {
      LINALG::SplitVector(*fa,*gslipdofs_,fsl,*gstdofs,fst);
    }
  }

  // fm: add alphaf * old contact forces (t_n)
  RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mold_->Multiply(true,*zold_,*tempvecm);
  fm->Update(alphaf_,*tempvecm,1.0);

  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  // fst: mutliply with linstickLM
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> fstmod = rcp(new Epetra_Vector(*gstickt));
  RCP<LINALG::SparseMatrix> temp1 = LINALG::Multiply(*linstickLM_,false,*invdst,false,true);

  if(gstickdofs->NumGlobalElements())
  {
    temp1->Multiply(false,*fst,*fstmod);
  }

  // fsl: mutliply with linslipLM
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> fslmod = rcp(new Epetra_Vector(*gslipt_));
  RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*linslipLM_,false,*invdsl,false,true);

  if(gslipdofs_->NumGlobalElements())
  {
    temp->Multiply(false,*fsl,*fslmod);
  }

  // gactive: nothing to do

 #ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> deriv = rcp(new LINALG::SparseMatrix(*gactiven_,81));
    deriv->Add(*nmatrix_,false,1.0,1.0);
    deriv->Add(*smatrix_,false,1.0,1.0);
    deriv->Add(*nmhata,false,-1.0,1.0);
    deriv->Complete(*gsmdofs,*gactiven_);
    cout << *deriv << endl;
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives
  for (int i=0; i<(int)interface_.size();++i)
  {
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
  // FD check of stick condition
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      RCP<LINALG::SparseMatrix> deriv1 = rcp(new LINALG::SparseMatrix(*gactivet_,81));
      RCP<LINALG::SparseMatrix> deriv2 = rcp(new LINALG::SparseMatrix(*gactivet_,81));

      deriv1->Add(*linstickLM_,false,1.0,1.0);
      deriv1->Complete(*gsmdofs,*gactivet_);

      deriv2->Add(*linstickDIS_,false,1.0,1.0);
      deriv2->Complete(*gsmdofs,*gactivet_);

      cout << *deriv1 << endl;
      cout << *deriv2 << endl;

      interface_[i]->FDCheckStickDeriv();
    }
  }
#endif // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
  // FD check of tresca slip condition
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> deriv1 = rcp(new LINALG::SparseMatrix(*gactivet_,81));
    RCP<LINALG::SparseMatrix> deriv2 = rcp(new LINALG::SparseMatrix(*gactivet_,81));

    deriv1->Add(*linslipLM_,false,1.0,1.0);
    deriv1->Complete(*gsmdofs,*gslipt_);

    deriv2->Add(*linslipDIS_,false,1.0,1.0);
    deriv2->Complete(*gsmdofs,*gslipt_);

    cout << *deriv1 << endl;
    cout << *deriv2 << endl;

    interface_[i]->FDCheckSlipDeriv();
  }
}
#endif // #ifdef CONTACTFDSLIP

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81, true, false, kteff->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);

  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmimod,false,1.0,1.0);
  kteffnew->Add(*kmamod,false,1.0,1.0);

  // add i submatrices to kteffnew
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kin,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kim,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kii,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kia,false,1.0,1.0);

  // add matrices n and nmhata to kteffnew
  // this is only done for the "NO full linearization" case
  if (!fulllin)
  {
    if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmatrix_,false,1.0,1.0);
    if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmhata,false,-1.0,1.0);
  }

  // add full linearization terms to kteffnew
 if (fulllin)
  {
   if (gactiven_->NumGlobalElements()) kteffnew->Add(*smatrix_,false,-1.0,1.0);
  }

  // add terms of linearization of sick condition to kteffnew
  if (gstickt->NumGlobalElements()) kteffnew->Add(*linstickDIS_,false,-1.0,1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (gslipt_->NumGlobalElements())
  {
    kteffnew->Add(*linslipDIS_,false,-1.0,+1.0);

    RCP<Epetra_Vector> linslipRHSexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*linslipRHS_,*linslipRHSexp);
    feffnew->Update(-1.0,*linslipRHSexp,1.0);
  }

  // add terms of linearization feffnew
  // this is done also for evalutating the relative velocity with material
  // velocities
   if (gstickt->NumGlobalElements())
   {
      RCP<Epetra_Vector> linstickRHSexp = rcp(new Epetra_Vector(*problemrowmap_));
      LINALG::Export(*linstickRHS_,*linstickRHSexp);
      feffnew->Update(-1.0,*linstickRHSexp,1.0);
   }
   // add a submatrices to kteffnew
   if (gstickt->NumGlobalElements()) kteffnew->Add(*kstnmod,false,1.0,1.0);
   if (gstickt->NumGlobalElements()) kteffnew->Add(*kstmmod,false,1.0,1.0);
   if (gstickt->NumGlobalElements()) kteffnew->Add(*kstimod,false,1.0,1.0);
   if (gstickt->NumGlobalElements()) kteffnew->Add(*kstslmod,false,1.0,1.0);
   if (gstickt->NumGlobalElements()) kteffnew->Add(*kststmod,false,1.0,1.0);

  // add a submatrices to kteffnew
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslnmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslmmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslimod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslslmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslstmod,false,1.0,1.0);

  // FillComplete kteffnew (square)
  kteffnew->Complete();

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add i and sl subvector to feffnew
  RCP<Epetra_Vector> fiexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fi,*fiexp);
  if (gidofs->NumGlobalElements()) feffnew->Update(1.0,*fiexp,1.0);

  // add a subvector to feffnew
  RCP<Epetra_Vector> fstmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fstmod,*fstmodexp);
  if (gstickdofs->NumGlobalElements())feffnew->Update(1.0,*fstmodexp,+1.0);

  // add a subvector to feffnew
  RCP<Epetra_Vector> fslmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fslmod,*fslmodexp);
  if (gslipnodes_->NumGlobalElements())feffnew->Update(1.0,*fslmodexp,1.0);

  // add weighted gap vector to feffnew, if existing
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*gact,*gexp);
  if (gact->GlobalLength()) feffnew->Update(1.0,*gexp,1.0);

  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::EvaluateContact(RCP<LINALG::SparseMatrix> kteff,
                                                RCP<Epetra_Vector> feff)
{
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...

  // input parameters
  bool fulllin = Teuchos::getIntegralValue<int>(Params(),"FULL_LINEARIZATION");

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }

  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /* and global matrix s with normal derivatives of active nodes        */
  /**********************************************************************/
  // here and for the splitting later, we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleP(*pmatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
  }

  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);

  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofs,*gactiven_);

  // FillComplete() global matrix P
  // (actually gsdofrowmap_ is in general sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  pmatrix_->Complete(*gsmdofs,*gactivet_);

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofs,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofs,*gmdofrowmap_);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);

  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);

  // we want to split fsm into 2 groups s,m
  fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  if (!gsdofrowmap_->NumGlobalElements())
    *fm = *fsm;
  else if (!gmdofrowmap_->NumGlobalElements())
    *fs = *fsm;
  else
  {
    LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  }

  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;

  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*gactivedofs_));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    *fa = *fs;
  else if (!gactivedofs_->NumGlobalElements())
    *fi = *fs;
  else
  {
    LINALG::SplitVector(*fs,*gactivedofs_,fa,*gidofs,fi);
  }

  /**********************************************************************/
  /* Isolate active part from mhat, invd and dold                       */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> mhata;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);
  mhata_=mhata;

  RCP<LINALG::SparseMatrix> invda;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));

  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);

  /**********************************************************************/
  /* Split LinD and LinM into blocks                                    */
  /**********************************************************************/
  // we want to split lindmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> lindai, lindaa, lindam, lindas;

  // we want to split linmmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> linmmi, linmma, linmmm, linmms;

  if (fulllin)
  {
    // do the splitting
    LINALG::SplitMatrix2x2(lindmatrix_,gactivedofs_,gidofs,gmdofrowmap_,gsdofrowmap_,lindam,lindas,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(lindas,gactivedofs_,tempmap,gactivedofs_,gidofs,lindaa,lindai,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmmatrix_,gmdofrowmap_,tempmap,gmdofrowmap_,gsdofrowmap_,linmmm,linmms,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,linmma,linmmi,tempmtx1,tempmtx2);

    // modification of kai, kaa, kam
    // (this has to be done first as they are needed below)
    // (note, that kai, kaa, kam have to be UNcompleted again first!!!)
    kai->UnComplete();
    kaa->UnComplete();
    kam->UnComplete();
    kai->Add(*lindai,false,1.0-alphaf_,1.0);
    kaa->Add(*lindaa,false,1.0-alphaf_,1.0);
    kam->Add(*lindam,false,1.0-alphaf_,1.0);
    kai->Complete(*gidofs,*gactivedofs_);
    kaa->Complete();
    kam->Complete(*gmdofrowmap_,*gactivedofs_);
  }

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhata,true,*kan,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhata,true,*kam,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  if (fulllin) kmmmod->Add(*linmmm,false,1.0-alphaf_,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod = LINALG::Multiply(*mhata,true,*kai,false,false);
  kmimod->Add(*kmi,false,1.0,1.0);
  if (fulllin) kmimod->Add(*linmmi,false,1.0-alphaf_,1.0);
  kmimod->Complete(kmi->DomainMap(),kmi->RowMap());

  // kma: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod = LINALG::Multiply(*mhata,true,*kaa,false,false);
  kmamod->Add(*kma,false,1.0,1.0);
  if (fulllin) kmamod->Add(*linmma,false,1.0-alphaf_,1.0);
  kmamod->Complete(kma->DomainMap(),kma->RowMap());

  // kin: nothing to do

  // kim: nothing to do

  // kii: nothing to do

  // kia: nothing to do

  // n*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata,false,true);

  // nmatrix: nothing to do

  // kan: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kanmod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kanmod = LINALG::Multiply(*kanmod,false,*kan,false,true);

  // kam: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kammod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kammod = LINALG::Multiply(*kammod,false,*kam,false,true);

  // kai: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaimod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaimod = LINALG::Multiply(*kaimod,false,*kai,false,true);

  // kaa: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaamod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaamod = LINALG::Multiply(*kaamod,false,*kaa,false,true);

  // t*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> tmhata = LINALG::Multiply(*tmatrix_,false,*mhata,false,true);

  // fn: nothing to do

  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fi->Update(-alphaf_,*tempveci,1.0);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    fa->Update(-alphaf_,*tempvec,1.0);
  }

  // fm: add alphaf * old contact forces (t_n)

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (IsSelfContact())
  {
    RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
    RCP<Epetra_Vector> tempvecm2  = rcp(new Epetra_Vector(mold_->DomainMap()));
    RCP<Epetra_Vector> zoldexp  = rcp(new Epetra_Vector(mold_->RowMap()));
    if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    mold_->Multiply(true,*zoldexp,*tempvecm2);
    if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*tempvecm2,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->Multiply(true,*zold_,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }

  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  // fa: mutliply with tmatrix
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> famod = rcp(new Epetra_Vector(*gactivet_));
  RCP<LINALG::SparseMatrix> tinvda = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  tinvda->Multiply(false,*fa,*famod);

  // gactive: nothing to do

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    cout << *smatrix_ << endl;
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i=0; i<(int)interface_.size();++i)
  {
    cout << *pmatrix_ << endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81, true, false, kteff->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);

  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmimod,false,1.0,1.0);
  kteffnew->Add(*kmamod,false,1.0,1.0);

  // add i submatrices to kteffnew
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kin,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kim,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kii,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kia,false,1.0,1.0);

  // add matrices n and nmhata to kteffnew
  // this is only done for the "NO full linearization" case
  if (!fulllin)
  {
    if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmatrix_,false,1.0,1.0);
    if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmhata,false,-1.0,1.0);
  }

  // add full linearization terms to kteffnew
  if (fulllin)
  {
   if (gactiven_->NumGlobalElements()) kteffnew->Add(*smatrix_,false,-1.0,1.0);
   if (gactivet_->NumGlobalElements()) kteffnew->Add(*pmatrix_,false,-1.0,1.0);
  }

  // add a submatrices to kteffnew
  if (gactivet_->NumGlobalElements()) kteffnew->Add(*kanmod,false,1.0,1.0);
  if (gactivet_->NumGlobalElements()) kteffnew->Add(*kammod,false,1.0,1.0);
  if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaimod,false,1.0,1.0);
  if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaamod,false,1.0,1.0);

  // FillComplete kteffnew (square)
  kteffnew->Complete();

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fi,*fiexp);
  if (gidofs->NumGlobalElements()) feffnew->Update(1.0,*fiexp,1.0);

  // add weighted gap vector to feffnew, if existing
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*gact,*gexp);
  if (gact->GlobalLength()) feffnew->Update(1.0,*gexp,1.0);

  // add a subvector to feffnew
  RCP<Epetra_Vector> famodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*famod,*famodexp);
  if (gactivenodes_->NumGlobalElements())feffnew->Update(1.0,*famodexp,1.0);

  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;

  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::Recover(RCP<Epetra_Vector> disi)
{
  // extract slave displacements from disi
  RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
  if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);

  // extract other displacements from disi
  RCP<Epetra_Vector> disin = rcp(new Epetra_Vector(*gndofrowmap_));
  if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);

  //  // recover incremental jump (for active set)
  //  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  //  mhatmatrix_->Multiply(false,*disim,*incrjump_);
  //  incrjump_->Update(1.0,*disis,-1.0);
  //
  //  // friction
  //  // sum up incremental jumps from active set nodes
  //  jump_->Update(1.0,*incrjump_,1.0);
  //  // friction
  //  // store updaded jumps to nodes
  //  StoreNodalQuantities(ManagerBase::jump);

  /**********************************************************************/
  /* Update Lagrange multipliers z_n+1                                  */
  /**********************************************************************/

  // for self contact, slave and master sets may have changed,
  // thus we have to export the products Dold * zold and Mold^T * zold to fit
  if (IsSelfContact())
  {
    // approximate update
    //z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    //invd_->Multiply(false,*fs_,*z_);

    // full update
    z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    z_->Update(1.0,*fs_,0.0);
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gsdofrowmap_));
    kss_->Multiply(false,*disis,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksm_->Multiply(false,*disim,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksn_->Multiply(false,*disin,*mod);
    z_->Update(-1.0,*mod,1.0);
    RCP<Epetra_Vector> mod2 = rcp(new Epetra_Vector((dold_->RowMap())));
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*mod2);
    RCP<Epetra_Vector> mod3 = rcp(new Epetra_Vector((dold_->RowMap())));
    dold_->Multiply(false,*mod2,*mod3);
    RCP<Epetra_Vector> mod4 = rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*mod3,*mod4);
    z_->Update(-alphaf_,*mod4,1.0);
    RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
    invd_->Multiply(false,*zcopy,*z_);
    z_->Scale(1/(1-alphaf_));
  }
  else
  {
    // approximate update
    //invd_->Multiply(false,*fs_,*z_);

    // full update
    z_->Update(1.0,*fs_,0.0);
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gsdofrowmap_));
    kss_->Multiply(false,*disis,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksm_->Multiply(false,*disim,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksn_->Multiply(false,*disin,*mod);
    z_->Update(-1.0,*mod,1.0);
    dold_->Multiply(false,*zold_,*mod);
    z_->Update(-alphaf_,*mod,1.0);
    RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
    invd_->Multiply(false,*zcopy,*z_);
    z_->Scale(1/(1-alphaf_));
  }

  // store updated LM into nodes
  StoreNodalQuantities(AbstractStrategy::lmupdate);

  /*
   // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
   #ifdef DEBUG
   //debugging (check for z_i = 0)
   RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
   if (gidofs->NumGlobalElements())
   {
   RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
   LINALG::Export(*z_,*zinactive);
   cout << *zinactive << endl;
   }

   bool fulllin = Teuchos::getIntegralValue<int>(Params(),"FULL_LINEARIZATION");

   // debugging (check for N*[d_a] = g_a and T*z_a = 0)
   if (gactivedofs_->NumGlobalElements())
   {
   RCP<Epetra_Vector> activejump = rcp(new Epetra_Vector(*gactivedofs_));
   RCP<Epetra_Vector> gtest = rcp(new Epetra_Vector(*gactiven_));
   RCP<Epetra_Vector> gtest2 = rcp(new Epetra_Vector(*gactiven_));
   LINALG::Export(*incrjump_,*activejump);
   nmatrix_->Multiply(false,*activejump,*gtest);

   RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
   RCP<Epetra_Vector> disism = rcp(new Epetra_Vector(*gsmdofs));
   LINALG::Export(*disi,*disism);
   if (fulllin)
   {
   smatrix_->Multiply(false,*disism,*gtest2);
   gtest->Update(1.0,*gtest2,1.0);
   }
   cout << *gtest << endl << *g_ << endl;

   RCP<Epetra_Vector> zactive = rcp(new Epetra_Vector(*gactivedofs_));
   RCP<Epetra_Vector> zerotest = rcp(new Epetra_Vector(*gactivet_));
   RCP<Epetra_Vector> zerotest2 = rcp(new Epetra_Vector(*gactivet_));
   LINALG::Export(*z_,*zactive);
   tmatrix_->Multiply(false,*zactive,*zerotest);
   if (fulllin)
   {
   pmatrix_->Multiply(false,*disis,*zerotest2);
   zerotest->Update(1.0,*zerotest2,1.0);
   }
   cout << *zerotest << endl;
   }
   #endif // #ifdef DEBUG
   // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
   */

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence               popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::UpdateActiveSet()
{
  // get input parameter ctype
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(Params(),"CONTACT");
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

  // assume that active set has converged and check for opposite
  activesetconv_=true;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);

      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }

      // friction

      double tz = 0.0;
      double tjump = 0.0;

      if(ftype == INPAR::CONTACT::friction_tresca || ftype == INPAR::CONTACT::friction_coulomb)
      {
        // compute tangential part of Lagrange multiplier
        tz = cnode->txi()[0]*cnode->lm()[0] + cnode->txi()[1]*cnode->lm()[1];

        // compute tangential part of jump
        tjump = cnode->txi()[0]*cnode->jump()[0] + cnode->txi()[1]*cnode->jump()[1];
      }

      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << endl;

        // check for penetration
        if (wgap < 0)
        {
          cnode->Active() = true;
          activesetconv_ = false;
#ifdef CONTACTFRICTIONLESSFIRST
       if (cnode->ActiveOld()==false) cnode->Slip() = true;
#endif
        }
      }

      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << endl;

        // check for tensile contact forces
        if (nz <= 0) // no averaging of Lagrange multipliers
        //if (0.5*nz+0.5*nzold <= 0) // averaging of Lagrange multipliers
        {
          if (ctype != INPAR::CONTACT::contact_meshtying)
          {
            cnode->Active() = false;

            // friction
            if(ftype == INPAR::CONTACT::friction_tresca || ftype == INPAR::CONTACT::friction_coulomb)
            {
              cnode->Slip() = false;
            }
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        }

        // friction
        else
        {
          // friction tresca
          if(ftype == INPAR::CONTACT::friction_tresca)
          {
            double frbound = Params().get<double>("FRBOUND");
            double ct = Params().get<double>("SEMI_SMOOTH_CT");

            if(cnode->Slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if(abs(tz+ct*tjump)-frbound > 0) {}
                // do nothing (slip was correct)
              else
              {
#ifdef CONTACTFRICTIONLESSFIRST
                if(cnode->ActiveOld()==false)
                {}
                else
                {
                 cnode->Slip() = false;
                 activesetconv_ = false;
                }
#else
                cnode->Slip() = false;
                activesetconv_ = false;
#endif
              }
            }
          } // if(ftype == INPAR::CONTACT::friction_tresca)

          // friction coulomb
          if(ftype == INPAR::CONTACT::friction_coulomb)
          {
            double frcoeff = Params().get<double>("FRCOEFF");
            double ct = Params().get<double>("SEMI_SMOOTH_CT");

            if(cnode->Slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frcoeff*nz <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if(abs(tz+ct*tjump)-frcoeff*nz > 0) {}
                // do nothing (slip was correct)
              else
              {
#ifdef CONTACTFRICTIONLESSFIRST
                if(cnode->ActiveOld()==false)
                {}
                else
                {
                 cnode->Slip() = false;
                 activesetconv_ = false;
                }
#else
                cnode->Slip() = false;
                activesetconv_ = false;
#endif
              }
            }
          } // if(ftype == INPAR::CONTACT::friction_coulomb)
        } // if (nz <= 0)
      } // if (cnode->Active()==false)
    } // loop over all slave nodes
  } // loop over all interfaces

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }

  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=null) zigzagthree_  = rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=null) zigzagtwo_    = rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=null) zigzagone_ = rcp(new Epetra_Map(*gactivenodes_));

  // update zig-zagging history for slip nodes (shift by one)
  if (zigzagsliptwo_!=null) zigzagslipthree_  = rcp(new Epetra_Map(*zigzagsliptwo_));
  if (zigzagslipone_!=null) zigzagsliptwo_    = rcp(new Epetra_Map(*zigzagslipone_));
  if (gslipnodes_!=null) zigzagslipone_ = rcp(new Epetra_Map(*gslipnodes_));


  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
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
  if(ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)
  {
    if (ActiveSetSteps()>2)
    {
      if (zigzagtwo_!=null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_) and zigzagsliptwo_->SameAs(*gslipnodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (Comm().MyPID()==0)
            cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << endl;
        }
      }

      if (zigzagthree_!=null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_) and zigzagslipthree_->SameAs(*gslipnodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (Comm().MyPID()==0)
            cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << endl;
        }
      }
    }
  } // if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)

  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = null;
    zigzagtwo_  = null;
    zigzagthree_= null;
  }

  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET ITERATION " << ActiveSetSteps()-1
         << " NOT CONVERGED - REPEAT TIME STEP................." << endl;
  else if (Comm().MyPID()==0 && activesetconv_==true)
    cout << "ACTIVE SET CONVERGED IN " << ActiveSetSteps()-zigzagging
         << " STEP(S)................." << endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::LagrangeStrategy::UpdateActiveSetSemiSmooth()
{
  // FIXME: Here we do not consider zig-zagging yet!

  // get input parameter ctype
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(Params(),"CONTACT");
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weigh the two violations via cn!
  double cn = Params().get<double>("SEMI_SMOOTH_CN");

  // assume that active set has converged and check for opposite
  activesetconv_=true;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);

      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }

      // friction
      double ct = Params().get<double>("SEMI_SMOOTH_CT");
      vector<double> tz (Dim()-1,0);
      vector<double> tjump (Dim()-1,0);
      double euclidean = 0.0;

      if(ftype == INPAR::CONTACT::friction_tresca || ftype == INPAR::CONTACT::friction_coulomb)
       {
       	// compute tangential parts and of Lagrange multiplier and incremental jumps
       	for (int i=0;i<Dim();++i)
       	{
       		tz[0] += cnode->txi()[i]*cnode->lm()[i];
           if(Dim()==3) tz[1] += cnode->teta()[i]*cnode->lm()[i];

           tjump[0] += cnode->txi()[i]*cnode->jump()[i];
           if(Dim()==3) tjump[1] += cnode->teta()[i]*cnode->jump()[i];
       	}

       	// evaluate euclidean norm |tz+ct.tjump|
       	vector<double> sum (Dim()-1,0);
       	sum[0] = tz[0]+ct*tjump[0];
       	if (Dim()==3) sum[1] = tz[1]+ct*tjump[1];
       	if (Dim()==2) euclidean = abs(sum[0]);
       	if (Dim()==3) euclidean = sqrt(sum[0]*sum[0]+sum[1]*sum[1]);
       }

      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << endl;

        // check for penetration and/or tensile contact forces
        if (nz - cn*wgap > 0)
        {
          cnode->Active() = true;

          //friction
          if(ftype == INPAR::CONTACT::friction_tresca || ftype == INPAR::CONTACT::friction_coulomb)
          {
          	// nodes coming into contact
            cnode->Slip() = true;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) cnode->Slip() = true;
#endif
          }

          activesetconv_ = false;
        }
      }

      // check nodes of active set ***************************************
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << endl;

        // check for tensile contact forces and/or penetration
        if (nz - cn*wgap <= 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
        {
          if (ctype != INPAR::CONTACT::contact_meshtying)
          {
            cnode->Active() = false;

            // friction
            if(ftype == INPAR::CONTACT::friction_tresca || ftype == INPAR::CONTACT::friction_coulomb)
            {
              cnode->Slip() = false;
            }
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        }

        // friction
        else
        {
          // friction tresca
          if(ftype == INPAR::CONTACT::friction_tresca)
          {
            double frbound = Params().get<double>("FRBOUND");

            if(cnode->Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if(euclidean-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (euclidean)-frbound > 0
              if(euclidean-frbound > 0) {}
               // do nothing (slip was correct)
              else
              {
#ifdef CONTACTFRICTIONLESSFIRST
                if(cnode->ActiveOld()==false)
                {}
                else
                {
                 cnode->Slip() = false;
                 activesetconv_ = false;
                }
#else
                cnode->Slip() = false;
                activesetconv_ = false;
#endif
              }
            }
          } // if (fytpe=="tresca")

          // friction coulomb
          if(ftype == INPAR::CONTACT::friction_coulomb)
          {
            double frcoeff = Params().get<double>("FRCOEFF");
            if(cnode->Slip() == false)
            {
              // check (euclidean)-frbound <= 0
#ifdef CONTACTCOMPHUEBER
            if(euclidean-frcoeff*(nz-cn*wgap) <= 0) {}
#else
            if(euclidean-frcoeff*nz <= 0) {}
#endif
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (euclidean)-frbound > 0
#ifdef CONTACTCOMPHUEBER
              if(euclidean-frcoeff*(nz-cn*wgap) > 0) {}
#else
              if(euclidean-frcoeff*nz > 0) {}
#endif
              // do nothing (slip was correct)
              else
              {
#ifdef CONTACTFRICTIONLESSFIRST
                if(cnode->ActiveOld()==false)
                {}
                else
                {
                 cnode->Slip() = false;
                 activesetconv_ = false;
                }
#else
                cnode->Slip() = false;
                activesetconv_ = false;
#endif
              }
            }
          } // if(ftype == INPAR::CONTACT::friction_coulomb)
        } // if (nz - cn*wgap <= 0)
      } // if (cnode->Active()==false)
    } // loop over all slave nodes
  } // loop over all interfaces

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
  }

  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET HAS CHANGED... CHANGE No. " << ActiveSetSteps()-1 << endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;

  return;
}

#endif // CCADISCRET
