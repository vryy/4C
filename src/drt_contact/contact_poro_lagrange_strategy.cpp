/*!----------------------------------------------------------------------
\file contact_poro_lagrange_strategy.cpp

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
Maintainer: Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "contact_poro_lagrange_strategy.H"
#include "contact_interface.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              ager 08/14|
 *----------------------------------------------------------------------*/
CONTACT::PoroLagrangeStrategy::PoroLagrangeStrategy(DRT::Discretization& probdiscret,
                                                Teuchos::ParameterList params,
                                                std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
                                                int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof):
CoLagrangeStrategy(probdiscret,params,interface,dim,comm,alphaf,maxdof),
no_penetration_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 | Activate No-Penetration for the contact surface (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SetupNoPenetrationCondition()
{
  no_penetration_ = true;
  lambda_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
}

/*----------------------------------------------------------------------*
 | initialize global poro contact variables                             |
 |                            for next Newton step (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroInitialize(ADAPTER::Coupling& coupfs, Teuchos::RCP<const Epetra_Map> fluiddofs)
{
    //  (1)                                                          //
   //      Get required fluid maps from structural maps             //
  //
  fgsdofrowmap_ = coupfs.MasterToSlaveMap(gsdofrowmap_);
  fgmdofrowmap_ = coupfs.MasterToSlaveMap(gmdofrowmap_);
  fgsmdofrowmap_ = coupfs.MasterToSlaveMap(gsmdofrowmap_);
  fgndofrowmap_  = LINALG::SplitMap(*fluiddofs,*fgsmdofrowmap_); //Not equal to transforming gndofrowmap_ (pressure dofs missing!)
  fgactivedofs_ = coupfs.MasterToSlaveMap(gactivedofs_);
  falldofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(*fluiddofs));
  fgactiven_ = coupfs.MasterToSlaveMap(gactiven_);
  fgactivet_ = coupfs.MasterToSlaveMap(gactivet_);

    //  (2)                                                          //
   //      Initialize Matrices                                      //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  {
    // (re)setup global nCoup Vector
    NCoup_ = Teuchos::rcp(new Epetra_Vector(*gactiven_));

    // (re)setup global linearisation matrices of nCoup
    NCoup_lindisp_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_, 10));
    NCoup_linvel_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_, 10));

    // (re)setup global tangential and lin(tangential)*lambda matrices
    Tangential_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 10));
    linTangentiallambda_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 10));

    // (re)setup global lin of D & M * lambda - Matrix
    porolindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
    porolinmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  }

    //  (3)                                                          //
   //      Assemble Matrices                                        //
  //                                                               //
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
    {
     interface_[i]->AssembleNCoup(*NCoup_);

     interface_[i]->AssembleNCoupLin(*NCoup_lindisp_);
     interface_[i]->AssembleNCoupLin(*NCoup_linvel_,true);

     interface_[i]->AssembleT(*Tangential_);
     interface_[i]->AssembleP(*linTangentiallambda_,true); //use lambda(n +1) for tangential condition!!!

     interface_[i]->AssembleLinDM(*porolindmatrix_,*porolinmmatrix_);
    }
  }

    //  (4)                                                          //
   //      Complete Matrices                                        //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  {
    NCoup_lindisp_->Complete(*gsmdofrowmap_, *gactiven_);
    NCoup_linvel_->Complete(*gsdofrowmap_, *gactiven_);

    Tangential_->Complete(*gactivedofs_,*gactivet_ );
    linTangentiallambda_->Complete(*gsmdofrowmap_, *gactivet_);

    porolindmatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
    porolinmmatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);
  }

    //  (5)                                                          //
   //      Reset Matrix Transform Objects                           //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  {
  linncoupveltransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  linncoupdisptransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  tanginvtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  lintangentlambdatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  porolindmatrixtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  mhataamtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  dhattransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  doldtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  moldtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  invDatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  }

    //  (6)                                                          //
   //      Transform Matrices from structural dofs to fluid dofs    //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
    //eventually a coupling object just on the mortar interface would make sense!!! ChrAg
  {
    // transform matrices coming from contact to fluid maps, as they are all in structure maps!
    //
    // A generell problem here is that we would need to update coupling objects in everey newton step if the active
    // set changes. To avoid this, a 'bigger' coupling object is used - but therefore now the Row-Maps of the
    // created Sparse Matrixes are to big! --- Leads to problems for Matrix - Multiplications where the Row - Map is used!
    //
    // At the moment this is solved by using SplitMatrix2x2 to cut out just the relevant part of the matrix, but this shouldn't
    // be the final solution!
    //
    //************************************************************************************************
    //
    Teuchos::RCP<Epetra_Vector> tmpfullncoup = Teuchos::rcp(new Epetra_Vector(*coupfs.MasterDofMap()));
    LINALG::Export(*NCoup_,*tmpfullncoup);
    tmpfullncoup = coupfs.MasterToSlave(tmpfullncoup);
    fNCoup_ = Teuchos::rcp(new Epetra_Vector(*fgactiven_));
    LINALG::Export(*tmpfullncoup,*fNCoup_);
    //
    //************************************************************************************************
    //
    fdoldtransp_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*doldtransform_)(*dold_->Transpose(),
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *fdoldtransp_,
                           false);
    fdoldtransp_->Complete(dold_->DomainMap(),*fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    fmoldtransp_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*moldtransform_)(*mold_->Transpose(),
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *fmoldtransp_,
                           false);
    fmoldtransp_->Complete(mold_->DomainMap(),*fgmdofrowmap_);
    //
    //************************************************************************************************
    //
    fporolindmatrix_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*porolindmatrixtransform_)(*porolindmatrix_,
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *fporolindmatrix_,
                           false);
    fporolindmatrix_->Complete(porolindmatrix_->DomainMap(),*fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    //porolinmmatrixtransform_ is missing as atm just onesided porocontact is done!!!
    //
    //************************************************************************************************
    //
    fmhataam_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*mhataamtransform_)(*mhataam_,
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *fmhataam_,
                           false);
    fmhataam_->Complete(mhataam_->DomainMap(),*fgactivedofs_);
    //
    //************************************************************************************************
    //
    fdhat_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*dhattransform_)(*dhat_,
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *fdhat_,
                           false);
    fdhat_->Complete(dhat_->DomainMap(),*fgactivedofs_);
    //
    //************************************************************************************************
    //
    if (gactivedofs_->NumGlobalElements())
    {
      Teuchos::RCP<LINALG::SparseMatrix> tanginvD  = LINALG::MLMultiply(*Tangential_,false,*invda_,true,false,false,true);
      Teuchos::RCP<LINALG::SparseMatrix> tmpftanginvD = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
      (*tanginvtransform_)(*tanginvD,
                          1.0,
                          ADAPTER::CouplingMasterConverter(coupfs),
                          ADAPTER::CouplingMasterConverter(coupfs),
                          *tmpftanginvD,
                          false,
                          false);
     tmpftanginvD->Complete(*fgactivedofs_,*falldofrowmap_);
      //better solution to get maps as wanted? -- for this matrix map as important as there will be a matrix-matrix multiplication

     Teuchos::RCP<Epetra_Map> restfgactivet_, restfgactivedofs;
     Teuchos::RCP<LINALG::SparseMatrix> tmpm1, tmpm2, tmpm3;

      //This should just be a temporary solution to change the row map of the matrix ...
      LINALG::SplitMatrix2x2(tmpftanginvD,fgactivet_,restfgactivet_,fgactivedofs_,restfgactivedofs,ftanginvD_,tmpm1,tmpm2,tmpm3);

#if(0)
//Some solutions that do not work!

    //    ftanginvD_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*fgactivet_));
    //  ftanginvD_->Assign(View,*tmpftanginvD);

   //   ftanginvD_->Add(*tmpftanginvD,false,1.0,1.0);
   //   ftanginvD_->Complete(*fgactivedofs_,*fgactivet_);

      //  int err = ftanginvD_->EpetraMatrix()->ReplaceRowMap(*fgactivet_);
       // if (err) dserror("RRM STOPPED .... with err! %d",err);
#endif
    //
    //************************************************************************************************
    //
      fNCoup_linvel_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_,108,false));

      (*linncoupveltransform_)(*NCoup_linvel_,
                              1.0,
                              ADAPTER::CouplingMasterConverter(coupfs),
                              ADAPTER::CouplingMasterConverter(coupfs),
                              *fNCoup_linvel_,
                              false,
                              false);
      fNCoup_linvel_->Complete(*falldofrowmap_,*fgactiven_);
    //
    //************************************************************************************************
    //
      fNCoup_lindisp_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_,81,false));
      (*linncoupdisptransform_)(*NCoup_lindisp_,
                              1.0,
                              ADAPTER::CouplingMasterConverter(coupfs),
                              *fNCoup_lindisp_,
                              false);
      fNCoup_lindisp_->Complete(*ProblemDofs(),*fgactiven_);
    //
    //************************************************************************************************
    //
      flinTangentiallambda_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_,81,false));
      (*lintangentlambdatransform_)(*linTangentiallambda_,
                                    1.0,
                                    ADAPTER::CouplingMasterConverter(coupfs),
                                    *flinTangentiallambda_,
                                    false);
      flinTangentiallambda_->Complete(*gsdofrowmap_,*fgactivedofs_);

#if(0) //just in case
     // only for parallel redistribution case
//     if (parredist) //care about that if its in contact ... ChrAg Todo
//     {
//       NCoup_lindisp_    = MORTAR::MatrixRowColTransform(NCoup_lindisp_, interface_[i]->NormalDofs(), masterdofrowmap_);
//       NCoup_linvel_     = MORTAR::MatrixRowColTransform(NCoup_linvel_, interface_[i]->NormalDofs(), slavedofrowmap_);
//     }

     // only for parallel redistribution case
//     if (parredist) //care about that if its in contact ... ChrAg Todo
//     {
//       Tangential_    = MORTAR::MatrixRowColTransform(Tangential_, interface_[i]->TangentialDofs(), slavedofrowmap_);
//       linTangentiallambda_     = MORTAR::MatrixRowColTransform(linTangentiallambda_, interface_[i]->TangentialDofs(), slavedofrowmap_);
//     }


   //  Assemble lin of D & M - Matrix --
   //  Teuchos::RCP<LINALG::SparseMatrix> lindglobal = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX)); //ChrAg
   //  Teuchos::RCP<LINALG::SparseMatrix> linmglobal = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX)); //ChrAg






     //porolindmatrix_ = lindglobal;
     //porolinmmatrix_ = linmglobal;
#endif
  //
  //************************************************************************************************
  //
    finvda_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*falldofrowmap_));
    (*invDatransform_)(*invda_,
                           1.0,
                           ADAPTER::CouplingMasterConverter(coupfs),
                           *finvda_,
                           false);
    finvda_->Complete(invda_->DomainMap(),*fgsdofrowmap_);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices (public)          ager 07/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluatePoroContact(Teuchos::RCP<LINALG::SparseMatrix>& kteff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
      kteff->Complete();

  // system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    /**********************************************************************/
    /* (3) Split kteff into 3x3 matrix blocks                             */
    /**********************************************************************/

    // we want to split k into 3 groups s,m,n = 9 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

    // temporarily we need the blocks ksmsm, ksmn, knsm
    // (FIXME: because a direct SplitMatrix3x3 is still missing!)
    Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Epetra_Map> ftempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

    // split into slave/master part + structure part
    Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);

    if (ParRedist())  //TODO Check if how to modifiy
    {
      dserror("ParRedist(): CHECK ME!");
      // split and transform to redistributed maps
      LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
      ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
      ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
      knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
    }
    else
    {
      // only split, no need to transform
      LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,fgsmdofrowmap_,fgndofrowmap_,ksmsm,ksmn,knsm,knn);
    }

    // further splits into slave part + master part
    LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,fgsdofrowmap_,fgmdofrowmap_,kss,ksm,kms,kmm);
    LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,fgndofrowmap_,ftempmap,ksn,tempmtx1,kmn,tempmtx2);
    LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,fgsdofrowmap_,fgmdofrowmap_,kns,knm,tempmtx1,tempmtx2);


    // abbreviations for slave  and master set
    int sset = fgsdofrowmap_->NumGlobalElements();
    int mset = fgmdofrowmap_->NumGlobalElements();  //eval if there is an mset on the fluid side -- not for onsided poro contact!


    // store some stuff for static condensation of LM
    csfsn_  = ksn;
    csfsm_  = ksm;
    csfss_  = kss;

    /**********************************************************************/
    /* (5) Split slave quantities into active / inactive                  */
    /**********************************************************************/

    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;
    Teuchos::RCP<Epetra_Map> fgidofs;

    // do the splitting
    LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,fgactivedofs_,fgidofs,kaa,kai,kia,kii);
    LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,fgndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
    LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,fgmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
    LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,fgactivedofs_,fgidofs,kma,kmi,tempmtx1,tempmtx2);

    // abbreviations for master, active and inactive set
    int aset = fgactivedofs_->NumGlobalElements();
    int iset = fgidofs->NumGlobalElements();

    /**********************************************************************/
    /* (7) Build the final K blocks                                       */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // knn: nothing to do

    // knm: nothing to do

    // kns: nothing to do

    //---------------------------------------------------------- SECOND LINE
    // kmn: add T(mhataam)*kan
    Teuchos::RCP<LINALG::SparseMatrix> kmnmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmnmod->Add(*kmn,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhataam_,true,*kan,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

    // kmm: add T(mhataam)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kmmmod;
    if (mset)
    {
      kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmmmod->Add(*kmm,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhataam_,true,*kam,false,false,false,true);
      kmmmod->Add(*kmmadd,false,1.0,1.0);
      kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
    }

    // kmi: add T(mhataam)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmimod->Add(*kmi,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhataam_,true,*kai,false,false,false,true);
      kmimod->Add(*kmiadd,false,1.0,1.0);
      kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
    }

    // kma: add T(mhataam)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmamod->Add(*kma,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhataam_,true,*kaa,false,false,false,true);
      kmamod->Add(*kmaadd,false,1.0,1.0);
      kmamod->Complete(kma->DomainMap(),kma->RowMap());
    }

    //----------------------------------------------------------- THIRD LINE
    //------------------- FOR 3D QUADRATIC CASE ----------------------------

    //--- For using non diagonal D-Matrix, it should be checked if this assumtion isn't somewhere else!!!

    // kin: subtract T(dhat)*kan --
    Teuchos::RCP<LINALG::SparseMatrix> kinmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kinmod->Add(*kin,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kinadd = LINALG::MLMultiply(*dhat_,true,*kan,false,false,false,true);
    kinmod->Add(*kinadd,false,-1.0,1.0);
    kinmod->Complete(kin->DomainMap(),kin->RowMap());

    // kim: subtract T(dhat)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kimmod;
    if (mset)
    {
      kimmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kimmod->Add(*kim,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kimadd = LINALG::MLMultiply(*dhat_,true,*kam,false,false,false,true);
      kimmod->Add(*kimadd,false,-1.0,1.0);
      kimmod->Complete(kim->DomainMap(),kim->RowMap());
    }

    // kii: subtract T(dhat)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kiimod;
    if (iset)
    {
      kiimod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiimod->Add(*kii,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kiiadd = LINALG::MLMultiply(*dhat_,true,*kai,false,false,false,true);
      kiimod->Add(*kiiadd,false,-1.0,1.0);
      kiimod->Complete(kii->DomainMap(),kii->RowMap());
    }
    // kia: subtract T(dhat)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kiamod;
    if (iset && aset)
    {
      kiamod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiamod->Add(*kia,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kiaadd = LINALG::MLMultiply(*dhat_,true,*kaa,false,false,false,true);
      kiamod->Add(*kiaadd,false,-1.0,1.0);
      kiamod->Complete(kia->DomainMap(),kia->RowMap());
    }

    //---------------------------------------------------------- FOURTH LINE
    // nothing to do

    //----------------------------------------------------------- FIFTH LINE
    // kan: multiply tmatrix with invda and kan
    Teuchos::RCP<LINALG::SparseMatrix> kanmod;
    if (aset)
    {
      kanmod = LINALG::MLMultiply(*tmatrix_,false,*invda_,true,false,false,true);
      kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
    }
    // kam: multiply tmatrix with invda and kam
    Teuchos::RCP<LINALG::SparseMatrix> kammod;
    if (aset && mset)
    {
      kammod = LINALG::MLMultiply(*tmatrix_,false,*invda_,true,false,false,true);
      kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
    }
    // kai: multiply tmatrix with invda and kai
    Teuchos::RCP<LINALG::SparseMatrix> kaimod;
    if (aset && iset)
    {
      kaimod = LINALG::MLMultiply(*tmatrix_,false,*invda_,true,false,false,true);
      kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
    }
    // kaa: multiply tmatrix with invda and kaa
    Teuchos::RCP<LINALG::SparseMatrix> kaamod;
    if (aset)
    {
      kaamod = LINALG::MLMultiply(*tmatrix_,false,*invda_,true,false,false,true);
      kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
    }
    /********************************************************************/
    /* (9) Transform the final K blocks                                 */
    /********************************************************************/
    // The row maps of all individual matrix blocks are transformed to
    // the parallel layout of the underlying problem discretization.
    // Of course, this is only necessary in the parallel redistribution
    // case, where the contact interfaces have been redistributed
    // independently of the underlying problem discretization.

    if (ParRedist()) //check what to do
    {
      //----------------------------------------------------------- FIRST LINE
      // nothing to do (ndof-map independent of redistribution)

      //---------------------------------------------------------- SECOND LINE
      kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
      kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
      if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
      if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);

      //----------------------------------------------------------- THIRD LINE
      if (iset)
      {
        kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
        kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
        kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
        if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
      }

      //---------------------------------------------------------- FOURTH LINE
      if (aset) smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);

      //----------------------------------------------------------- FIFTH LINE
      if (aset)
      {
        kanmod = MORTAR::MatrixRowTransform(kanmod,pgsdofrowmap_);
        kammod = MORTAR::MatrixRowTransform(kammod,pgsdofrowmap_);
        kaamod = MORTAR::MatrixRowTransform(kaamod,pgsdofrowmap_);
        if (iset) kaimod = MORTAR::MatrixRowTransform(kaimod,pgsdofrowmap_);
        pmatrix_ = MORTAR::MatrixRowTransform(pmatrix_,pgsdofrowmap_);
      }
    }

    /**********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                  */
    /**********************************************************************/

    Teuchos::RCP<LINALG::SparseMatrix> kteffnew = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(),81,true,false,kteffmatrix->GetMatrixtype()));
   // Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*ProblemDofs());

    //----------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->Add(*knn,false,1.0,1.0);
    if (mset) kteffnew->Add(*knm,false,1.0,1.0);
    if (sset) kteffnew->Add(*kns,false,1.0,1.0);
    //---------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->Add(*kmnmod,false,1.0,1.0);
    if (mset) kteffnew->Add(*kmmmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kmimod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);
    //----------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->Add(*kinmod,false,1.0,1.0);
    if (iset && mset) kteffnew->Add(*kimmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kiimod,false,1.0,1.0);
    if (iset && aset) kteffnew->Add(*kiamod,false,1.0,1.0);

    //---------------------------------------------------------- FOURTH LINE
    // add a submatrices to kteffnew
    //if (aset) kteffnew->Add(*smatrix_,false,1.0,1.0);

    //----------------------------------------------------------- FIFTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->Add(*kanmod,false,1.0,1.0);
    if (aset && mset) kteffnew->Add(*kammod,false,1.0,1.0);
    if (aset && iset) kteffnew->Add(*kaimod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kaamod,false,1.0,1.0);
   // if (aset) kteffnew->Add(*pmatrix_,false,-1.0,1.0);



    // FillComplete kteffnew (square)
    kteffnew->Complete(*falldofrowmap_,*gdisprowmap_);



    // finally do the replacement
    kteff = kteffnew;
  }
  else
  {
    dserror("Trying to use not condensed PoroContact --- Feel Free to implement!");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration          |
 |  condition on contact surface (public)                     ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluatePoroNoPenContact(Teuchos::RCP<LINALG::SparseMatrix>& k_fseff,
                                                             Teuchos::RCP<LINALG::SparseMatrix>& Feff,
                                                             Teuchos::RCP<Epetra_Vector>& feff)
{

    // check if contact contributions are present,
    // if not we can skip this routine to speed things up
    if (!no_penetration_ || (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())) return;

    // complete stiffness matrix
    // (this is a prerequisite for the Split2x2 methods to be called later)
    k_fseff->Complete();
    Feff->Complete();

       // shape function
    INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");

    //**********************************************************************
    //**********************************************************************
    // CASE A: CONDENSED SYSTEM (DUAL)
    //**********************************************************************
    //**********************************************************************

    // double-check if this is a dual LM system
      if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
        dserror("Condensation only for dual LM");

      /**********************************************************************/
      /* (2) Add contact stiffness terms to kteff                           */
      /**********************************************************************/

      // transform if necessary
      if (ParRedist())
      {
        dserror("CHECK ME!!!");
        lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
        linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
      }

      k_fseff->UnComplete();
      k_fseff->Add(*fporolindmatrix_,false,(1.0-alphaf_)*1.0,1.0);
      //k_fseff->Add(*fporolinmmatrix_,false,1.0-alphaf_,1.0); will be needed for twosided poro contact
      k_fseff->Complete(*ProblemDofs(), *falldofrowmap_); //gets bigger because of linearisation w.r.t. to pure structural displacements!
      /**********************************************************************/
      /* (3) Split k_fseff and Feff into 3x3 matrix blocks                             */
      /**********************************************************************/

      // we want to split k into 3 groups s,m,n = 9 blocks
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_ss, k_fs_sm, k_fs_sn, k_fs_ms, k_fs_mm, k_fs_mn, k_fs_ns, k_fs_nm, k_fs_nn;
      Teuchos::RCP<LINALG::SparseMatrix> F_ss, F_sm, F_sn, F_ms, F_mm, F_mn, F_ns, F_nm, F_nn;

      // temporarily we need the blocks ksmsm, ksmn, knsm
      // (FIXME: because a direct SplitMatrix3x3 is still missing!)
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_smsm, k_fs_smn, k_fs_nsm;
      Teuchos::RCP<LINALG::SparseMatrix> F_smsm, F_smn, F_nsm;

      // some temporary Teuchos::RCPs
      Teuchos::RCP<Epetra_Map> tempmap;
      Teuchos::RCP<Epetra_Map> ftempmap;
      Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
      Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
      Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

      // split into slave/master part + structure part
      Teuchos::RCP<LINALG::SparseMatrix> k_fseffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fseff);
      Teuchos::RCP<LINALG::SparseMatrix> Feffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(Feff);
      if (ParRedist())
      {
        dserror("CHECK ME!");
        // split and transform to redistributed maps
//        LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
//        ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
//        ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
//        knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
      }
      else
      {
        // only split, no need to transform
        LINALG::SplitMatrix2x2(k_fseffmatrix,fgsmdofrowmap_,fgndofrowmap_,gsmdofrowmap_,gndofrowmap_,k_fs_smsm,k_fs_smn,k_fs_nsm,k_fs_nn);
        LINALG::SplitMatrix2x2(Feffmatrix,fgsmdofrowmap_,fgndofrowmap_,fgsmdofrowmap_,fgndofrowmap_,F_smsm,F_smn,F_nsm,F_nn);
      }

      // further splits into slave part + master part
      LINALG::SplitMatrix2x2(k_fs_smsm,fgsdofrowmap_,fgmdofrowmap_,gsdofrowmap_,gmdofrowmap_,k_fs_ss,k_fs_sm,k_fs_ms,k_fs_mm);
      LINALG::SplitMatrix2x2(k_fs_smn,fgsdofrowmap_,fgmdofrowmap_,gndofrowmap_,tempmap,k_fs_sn,tempmtx1,k_fs_mn,tempmtx2);
      LINALG::SplitMatrix2x2(k_fs_nsm,fgndofrowmap_,ftempmap,gsdofrowmap_,gmdofrowmap_,k_fs_ns,k_fs_nm,tempmtx1,tempmtx2);

      LINALG::SplitMatrix2x2(F_smsm,fgsdofrowmap_,fgmdofrowmap_,fgsdofrowmap_,fgmdofrowmap_,F_ss,F_sm,F_ms,F_mm);
      LINALG::SplitMatrix2x2(F_smn,fgsdofrowmap_,fgmdofrowmap_,fgndofrowmap_,ftempmap,F_sn,tempmtx1,F_mn,tempmtx2);
      LINALG::SplitMatrix2x2(F_nsm,fgndofrowmap_,ftempmap,fgsdofrowmap_,fgmdofrowmap_,F_ns,F_nm,tempmtx1,tempmtx2);

      /**********************************************************************/
      /* (4) Split feff into 3 subvectors                                   */
      /**********************************************************************/

      // we want to split f into 3 groups s.m,n
      Teuchos::RCP<Epetra_Vector> fs, fm, fn;

      // temporarily we need the group sm
      Teuchos::RCP<Epetra_Vector> fsm;

      // do the vector splitting smn -> sm+n
      if (ParRedist())
      {
        dserror("CHECK ME!");
        // split and transform to redistributed maps
        LINALG::SplitVector(*ProblemDofs(),*feff,pgsmdofrowmap_,fsm,gndofrowmap_,fn);
        Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
        LINALG::Export(*fsm,*fsmtemp);
        fsm = fsmtemp;
      }
      else
      {
        // only split, no need to transform
        LINALG::SplitVector(*falldofrowmap_,*feff,fgsmdofrowmap_,fsm,fgndofrowmap_,fn);
      }

      // abbreviations for slave  and master set
      int sset = fgsdofrowmap_->NumGlobalElements();
      int mset = fgmdofrowmap_->NumGlobalElements();

      // we want to split fsm into 2 groups s,m
      fs = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));
      fm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));

      // do the vector splitting sm -> s+m
      LINALG::SplitVector(*fgsmdofrowmap_,*fsm,fgsdofrowmap_,fs,fgmdofrowmap_,fm);

      // store some stuff for static condensation of poro no pen. LM

      ffs_   = fs;

      cfssn_ = k_fs_sn;
      cfssm_ = k_fs_sm;
      cfsss_ = k_fs_ss;

      fsn_ = F_sn;
      fsm_ = F_sm;
      fss_ = F_ss;

      /**********************************************************************/
      /* (5) Split slave quantities into active / inactive                  */
      /**********************************************************************/

      // we want to split kssmod into 2 groups a,i = 4 blocks
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_aa, k_fs_ai, k_fs_ia, k_fs_ii;
      Teuchos::RCP<LINALG::SparseMatrix> F_aa, F_ai, F_ia, F_ii;

      // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_an, k_fs_in, k_fs_am, k_fs_im, k_fs_ma, k_fs_mi;
      Teuchos::RCP<LINALG::SparseMatrix> F_an, F_in, F_am, F_im, F_ma, F_mi;

      // we will get the i rowmap as a by-product
      Teuchos::RCP<Epetra_Map> gidofs;
      Teuchos::RCP<Epetra_Map> fgidofs;

      // do the splitting
      LINALG::SplitMatrix2x2(k_fs_ss,fgactivedofs_,fgidofs,gactivedofs_,gidofs,k_fs_aa,k_fs_ai,k_fs_ia,k_fs_ii);
      LINALG::SplitMatrix2x2(k_fs_sn,fgactivedofs_,fgidofs,gndofrowmap_,tempmap,k_fs_an,tempmtx1,k_fs_in,tempmtx2);
      LINALG::SplitMatrix2x2(k_fs_sm,fgactivedofs_,fgidofs,gmdofrowmap_,tempmap,k_fs_am,tempmtx1,k_fs_im,tempmtx2);
      LINALG::SplitMatrix2x2(k_fs_ms,fgmdofrowmap_,tempmap,gactivedofs_,gidofs,k_fs_ma,k_fs_mi,tempmtx1,tempmtx2);

      LINALG::SplitMatrix2x2(F_ss,fgactivedofs_,fgidofs,fgactivedofs_,fgidofs,F_aa,F_ai,F_ia,F_ii);
      LINALG::SplitMatrix2x2(F_sn,fgactivedofs_,fgidofs,fgndofrowmap_,tempmap,F_an,tempmtx1,F_in,tempmtx2);
      LINALG::SplitMatrix2x2(F_sm,fgactivedofs_,fgidofs,fgmdofrowmap_,tempmap,F_am,tempmtx1,F_im,tempmtx2);
      LINALG::SplitMatrix2x2(F_ms,fgmdofrowmap_,tempmap,fgactivedofs_,fgidofs,F_ma,F_mi,tempmtx1,tempmtx2);

      // abbreviations for active and inactive set
      int aset = fgactivedofs_->NumGlobalElements();
      int iset = fgidofs->NumGlobalElements();

      // we want to split fsmod into 2 groups a,i
      Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
      Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*fgidofs));

      // do the vector splitting s -> a+i
      LINALG::SplitVector(*fgsdofrowmap_,*fs,fgactivedofs_,fa,fgidofs,fi);

      /**********************************************************************/
      /* (7) Build the final K blocks                                       */
      /* where K stands for k_fs and F!!!                                   */
      /**********************************************************************/

      //----------------------------------------------------------- FIRST LINE
      // knn: nothing to do

      // knm: nothing to do

      // kns: nothing to do
      //---------------------------------------------------------- SECOND LINE --- Will just exist when starting with two sided poro contact!!!
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_mnmod;
      Teuchos::RCP<LINALG::SparseMatrix> F_mnmod;

      Teuchos::RCP<LINALG::SparseMatrix> k_fs_mmmod;
      Teuchos::RCP<LINALG::SparseMatrix> F_mmmod;

      Teuchos::RCP<LINALG::SparseMatrix> k_fs_mimod;
      Teuchos::RCP<LINALG::SparseMatrix> F_mimod;

      Teuchos::RCP<LINALG::SparseMatrix> k_fs_mamod;
      Teuchos::RCP<LINALG::SparseMatrix> F_mamod;
      if (mset)
      {
        // kmn: add T(mhataam)*kan
        k_fs_mnmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
        k_fs_mnmod->Add(*k_fs_mn,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> k_fs_mnadd = LINALG::MLMultiply(*fmhataam_,true,*k_fs_an,false,false,false,true);
        k_fs_mnmod->Add(*k_fs_mnadd,false,1.0,1.0);
        k_fs_mnmod->Complete(k_fs_mn->DomainMap(),k_fs_mn->RowMap());

        F_mnmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
        F_mnmod->Add(*F_mn,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> F_mnadd = LINALG::MLMultiply(*fmhataam_,true,*F_an,false,false,false,true);
        F_mnmod->Add(*k_fs_mnadd,false,1.0,1.0);
        F_mnmod->Complete(F_mn->DomainMap(),F_mn->RowMap());

      // kmm: add T(mhataam)*kam
        k_fs_mmmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
        k_fs_mmmod->Add(*k_fs_mm,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> k_fs_mmadd = LINALG::MLMultiply(*fmhataam_,true,*k_fs_am,false,false,false,true);
        k_fs_mmmod->Add(*k_fs_mmadd,false,1.0,1.0);
        k_fs_mmmod->Complete(k_fs_mm->DomainMap(),k_fs_mm->RowMap());

        F_mmmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
        F_mmmod->Add(*F_mm,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> F_mmadd = LINALG::MLMultiply(*fmhataam_,true,*F_am,false,false,false,true);
        F_mmmod->Add(*F_mmadd,false,1.0,1.0);
        F_mmmod->Complete(F_mm->DomainMap(),F_mm->RowMap());

        // kmi: add T(mhataam)*kai
        if (iset)
        {
          k_fs_mimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
          k_fs_mmmod->Add(*k_fs_mm,false,1.0,1.0);
          Teuchos::RCP<LINALG::SparseMatrix> k_fs_mmadd = LINALG::MLMultiply(*fmhataam_,true,*k_fs_am,false,false,false,true);
          k_fs_mmmod->Add(*k_fs_mmadd,false,1.0,1.0);
          k_fs_mmmod->Complete(k_fs_mm->DomainMap(),k_fs_mm->RowMap());

          F_mimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
          F_mmmod->Add(*F_mm,false,1.0,1.0);
          Teuchos::RCP<LINALG::SparseMatrix> F_mmadd = LINALG::MLMultiply(*fmhataam_,true,*F_am,false,false,false,true);
          F_mmmod->Add(*k_fs_mmadd,false,1.0,1.0);
          F_mmmod->Complete(F_mm->DomainMap(),F_mm->RowMap());
        }

        // kma: add T(mhataam)*kaa
        if (aset)
        {
          k_fs_mamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
          k_fs_mamod->Add(*k_fs_ma,false,1.0,1.0);
          Teuchos::RCP<LINALG::SparseMatrix> k_fs_maadd = LINALG::MLMultiply(*fmhataam_,true,*k_fs_aa,false,false,false,true);
          k_fs_mamod->Add(*k_fs_maadd,false,1.0,1.0);
          k_fs_mamod->Complete(k_fs_ma->DomainMap(),k_fs_ma->RowMap());

          F_mamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_,100));
          F_mamod->Add(*F_ma,false,1.0,1.0);
          Teuchos::RCP<LINALG::SparseMatrix> F_maadd = LINALG::MLMultiply(*fmhataam_,true,*F_aa,false,false,false,true);
          F_mamod->Add(*F_maadd,false,1.0,1.0);
          F_mamod->Complete(F_ma->DomainMap(),F_ma->RowMap());
        }
      }

      //----------------------------------------------------------- THIRD LINE
      //------------------- FOR 3D QUADRATIC CASE ----------------------------
      // kin: subtract T(dhat)*kan --
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_inmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
      k_fs_inmod->Add(*k_fs_in,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_inadd = LINALG::MLMultiply(*fdhat_,true,*k_fs_an,false,false,false,true);
      k_fs_inmod->Add(*k_fs_inadd,false,-1.0,1.0);
      k_fs_inmod->Complete(k_fs_in->DomainMap(),k_fs_in->RowMap());

      Teuchos::RCP<LINALG::SparseMatrix> F_inmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
      F_inmod->Add(*F_in,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> F_inadd = LINALG::MLMultiply(*fdhat_,true,*F_an,false,false,false,true);
      F_inmod->Add(*F_inadd,false,-1.0,1.0);
      F_inmod->Complete(F_in->DomainMap(),F_in->RowMap());

      // kim: subtract T(dhat)*kam
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_immod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
      k_fs_immod->Add(*k_fs_im,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_imadd = LINALG::MLMultiply(*fdhat_,true,*k_fs_am,false,false,false,true);
      k_fs_immod->Add(*k_fs_imadd,false,-1.0,1.0);
      k_fs_immod->Complete(k_fs_im->DomainMap(),k_fs_im->RowMap());

      Teuchos::RCP<LINALG::SparseMatrix> F_immod;
      if (mset)
      {
        F_immod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
        F_immod->Add(*F_im,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> F_imadd = LINALG::MLMultiply(*fdhat_,true,*F_am,false,false,false,true);
        F_immod->Add(*F_imadd,false,-1.0,1.0);
        F_immod->Complete(F_im->DomainMap(),F_im->RowMap());
      }

      // kii: subtract T(dhat)*kai
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_iimod;
      Teuchos::RCP<LINALG::SparseMatrix> F_iimod;
      if (iset)
      {
        k_fs_iimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
        k_fs_iimod->Add(*k_fs_ii,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> k_fs_iiadd = LINALG::MLMultiply(*fdhat_,true,*k_fs_ai,false,false,false,true);
        k_fs_iimod->Add(*k_fs_iiadd,false,-1.0,1.0);
        k_fs_iimod->Complete(k_fs_ii->DomainMap(),k_fs_ii->RowMap());

        F_iimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
        F_iimod->Add(*F_ii,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> F_iiadd = LINALG::MLMultiply(*fdhat_,true,*F_ai,false,false,false,true);
        F_iimod->Add(*F_iiadd,false,-1.0,1.0);
        F_iimod->Complete(F_ii->DomainMap(),F_ii->RowMap());
      }

      // kia: subtract T(dhat)*kaa
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_iamod;
      Teuchos::RCP<LINALG::SparseMatrix> F_iamod;
      if (iset && aset)
      {
        k_fs_iamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
        k_fs_iamod->Add(*k_fs_ia,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> k_fs_iaadd = LINALG::MLMultiply(*fdhat_,true,*k_fs_aa,false,false,false,true);
        k_fs_iamod->Add(*k_fs_iaadd,false,-1.0,1.0);
        k_fs_iamod->Complete(k_fs_ia->DomainMap(),k_fs_ia->RowMap());

        F_iamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs,100));
        F_iamod->Add(*F_ia,false,1.0,1.0);
        Teuchos::RCP<LINALG::SparseMatrix> F_iaadd = LINALG::MLMultiply(*fdhat_,true,*F_aa,false,false,false,true);
        F_iamod->Add(*F_iaadd,false,-1.0,1.0);
        F_iamod->Complete(F_ia->DomainMap(),F_ia->RowMap());
      }

      //---------------------------------------------------------- FOURTH LINE
      // nothing to do
      //----------------------------------------------------------- FIFTH LINE
      // kan: multiply tmatrix with invda and kan
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_anmod;
      Teuchos::RCP<LINALG::SparseMatrix> F_anmod;
      if (aset)
      {
        k_fs_anmod = LINALG::MLMultiply(*ftanginvD_,false,*k_fs_an,false,false,false,true);
        F_anmod = LINALG::MLMultiply(*ftanginvD_,false,*F_an,false,false,false,true);
      }

      // kam: multiply tmatrix with invda and kam
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_ammod;
      Teuchos::RCP<LINALG::SparseMatrix> F_ammod;
      if (aset)
      {
        k_fs_ammod = LINALG::MLMultiply(*ftanginvD_,false,*k_fs_am,false,false,false,true);
        if (mset)
          F_ammod = LINALG::MLMultiply(*ftanginvD_,false,*F_am,false,false,false,true);
      }
      // kai: multiply tmatrix with invda and kai
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_aimod;
      Teuchos::RCP<LINALG::SparseMatrix> F_aimod;
      if (aset && iset)
      {
        k_fs_aimod = LINALG::MLMultiply(*ftanginvD_,false,*k_fs_ai,false,false,false,true);
        F_aimod = LINALG::MLMultiply(*ftanginvD_,false,*F_ai,false,false,false,true);
      }
      // kaa: multiply tmatrix with invda and kaa
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_aamod;
      Teuchos::RCP<LINALG::SparseMatrix> F_aamod;
      if (aset)
      {
        k_fs_aamod = LINALG::MLMultiply(*ftanginvD_,false,*k_fs_aa,false,false,false,true);
        F_aamod = LINALG::MLMultiply(*ftanginvD_,false,*F_aa,false,false,false,true);
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
      if (IsSelfContact())
      {
        dserror("CHECK ME!");
//        Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
//        Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
//        Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(mold_->RowMap()));
//        if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
//        mold_->Multiply(true,*zoldexp,*tempvecm2);
//        if (mset) LINALG::Export(*tempvecm2,*tempvecm);
//        fm->Update(alphaf_,*tempvecm,1.0);
      }
      // if there is no self contact everything is ok
      else
      {
        Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
        fmoldtransp_->Multiply(false,*lambdaold_,*tempvecm);
        fm->Update(alphaf_,*tempvecm,1.0);
      }

      // fs: prepare alphaf * old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

      // for self contact, slave and master sets may have changed,
      // thus we have to export the product Dold^T * zold to fit
      if (IsSelfContact())
      {
        dserror("CHECK ME!");
//        Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
//        Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
//        if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
//        dold_->Multiply(true,*zoldexp,*tempvec);
//        if (sset) LINALG::Export(*tempvec,*fsadd);
      }
      // if there is no self contact everything is ok
      else
      {
        fdoldtransp_->Multiply(false,*lambdaold_,*fsadd);
      }

      // fa: subtract alphaf * old contact forces (t_n)
      if (aset)
      {
        Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
        LINALG::Export(*fsadd,*faadd);

        fa->Update(-alphaf_,*faadd,1.0);
      }

      // fm: add T(mhat)*fa
      Teuchos::RCP<Epetra_Vector> fmmod;
      if (mset)
      {
        fmmod = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
        if (aset) fmhataam_->Multiply(true,*fa,*fmmod);
        fmmod->Update(1.0,*fm,1.0);
      }

      //----------------------------------------------------------- THIRD LINE
      // fi: subtract alphaf * old contact forces (t_n)
      if (iset)
      {
        Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*fgidofs));
        LINALG::Export(*fsadd,*fiadd);
        fi->Update(-alphaf_,*fiadd,1.0);
      }

      // fi: add T(dhat)*fa
      Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*fgidofs));
      if (aset) fdhat_->Multiply(true,*fa,*fimod);
      fimod->Update(1.0,*fi,-1.0);

      //---------------------------------------------------------- FOURTH LINE
      // gactive: nothing to do
      //----------------------------------------------------------- FIFTH LINE
      // fa: multiply tmatrix with invda and fa
      Teuchos::RCP<Epetra_Vector> famod;
      if (aset)
      {
        famod = Teuchos::rcp(new Epetra_Vector(*fgactivet_));
        ftanginvD_->Multiply(false,*fa,*famod);
      }
      /********************************************************************/
      /* (9) Transform the final K blocks                                 */
      /********************************************************************/
      // The row maps of all individual matrix blocks are transformed to
      // the parallel layout of the underlying problem discretization.
      // Of course, this is only necessary in the parallel redistribution
      // case, where the contact interfaces have been redistributed
      // independently of the underlying problem discretization.

      if (ParRedist())
      {
        dserror("CHECK ME!");
//        //----------------------------------------------------------- FIRST LINE
//        // nothing to do (ndof-map independent of redistribution)
//
//        //---------------------------------------------------------- SECOND LINE
//        kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
//        kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
//        if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
//        if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);
//
//        //----------------------------------------------------------- THIRD LINE
//        if (iset)
//        {
//          kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
//          kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
//          kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
//          if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
//        }
//
//        //---------------------------------------------------------- FOURTH LINE
//        if (aset) smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);
//
//        //----------------------------------------------------------- FIFTH LINE
//        if (aset)
//        {
//          kanmod = MORTAR::MatrixRowTransform(kanmod,pgsdofrowmap_);
//          kammod = MORTAR::MatrixRowTransform(kammod,pgsdofrowmap_);
//          kaamod = MORTAR::MatrixRowTransform(kaamod,pgsdofrowmap_);
//          if (iset) kaimod = MORTAR::MatrixRowTransform(kaimod,pgsdofrowmap_);
//          pmatrix_ = MORTAR::MatrixRowTransform(pmatrix_,pgsdofrowmap_);
//        }
      }
      /**********************************************************************/
      /* (10) Global setup of kteffnew (including contact)                  */
      /**********************************************************************/

      Teuchos::RCP<LINALG::SparseMatrix> k_fs_effnew = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_,81,true,false,k_fseff->GetMatrixtype()));
      Teuchos::RCP<LINALG::SparseMatrix> F_effnew = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_,108,true,false,Feff->GetMatrixtype()));
      Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*falldofrowmap_);

      //----------------------------------------------------------- FIRST LINE
      // add n submatrices to kteffnew
      k_fs_effnew->Add(*k_fs_nn,false,1.0,1.0);
      k_fs_effnew->Add(*k_fs_nm,false,1.0,1.0);
      if (sset) k_fs_effnew->Add(*k_fs_ns,false,1.0,1.0);

      F_effnew->Add(*F_nn,false,1.0,1.0);
      F_effnew->Add(*F_nm,false,1.0,1.0);
      if (sset)F_effnew->Add(*F_ns,false,1.0,1.0);

      //---------------------------------------------------------- SECOND LINE
      // add m submatrices to kteffnew
      if (mset)
      {
        k_fs_effnew->Add(*k_fs_mnmod,false,1.0,1.0);
        k_fs_effnew->Add(*k_fs_mmmod,false,1.0,1.0);
        if (iset) k_fs_effnew->Add(*k_fs_mimod,false,1.0,1.0);
        if (aset) k_fs_effnew->Add(*k_fs_mamod,false,1.0,1.0);

        F_effnew->Add(*k_fs_mnmod,false,1.0,1.0);
        F_effnew->Add(*k_fs_mmmod,false,1.0,1.0);
        if (iset) F_effnew->Add(*F_mimod,false,1.0,1.0);
        if (aset) F_effnew->Add(*F_mamod,false,1.0,1.0);
      }
      //----------------------------------------------------------- THIRD LINE
      // add i submatrices to kteffnew
      if (iset) k_fs_effnew->Add(*k_fs_inmod,false,1.0,1.0);
      if (iset) k_fs_effnew->Add(*k_fs_immod,false,1.0,1.0);
      if (iset) k_fs_effnew->Add(*k_fs_iimod,false,1.0,1.0);
      if (iset && aset) k_fs_effnew->Add(*k_fs_iamod,false,1.0,1.0);

      if (iset) F_effnew->Add(*F_inmod,false,1.0,1.0);
      if (iset && mset) F_effnew->Add(*F_immod,false,1.0,1.0);
      if (iset) F_effnew->Add(*F_iimod,false,1.0,1.0);
      if (iset && aset) F_effnew->Add(*F_iamod,false,1.0,1.0);
      //---------------------------------------------------------- FOURTH LINE
      // add a submatrices to kteffnew
      if (aset) k_fs_effnew->Add(*fNCoup_lindisp_,false,1.0,1.0);
      if (aset) F_effnew->Add(*fNCoup_linvel_,false,1.0,1.0);
      //----------------------------------------------------------- FIFTH LINE
      // add a submatrices to kteffnew
      if (aset) k_fs_effnew->Add(*k_fs_anmod,false,-1.0,1.0);
      if (aset) k_fs_effnew->Add(*k_fs_ammod,false,-1.0,1.0);
      if (aset && iset) k_fs_effnew->Add(*k_fs_aimod,false,-1.0,1.0);
      if (aset) k_fs_effnew->Add(*k_fs_aamod,false,-1.0,1.0);
      if (aset) k_fs_effnew->Add(*flinTangentiallambda_,false,1.0,1.0);

      if (aset) F_effnew->Add(*F_anmod,false,-1.0,1.0);
      if (aset && mset) F_effnew->Add(*F_ammod,false,-1.0,1.0);
      if (aset && iset) F_effnew->Add(*F_aimod,false,-1.0,1.0);
      if (aset) F_effnew->Add(*F_aamod,false,-1.0,1.0);

      // FillComplete kteffnew (square)
      k_fs_effnew->Complete(*ProblemDofs(), *falldofrowmap_);
      F_effnew->Complete(*falldofrowmap_,*falldofrowmap_);
      /**********************************************************************/
      /* (11) Global setup of feffnew (including contact)                   */
      /**********************************************************************/

      //----------------------------------------------------------- FIRST LINE
      // add n subvector to feffnew
      Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));

       LINALG::Export(*fn,*fnexp);

      feffnew->Update(1.0,*fnexp,1.0);
      //---------------------------------------------------------- SECOND LINE
      // add m subvector to feffnew
      if (mset)
      {
        Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
        LINALG::Export(*fmmod,*fmmodexp);
        feffnew->Update(1.0,*fmmodexp,1.0);
      }
      //----------------------------------------------------------- THIRD LINE
      // add i subvector to feffnew
      Teuchos::RCP<Epetra_Vector> fimodexp;
      if (iset)
      {
        fimodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
        LINALG::Export(*fimod,*fimodexp);
        feffnew->Update(1.0,*fimodexp,1.0);
      }

      //---------------------------------------------------------- FOURTH LINE
      // add weighted nCoup vector to feffnew, if existing
      Teuchos::RCP<Epetra_Vector> nCoupexp;
      if (aset)
      {
        nCoupexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
        LINALG::Export(*fNCoup_,*nCoupexp);
        feffnew->Update(-1.0,*nCoupexp,1.0);
      }

      //----------------------------------------------------------- FIFTH LINE
      // add a subvector to feffnew
      Teuchos::RCP<Epetra_Vector> famodexp;
      if (aset)
      {
        famodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
        LINALG::Export(*famod,*famodexp);
        feffnew->Update(-1.0,*famodexp,1.0);
      }

      // finally do the replacement
      k_fseff = k_fs_effnew;
      Feff = F_effnew;
      feff = feffnew;
}


/*----------------------------------------------------------------------*
 | Poro Recovery method for contact LM                        ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::Recover(Teuchos::RCP<Epetra_Vector> disi,Teuchos::RCP<Epetra_Vector> veli)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  CoLagrangeStrategy::Recover(disi); //Update Structural Part! --> Here just Part from Coupling Matrix will be added!

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");
    // extract slave velocity from veli
    Teuchos::RCP<Epetra_Vector> velis = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));
    if (fgsdofrowmap_->NumGlobalElements()) LINALG::Export(*veli, *velis);
    // extract master velocity from veli
    Teuchos::RCP<Epetra_Vector> velim = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    if (fgmdofrowmap_->NumGlobalElements()) LINALG::Export(*veli, *velim);

    // extract other velocity from veli
    Teuchos::RCP<Epetra_Vector> velin = Teuchos::rcp(new Epetra_Vector(*fgndofrowmap_));
    if (fgndofrowmap_->NumGlobalElements()) LINALG::Export(*veli,*velin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(invd_,gactivedofs_,tempmap,gactivedofs_,tempmap,invda,tempmtx1,tempmtx2,tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> invdmod = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    invdmod->Add(*invda,false,1.0,1.0);
    invdmod->Complete();

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/
    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold * zold and Mold^T * zold to fit
    if (IsSelfContact()) //is not considered yet!
    {
      dserror("Trying to make poro selfcontact ... Check if this makes any sense!!!");
      // approximate update
      //z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      //invdmod->Multiply(false,*fs_,*z_);

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
//      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
//      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*mod2);
//      Teuchos::RCP<Epetra_Vector> mod3 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
//      dold_->Multiply(true,*mod2,*mod3);
//      Teuchos::RCP<Epetra_Vector> mod4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
//      if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*mod3,*mod4);
//      z_->Update(-alphaf_,*mod4,1.0);
//      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
//      invdmod->Multiply(true,*zcopy,*z_);
//      z_->Scale(1/(1-alphaf_));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> zfluid = Teuchos::rcp(new Epetra_Vector(z_->Map(),true));

      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      csfss_->Multiply(false,*velis,*mod);
      zfluid->Update(-1.0,*mod,0.0);
      csfsm_->Multiply(false,*velim,*mod);
      zfluid->Update(-1.0,*mod,1.0);
      csfsn_->Multiply(false,*velin,*mod);
      zfluid->Update(-1.0,*mod,1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*zfluid));
      invdmod->Multiply(true,*zcopy,*zfluid);
      zfluid->Scale(1/(1-alphaf_));

      z_->Update(1.0,*zfluid,1.0); // Add FluidCoupling Contribution to LM!!!
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
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate); //For Poro Contact done twice: already in structural contact --> not wanted

  return;
}

/*----------------------------------------------------------------------*
 | Poro Recovery method for no penetration LM                 ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::RecoverPoroNoPen(Teuchos::RCP<Epetra_Vector> disi,Teuchos::RCP<Epetra_Vector> veli)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (no_penetration_ && !IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);

    // extract slave velocity from veli
    Teuchos::RCP<Epetra_Vector> velis = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));
    if (fgsdofrowmap_->NumGlobalElements()) LINALG::Export(*veli, *velis);

    // extract master velocity from veli
    Teuchos::RCP<Epetra_Vector> velim = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    if (fgmdofrowmap_->NumGlobalElements()) LINALG::Export(*veli, *velim);

    // extract other velocity from veli
    Teuchos::RCP<Epetra_Vector> velin = Teuchos::rcp(new Epetra_Vector(*fgndofrowmap_));
    if (fgndofrowmap_->NumGlobalElements()) LINALG::Export(*veli,*velin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> finvda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(finvda_,fgactivedofs_,tempmap,gactivedofs_,tempmap,finvda,tempmtx1,tempmtx2,tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> finvdmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgsdofrowmap_,10));
    finvdmod->Add(*finvda,false,1.0,1.0);
    finvdmod->Complete(*gsdofrowmap_,*fgsdofrowmap_);

    /**********************************************************************/
    /* Update Lagrange multipliers lambda_n+1                                  */
    /**********************************************************************/
    {
      Teuchos::RCP<Epetra_Vector> flambda = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_,true));

      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

      cfssn_->Multiply(false,*disin,*mod);
      flambda->Update(-1.0,*mod,0.0);
      cfssm_->Multiply(false,*disim,*mod);
      flambda->Update(-1.0,*mod,1.0);
      cfsss_->Multiply(false,*disis,*mod);
      flambda->Update(-1.0,*mod,1.0);

      fss_->Multiply(false,*velis,*mod);
      flambda->Update(-1.0,*mod,1.0);
      fsm_->Multiply(false,*velim,*mod);
      flambda->Update(-1.0,*mod,1.0);
      fsn_->Multiply(false,*velin,*mod);
      flambda->Update(-1.0,*mod,1.0);

      flambda->Update(1.0,*ffs_,1.0);



      fdoldtransp_->Multiply(false,*lambdaold_,*mod);
      flambda->Update(-alphaf_,*mod,1.0);

      Teuchos::RCP<Epetra_Vector> lambdacopy = Teuchos::rcp(new Epetra_Vector(*flambda));

      finvdmod->Multiply(true,*lambdacopy,*lambda_); // should be lambda_ at the ende!!!

     //lambda_->Scale((1-alphaf_)); //-- is already scaled by this factor!!!
    }
  }
   // store updated LM into nodes
  SetState("lm",lambda_);

  return;
}

/*------------------------------------------------------------------------------*
 |  Additional update and output poro contact at end of time step     ager 08/14|
 *-----------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::UpdatePoroContact()
{
  if(no_penetration_)
  {
    //std::cout << "print lambda: " << *lambda_ << std::endl;
    lambdaold_->Update(1.0,*lambda_,0.0);
  }
  return;
}


/*------------------------------------------------------------------------*
 | Save mortar coupling matrices for evaluation of poro contact! ager 08/14|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SaveCouplingMatrices(Teuchos::RCP<LINALG::SparseMatrix> dhat,
                Teuchos::RCP<LINALG::SparseMatrix> mhataam, Teuchos::RCP<LINALG::SparseMatrix> invda)
{
  dhat_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*dhat));
  mhataam_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*mhataam));
  invda_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*invda));
}

/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                          ager 08/14|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec)
{
  if (statename=="displacement" || statename=="olddisplacement")
  {
    CONTACT::CoAbstractStrategy::SetState(statename,vec);
    vec->Comm().Barrier();
    return;
  }

  if (statename=="fvelocity" || statename=="svelocity" || statename=="lm")
  {
    //set state on interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {

      interface_[i]->SetState(statename, vec);
      DRT::Discretization& idiscret_ = interface_[i]->Discret();

      if (statename=="fvelocity" or statename=="svelocity")
       {
         bool isfvel = (statename=="fvelocity");
         // alternative method to get vec to full overlap
         Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*idiscret_.DofColMap(),true));

         LINALG::Export(*vec,*global);

         // loop over all nodes to set current velocity
         // (use fully overlapping column map)
         for (int i=0;i<idiscret_.NumMyColNodes();++i)
         {
           CONTACT::CoNode* node = static_cast<CONTACT::CoNode*>(idiscret_.lColNode(i));
           const int numdof = node->NumDof();
           std::vector<double> myvel(numdof);
           std::vector<int> lm(numdof);

           for (int j=0;j<numdof;++j)
             lm[j]=node->Dofs()[j];
           DRT::UTILS::ExtractMyValues(*global,myvel,lm);

           // add myvel[2]=0 for 2D problems
           if (myvel.size()<3)
             myvel.resize(3);
           // set current configuration
           for (int j=0;j<3;++j)
           {
             if (isfvel)
               node->CoPoroData().fvel()[j] = myvel[j];
             else
               node->CoPoroData().svel()[j] = myvel[j];
           }
         }
       }
       if (statename=="lm")
       {
         // alternative method to get vec to full overlap
         Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*idiscret_.DofColMap(),true));
         LINALG::Export(*vec,*global);

         // loop over all nodes to set current velocity
         // (use fully overlapping column map)
         for (int i=0;i<idiscret_.NumMyColNodes();++i)
         {
           CONTACT::CoNode* node = static_cast<CONTACT::CoNode*>(idiscret_.lColNode(i));

           const int numdof = node->NumDof();
           std::vector<double> mylm(numdof);
           std::vector<int> lm(numdof);

           for (int j=0;j<numdof;++j)
             lm[j]=node->Dofs()[j];

           DRT::UTILS::ExtractMyValues(*global,mylm,lm);

           // add myvel[2]=0 for 2D problems
           if (mylm.size()<3)
             mylm.resize(3);
           // set current configuration
           for (int j=0;j<3;++j)
           {
             node->CoPoroData().poroLM()[j] = mylm[j];
           }
         }
       }
    }
  }
  return;
}
