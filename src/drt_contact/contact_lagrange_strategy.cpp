/*---------------------------------------------------------------------*/
/*!
\file contact_lagrange_strategy.cpp

\brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*---------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "contact_lagrange_strategy.H"
#include "contact_interface.H"
#include "contact_defines.H"
#include "contact_paramsinterface.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"

///*----------------------------------------------------------------------*
// | ctor (public)                                              popp 05/09|
// *----------------------------------------------------------------------*/
//CONTACT::CoLagrangeStrategy::CoLagrangeStrategy(
//    const Epetra_Map* DofRowMap,
//    const Epetra_Map* NodeRowMap,
//    Teuchos::ParameterList params,
//    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
//    int dim,
//    Teuchos::RCP<const Epetra_Comm> comm,
//    double alphaf,
//    int maxdof)
//    : CoAbstractStrategy(Teuchos::rcp(new CONTACT::AbstractStratDataContainer()),
//        DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof),
//      activesetssconv_(false),
//      activesetconv_(false),
//      activesetsteps_(1)
//{
//  // empty constructor body
//  return;
//}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoLagrangeStrategy::CoLagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
    int dim,
    Teuchos::RCP<const Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoAbstractStrategy(data_ptr,DofRowMap,NodeRowMap,params,dim,comm,alphaf,maxdof),
      interface_( interface ),
      evalForceCalled_(false),
      activesetssconv_(false),
      activesetconv_(false),
      activesetsteps_(1)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::Initialize()
{
  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    // (re)setup global tangent matrix
    if (!friction_)
      tmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,3));

    // (re)setup global normal matrix
    if (regularized_)
      nmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,3));

    // (re)setup global matrix containing gap derivatives
    smatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,3));

    // inactive rhs for the saddle point problem
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_, *gactivedofs_);
    inactiverhs_ = LINALG::CreateVector(*gidofs, true);

    // further terms depend on friction case
    // (re)setup global matrix containing "no-friction"-derivatives
    if (!friction_)
    {
      // tangential rhs
      tangrhs_ = LINALG::CreateVector(*gactivedofs_, true);
      tderivmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,3));
      if (regularized_)
        nderivmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,3));
    }
    // (re)setup of global friction
    else
    {
      // here the calculation of gstickt is necessary
      Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);
      linstickLM_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickdofs,3));
      linstickDIS_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickdofs,3));
      linstickRHS_ = LINALG::CreateVector(*gstickdofs,true);

      linslipLM_ = Teuchos::rcp(new LINALG::SparseMatrix(*gslipdofs_,3));
      linslipDIS_ = Teuchos::rcp(new LINALG::SparseMatrix(*gslipdofs_,3));
      linslipRHS_ = LINALG::CreateVector(*gslipdofs_,true);
    }
  }

  else
  {
    // (re)setup global tangent matrix
    if (!friction_)
      tmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,3));

    // (re)setup global normal matrix
    if (regularized_)
      nmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3));

    // (re)setup global matrix containing gap derivatives
    smatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3));

    // inactive rhs for the saddle point problem
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_, *gactivedofs_);
    inactiverhs_ = LINALG::CreateVector(*gidofs, true);

    // further terms depend on friction case
    // (re)setup global matrix containing "no-friction"-derivatives
    if (!friction_)
    {
      // tangential rhs
      tangrhs_ = LINALG::CreateVector(*gactivet_, true);
      tderivmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,3));
      if (regularized_)
        nderivmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3));
    }
    // (re)setup of global friction
    else
    {
      // here the calculation of gstickt is necessary
      Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
      linstickLM_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickt,3));
      linstickDIS_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickt,3));
      linstickRHS_ = LINALG::CreateVector(*gstickt,true);

      linslipLM_ = Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));
      linslipDIS_ = Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));
      linslipRHS_ = LINALG::CreateVector(*gslipt_,true);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate frictional contact (public)                   gitterle 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                   Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
      LINALG::Export(*g_,*gact);
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
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
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
    interface_[i]->AssembleLinStick(*linstickLM_,*linstickDIS_,*linstickRHS_);
    interface_[i]->AssembleLinSlip(*linslipLM_,*linslipDIS_,*linslipRHS_);
    if (SystemType() != INPAR::CONTACT::system_condensed)
      interface_[i]->AssembleInactiverhs(*inactiverhs_);
  }
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    smatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
  }
  else
  {
    // FillComplete() global matrix S
    smatrix_->Complete(*gsmdofrowmap_,*gactiven_);
  }

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);

  // FillComplete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);

  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    linstickLM_->Complete(*gstickdofs,*gstickdofs);
    linstickDIS_->Complete(*gsmdofrowmap_,*gstickdofs);
    linslipLM_->Complete(*gslipdofs_,*gslipdofs_);
    linslipDIS_->Complete(*gsmdofrowmap_,*gslipdofs_);
  }
  else
  {
    linstickLM_->Complete(*gstickdofs,*gstickt);
    linstickDIS_->Complete(*gsmdofrowmap_,*gstickt);

    // FillComplete global Matrix linslipLM_ and linslipDIS_
    linslipLM_->Complete(*gslipdofs_,*gslipt_);
    linslipDIS_->Complete(*gsmdofrowmap_,*gslipt_);
  }
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify lindmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*invtrafo_,true,*lindmatrix_,false,false,false,true);
    lindmatrix_   = temp1;
  }

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype_ == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    /********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                    */
    /********************************************************************/
    Teuchos::RCP<LINALG::SparseMatrix> invd = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix_));
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
    int err = 0;

    // extract diagonal of invd into diag
    invd->ExtractDiagonalCopy(*diag);

    // set zero diagonal values to dummy 1.0
    for (int i=0;i<diag->MyLength();++i)
      if ((*diag)[i]==0.0) (*diag)[i]=1.0;

    // scalar inversion of diagonal values
    err = diag->Reciprocal(*diag);
    if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      Teuchos::RCP<Epetra_Vector> lmDBC=LINALG::CreateVector(*gsdofrowmap_,true);
      LINALG::Export(*pgsdirichtoggle_,*lmDBC);
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*gsdofrowmap_,true);
      tmp->Multiply(1.,*diag,*lmDBC,0.);
      diag->Update(-1.,*tmp,1.);
    }

    // re-insert inverted diagonal into invd
    err = invd->ReplaceDiagonalValues(*diag);
    // we cannot use this check, as we deliberately replaced zero entries
    //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = LINALG::MLMultiply(*invd,false,*mmatrix_,false,false,false,true);

    /********************************************************************/
    /* (2) Add contact stiffness terms to kteff                         */
    /********************************************************************/

    // transform if necessary
    if (ParRedist())
    {
      lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
      linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
    }

    kteff->UnComplete();
    kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
    kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
    kteff->Complete();

    /********************************************************************/
    /* (3) Split kteff into 3x3 matrix blocks                           */
    /********************************************************************/

    // we want to split k into 3 groups s,m,n = 9 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

    // temporarily we need the blocks ksmsm, ksmn, knsm
    // (FIXME: because a direct SplitMatrix3x3 is still missing!)
    Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

    // split into slave/master part + structure part
    Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);
    if (ParRedist())
    {
      // split and transform to redistributed maps
      LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
      ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
      ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
      knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
    }
    else
    {
      // only split, no need to transform
      LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    }

    // further splits into slave part + master part
    LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
    LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
    LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

    /********************************************************************/
    /* (4) Split feff into 3 subvectors                                 */
    /********************************************************************/

    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fs, fm, fn;

    // temporarily we need the group sm
    Teuchos::RCP<Epetra_Vector> fsm;

    // do the vector splitting smn -> sm+n
    if (ParRedist())
    {
      // split and transform to redistributed maps
      LINALG::SplitVector(*ProblemDofs(),*feff,pgsmdofrowmap_,fsm,gndofrowmap_,fn);
      Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
      LINALG::Export(*fsm,*fsmtemp);
      fsm = fsmtemp;
    }
    else
    {
      // only split, no need to transform
      LINALG::SplitVector(*ProblemDofs(),*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);
    }

    // abbreviations for slave and master set
    const int sset = gsdofrowmap_->NumGlobalElements();
    const int mset = gmdofrowmap_->NumGlobalElements();

    // we want to split fsm into 2 groups s,m
    fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

    // do the vector splitting sm -> s+m
    LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

    // store some stuff for static condensation of LM
    fs_   = fs;
    invd_ = invd;
    ksn_  = ksn;
    ksm_  = ksm;
    kss_  = kss;

    //--------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //--------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    // D^(-1)    ---->   T * D^(-1)
    // \hat{M}   ---->   T * \hat{M}
    //--------------------------------------------------------------------
    if (Dualquadslave3d())
    {
      // modify dmatrix_, invd_ and mhatmatrix_
      Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
      Teuchos::RCP<LINALG::SparseMatrix> temp3 = LINALG::MLMultiply(*trafo_,false,*invd_,false,false,false,true);
      Teuchos::RCP<LINALG::SparseMatrix> temp4 = LINALG::MLMultiply(*trafo_,false,*mhatmatrix_,false,false,false,true);
      dmatrix_    = temp2;
      invd_       = temp3;
      mhatmatrix_ = temp4;
    }

    /********************************************************************/
    /* (5) Split slave quantities into active / inactive, stick / slip  */
    /********************************************************************/
    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;

    // do the splitting
    LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
    LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
    LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
    LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

    // we want to split kaa into 2 groups sl,st = 4 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

    // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
    Teuchos::RCP<LINALG::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> temp1map;
    Teuchos::RCP<LINALG::SparseMatrix> temp1mtx4,temp1mtx5;

    // we will get the stick rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gstdofs;

    LINALG::SplitMatrix2x2(kaa,gactivedofs_,gidofs,gstdofs,gslipdofs_,kast,kasl,temp1mtx4,temp1mtx5);

    // abbreviations for active and inactive set, stick and slip set
    const int aset = gactivedofs_->NumGlobalElements();
    const int iset = gidofs->NumGlobalElements();
    const int stickset = gstdofs->NumGlobalElements();
    const int slipset = gslipdofs_->NumGlobalElements();

    // we want to split fs into 2 groups a,i
    Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // do the vector splitting s -> a+i
    LINALG::SplitVector(*gsdofrowmap_,*fs,gactivedofs_,fa,gidofs,fi);

    // we want to split fa into 2 groups sl,st
    Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
    Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

    // do the vector splitting a -> sl+st
    if(aset)
      LINALG::SplitVector(*gactivedofs_,*fa,gslipdofs_,fsl,gstdofs,fst);

    /********************************************************************/
    /* (6) Isolate necessary parts from invd and mhatmatrix             */
    /********************************************************************/

    // active, stick and slip part of invd
    Teuchos::RCP<LINALG::SparseMatrix> invda, invdsl, invdst;
    LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
    LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,tempmtx3);
    LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gstdofs,gslipdofs_,invdst,tempmtx1,tempmtx2,tempmtx3);

    // coupling part of dmatrix (only nonzero for 3D quadratic case!)
    Teuchos::RCP<LINALG::SparseMatrix> dai;
    LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,tempmtx1,dai,tempmtx2,tempmtx3);

     // do the multiplication dhat = invda * dai
    Teuchos::RCP<LINALG::SparseMatrix> dhat = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
    if (aset && iset) dhat = LINALG::MLMultiply(*invda,false,*dai,false,false,false,true);
    dhat->Complete(*gidofs,*gactivedofs_);

    // active part of mmatrix
    Teuchos::RCP<LINALG::SparseMatrix> mmatrixa;
    LINALG::SplitMatrix2x2(mmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mmatrixa,tempmtx1,tempmtx2,tempmtx3);

    // do the multiplication mhataam = invda * mmatrixa
    // (this is only different from mhata for 3D quadratic case!)
    Teuchos::RCP<LINALG::SparseMatrix> mhataam = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
    if (aset) mhataam = LINALG::MLMultiply(*invda,false,*mmatrixa,false,false,false,true);
    mhataam->Complete(*gmdofrowmap_,*gactivedofs_);

    // for the case without full linearization, we still need the
    // "classical" active part of mhat, which is isolated here
    Teuchos::RCP<LINALG::SparseMatrix> mhata;
    LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);

    // scaling of invd and dai
    invda->Scale(1/(1-alphaf_));
    invdsl->Scale(1/(1-alphaf_));
    invdst->Scale(1/(1-alphaf_));
    dai->Scale(1-alphaf_);

    /********************************************************************/
    /* (7) Build the final K blocks                                     */
    /********************************************************************/

    //--------------------------------------------------------- FIRST LINE
    // knn: nothing to do

    // knm: nothing to do

    // kns: nothing to do

    //-------------------------------------------------------- SECOND LINE
    // kmn: add T(mhataam)*kan
    Teuchos::RCP<LINALG::SparseMatrix> kmnmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmnmod->Add(*kmn,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhataam,true,*kan,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

    // kmm: add T(mhataam)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmmmod->Add(*kmm,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhataam,true,*kam,false,false,false,true);
    kmmmod->Add(*kmmadd,false,1.0,1.0);
    kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

    // kmi: add T(mhataam)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmimod->Add(*kmi,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhataam,true,*kai,false,false,false,true);
      kmimod->Add(*kmiadd,false,1.0,1.0);
      kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
    }

    // kma: add T(mhataam)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmamod->Add(*kma,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhataam,true,*kaa,false,false,false,true);
      kmamod->Add(*kmaadd,false,1.0,1.0);
      kmamod->Complete(kma->DomainMap(),kma->RowMap());
    }

    //--------------------------------------------------------- THIRD LINE
    // kin: subtract T(dhat)*kan
    Teuchos::RCP<LINALG::SparseMatrix> kinmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kinmod->Add(*kin,false,1.0,1.0);
    if (aset && iset)
    {
      Teuchos::RCP<LINALG::SparseMatrix> kinadd = LINALG::MLMultiply(*dhat,true,*kan,false,false,false,true);
      kinmod->Add(*kinadd,false,-1.0,1.0);
    }
    kinmod->Complete(kin->DomainMap(),kin->RowMap());

    // kim: subtract T(dhat)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kimmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kimmod->Add(*kim,false,1.0,1.0);
    if (aset && iset)
    {
      Teuchos::RCP<LINALG::SparseMatrix> kimadd = LINALG::MLMultiply(*dhat,true,*kam,false,false,false,true);
      kimmod->Add(*kimadd,false,-1.0,1.0);
    }
    kimmod->Complete(kim->DomainMap(),kim->RowMap());

    // kii: subtract T(dhat)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kiimod;
    if (iset)
    {
      kiimod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiimod->Add(*kii,false,1.0,1.0);
      if (aset)
      {
        Teuchos::RCP<LINALG::SparseMatrix> kiiadd = LINALG::MLMultiply(*dhat,true,*kai,false,false,false,true);
        kiimod->Add(*kiiadd,false,-1.0,1.0);
      }
      kiimod->Complete(kii->DomainMap(),kii->RowMap());
    }

    // kia: subtract T(dhat)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kiamod;
    if (iset && aset)
    {
      kiamod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiamod->Add(*kia,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kiaadd = LINALG::MLMultiply(*dhat,true,*kaa,false,false,false,true);
      kiamod->Add(*kiaadd,false,-1.0,1.0);
      kiamod->Complete(kia->DomainMap(),kia->RowMap());
    }

    //-------------------------------------------------------- FOURTH LINE

    //--------------------------------------------------------- FIFTH LINE
    // blocks for complementary conditions (stick nodes)

    // kstn: multiply with linstickLM
    Teuchos::RCP<LINALG::SparseMatrix> kstnmod;
    if (stickset)
    {
      kstnmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      kstnmod = LINALG::MLMultiply(*kstnmod,false,*kan,false,false,false,true);
    }

    // kstm: multiply with linstickLM
    Teuchos::RCP<LINALG::SparseMatrix> kstmmod;
    if(stickset)
    {
      kstmmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      kstmmod = LINALG::MLMultiply(*kstmmod,false,*kam,false,false,false,true);
    }

    // ksti: multiply with linstickLM
    Teuchos::RCP<LINALG::SparseMatrix> kstimod;
    if(stickset && iset)
    {
      kstimod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      kstimod = LINALG::MLMultiply(*kstimod,false,*kai,false,false,false,true);
    }

    // kstsl: multiply with linstickLM
    Teuchos::RCP<LINALG::SparseMatrix> kstslmod;
    if(stickset && slipset)
    {
      kstslmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      kstslmod = LINALG::MLMultiply(*kstslmod,false,*kasl,false,false,false,true);
    }

    // kststmod: multiply with linstickLM
    Teuchos::RCP<LINALG::SparseMatrix> kststmod;
    if (stickset)
    {
      kststmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      kststmod = LINALG::MLMultiply(*kststmod,false,*kast,false,false,false,true);
    }

    //--------------------------------------------------------- SIXTH LINE
    // blocks for complementary conditions (slip nodes)

    // ksln: multiply with linslipLM
    Teuchos::RCP<LINALG::SparseMatrix> kslnmod;
    if(slipset)
    {
      kslnmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      kslnmod = LINALG::MLMultiply(*kslnmod,false,*kan,false,false,false,true);
    }

    // kslm: multiply with linslipLM
    Teuchos::RCP<LINALG::SparseMatrix> kslmmod;
    if(slipset)
    {
      kslmmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      kslmmod = LINALG::MLMultiply(*kslmmod,false,*kam,false,false,false,true);
    }

    // ksli: multiply with linslipLM
    Teuchos::RCP<LINALG::SparseMatrix> kslimod;
    if (slipset && iset)
    {
      kslimod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      kslimod = LINALG::MLMultiply(*kslimod,false,*kai,false,false,false,true);
    }

    // kslsl: multiply with linslipLM
    Teuchos::RCP<LINALG::SparseMatrix> kslslmod;
    if(slipset)
    {
      kslslmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      kslslmod = LINALG::MLMultiply(*kslslmod,false,*kasl,false,false,false,true);
    }

    // slstmod: multiply with linslipLM
    Teuchos::RCP<LINALG::SparseMatrix> kslstmod;
    if (slipset && stickset)
    {
      kslstmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      kslstmod = LINALG::MLMultiply(*kslstmod,false,*kast,false,false,false,true);
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
    if (IsSelfContact())
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(mold_->RowMap()));
      if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      mold_->Multiply(true,*zoldexp,*tempvecm2);
      if (mset) LINALG::Export(*tempvecm2,*tempvecm);
      fm->Update(alphaf_,*tempvecm,1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->Multiply(true,*zold_,*tempvecm);
      fm->Update(alphaf_,*tempvecm,1.0);
    }

    // fs: prepare alphaf * old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Dold^T * zold to fit
    if (IsSelfContact())
    {
      Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      dold_->Multiply(true,*zoldexp,*tempvec);
      if (sset) LINALG::Export(*tempvec,*fsadd);
    }
    // if there is no self contact everything is ok
    else
    {
      dold_->Multiply(true,*zold_,*fsadd);
    }

    // fa: subtract alphaf * old contact forces (t_n)
    if (aset)
    {
      Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      LINALG::Export(*fsadd,*faadd);
      fa->Update(-alphaf_,*faadd,1.0);
    }

    // fm: add T(mhat)*fa
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (aset) mhataam->Multiply(true,*fa,*fmmod);
    fmmod->Update(1.0,*fm,1.0);

    //--------------------------------------------------------- THIRD LINE
    // fi: subtract alphaf * old contact forces (t_n)
    if (iset)
    {
      Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
      LINALG::Export(*fsadd,*fiadd);
      fi->Update(-alphaf_,*fiadd,1.0);
    }

    // fi: add T(dhat)*fa
    Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
    if (aset && iset) dhat->Multiply(true,*fa,*fimod);
    fimod->Update(1.0,*fi,-1.0);

    //-------------------------------------------------------- FOURTH LINE

    //--------------------------------------------------------- FIFTH LINE
    Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);  // get global stick dofs

    // split the lagrange multiplier vector in stick and slip part
    Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
    Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
    Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

    LINALG::SplitVector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
    LINALG::SplitVector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
    Teuchos::RCP<Epetra_Vector> tempvec1;

    // fst: mutliply with linstickLM
    Teuchos::RCP<Epetra_Vector> fstmod;
    if (stickset)
    {
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        fstmod = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
      else
        fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
      Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
      temp1->Multiply(false,*fa,*fstmod);

      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
      else
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

      linstickLM_->Multiply(false, *zst, *tempvec1);
      fstmod->Update(-1.0,*tempvec1,1.0);
    }

    //--------------------------------------------------------- SIXTH LINE
    // fsl: mutliply with linslipLM
    Teuchos::RCP<Epetra_Vector> fslmod;
    Teuchos::RCP<Epetra_Vector> fslwmod;

    if (slipset)
    {
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        fslmod = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
      else
        fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
      Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
      temp->Multiply(false,*fa,*fslmod);

      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
      else
        tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

      linslipLM_->Multiply(false, *zsl, *tempvec1);

      fslmod->Update(-1.0,*tempvec1,1.0);
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
      if (aset)
      {
        smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);
      }

      //----------------------------------------------------------- FIFTH LINE
      if (stickset)
      {
        kstnmod = MORTAR::MatrixRowTransform(kstnmod,pgsdofrowmap_);
        kstmmod = MORTAR::MatrixRowTransform(kstmmod,pgsdofrowmap_);
        if (iset) kstimod = MORTAR::MatrixRowTransform(kstimod,pgsdofrowmap_);
        if (slipset) kstslmod = MORTAR::MatrixRowTransform(kstslmod,pgsdofrowmap_);
        kststmod = MORTAR::MatrixRowTransform(kststmod,pgsdofrowmap_);
        linstickDIS_ = MORTAR::MatrixRowTransform(linstickDIS_,pgsdofrowmap_);
      }

      //----------------------------------------------------------- SIXTH LINE
      if (slipset)
      {
        kslnmod = MORTAR::MatrixRowTransform(kslnmod,pgsdofrowmap_);
        kslmmod = MORTAR::MatrixRowTransform(kslmmod,pgsdofrowmap_);
        if (iset) kslimod = MORTAR::MatrixRowTransform(kslimod,pgsdofrowmap_);
        if (stickset) kslstmod = MORTAR::MatrixRowTransform(kslstmod,pgsdofrowmap_);
        kslslmod = MORTAR::MatrixRowTransform(kslslmod,pgsdofrowmap_);
        linslipDIS_ = MORTAR::MatrixRowTransform(linslipDIS_,pgsdofrowmap_);
      }
    }

    /********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                */
    /********************************************************************/

    Teuchos::RCP<LINALG::SparseMatrix> kteffnew = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(),81,true,false,kteffmatrix->GetMatrixtype()));
    Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*ProblemDofs());

    //--------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->Add(*knn,false,1.0,1.0);
    kteffnew->Add(*knm,false,1.0,1.0);
    if (sset) kteffnew->Add(*kns,false,1.0,1.0);

    //-------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->Add(*kmnmod,false,1.0,1.0);
    kteffnew->Add(*kmmmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kmimod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);

    //--------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->Add(*kinmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kimmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kiimod,false,1.0,1.0);
    if (iset && aset) kteffnew->Add(*kiamod,false,1.0,1.0);

    //-------------------------------------------------------- FOURTH LINE

    // add a submatrices to kteffnew
    if (aset) kteffnew->Add(*smatrix_,false,1.0,1.0);

    //--------------------------------------------------------- FIFTH LINE
    // add st submatrices to kteffnew
    if (stickset) kteffnew->Add(*kstnmod,false,1.0,1.0);
    if (stickset) kteffnew->Add(*kstmmod,false,1.0,1.0);
    if (stickset && iset) kteffnew->Add(*kstimod,false,1.0,1.0);
    if (stickset && slipset) kteffnew->Add(*kstslmod,false,1.0,1.0);
    if (stickset) kteffnew->Add(*kststmod,false,1.0,1.0);

    // add terms of linearization of sick condition to kteffnew
    if (stickset) kteffnew->Add(*linstickDIS_,false,-1.0,1.0);

    //--------------------------------------------------------- SIXTH LINE
    // add sl submatrices to kteffnew
    if (slipset) kteffnew->Add(*kslnmod,false,1.0,1.0);
    if (slipset) kteffnew->Add(*kslmmod,false,1.0,1.0);
    if (slipset && iset) kteffnew->Add(*kslimod,false,1.0,1.0);
    if (slipset) kteffnew->Add(*kslslmod,false,1.0,1.0);
    if (slipset && stickset) kteffnew->Add(*kslstmod,false,1.0,1.0);

    // add terms of linearization of slip condition to kteffnew and feffnew
    if (slipset) kteffnew->Add(*linslipDIS_,false,-1.0,+1.0);

    // FillComplete kteffnew (square)
    kteffnew->Complete();

    /********************************************************************/
    /* (11) Global setup of feffnew (including contact)                 */
    /********************************************************************/

    //--------------------------------------------------------- FIRST LINE
    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fn,*fnexp);
    feffnew->Update(1.0,*fnexp,1.0);

    //-------------------------------------------------------- SECOND LINE
    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fmmod,*fmmodexp);
    feffnew->Update(1.0,*fmmodexp,1.0);

    //--------------------------------------------------------- THIRD LINE
    // add i subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fimodexp;
    if (iset)
    {
      fimodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fimod,*fimodexp);
      feffnew->Update(1.0,*fimodexp,1.0);
    }

    //-------------------------------------------------------- FOURTH LINE
    // add weighted gap vector to feffnew, if existing
    Teuchos::RCP<Epetra_Vector> gexp;
    Teuchos::RCP<Epetra_Vector> fwexp;
    Teuchos::RCP<Epetra_Vector> fgmodexp;

    if (aset)
    {
      gexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*gact,*gexp);
      feffnew->Update(-1.0,*gexp,1.0);
    }

    //--------------------------------------------------------- FIFTH LINE
    // add st subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fstmodexp;
    if (stickset)
    {
      fstmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fstmod,*fstmodexp);
      feffnew->Update(1.0,*fstmodexp,+1.0);
    }

    // add terms of linearization feffnew
     if (stickset)
     {
        Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
        LINALG::Export(*linstickRHS_,*linstickRHSexp);
        feffnew->Update(-1.0,*linstickRHSexp,1.0);
     }

    //--------------------------------------------------------- SIXTH LINE

    // add a subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fslmodexp;
    Teuchos::RCP<Epetra_Vector> fwslexp;
    Teuchos::RCP<Epetra_Vector> fslwmodexp;


    if (slipset)
    {
      fslmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fslmod,*fslmodexp);
      feffnew->Update(1.0,*fslmodexp,1.0);
    }

    if (slipset)
    {
      Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*linslipRHS_,*linslipRHSexp);
      feffnew->Update(-1.0,*linslipRHSexp,1.0);
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
    if (Dualquadslave3d())
    {
      // modify dmatrix_
      Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
      dmatrix_    = temp2;
    }

    // transform if necessary
    if (ParRedist())
    {
      lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
      linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
    }

    // add contact stiffness
    kteff->UnComplete();
    kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
    kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
    kteff->Complete();

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold^T * zold / D^T * z to fit
    // thus we have to export the products Mold^T * zold / M^T * z to fit
    if (IsSelfContact())
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecd  = Teuchos::rcp(new Epetra_Vector(dmatrix_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zexp  = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
      if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
      dmatrix_->Multiply(true,*zexp,*tempvecd);
      LINALG::Export(*tempvecd,*fsexp);
      feff->Update(-(1.0-alphaf_),*fsexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecm  = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));
      mmatrix_->Multiply(true,*zexp,*tempvecm);
      LINALG::Export(*tempvecm,*fmexp);
      feff->Update(1.0-alphaf_,*fmexp,1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecdold  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      dold_->Multiply(true,*zoldexp,*tempvecdold);
      LINALG::Export(*tempvecdold,*fsoldexp);
      feff->Update(-alphaf_,*fsoldexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecmold  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
      mold_->Multiply(true,*zoldexp,*tempvecmold);
      LINALG::Export(*tempvecmold,*fmoldexp);
      feff->Update(alphaf_,*fmoldexp,1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dmatrix_->Multiply(true,*z_,*fs);
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fs,*fsexp);
      feff->Update(-(1.0-alphaf_),*fsexp,1.0);

      Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mmatrix_->Multiply(true,*z_,*fm);
      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fm,*fmexp);
      feff->Update(1.0-alphaf_,*fmexp,1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dold_->Multiply(true,*zold_,*fsold);
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fsold,*fsoldexp);
      feff->Update(-alphaf_,*fsoldexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->Multiply(true,*zold_,*fmold);
      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fmold,*fmoldexp);
      feff->Update(alphaf_,*fmoldexp,1.0);
    }
  }

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSLIPINCR
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckSlipIncrDerivTXI();
    if (Dim()==3)
      interface_[i]->FDCheckSlipIncrDerivTETA();

  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
    // FD check of stick condition
    for (int i=0; i<(int)interface_.size(); ++i)
    {
//      Teuchos::RCP<LINALG::SparseMatrix> deriv1 = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,81));
//      Teuchos::RCP<LINALG::SparseMatrix> deriv2 = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,81));
//
//      deriv1->Add(*linstickLM_,false,1.0,1.0);
//      deriv1->Complete(*gsmdofrowmap_,*gactivet_);
//
//      deriv2->Add(*linstickDIS_,false,1.0,1.0);
//      deriv2->Complete(*gsmdofrowmap_,*gactivet_);
//
//      std::cout << "DERIV 1 *********** "<< *deriv1 << std::endl;
//      std::cout << "DERIV 2 *********** "<< *deriv2 << std::endl;

      interface_[i]->FDCheckStickDeriv(*linstickLM_,*linstickDIS_);
    }
  }
#endif // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
    // FD check of slip condition
    for (int i=0; i<(int)interface_.size(); ++i)
    {
//      Teuchos::RCP<LINALG::SparseMatrix> deriv1 = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,81));
//      Teuchos::RCP<LINALG::SparseMatrix> deriv2 = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,81));
//
//      deriv1->Add(*linslipLM_,false,1.0,1.0);
//      deriv1->Complete(*gsmdofrowmap_,*gslipt_);
//
//      deriv2->Add(*linslipDIS_,false,1.0,1.0);
//      deriv2->Complete(*gsmdofrowmap_,*gslipt_);
//
//      std::cout << *deriv1 << std::endl;
//      std::cout << *deriv2 << std::endl;

      interface_[i]->FDCheckSlipDeriv(*linslipLM_,*linslipDIS_);
    }
  }
#endif // #ifdef CONTACTFDSLIP

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                  Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
    }
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
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
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleTN(tmatrix_,nmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleTNderiv(tderivmatrix_,nderivmatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);

    if (SystemType() != INPAR::CONTACT::system_condensed)
    {
      interface_[i]->AssembleInactiverhs(*inactiverhs_);
      interface_[i]->AssembleTangrhs(*tangrhs_);
    }
  }
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    // FillComplete() global matrix T
    tmatrix_->Complete(*gactivedofs_,*gactivedofs_);

    // FillComplete() global matrix N
    if (nmatrix_ != Teuchos::null)
      nmatrix_->Complete(*gactivedofs_,*gactivedofs_);
    smatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
    tderivmatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
    if (nderivmatrix_ != Teuchos::null)
      nderivmatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
  }
  else
  {
    // FillComplete() global matrix T
    tmatrix_->Complete(*gactivedofs_,*gactivet_);

    // FillComplete() global matrix N
    if (nmatrix_ != Teuchos::null)
      nmatrix_->Complete(*gactivedofs_,*gactiven_);

    // FillComplete() global matrix S
    smatrix_->Complete(*gsmdofrowmap_,*gactiven_);

    // FillComplete() global matrix Tderiv
    // (actually gsdofrowmap_ is in general sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    tderivmatrix_->Complete(*gsmdofrowmap_,*gactivet_);

    // FillComplete() global matrix Nderiv
    if (nderivmatrix_ != Teuchos::null)
      nderivmatrix_->Complete(*gsmdofrowmap_,*gactiven_);
  }

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify lindmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*invtrafo_,true,*lindmatrix_,false,false,false,true);
    lindmatrix_   = temp1;
  }

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"LM_SHAPEFCN");

  // do reagularization scaling
  if(regularized_)
    EvaluateRegularizationScaling(gact);

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype_ == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    /**********************************************************************/
    /* (1) Multiply Mortar matrices: m^ = inv(d) * m                      */
    /**********************************************************************/

    Teuchos::RCP<LINALG::SparseMatrix> invd = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix_));
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
    int err = 0;

    // extract diagonal of invd into diag
    invd->ExtractDiagonalCopy(*diag);

    // set zero diagonal values to dummy 1.0
    for (int i=0;i<diag->MyLength();++i)
      if ((*diag)[i]==0.0) (*diag)[i]=1.0;

    // scalar inversion of diagonal values
    err = diag->Reciprocal(*diag);
    if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

    Teuchos::RCP<Epetra_Vector> lmDBC=LINALG::CreateVector(*gsdofrowmap_,true);
    LINALG::Export(*pgsdirichtoggle_,*lmDBC);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*gsdofrowmap_,true);
    tmp->Multiply(1.,*diag,*lmDBC,0.);
    diag->Update(-1.,*tmp,1.);

    // re-insert inverted diagonal into invd
    err = invd->ReplaceDiagonalValues(*diag);
    // we cannot use this check, as we deliberately replaced zero entries
    //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

    // do the multiplication mhat = inv(D) * M
    mhatmatrix_ = LINALG::MLMultiply(*invd,false,*mmatrix_,false,false,false,true);

    /**********************************************************************/
    /* (2) Add contact stiffness terms to kteff                           */
    /**********************************************************************/

    // transform if necessary
    if (ParRedist())
    {
      lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
      linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
    }

    kteff->UnComplete();
    kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
    kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
    kteff->Complete();

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
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

    // split into slave/master part + structure part
    Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);
    if (ParRedist())
    {
      // split and transform to redistributed maps
      LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
      ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
      ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
      knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
    }
    else
    {
      // only split, no need to transform
      LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    }

    // further splits into slave part + master part
    LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
    LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
    LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

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
      // split and transform to redistributed maps
      LINALG::SplitVector(*ProblemDofs(),*feff,pgsmdofrowmap_,fsm,gndofrowmap_,fn);
      Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
      LINALG::Export(*fsm,*fsmtemp);
      fsm = fsmtemp;
    }
    else
    {
      // only split, no need to transform
      LINALG::SplitVector(*ProblemDofs(),*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);
    }

    // abbreviations for slave  and master set
    const int sset = gsdofrowmap_->NumGlobalElements();
    const int mset = gmdofrowmap_->NumGlobalElements();

    // we want to split fsm into 2 groups s,m
    fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

    // do the vector splitting sm -> s+m
    LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

    // store some stuff for static condensation of LM
    fs_   = fs;
    invd_ = invd;
    ksn_  = ksn;
    ksm_  = ksm;
    kss_  = kss;

    //----------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //----------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // D         ---->   D * T^(-1)
    // D^(-1)    ---->   T * D^(-1)
    // \hat{M}   ---->   T * \hat{M}
    //----------------------------------------------------------------------
    if (Dualquadslave3d())
    {
      // modify dmatrix_, invd_ and mhatmatrix_
      Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
      Teuchos::RCP<LINALG::SparseMatrix> temp3 = LINALG::MLMultiply(*trafo_,false,*invd_,false,false,false,true);
      Teuchos::RCP<LINALG::SparseMatrix> temp4 = LINALG::MLMultiply(*trafo_,false,*mhatmatrix_,false,false,false,true);
      dmatrix_    = temp2;
      invd_       = temp3;
      mhatmatrix_ = temp4;
    }

    /**********************************************************************/
    /* (5) Split slave quantities into active / inactive                  */
    /**********************************************************************/

    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

    // we will get the i rowmap as a by-product
    Teuchos::RCP<Epetra_Map> gidofs;

    // do the splitting
    LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
    LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
    LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
    LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

    // abbreviations for active and inactive set
    const int aset = gactivedofs_->NumGlobalElements();
    const int iset = gidofs->NumGlobalElements();

    // we want to split fsmod into 2 groups a,i
    Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // do the vector splitting s -> a+i
    LINALG::SplitVector(*gsdofrowmap_,*fs,gactivedofs_,fa,gidofs,fi);

    /**********************************************************************/
    /* (6) Isolate necessary parts from invd and mhatmatrix               */
    /**********************************************************************/

    // active part of invd
    Teuchos::RCP<LINALG::SparseMatrix> invda;
    LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);

    // coupling part of dmatrix (only nonzero for 3D quadratic case!)
    Teuchos::RCP<LINALG::SparseMatrix> dai;
    LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,tempmtx1,dai,tempmtx2,tempmtx3);

     // do the multiplication dhat = invda * dai
    Teuchos::RCP<LINALG::SparseMatrix> dhat = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
    if (aset && iset) dhat = LINALG::MLMultiply(*invda,false,*dai,false,false,false,true);
    dhat->Complete(*gidofs,*gactivedofs_);

    // active part of mmatrix
    Teuchos::RCP<LINALG::SparseMatrix> mmatrixa;
    LINALG::SplitMatrix2x2(mmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mmatrixa,tempmtx1,tempmtx2,tempmtx3);

    // do the multiplication mhataam = invda * mmatrixa
    // (this is only different from mhata for 3D quadratic case!)
    Teuchos::RCP<LINALG::SparseMatrix> mhataam = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
    if (aset) mhataam = LINALG::MLMultiply(*invda,false,*mmatrixa,false,false,false,true);
    mhataam->Complete(*gmdofrowmap_,*gactivedofs_);

    // for the case without full linearization, we still need the
    // "classical" active part of mhat, which is isolated here
    Teuchos::RCP<LINALG::SparseMatrix> mhata;
    LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);

    // scaling of invd and dai
    invda->Scale(1/(1-alphaf_));
    dai->Scale(1-alphaf_);

    SaveCouplingMatrices(dhat, mhataam, invda);
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
    Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhataam,true,*kan,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

    // kmm: add T(mhataam)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmmmod->Add(*kmm,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhataam,true,*kam,false,false,false,true);
    kmmmod->Add(*kmmadd,false,1.0,1.0);
    kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

    // kmi: add T(mhataam)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmimod->Add(*kmi,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhataam,true,*kai,false,false,false,true);
      kmimod->Add(*kmiadd,false,1.0,1.0);
      kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
    }

    // kma: add T(mhataam)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
      kmamod->Add(*kma,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhataam,true,*kaa,false,false,false,true);
      kmamod->Add(*kmaadd,false,1.0,1.0);
      kmamod->Complete(kma->DomainMap(),kma->RowMap());
    }

    //----------------------------------------------------------- THIRD LINE
    //------------------- FOR 3D QUADRATIC CASE ----------------------------
    // kin: subtract T(dhat)*kan --
    Teuchos::RCP<LINALG::SparseMatrix> kinmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kinmod->Add(*kin,false,1.0,1.0);
    if (aset && iset)
    {
      Teuchos::RCP<LINALG::SparseMatrix> kinadd = LINALG::MLMultiply(*dhat,true,*kan,false,false,false,true);
      kinmod->Add(*kinadd,false,-1.0,1.0);
    }
    kinmod->Complete(kin->DomainMap(),kin->RowMap());

    // kim: subtract T(dhat)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kimmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kimmod->Add(*kim,false,1.0,1.0);
    if (aset && iset)
    {
      Teuchos::RCP<LINALG::SparseMatrix> kimadd = LINALG::MLMultiply(*dhat,true,*kam,false,false,false,true);
      kimmod->Add(*kimadd,false,-1.0,1.0);
    }
    kimmod->Complete(kim->DomainMap(),kim->RowMap());

    // kii: subtract T(dhat)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kiimod;
    if (iset)
    {
      kiimod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiimod->Add(*kii,false,1.0,1.0);
      if(aset)
      {
        Teuchos::RCP<LINALG::SparseMatrix> kiiadd = LINALG::MLMultiply(*dhat,true,*kai,false,false,false,true);
        kiimod->Add(*kiiadd,false,-1.0,1.0);
      }
      kiimod->Complete(kii->DomainMap(),kii->RowMap());
    }

    // kia: subtract T(dhat)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kiamod;
    if (iset && aset)
    {
      kiamod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
      kiamod->Add(*kia,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kiaadd = LINALG::MLMultiply(*dhat,true,*kaa,false,false,false,true);
      kiamod->Add(*kiaadd,false,-1.0,1.0);
      kiamod->Complete(kia->DomainMap(),kia->RowMap());
    }

    //---------------------------------------------------------- FOURTH LINE
    // nothing to do - execpt its regularized contact see --> DoRegularizationScaling()

    //----------------------------------------------------------- FIFTH LINE
    // kan: multiply tmatrix with invda and kan
    Teuchos::RCP<LINALG::SparseMatrix> kanmod;
    if (aset)
    {
      kanmod = LINALG::MLMultiply(*tmatrix_,false,*invda,true,false,false,true);
      kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
    }

    // kam: multiply tmatrix with invda and kam
    Teuchos::RCP<LINALG::SparseMatrix> kammod;
    if (aset)
    {
      kammod = LINALG::MLMultiply(*tmatrix_,false,*invda,true,false,false,true);
      kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
    }

    // kai: multiply tmatrix with invda and kai
    Teuchos::RCP<LINALG::SparseMatrix> kaimod;
    if (aset && iset)
    {
      kaimod = LINALG::MLMultiply(*tmatrix_,false,*invda,true,false,false,true);
      kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
    }

    // kaa: multiply tmatrix with invda and kaa
    Teuchos::RCP<LINALG::SparseMatrix> kaamod;
    if (aset)
    {
      kaamod = LINALG::MLMultiply(*tmatrix_,false,*invda,true,false,false,true);
      kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
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
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(mold_->RowMap()));
      if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      mold_->Multiply(true,*zoldexp,*tempvecm2);
      if (mset) LINALG::Export(*tempvecm2,*tempvecm);
      fm->Update(alphaf_,*tempvecm,1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->Multiply(true,*zold_,*tempvecm);
      fm->Update(alphaf_,*tempvecm,1.0);
    }

    // fs: prepare alphaf * old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

    // for self contact, slave and master sets may have changed,
    // thus we have to export the product Dold^T * zold to fit
    if (IsSelfContact())
    {
      Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      dold_->Multiply(true,*zoldexp,*tempvec);
      if (sset) LINALG::Export(*tempvec,*fsadd);
    }
    // if there is no self contact everything is ok
    else
    {
      dold_->Multiply(true,*zold_,*fsadd);
    }

    // fa: subtract alphaf * old contact forces (t_n)
    if (aset)
    {
      Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      LINALG::Export(*fsadd,*faadd);
      fa->Update(-alphaf_,*faadd,1.0);
    }

    // fm: add T(mhat)*fa
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (aset) mhataam->Multiply(true,*fa,*fmmod);
    fmmod->Update(1.0,*fm,1.0);

    //----------------------------------------------------------- THIRD LINE
    // fi: subtract alphaf * old contact forces (t_n)
    if (iset)
    {
      Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
      LINALG::Export(*fsadd,*fiadd);
      fi->Update(-alphaf_,*fiadd,1.0);
    }

    // fi: add T(dhat)*fa
    Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
    if (iset && aset) dhat->Multiply(true,*fa,*fimod);
    fimod->Update(1.0,*fi,-1.0);

    //---------------------------------------------------------- FOURTH LINE
    // gactive: nothing to do  - execpt its regularized contact see --> DoRegularizationScaling()

    //----------------------------------------------------------- FIFTH LINE
    // fa: multiply tmatrix with invda and fa
    Teuchos::RCP<Epetra_Vector> famod;
    Teuchos::RCP<LINALG::SparseMatrix> tinvda;
    if (aset)
    {
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        famod = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    else
      famod = Teuchos::rcp(new Epetra_Vector(*gactivet_));

    tinvda = LINALG::MLMultiply(*tmatrix_,false,*invda,true,false,false,true);
      tinvda->Multiply(false,*fa,*famod);
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
        tderivmatrix_ = MORTAR::MatrixRowTransform(tderivmatrix_,pgsdofrowmap_);
      }
    }

    /**********************************************************************/
    /* (10) Global setup of kteffnew (including contact)                  */
    /**********************************************************************/

    Teuchos::RCP<LINALG::SparseMatrix> kteffnew = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(),81,true,false,kteffmatrix->GetMatrixtype()));
    Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*ProblemDofs());

    //Add all additional contributions from regularized contact to kteffnew and feffnew!
    if (regularized_)
      DoRegularizationScaling(aset, iset, invda, kan, kam, kai, kaa, fa, kteffnew, feffnew);

    //----------------------------------------------------------- FIRST LINE
    // add n submatrices to kteffnew
    kteffnew->Add(*knn,false,1.0,1.0);
    kteffnew->Add(*knm,false,1.0,1.0);
    if (sset) kteffnew->Add(*kns,false,1.0,1.0);

    //---------------------------------------------------------- SECOND LINE
    // add m submatrices to kteffnew
    kteffnew->Add(*kmnmod,false,1.0,1.0);
    kteffnew->Add(*kmmmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kmimod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);

    //----------------------------------------------------------- THIRD LINE
    // add i submatrices to kteffnew
    if (iset) kteffnew->Add(*kinmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kimmod,false,1.0,1.0);
    if (iset) kteffnew->Add(*kiimod,false,1.0,1.0);
    if (iset && aset) kteffnew->Add(*kiamod,false,1.0,1.0);

    //---------------------------------------------------------- FOURTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->Add(*smatrix_,false,1.0,1.0);

    //----------------------------------------------------------- FIFTH LINE
    // add a submatrices to kteffnew
    if (aset) kteffnew->Add(*kanmod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kammod,false,1.0,1.0);
    if (aset && iset) kteffnew->Add(*kaimod,false,1.0,1.0);
    if (aset) kteffnew->Add(*kaamod,false,1.0,1.0);
    if (aset) kteffnew->Add(*tderivmatrix_,false,-1.0,1.0);

    // FillComplete kteffnew (square)
    kteffnew->Complete();

    /**********************************************************************/
    /* (11) Global setup of feffnew (including contact)                   */
    /**********************************************************************/

    //----------------------------------------------------------- FIRST LINE
    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fn,*fnexp);
    feffnew->Update(1.0,*fnexp,1.0);

    //---------------------------------------------------------- SECOND LINE
    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fmmod,*fmmodexp);
    feffnew->Update(1.0,*fmmodexp,1.0);

    //----------------------------------------------------------- THIRD LINE
    // add i subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fimodexp;
    if (iset)
    {
      fimodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fimod,*fimodexp);
      feffnew->Update(1.0,*fimodexp,1.0);
    }

    //---------------------------------------------------------- FOURTH LINE
        // add weighted gap vector to feffnew, if existing
    Teuchos::RCP<Epetra_Vector> gexp;
    if (aset)
    {
      gexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*gact,*gexp);
      feffnew->Update(-1.0,*gexp,1.0);
    }

    //----------------------------------------------------------- FIFTH LINE
    // add a subvector to feffnew
    Teuchos::RCP<Epetra_Vector> famodexp;
    if (aset)
    {
      famodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*famod,*famodexp);
      feffnew->Update(1.0,*famodexp,1.0);
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
    if (Dualquadslave3d())
    {
      // modify dmatrix_
      Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
      dmatrix_    = temp2;
    }

    // transform if necessary
    if (ParRedist())
    {
      lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
      linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
    }

    // add contact stiffness
    kteff->UnComplete();
    kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
    kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
    kteff->Complete();

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold^T * zold / D^T * z to fit
    // thus we have to export the products Mold^T * zold / M^T * z to fit
    if (IsSelfContact())
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecd  = Teuchos::rcp(new Epetra_Vector(dmatrix_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zexp  = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
      if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
      dmatrix_->Multiply(true,*zexp,*tempvecd);
      LINALG::Export(*tempvecd,*fsexp);
      feff->Update(-(1.0-alphaf_),*fsexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecm  = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));
      mmatrix_->Multiply(true,*zexp,*tempvecm);
      LINALG::Export(*tempvecm,*fmexp);
      feff->Update(1.0-alphaf_,*fmexp,1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecdold  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
      Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
      dold_->Multiply(true,*zoldexp,*tempvecdold);
      LINALG::Export(*tempvecdold,*fsoldexp);
      feff->Update(-alphaf_,*fsoldexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      Teuchos::RCP<Epetra_Vector> tempvecmold  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
      mold_->Multiply(true,*zoldexp,*tempvecmold);
      LINALG::Export(*tempvecmold,*fmoldexp);
      feff->Update(alphaf_,*fmoldexp,1.0);
    }
    // if there is no self contact everything is ok
    else
    {
      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dmatrix_->Multiply(true,*z_,*fs);
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fs,*fsexp);
      feff->Update(-(1.0-alphaf_),*fsexp,1.0);

      Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mmatrix_->Multiply(true,*z_,*fm);
      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fm,*fmexp);
      feff->Update(1.0-alphaf_,*fmexp,1.0);

      // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
      CheckConservationLaws(*fs,*fm);
#endif

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dold_->Multiply(true,*zold_,*fsold);
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fsold,*fsoldexp);
      feff->Update(-alphaf_,*fsoldexp,1.0);

      Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      mold_->Multiply(true,*zold_,*fmold);
      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fmold,*fmoldexp);
      feff->Update(alphaf_,*fmoldexp,1.0);
    }
  }

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDALPHA
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckAlphaDeriv();
  }
#endif // #ifdef CONTACTFDGAP
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i=0; i<(int)interface_.size();++i)
  {
    std::cout << *tderivmatrix_ << std::endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM

  return;
}

/*----------------------------------------------------------------------*
 | Setup 2x2 saddle point system for contact problems      wiesner 11/14|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
                                        Teuchos::RCP<Epetra_Vector> fd,
                                        Teuchos::RCP<Epetra_Vector> sold,
                                        Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
                                        int numiter,
                                        Teuchos::RCP<Epetra_Operator>& blockMat,
                                        Teuchos::RCP<Epetra_Vector>& blocksol,
                                        Teuchos::RCP<Epetra_Vector>& blockrhs)
{
  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->FullMap())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp,*dirichtoggle);

  //**********************************************************************
  // prepare saddle point system
  //**********************************************************************
  // the standard stiffness matrix
  Teuchos::RCP<LINALG::SparseMatrix> stiffmt = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kdd);

  // initialize merged system (matrix, rhs, sol)
  Teuchos::RCP<Epetra_Map>           mergedmap   = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
  Teuchos::RCP<LINALG::SparseMatrix> mergedmt    = Teuchos::null;
  Teuchos::RCP<Epetra_Vector>        mergedrhs   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedsol   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedzeros = LINALG::CreateVector(*mergedmap);

  // initialize constraint matrices
  Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,true));
  Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
  Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));

  // initialize transformed constraint matrices
  Teuchos::RCP<LINALG::SparseMatrix> trkdz, trkzd, trkzz;

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // build constraint matrix kdz
  kdz->Add(*dmatrix_,true,1.0-alphaf_,1.0);
  kdz->Add(*mmatrix_,true,-(1.0-alphaf_),1.0);
  kdz->Complete(*gsdofrowmap_,*gdisprowmap_);

  // *** CASE 1: FRICTIONLESS CONTACT ************************************
  if (!friction_)
  {
    // build constraint matrix kzd
    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      if (gactivedofs_->NumGlobalElements())
      {
        kzd->Add(*smatrix_,false,1.0,1.0);
        kzd->Add(*tderivmatrix_,false,1.0,1.0);
      }
    }
    else
    {
      if (gactiven_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
      if (gactivet_->NumGlobalElements()) kzd->Add(*tderivmatrix_,false,1.0,1.0);
    }
    kzd->Complete(*gdisprowmap_,*gsdofrowmap_);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kzz->Add(*tmatrix_,false,1.0,1.0);
    kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);
  }

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 2: FRICTIONAL CONTACT **************************************
  else
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
    Teuchos::RCP<Epetra_Map> gstickdofs=LINALG::SplitMap(*gactivedofs_,*gslipdofs_);

    // build constraint matrix kzd
    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      if (gactivedofs_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
      if (gstickdofs->NumGlobalElements()) kzd->Add(*linstickDIS_,false,1.0,1.0);
      if (gslipdofs_->NumGlobalElements()) kzd->Add(*linslipDIS_,false,1.0,1.0);
    }
    else
    {
      if (gactiven_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
      if (gstickt->NumGlobalElements()) kzd->Add(*linstickDIS_,false,1.0,1.0);
      if (gslipt_->NumGlobalElements()) kzd->Add(*linslipDIS_,false,1.0,1.0);
    }
    kzd->Complete(*gdisprowmap_,*gsdofrowmap_);;

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
      if (gstickdofs->NumGlobalElements()) kzz->Add(*linstickLM_,false,1.0,1.0);
      if (gslipdofs_->NumGlobalElements()) kzz->Add(*linslipLM_,false,1.0,1.0);
    }
    else
    {
      if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
      if (gstickt->NumGlobalElements()) kzz->Add(*linstickLM_,false,1.0,1.0);
      if (gslipt_->NumGlobalElements()) kzz->Add(*linslipLM_,false,1.0,1.0);
    }
    kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);
  }

  //**********************************************************************
  // transform matrix blocks to lm-dof-map
  //**********************************************************************

  // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
  trkzd = MORTAR::MatrixRowTransformGIDs(kzd,glmdofrowmap_);

  // transform parallel column distribution of constraint matrix kzd
  // (only necessary in the parallel redistribution case)
  if (ParRedist()) trkzd = MORTAR::MatrixColTransform(trkzd,ProblemDofs());

  // transform constraint matrix kzz to lmdofmap (MatrixRowColTransform)
  trkzz = MORTAR::MatrixRowColTransformGIDs(kzz,glmdofrowmap_,glmdofrowmap_);

  // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
  trkdz = MORTAR::MatrixColTransformGIDs(kdz,glmdofrowmap_);

  // transform parallel row distribution of constraint matrix kdz
  // (only necessary in the parallel redistribution case)
  if (ParRedist()) trkdz = MORTAR::MatrixRowTransform(trkdz,ProblemDofs());

  //**********************************************************************
  // build and solve saddle point system
  //**********************************************************************
  if (SystemType()==INPAR::CONTACT::system_saddlepoint)
  {
    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*dirichtoggle,*dirichtoggleexp);
    Teuchos::RCP<Epetra_Vector> lmDBC = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*dirichtoggle,*lmDBC);
    lmDBC->ReplaceMap(*glmdofrowmap_);
    Teuchos::RCP<Epetra_Vector> lmDBCexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*lmDBC,*lmDBCexp);
    dirichtoggleexp->Update(1.,*lmDBCexp,1.);
    trkzd->ApplyDirichlet(lmDBC,false);

    trkzz->Complete();
    trkzz->ApplyDirichlet(lmDBC,true);

    // apply Dirichlet conditions to (0,0) and (0,1) blocks
    Teuchos::RCP<Epetra_Vector> zeros   = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));
    Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
    LINALG::ApplyDirichlettoSystem(stiffmt,sold,rhscopy,zeros,dirichtoggle);
    trkdz->ApplyDirichlet(dirichtoggle,false);

    // row map (equals domain map) extractor
    LINALG::MapExtractor rowmapext(*mergedmap,glmdofrowmap_,ProblemDofs());
    LINALG::MapExtractor dommapext(*mergedmap,glmdofrowmap_,ProblemDofs());

    // set a helper flag for the CheapSIMPLE preconditioner (used to detect, if Teuchos::nullspace has to be set explicitly)
    // do we need this? if we set the Teuchos::nullspace when the solver is constructed?
    //solver.Params().set<bool>("CONTACT",true); // for simpler precond

    // build block matrix for SIMPLER
    blockMat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> >(blockMat);
    mat->Assign(0,0,LINALG::View,*stiffmt);
    mat->Assign(0,1,LINALG::View,*trkdz);
    mat->Assign(1,0,LINALG::View,*trkzd);
    mat->Assign(1,1,LINALG::View,*trkzz);
    mat->Complete();

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*fd,*fresmexp);
    mergedrhs->Update(1.0,*fresmexp,1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*constrrhs_,*constrexp);
    mergedrhs->Update(1.0,*constrexp,1.0);

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    LINALG::ApplyDirichlettoSystem(mergedsol,mergedrhs,mergedzeros,dirichtoggleexp);

    blocksol = mergedsol;
    blockrhs = mergedrhs;
    return;
  }

  //**********************************************************************
  // invalid system types
  //**********************************************************************
  else dserror("ERROR: Invalid system type in SaddlePointSolve");

  return;
}

/*------------------------------------------------------------------------*
 | Update internal member variables after saddle point solve wiesner 11/14|
 *------------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol)
{
  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
  LINALG::MapExtractor mapext(*mergedmap,ProblemDofs(),glmdofrowmap_);
  mapext.ExtractCondVector(blocksol,sold);
  mapext.ExtractOtherVector(blocksol,sollm);
  sollm->ReplaceMap(*gsdofrowmap_);

  if (IsSelfContact())
  // for self contact, slave and master sets may have changed,
  // thus we have to reinitialize the LM vector map
  {
    zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
    LINALG::Export(*z_, *zincr_);                     // change the map of z_
    z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
    zincr_->Update(1.0, *sollm, 0.0);                 // save sollm in zincr_
    z_->Update(1.0, *zincr_, 1.0);                    // update z_
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
void CONTACT::CoLagrangeStrategy::EvalConstrRHS()
{
  if (SystemType()==INPAR::CONTACT::system_condensed)
    return;

  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    constrrhs_          = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) * z^(i)_(n+1)

  // export weighted gap vector
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
    }
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gactiven_->NumGlobalElements())
    {
      LINALG::Export(*g_,*gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  Teuchos::RCP<Epetra_Vector> gact_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*gact,*gact_exp);

  constrrhs->Update(-1.0,*gact_exp,1.0);

  // export inactive rhs
  Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*inactiverhs_, *inactiverhsexp);

  // build constraint rhs (1)
  constrrhs->Update(1.0, *inactiverhsexp, 1.0);

  // *** CASE 1: FRICTIONLESS CONTACT *******************************************************
  if (!friction_)
  {
    // export tangential rhs
    Teuchos::RCP<Epetra_Vector> tangrhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*tangrhs_, *tangrhs_exp);

    // build constraint rhs (2)
    constrrhs->Update(1.0, *tangrhs_exp, 1.0);

  }
  // *** CASE 2: FRICTIONAL CONTACT *******************************************************
  else
  {
    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linstickRHS_,*stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linslipRHS_,*slipexp);

    // build constraint rhs
    constrrhs->Update(1.0,*stickexp,1.0);
    constrrhs->Update(1.0,*slipexp,1.0);
  }

  constrrhs->ReplaceMap(*glmdofrowmap_);

  // export and set constraint rhs vector
  if (ParRedist())
  {
    constrrhs_ = Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(false)));
    LINALG::Export(*constrrhs,*constrrhs_);
  }
  else
    constrrhs_ = constrrhs;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::EvalForce(
    CONTACT::ParamsInterface& cparams)
{
  //---------------------------------------------------------------
  // For selfcontact the master/slave sets are updated within the -
  // contact search, see SelfBinaryTree.                          -
  // Therefore, we have to initialize the mortar matrices after   -
  // interface evaluations.                                       -
  //---------------------------------------------------------------
  if (IsSelfContact())
  {
    InitEvalInterface(); // evaluate mortar terms (integrate...)
    InitMortar();        // initialize mortar matrices and vectors
    AssembleMortar();    // assemble mortar terms into global matrices
  }
  else
  {
    InitMortar();        // initialize mortar matrices and vectors
    InitEvalInterface(); // evaluate mortar terms (integrate...)
    AssembleMortar();    // assemble mortar terms into global matrices
  }

  // evaluate relative movement for friction
  if (cparams.IsPredictor())
    EvaluateRelMovPredict();
  else
    EvaluateRelMov();

  // update active set
  UpdateActiveSetSemiSmooth();

  // apply contact forces and stiffness
  Initialize();         // init lin-matrices

  // --- Assemble the ride hand side terms -------------------------------
  AssembleContactRHS();

  EvalStrContactRHS();  // evaluate the structure/displacement rhs
  EvalConstrRHS();      // evaluate the constraint rhs (saddle-point system only)

  if(constrrhs_ != Teuchos::null)
    constrrhs_->Scale(-1.0);// scale with -1.0 --> when old structure is deleted change this!!!

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::AssembleContactRHS()
{

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (SystemType() != INPAR::CONTACT::system_condensed)
    {
      interface_[i]->AssembleInactiverhs(*inactiverhs_);
      interface_[i]->AssembleTangrhs(*tangrhs_);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    strcontactrhs_ = Teuchos::null;
    return;
  }

  strcontactrhs_ =
      Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (IsSelfContact())
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecd  = Teuchos::rcp(new Epetra_Vector(dmatrix_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zexp  = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
    dmatrix_->Multiply(true,*zexp,*tempvecd);
    LINALG::Export(*tempvecd,*fsexp);
    strcontactrhs_->Update(-1.0,*fsexp,1.0);

    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecm  = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));
    mmatrix_->Multiply(true,*zexp,*tempvecm);
    LINALG::Export(*tempvecm,*fmexp);
    strcontactrhs_->Update(1.0,*fmexp,1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*z_,*fs);
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fs,*fsexp);
    strcontactrhs_->Update(-1.0,*fsexp,1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*z_,*fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fm,*fmexp);
    strcontactrhs_->Update(1.0,*fmexp,1.0);

    // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
    CheckConservationLaws(*fs,*fm);
#endif
  }
  strcontactrhs_->Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------*
 | set force evaluation flag before evaluation step          farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::PreEvaluate(
    CONTACT::ParamsInterface& cparams)
{
  const enum MORTAR::ActionType& act = cparams.GetActionType();

  switch(act)
  {
  // -------------------------------------------------------------------
  // reset force evaluation flag for predictor step
  // -------------------------------------------------------------------
    case MORTAR::eval_force_stiff:
    {
      if(cparams.IsPredictor())
        evalForceCalled_ = false;
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
void CONTACT::CoLagrangeStrategy::PostEvaluate(
    CONTACT::ParamsInterface& cparams)
{
  const enum MORTAR::ActionType& act = cparams.GetActionType();

  switch(act)
  {
    // -------------------------------------------------------------------
    // set flag to false after force stiff evaluation
    // -------------------------------------------------------------------
    case MORTAR::eval_force_stiff:
    {
      evalForceCalled_ = false;
      break;
    }
    // -------------------------------------------------------------------
    // set flag for force evaluation to true
    // -------------------------------------------------------------------
    case MORTAR::eval_force:
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
void CONTACT::CoLagrangeStrategy::EvalForceStiff(
    CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine if not done before
  if(!evalForceCalled_)
    EvalForce(cparams);

  // --- Assemble stiffness matrix ---------------------------------------
  AssembleContactStiff();

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::AssembleContactStiff()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
    }
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
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
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleTN(tmatrix_,nmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleTNderiv(tderivmatrix_,nderivmatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
  }
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    // FillComplete() global matrix T
    tmatrix_->Complete(*gactivedofs_,*gactivedofs_);

    // FillComplete() global matrix N
    if (nmatrix_ != Teuchos::null)
      nmatrix_->Complete(*gactivedofs_,*gactivedofs_);
    smatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
    tderivmatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
    if (nderivmatrix_ != Teuchos::null)
      nderivmatrix_->Complete(*gsmdofrowmap_,*gactivedofs_);
  }
  else
  {
    // FillComplete() global matrix T
    tmatrix_->Complete(*gactivedofs_,*gactivet_);

    // FillComplete() global matrix N
    if (nmatrix_ != Teuchos::null)
      nmatrix_->Complete(*gactivedofs_,*gactiven_);

    // FillComplete() global matrix S
    smatrix_->Complete(*gsmdofrowmap_,*gactiven_);

    // FillComplete() global matrix Tderiv
    // (actually gsdofrowmap_ is in general sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    tderivmatrix_->Complete(*gsmdofrowmap_,*gactivet_);

    // FillComplete() global matrix Nderiv
    if (nderivmatrix_ != Teuchos::null)
      nderivmatrix_->Complete(*gsmdofrowmap_,*gactiven_);
  }

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify lindmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*invtrafo_,true,*lindmatrix_,false,false,false,true);
    lindmatrix_   = temp1;
  }

  // do reagularization scaling
  if(regularized_)
    EvaluateRegularizationScaling(gact);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoLagrangeStrategy::GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new LINALG::SparseMatrix(SlMaDoFRowMap(true),100,false,true));

      // build matrix kdd
      mat_ptr->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
      mat_ptr->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
      mat_ptr->Complete();

      // transform parallel row/column distribution of matrix kdd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            SlMaDoFRowMapPtr(false),SlMaDoFRowMapPtr(false));

      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      // build constraint matrix kdz
      Teuchos::RCP<LINALG::SparseMatrix> kdz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));

      kdz_ptr->Add(*dmatrix_,true,1.0-alphaf_,1.0);
      kdz_ptr->Add(*mmatrix_,true,-(1.0-alphaf_),1.0);
      kdz_ptr->Complete(*gsdofrowmap_,*gdisprowmap_);

      // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
      mat_ptr = MORTAR::MatrixColTransformGIDs(kdz_ptr,LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kdz
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            ProblemDofs(),LMDoFRowMapPtr(false));

      break;
    }
    case DRT::UTILS::block_lm_displ:
    {
      // build constraint matrix kzd
      Teuchos::RCP<LINALG::SparseMatrix> kzd_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));


      // build constraint matrix kzd
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        if (gactivedofs_->NumGlobalElements())
        {
          kzd_ptr->Add(*smatrix_,false,1.0,1.0);
          kzd_ptr->Add(*tderivmatrix_,false,1.0,1.0);
        }
      }
      else
      {
        if (gactiven_->NumGlobalElements())
          kzd_ptr->Add(*smatrix_,false,1.0,1.0);
        if (gactivet_->NumGlobalElements())
          kzd_ptr->Add(*tderivmatrix_,false,1.0,1.0);
      }
      kzd_ptr->Complete(*gdisprowmap_,*gsdofrowmap_);

      // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
      mat_ptr = MORTAR::MatrixRowTransformGIDs(kzd_ptr,LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kzd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            LMDoFRowMapPtr(false),ProblemDofs());

      break;
    }
    case DRT::UTILS::block_lm_lm:
    {
      // build constraint matrix kzz
      Teuchos::RCP<LINALG::SparseMatrix> kzz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));

      // build unity matrix for inactive dofs
      Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
      Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
      ones->PutScalar(1.0);
      Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
      onesdiag->Complete();

      // build constraint matrix kzz
      if (gidofs->NumGlobalElements())
        kzz_ptr->Add(*onesdiag,false,1.0,1.0);
      if (gactivet_->NumGlobalElements())
        kzz_ptr->Add(*tmatrix_,false,1.0,1.0);
      kzz_ptr->Complete(*gsdofrowmap_,*gsdofrowmap_);

      // transform constraint matrix kzz to lmdofmap
      mat_ptr = MORTAR::MatrixRowColTransformGIDs(kzz_ptr,
           LMDoFRowMapPtr(true),LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kzz
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            LMDoFRowMapPtr(false),LMDoFRowMapPtr(false));

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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::RecoverState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  Teuchos::RCP<Epetra_Vector> zdir_ptr =
      Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(true),true));
  LINALG::Export(dir,*zdir_ptr);
  // get the current step length
  const double stepLength = cparams.GetStepLength();
  // ---------------------------------------------------------------------
  // store the SCALED Lagrange multiplier increment in the contact
  // strategy
  // ---------------------------------------------------------------------
  zdir_ptr->ReplaceMap(zincr_->Map());
  zincr_->Scale(stepLength,*zdir_ptr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::CoLagrangeStrategy::GetRhsBlockPtr(
       const enum DRT::UTILS::VecBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      vec_ptr = strcontactrhs_;
      break;
    }
    case DRT::UTILS::block_constraint:
    {
      vec_ptr = constrrhs_;
      break;
    }
    default:
    {
      dserror("Unknown STR::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xnew)
{
  if (systype_ == INPAR::CONTACT::system_condensed)
  {
    dserror("ERROR: not yet implemented for condensed system!");
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> znew_ptr =
        Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(true),true));
    LINALG::Export(xnew,*znew_ptr);
    // ---------------------------------------------------------------------
    // Update the current lagrange multiplier
    // ---------------------------------------------------------------------
    znew_ptr->ReplaceMap(z_->Map());

    z_->Scale(1.0,*znew_ptr);

    // ---------------------------------------------------------------------
    // store the new Lagrange multiplier in the nodes
    // ---------------------------------------------------------------------
    StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::Recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (SystemType() == INPAR::CONTACT::system_condensed)
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
    if (IsSelfContact())
    {
      // approximate update
      //z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      //invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      z_->Update(1.0,*fs_,0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->Multiply(false,*disis,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksm_->Multiply(false,*disim,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksn_->Multiply(false,*disin,*mod);
      z_->Update(-1.0,*mod,1.0);
      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*mod2);
      Teuchos::RCP<Epetra_Vector> mod3 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
      dold_->Multiply(true,*mod2,*mod3);
      Teuchos::RCP<Epetra_Vector> mod4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*mod3,*mod4);
      z_->Update(-alphaf_,*mod4,1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->Multiply(true,*zcopy,*z_);
      z_->Scale(1/(1-alphaf_));
    }
    else
    {
      // approximate update
      //invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_->Update(1.0,*fs_,0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->Multiply(false,*disis,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksm_->Multiply(false,*disim,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksn_->Multiply(false,*disin,*mod);
      z_->Update(-1.0,*mod,1.0);
      dold_->Multiply(true,*zold_,*mod);
      z_->Update(-alphaf_,*mod,1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->Multiply(true,*zcopy,*z_);
      z_->Scale(1/(1-alphaf_));
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
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence               popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::UpdateActiveSet()
{
  // get input parameter ftype
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");

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
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];
        nzold += cnode->MoData().n()[k] * cnode->MoData().lmold()[k];
      }

      // friction
      double tz = 0.0;
      double tjump = 0.0;

      if (friction_)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute tangential part of Lagrange multiplier
        tz = frinode->CoData().txi()[0]*frinode->MoData().lm()[0] + frinode->CoData().txi()[1]*frinode->MoData().lm()[1];

        // compute tangential part of jump FIXME -- also the teta component should be considered
        if (DRT::INPUT::IntegralValue<int>(Params(),"GP_SLIP_INCR")==true)
          tjump = frinode->FriData().jump_var()[0];
        else
          tjump = frinode->CoData().txi()[0]*frinode->FriData().jump()[0] + frinode->CoData().txi()[1]*frinode->FriData().jump()[1];
      }

      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << std::endl;

        // check for penetration
        if (wgap < 0)
        {
          cnode->Active() = true;
          activesetconv_ = false;
        }
      }

      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << std::endl;

        // check for tensile contact forces
        if (nz <= 0) // no averaging of Lagrange multipliers
        //if (0.5*nz+0.5*nzold <= 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = false;
          activesetconv_ = false;

          // friction
          if (friction_) dynamic_cast<FriNode*>(cnode)->FriData().Slip() = false;
        }

        // only do something for friction
        else
        {
          // friction tresca
          if (ftype == INPAR::CONTACT::friction_tresca)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);
            double ct = interface_[i]->GetCtRef()[interface_[i]->GetCtRef().Map().LID(frinode->Id())];

            // CAREFUL: friction bound is now interface-local (popp 08/2012)
            double frbound = interface_[i]->IParams().get<double>("FRBOUND");

            if(frinode->FriData().Slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 frinode->FriData().Slip() = true;
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
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (ftype == INPAR::CONTACT::friction_tresca)

          // friction coulomb
          if (ftype == INPAR::CONTACT::friction_coulomb)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);
            double ct = interface_[i]->GetCtRef()[interface_[i]->GetCtRef().Map().LID(frinode->Id())];

            // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
            double frcoeff = interface_[i]->IParams().get<double>("FRCOEFF");

            if(frinode->FriData().Slip() == false)
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frcoeff*nz <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 frinode->FriData().Slip() = true;
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
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (ftype == INPAR::CONTACT::friction_coulomb)
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
    activesetsteps_ += 1;
  }

  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=Teuchos::null) zigzagthree_  = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=Teuchos::null) zigzagtwo_    = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=Teuchos::null) zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);

    if(friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
      gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
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
  if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)
  {
    // frictionless contact
    if (ActiveSetSteps()>2)
    {
      if (zigzagtwo_!=Teuchos::null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (Comm().MyPID()==0)
            std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << std::endl;
        }
      }

      if (zigzagthree_!=Teuchos::null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;

          // output to screen
          if (Comm().MyPID()==0)
            std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << std::endl;
        }
      }
    }
  } // if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)


  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = Teuchos::null;
    zigzagtwo_  = Teuchos::null;
    zigzagthree_= Teuchos::null;
  }

  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    std::cout << "ACTIVE SET ITERATION " << ActiveSetSteps()-1
         << " NOT CONVERGED - REPEAT TIME STEP................." << std::endl;
  else if (Comm().MyPID()==0 && activesetconv_==true)
    std::cout << "ACTIVE SET CONVERGED IN " << ActiveSetSteps()-zigzagging
         << " STEP(S)................." << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_=true;
    wasincontact_=true;
  }
  else
    isincontact_=false;

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::UpdateActiveSetSemiSmooth()
{
  // FIXME: Here we do not consider zig-zagging yet!
//  PrintActiveSet();

  // get out gof here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      // loop over all slave nodes on the current interface
      for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // The nested active set strategy cannot deal with the case of
        // active nodes that have no integration segments/cells attached,
        // as this leads to zero rows in D and M and thus to singular systems.
        // However, this case might possibly happen when slave nodes slide
        // over the edge of a master body within one fixed active set step.
        // (Remark: Semi-smooth Newton has no problems in this case, as it
        // updates the active set after EACH Newton step, see below, and thus
        // would always set the corresponding nodes to INACTIVE.)
        if (cnode->Active() && !cnode->HasSegment() && !cnode->IsOnBoundorCE())
          dserror("ERROR: Active node %i without any segment/cell attached",cnode->Id());
      }
    }
    return;
  }

  // get input parameter ftype
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");

//  // read weighting factor cn
//  // (this is necessary in semi-smooth Newton case, as the search for the
//  // active set is now part of the Newton iteration. Thus, we do not know
//  // the active / inactive status in advance and we can have a state in
//  // which both the condition znormal = 0 and wgap = 0 are violated. Here
//  // we have to weigh the two violations via cn!
//  const double ct_input = Params().get<double>("SEMI_SMOOTH_CT");

  // this is the complementarity parameter we use for the decision.
  // it might be scaled with a mesh-size dependent factor
  double cn = 0.;
  double ct = 0.;

  // do we use a mesh-size scaling for cn and ct?
  bool adaptive_cn = DRT::INPUT::IntegralValue<int>(Params(),"MESH_ADAPTIVE_CN");
  bool adaptive_ct = DRT::INPUT::IntegralValue<int>(Params(),"MESH_ADAPTIVE_CT");

  // do we apply a nodal scaling for better matrix condition
  bool scale = DRT::INPUT::IntegralValue<int>(scontact_,"LM_NODAL_SCALE");

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

      CoNode* cnode = dynamic_cast<CoNode*>(node);

      cn = interface_[i]->GetCnRef()[interface_[i]->GetCnRef().Map().LID(cnode->Id())];
      if(friction_)
        ct = interface_[i]->GetCtRef()[interface_[i]->GetCtRef().Map().LID(cnode->Id())];

      // get scaling factor
      double scalefac=1.;
      if (scale==true && cnode->MoData().GetScale() != 0.)
        scalefac = cnode->MoData().GetScale();

      // compute weighted gap
      double wgap = cnode->CoData().Getg()/scalefac;

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];
        nzold += cnode->MoData().n()[k] * cnode->MoData().lmold()[k];
      }

      // calculate mesh-size scaled version of cn

      if (adaptive_cn || adaptive_ct)
      {
        if (cnode->MoData().GetD().size()!=0) // only do that if there is a D matrix
        {
          // row sum of D matrix
          double sumd = 0.;
          GEN::pairedvector<int,double>::const_iterator p;
          for (p=cnode->MoData().GetD().begin(); p!=cnode->MoData().GetD().end(); p++)
            sumd += p->second;
          double mesh_h = pow(sumd,1./((double)Dim()-1.));
          if (adaptive_cn && mesh_h != 0.)
            cn /= mesh_h;
          if (adaptive_ct && mesh_h != 0.)
            ct /= mesh_h;
        }
      }

      // friction
      std::vector<double> tz (Dim()-1,0);
      std::vector<double> tjump (Dim()-1,0);
      double euclidean = 0.0;

      if (friction_)
      {
        // static cast
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute tangential parts and of Lagrange multiplier and incremental jumps
        for (int i=0;i<Dim();++i)
        {
          tz[0] += frinode->CoData().txi()[i]*frinode->MoData().lm()[i];
          if(Dim()==3) tz[1] += frinode->CoData().teta()[i]*frinode->MoData().lm()[i];

          if (DRT::INPUT::IntegralValue<int>(Params(),"GP_SLIP_INCR")==false)
          {
            tjump[0] += frinode->CoData().txi()[i]*frinode->FriData().jump()[i];
            if(Dim()==3) tjump[1] += frinode->CoData().teta()[i]*frinode->FriData().jump()[i];
          }
        }

        if (DRT::INPUT::IntegralValue<int>(Params(),"GP_SLIP_INCR")==true)
        {
          tjump[0] = frinode->FriData().jump_var()[0];
          if(Dim()==3) tjump[1] = frinode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm |tz+ct.tjump|
        std::vector<double> sum (Dim()-1,0);
        sum[0] = tz[0]+ct*tjump[0];
        if (Dim()==3) sum[1] = tz[1]+ct*tjump[1];
        if (Dim()==2) euclidean = abs(sum[0]);
        if (Dim()==3) euclidean = sqrt(sum[0]*sum[0]+sum[1]*sum[1]);
       }

      // adhesion
      double adhbound = 0.0;
      if(DRT::INPUT::IntegralValue<INPAR::CONTACT::AdhesionType>(Params(),"ADHESION") ==
          INPAR::CONTACT::adhesion_bound)
        adhbound = interface_[i]->IParams().get<double>("ADHESION_BOUND");

      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for penetration and/or tensile contact forces
        if (nz - cn*wgap > 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap > 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = true;
          activesetconv_ = false;

          // friction
          if (friction_)
          {
            // nodes coming into contact
            dynamic_cast<FriNode*>(cnode)->FriData().Slip() = true;
          }
        }
      }

      // check nodes of active set ***************************************
      else
      {
        //adhesion modification
        nz+=adhbound;

        // check for tensile contact forces and/or penetration
        if (nz - cn*wgap <= 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = false;
          activesetconv_ = false;

          // friction
          if (friction_) dynamic_cast<FriNode*>(cnode)->FriData().Slip() = false;
        }

        // only do something for friction
        else
        {
          // friction tresca
          if (ftype == INPAR::CONTACT::friction_tresca)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);

            // CAREFUL: friction bound is now interface-local (popp 08/2012)
            double frbound = interface_[i]->IParams().get<double>("FRBOUND");

            if(frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if(euclidean-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 frinode->FriData().Slip() = true;
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
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (fytpe=="tresca")

          // friction coulomb
          if (ftype == INPAR::CONTACT::friction_coulomb)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);

            // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
            double frcoeff = frinode->FrCoeff(interface_[i]->IParams().get<double>("FRCOEFF"));

            if(frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if(euclidean-frcoeff*(nz-cn*wgap) <= 1e-10) {}
              // do nothing (stick was correct)
              else
              {
                 frinode->FriData().Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (euclidean)-frbound > 0
              if(euclidean-frcoeff*(nz-cn*wgap) > -1e-10) {}
              // do nothing (slip was correct)
              else
              {
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (ftype == INPAR::CONTACT::friction_coulomb)
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
    activesetconv_ = false;
    activesetsteps_ += 1;
  }

  // only if it's a full Newton step...
  if (iterls_ <= 0)
  {
    // store the previous active set
    if (gactivenodes_ != Teuchos::null)
      gOldActiveSlaveNodes_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));
    else
      gOldActiveSlaveNodes_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
  }

  // also update special flag for semi-smooth Newton convergence
  activesetssconv_ = activesetconv_;

  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=Teuchos::null) zigzagthree_  = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=Teuchos::null) zigzagtwo_    = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=Teuchos::null) zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);

    if(friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
      gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
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
  if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)
  {
    // frictionless contact
    if (ActiveSetSteps()>2)
    {
      if (zigzagtwo_!=Teuchos::null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 1;
        }
      }

      if (zigzagthree_!=Teuchos::null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 2;
        }
      }
    }
  } // if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)

  // output to screen
  if (Comm().MyPID()==0)
  {
    if (zigzagging==1)
    {
      std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << std::endl;
    }
    else if (zigzagging==2)
    {
      std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << std::endl;
    }
    else
    {
      // do nothing, no zig-zagging
    }
  }

  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = Teuchos::null;
    zigzagtwo_  = Teuchos::null;
    zigzagthree_= Teuchos::null;
  }

  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    std::cout << "ACTIVE CONTACT SET HAS CHANGED... CHANGE No. " << ActiveSetSteps()-1 << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_=true;
    wasincontact_=true;
  }
  else
    isincontact_=false;

  return;
}

/*----------------------------------------------------------------------*
 |  Check conservation laws                              hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::CheckConservationLaws(const Epetra_Vector& fs,
                                                        const Epetra_Vector& fm)
{
  // change sign (adapt sign convention from the augmented Lagrange formulation)
  Teuchos::RCP<Epetra_Vector> tmpFm = Teuchos::rcp(new Epetra_Vector(fm));
  tmpFm->Scale(-1.0);
  /*-------------------------------*
   | LINEAR MOMENTUM CONSERVATION  |
   *-------------------------------*/
  double ssum = 0.0;
  double msum = 0.0;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Linear Momentum Conservation      <<" << std::endl;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Standard terms (lm)               <<" << std::endl;
  for (int i=0;i<fs.GlobalLength();++i) ssum+=fs[i];
  std::cout << "SLAVE:   " << std::setw(14) << ssum<< std::endl;
  for (int i=0;i<fm.GlobalLength();++i) msum+=(*tmpFm)[i];
  std::cout << "MASTER:  " << std::setw(14) << msum << std::endl;
  std::cout << "Balance: " << std::setw(14) << ssum+msum << std::endl;

  // write linear momentum conservation in xxx.lmom
  FILE* MyFile = NULL;
  std::ostringstream filename;
  const std::string filebase = "xxx";
  filename << filebase <<".lmom";
  MyFile = fopen(filename.str().c_str(), "at+");

  // store data
  if (MyFile)
  {
    fprintf(MyFile, "%g\t", ssum);
    fprintf(MyFile, "%g\t", msum);
    fprintf(MyFile, "%g\n", ssum+msum);
    fclose(MyFile);
  }
  else
    dserror("ERROR: File could not be opened.");

  /*-------------------------------*
   | ANGULAR MOMENTUM CONSERVATION |
   *-------------------------------*/
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Angular Momentum Conservation     <<" << std::endl;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    std::cout << ">>----- Interface " << std::setw(2) << i;
    std::cout << " ---------------------<<" << std::endl;
    std::cout << ">>      Standard terms (lm)               <<" << std::endl;
    interface_[i]->EvalResultantMoment(fs,*tmpFm);
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | do additional matrix manipulations for regularization scaling     ager 07/15|
 *----------------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::DoRegularizationScaling(bool aset,
                                                          bool iset,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& invda,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& kan,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& kam,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& kai,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& kaa,
                                                          Teuchos::RCP<Epetra_Vector>& fa,
                                                          Teuchos::RCP<LINALG::SparseMatrix>& kteffnew,
                                                          Teuchos::RCP<Epetra_Vector>& feffnew)
{
  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<LINALG::SparseMatrix> kanmod_n;
  if (aset)
  {
    kanmod_n = LINALG::MLMultiply(*nmatrix_,false,*invda,true,false,false,true);
    kanmod_n = LINALG::MLMultiply(*kanmod_n,false,*kan,false,false,false,true);
  }

  // kam: multiply tmatrix with invda and kam
  Teuchos::RCP<LINALG::SparseMatrix> kammod_n;
  if (aset)
  {
    kammod_n = LINALG::MLMultiply(*nmatrix_,false,*invda,true,false,false,true);
    kammod_n = LINALG::MLMultiply(*kammod_n,false,*kam,false,false,false,true);
  }

  // kai: multiply tmatrix with invda and kai
  Teuchos::RCP<LINALG::SparseMatrix> kaimod_n;
  if (aset && iset)
  {
    kaimod_n = LINALG::MLMultiply(*nmatrix_,false,*invda,true,false,false,true);
    kaimod_n = LINALG::MLMultiply(*kaimod_n,false,*kai,false,false,false,true);
  }

  // kaa: multiply tmatrix with invda and kaa
  Teuchos::RCP<LINALG::SparseMatrix> kaamod_n;
  if (aset)
  {
    kaamod_n = LINALG::MLMultiply(*nmatrix_,false,*invda,true,false,false,true);
    kaamod_n = LINALG::MLMultiply(*kaamod_n,false,*kaa,false,false,false,true);
  }

  /**********************************************************************/
  /* (8) Build the final f blocks                                       */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE

  Teuchos::RCP<Epetra_Vector> famod_n;
  Teuchos::RCP<LINALG::SparseMatrix> ninvda;
  if (aset)
  {
    famod_n = Teuchos::rcp(new Epetra_Vector(*gactiven_));
    ninvda = LINALG::MLMultiply(*nmatrix_,false,*invda,true,false,false,true);
    ninvda->Multiply(false,*fa,*famod_n);
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  if (ParRedist())
  {
    if (aset) nderivmatrix_ = MORTAR::MatrixRowTransform(nderivmatrix_,pgsdofrowmap_);
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  if (aset) kteffnew->Add(*nderivmatrix_,false,1.0,1.0);

  if (aset) kteffnew->Add(*kanmod_n,false,-1.0,1.0);
  if (aset) kteffnew->Add(*kammod_n,false,-1.0,1.0);
  if (aset && iset) kteffnew->Add(*kaimod_n,false,-1.0,1.0);
  if (aset) kteffnew->Add(*kaamod_n,false,-1.0,1.0);

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/
  //---------------------------------------------------------- FOURTH LINE
  Teuchos::RCP<Epetra_Vector> famodexp_n;
  if (aset)
  {
    famodexp_n = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*famod_n,*famodexp_n);
    feffnew->Update(-1.0,*famodexp_n,1.0);
  }
}

/*----------------------------------------------------------------------*
 | calculate regularization scaling and apply it to matrixes   ager 06/15|
 *----------------------------------------------------------------------*/
void CONTACT::CoLagrangeStrategy::EvaluateRegularizationScaling(Teuchos::RCP<Epetra_Vector> gact)
{
  /********************************************************************/
  /* (0) Get Nodal Gap Scaling                                        */
  /********************************************************************/
  //vector with all values alpha
  Teuchos::RCP<Epetra_Vector> scalevec = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));

  //vecor with (alpha + dalpha/dgap*scaling + d(1-alpha)/dgap)
  Teuchos::RCP<Epetra_Vector> derivscalevec = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));

  //vector with all values (1-alpha)
  Teuchos::RCP<Epetra_Vector> scalevec2 = Teuchos::rcp(new Epetra_Vector(*gactiven_, true));
  scalevec->PutScalar(1.0);

  {
    //setting of the parameters is still under investigation Ager Chr.
   // DRT::Problem* problem = DRT::Problem::Instance();
   // const Teuchos::ParameterList& structdyn   = problem->StructuralDynamicParams(); //just for now!

    double scaling = 1e6;//structdyn.get<double>("TOLRES"); //1e6
    //
    //sign of the term lambda*n (std is +1.0)?
    static double sign = 1.0;
    //scaling also with (1-alpha) of lambda*n term?
    static double scal1ma = 1.0;

    double epsilon = 1e-6;//structdyn.get<double>("TOLINCO");//1e-5;+

      // loop over all interfaces
      for (int i=0; i<(int)interface_.size(); ++i)
      {
        // loop over all slave nodes on the current interface
        for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
        {
          int gid = interface_[i]->SlaveRowNodes()->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);

          CoNode* cnode = dynamic_cast<CoNode*>(node);

          if (cnode->Active())
          {
            //normal lagrange multiplier! ... replace by multiplication
            double lmval = cnode->MoData().n()[0]*cnode->MoData().lm()[0] + cnode->MoData().n()[1]*cnode->MoData().lm()[1] + cnode->MoData().n()[2]*cnode->MoData().lm()[2];
            double alpha;
            double alphaderiv;

            std::cout << "cnode->CoData().Getg(): " << cnode->CoData().Getg() << "( eps = "<< epsilon << ", scal = " << scaling << " ) / lambdan: " << lmval << std::endl;


            if (cnode->CoData().Getg() >= 0)
            {
              alpha = 0.0;
              alphaderiv= 0.0;
            }
            else
            {
              double scaledgsum = (5*cnode->CoData().Getg()/epsilon+0.5);
              alpha = 0.5*(1-tanh(scaledgsum));
              alphaderiv= -0.5*5/(cosh(scaledgsum)*cosh(scaledgsum)*epsilon);

              alphaderiv= alphaderiv*(scaling*cnode->CoData().Getg()-sign*lmval);
            }

            int lid = gactiven_->LID(cnode->Dofs()[0]);
            if (lid != -1)
            {
              (*scalevec)[lid] = alpha;
              (*derivscalevec)[lid] = alphaderiv;
            }
            else
            {
              dserror("ERROR: lid = -1!");
            }

            std::cout << std::setprecision(32) <<"cnode->CoPoroData().GetnScale() for node " << cnode->Id() << " is: " << alpha << std::endl;
          }
        } // loop over all slave nodes
      } // loop over all interfaces

    //Evaluate final scaling vetors!
    {
      scalevec2->PutScalar(sign);
      scalevec2->Update(-sign*scal1ma,*scalevec,1.0); // 1-alpha

      scalevec->Scale(scaling);
      derivscalevec->Update(1.0,*scalevec,1.0);
    }
  }

  //scale all matrixes!!!
  int err = smatrix_->LeftScale(*derivscalevec); //scale with (alpha + dalpha/dgap*scaling + d(1-alpha)/dgap*lambda*n
  if (err) dserror("LeftScale failed!");

  err = nderivmatrix_->LeftScale(*scalevec2); //scale with 1.0-alpha
  if (err) dserror("LeftScale failed!");

  err = nmatrix_->LeftScale(*scalevec2); //scale with 1.0-alpha
  if (err) dserror("LeftScale failed!");

  Teuchos::RCP<Epetra_Vector> tmpgact = Teuchos::rcp(new Epetra_Vector(*gact)); //check later if this is required!!!
  gact->Multiply(1.0,*scalevec,*tmpgact,0.0);

  return;
}
