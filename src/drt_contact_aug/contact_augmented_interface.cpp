/*!----------------------------------------------------------------------
\file contact_augmented_interface.cpp

<pre>
Created on: Apr 16, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#include "contact_augmented_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_contact/contact_coupling3d.H"
#include "contact_augmented_integrator.H"
#include "../drt_mortar/mortar_binarytree.H"
#include "../drt_lib/drt_discret.H"
#include <Teuchos_Time.hpp>
#include <Epetra_Time.h>
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | ctor                                                  hiermeier 04/14|
 *----------------------------------------------------------------------*/
CONTACT::AugmentedInterface::AugmentedInterface(const int id,
                                                const Epetra_Comm& comm,
                                                const int dim,
                                                const Teuchos::ParameterList& icontact,
                                                bool selfcontact,
                                                INPAR::MORTAR::RedundantStorage redundant) :
CONTACT::CoInterface(id,comm,dim,icontact,selfcontact,redundant)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  initialize / reset augmented interface for contact   hiermeier 07/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // call initialization routine of the contact_interface
  CoInterface::Initialize();

  // *** Reset all augmented quantities *********************************
  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i=0;i<SlaveColNodesBound()->NumMyElements();++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // reset nodal weighted gap
    cnode->CoData().GetWGap() = 1.0e12;

    // reset nodal scaling factor
    cnode->CoData().GetKappa() = 1.0e12;
    cnode->CoData().GetAugA()  = 0.0;

    // reset variation of the weighted gap
    cnode->CoData().GetVarWGapSl().clear();
    cnode->CoData().GetVarWGapMa().clear();

    // reset kappaLin
    cnode->CoData().GetKappaLin().clear();
    cnode->CoData().GetAugALin().clear();

    // reset varWGapLin
    cnode->CoData().GetVarWGapLinSl().clear();
    cnode->CoData().GetVarWGapLinMa().clear();

    // reset linearization of the averaged weighted gap
    cnode->CoData().GetAWGapLin().clear();

    // reset linearization of the averaged weighted gap
    cnode->CoData().GetWGapLin().clear();
  }

  return;
}

/*----------------------------------------------------------------------*
 | Reduced evaluate of the contact interface             hiermeier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::RedEvaluate()
{
// interface needs to be complete
  if (!Filled() && Comm().MyPID()==0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(ele1);

    if (selement->ZeroSized())
      continue;

    /**************************************************************************
     *    Integrate all remaining quantities only over the slave interface    *
     **************************************************************************/
    // create a CONTACT integrator instance with correct NumGP and Dim
    // Augmented Lagrange integrator
    CONTACT::AugmentedIntegrator augIntegrator(IParams(),selement->Shape(),Comm(),AugActiveSlaveNodes());
    if (Dim()==2)
      augIntegrator.IntegrateDerivSlEle2D((*selement),Comm());
    else if (Dim()==3)
      augIntegrator.IntegrateDerivSlEle3D((*selement),Comm());
    else
      dserror("ERROR: RedEvaluate: Dim value has to be 2 or 3!");
    /**************************************************************************
     *                       !!! end integration !!!                          *
     **************************************************************************/
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global Dn and Mn matrices                    hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugDnMnMatrix(LINALG::SparseMatrix& augDnMatrix,
                                                        LINALG::SparseMatrix& augMnMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);

    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleAugDnMn: Node ownership inconsistency!");

    int rowId = augActiveSlaveNDofs_->GID(i);

    // typedef of the map iterator
    typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CII;
    typedef std::map<int,std::pair<int,double> >::const_iterator CI;

    /************************************* augDn-Matrix ******/
    GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->CoData().GetVarWGapSl();

    for(CII pp=varWGapSlMap.begin();pp!=varWGapSlMap.end();++pp)
    {
      int cDofId  = pp->first;
      double cval = -(pp->second).second;
      augDnMatrix.Assemble(cval,rowId,cDofId);
    }
    /************************************* augMn-Matrix ******/
    std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->CoData().GetVarWGapMa();

    for(CI pp=varWGapMaMap.begin();pp!=varWGapMaMap.end();++pp)
    {
      int cDofId  = pp->first;
      double cval = (pp->second).second;
      augMnMatrix.Assemble(cval,rowId,cDofId);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global dGLmLinMatrix                         hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmSlLinMatrix,
                                                        LINALG::SparseMatrix& dGLmMaLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    // get Lagrange multiplier in normal direction
    double lm_n = cnode->MoData().lm()[0];

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;
    typedef std::map<int,std::map<int,double> >::const_iterator CII;

    /* --------------------------- SLAVE SIDE ----------------------------------- */
    std::map<int,std::map<int,double> >& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl();
    // iteration over ALL slave Dof Ids
    for (CII p=varWGapLinSlMap.begin();p!=varWGapLinSlMap.end();++p)
    {
      int sRow = p->first;
      int col = 0;
      // *** linearization of varWGap w.r.t. displacements ***
      for (CI pp=varWGapLinSlMap[sRow].begin();pp!=varWGapLinSlMap[sRow].end();++pp)
      {
        double val = - (pp->second) * lm_n;
        col = pp->first;
        if (abs(val)>1.0e-12) dGLmSlLinMatrix.FEAssemble(val,sRow,col);
      }
    }
    /* --------------------------- MASTER SIDE ---------------------------------- */
    std::map<int,std::map<int,double> >& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa();
    // iteration over ALL master Dof Ids
    for (CII p=varWGapLinMaMap.begin();p!=varWGapLinMaMap.end();++p)
    {
      int mRow = p->first;
      int col = 0;
      // *** linearization of varWGap w.r.t. displacements ***
      for (CI pp=varWGapLinMaMap[mRow].begin();pp!=varWGapLinMaMap[mRow].end();++pp)
      {
        double val = (pp->second) * lm_n;
        col = pp->first;
        if (abs(val)>1.0e-12) dGLmMaLinMatrix.FEAssemble(val,mRow,col);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global dGGLinMatrix                         hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDGGLinMatrix(LINALG::SparseMatrix& dGGSlLinMatrix,
                                                       LINALG::SparseMatrix& dGGMaLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // get cn
  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

      // typedef of the map iterator
      typedef std::map<int,double>::const_iterator CI;
      typedef std::map<int,std::pair<int,double> >::const_iterator CII;
      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CIII;

      double aWGap  = cnode->CoData().GetWGap();
      double ckappa = cnode->CoData().GetKappa();
      if (aWGap!=1.0e12 and ckappa!=1.0e12)
        aWGap /= ckappa;
      else
        dserror("ERROR: GID % seems to be no active slave node!",gid);

      std::map<int,double>& aWGapLinMap = cnode->CoData().GetAWGapLin();

      /* --------------------------- SLAVE SIDE ----------------------------------- */
      GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->CoData().GetVarWGapSl();
      // iteration over ALL slave Dof Ids
      for (CIII p=varWGapSlMap.begin();p!=varWGapSlMap.end();++p)
      {
        int sRow = p->first;
        int col = 0;
        // *** linearization of varWGap w.r.t. displacements ***
        std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[sRow];
        for (CI pp=varWGapLinSlMap.begin();pp!=varWGapLinSlMap.end();++pp)
        {
          double val = cn*(pp->second)*aWGap;
          col = pp->first;
          if (abs(val)>1.0e-12) dGGSlLinMatrix.FEAssemble(val,sRow,col);
        }
        // *** linearization of the averaged weighted gap w.r.t. displacements ***
        for (CI pp=aWGapLinMap.begin();pp!=aWGapLinMap.end();++pp)
        {
          double val = cn*(p->second).second*(pp->second);
          col = pp->first;
          if (abs(val)>1.0e-12) dGGSlLinMatrix.FEAssemble(val,sRow,col);
        }
      }
      /* --------------------------- MASTER SIDE ---------------------------------- */
      std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->CoData().GetVarWGapMa();
      // iteration over ALL slave Dof Ids
      for (CII p=varWGapMaMap.begin();p!=varWGapMaMap.end();++p)
      {
        int mRow = p->first;
        int col = 0;
        // *** linearization of varWGap w.r.t. displacements ***
        std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[mRow];
        for (CI pp=varWGapLinMaMap.begin();pp!=varWGapLinMaMap.end();++pp)
        {
          double val = -cn*(pp->second)*aWGap;
          col = pp->first;
          if (abs(val)>1.0e-12) dGGMaLinMatrix.FEAssemble(val,mRow,col);
        }
        // *** linearization of the averaged weighted gap w.r.t. displacements ***
        for (CI pp=aWGapLinMap.begin();pp!=aWGapLinMap.end();++pp)
        {
          double val = -cn*(p->second).second*(pp->second);
          col = pp->first;
          if (abs(val)>1.0e-12) dGGMaLinMatrix.FEAssemble(val,mRow,col);
        }
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global augmented Lagrange multiplier vector  hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugLmVector(Epetra_Vector& augLmVec)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get cn
  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGLmrhs: Cannot find slave node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    double lm_n  = cnode->MoData().lm()[0];

    // calculate averaged weighted gap
    double aWGap = cnode->CoData().GetWGap();
    double kappa = cnode->CoData().GetKappa();
    aWGap /= kappa;

    Epetra_SerialDenseVector augLm(1);
    augLm[0] = lm_n-cn*aWGap;
    std::vector<int> rgid(1,augActiveSlaveNDofs_->GID(i));
    std::vector<int> rowner(1,cnode->Owner());

    LINALG::Assemble(augLmVec,augLm,rgid,rowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global averaged weighted gap vector          hiermeier 06/14|
 | Only necessary for the check of the conservation laws                |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAWGapRhs(Epetra_Vector& aWGapRhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleAWGapN: Node ownership inconsistency!");

    if (cnode->CoData().GetWGap()!=0.0)
    {
      double aWGap = cnode->CoData().GetWGap();
      double kappa = cnode->CoData().GetKappa();

      if (kappa!=1.0e12 and aWGap!=1.0e12)
      {
        aWGap = aWGap/kappa;

        Epetra_SerialDenseVector aWGapVec(1);
        std::vector<int> rGid(1,augActiveSlaveNDofs_->GID(i));
        std::vector<int> rOwner(1,cnode->Owner());

        aWGapVec[0] = aWGap;

        // Assemble vector
        LINALG::Assemble(aWGapRhs,aWGapVec,rGid,rOwner);
      }
      else
        dserror("ERROR: Kappa and/or the weighted gap should not be equal 1.0e12 for active nodes!");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global active normal constraint rhs          hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmNWGapRhs(Epetra_Vector& dLmNWGapRhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleAWGapN: Node ownership inconsistency!");

    double wGap = cnode->CoData().GetWGap();

    if (wGap!=1.0e12)
    {
      Epetra_SerialDenseVector rVal(1);
      std::vector<int> rGid(1,augActiveSlaveNDofs_->GID(i));
      std::vector<int> rOwner(1,cnode->Owner());

      rVal[0] = wGap;

      LINALG::Assemble(dLmNWGapRhs,rVal,rGid,rOwner);
    }
    else
      dserror("ERROR: The weighted gap should not be equal 1.0e12 for active nodes!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global dLmNWGapLinMatrix                     hiermeier 06/14|
 | Linearization w.r.t. the displ.                                      |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmNWGapLinMatrix(LINALG::SparseMatrix& dLmNWGapLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all active augmented slave nodes of the interface
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find node with gid %",gid);

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;

    double wGap = cnode->CoData().GetWGap();

    if (wGap==1.0e12)
      dserror("ERROR: The weighted gap should not be equal 1.0e12 for active nodes!");

    std::map<int,double>& wGapLinMap = cnode->CoData().GetWGapLin();
    double val = 0.0;

    // linearization of the weighted gap
    for (CI p=wGapLinMap.begin();p!=wGapLinMap.end();++p)
    {
      int rowId = augActiveSlaveNDofs_->GID(i);
      int colId = p->first;
      val = p->second;

      if (abs(val)>1.0e-12)
        dLmNWGapLinMatrix.Assemble(val,rowId,colId);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble augmented inactive right hand side          hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveRhs(Epetra_Vector& augInactiveRhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(*snoderowmap_, *augActiveSlaveNodes_);

  // get cn and invert it
  double cn_inv = IParams().get<double>("SEMI_SMOOTH_CN");
  cn_inv = 1/cn_inv;
  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleInactiverhs: Node ownership inconsistency!");

    double* lm  = cnode->MoData().lm();
    double augA = cnode->CoData().GetAugA();

    std::vector<int> rGid(Dim());
    std::vector<int> rOwner(Dim(),cnode->Owner());
    Epetra_SerialDenseVector rVal(Dim());

    for (int rDof=0; rDof<cnode->NumDof();++rDof)
    {
      rGid[rDof] = cnode->Dofs()[rDof];
      // normal direction
      if (rDof==0)
        rVal[rDof] = cn_inv * lm[rDof] * augA;
      // tangential direction
      else
        rVal[rDof] = ct_inv * lm[rDof] * augA;
    }

    LINALG::Assemble(augInactiveRhs,rVal,rGid,rOwner);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global dLmTLmT right hand side               hiermeier 06/14|
 | FRICTIONLESS TANGENTIAL CONSTRAINT                                   |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for(int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTRhs: Cannot find active slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // Get the Lagrange multiplier and txi of the current node
    double* lm  = cnode->MoData().lm();
    // Get weighted element area
    double augA  = cnode->CoData().GetAugA();

    std::vector<int> rGid(Dim()-1);
    std::vector<int> rOwner(Dim()-1,cnode->Owner());
    Epetra_SerialDenseVector rVal(Dim()-1);

    /*------------------------------------------------------------------*
     |* 2-D case *******************************************************|
     | (Dim()-1) = 1 and j = 0                                          |
     |                                                                  |
     |        i                      (Dim()-1)*i+j                      |
     |==================================================================|
     |        0             -->           0                             |
     |        1             -->           1                             |
     |                       :                                          |
     |                       :                                          |
     |                                                                  |
     |* 3-D case *******************************************************|
     | (Dim()-1) = 2 and j=0,1                                          |
     |                                                                  |
     |        i                      (Dim()-1)*i+j                      |
     |==================================================================|
     |        0             -->          0,1                            |
     |        1             -->          2,3                            |
     |                       :                                          |
     |                       :                                          |
     *------------------------------------------------------------------*/
    for (int j=0;j<(Dim()-1);++j)
    {
      rGid[j] = augActiveSlaveTDofs_->GID((Dim()-1)*i+j);
      rVal[j] = ct_inv*lm[j+1]*augA;
    }

    // Assemble
    LINALG::Assemble(dLmTLmTRhs,rVal,rGid,rOwner);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global DLmTLmTMatrix                         hiermeier 06/14|
 | Linearization matrix w.r.t. the LM                                   |
 | FRICTIONLESS TANGENTIAL CONSTRAINTS (active slave nodes)             |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTMatrix(LINALG::SparseMatrix& dLmTLmTMatrix)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for(int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTrhs: Cannot find slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // get weighted element area
    double augA = cnode->CoData().GetAugA();

    for (int j=0;j<(Dim()-1);++j)
    {
      int rowId = augActiveSlaveTDofs_->GID((Dim()-1)*i+j);
      int colId = rowId;
      /*-----------------------------------------------------------*
       | Attention:                                                |
       | D_{lm}[1/ct*(dlm_t * lm_t)] = D_{lm}[1/ct*[dlm_t*lm_t]]   |
       *-----------------------------------------------------------*/
      double val = ct_inv*augA;

      // Assemble
      if (abs(val)>1.0e-12) dLmTLmTMatrix.Assemble(val,rowId,colId);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global DLmTLmTLinMatrix                      hiermeier 05/14|
 | Linearization matrix w.r.t. the displ.                               |
 | FRICTIONLESS TANGENTIAL CONSTRAINTS (all slave nodes)                |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTLinMatrix(LINALG::SparseMatrix& dLmTLmTLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDGLmLin: Node ownership inconsistency!");

    double* lm  = cnode->MoData().lm();
    std::map<int,double>& augALinMap = cnode->CoData().GetAugALin();

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;

    /*---------------------------------------------------------*
     | Attention:                                              |
     | D_{d}[-1/ct*(dlm_t * lm_t)] = D_{d}[-1/ct*[dlm_t*lm_t]] |
     *---------------------------------------------------------*/
    for (int j=0;j<(Dim()-1);++j)
    {
      int rowId = augActiveSlaveTDofs_->GID((Dim()-1)*i+j);

      double tmp = ct_inv*lm[j+1];
      for (CI p=augALinMap.begin();p!=augALinMap.end();++p)
      {
        int colId = p->first;
        double val = tmp*(p->second);

        if (abs(val)>1.0e-12) dLmTLmTLinMatrix.Assemble(val,rowId,colId);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global AugInactiveMatrix                     hiermeier 06/14|
 | Linearization w.r.t. the LM                                          |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveMatrix(LINALG::SparseMatrix& augInactiveMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(*snoderowmap_, *augActiveSlaveNodes_);

  // get cn and invert it
  double cn_inv = IParams().get<double>("SEMI_SMOOTH_CN");
  cn_inv = 1/cn_inv;
  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix: Node ownership inconsistency!");

    double augA = cnode->CoData().GetAugA();

    for (int j=0; j<cnode->NumDof();++j)
    {
      int rowId = cnode->Dofs()[j];
      int colId = rowId;
      double val = 0.0;

      // normal part
      if (j==0)
        val = cn_inv*augA;
      // tangential part
      else
        val = ct_inv*augA;

      if (abs(val)>1.0e-12) augInactiveMatrix.Assemble(val,rowId,colId);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble global AugInactiveLinMatrix                  hiermeier 06/14|
 | Linearization w.r.t. the displ.                                      |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveLinMatrix(LINALG::SparseMatrix& augInactiveLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(*snoderowmap_, *augActiveSlaveNodes_);

  // get cn and invert it
  double cn_inv = IParams().get<double>("SEMI_SMOOTH_CN");
  cn_inv = 1/cn_inv;
  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix: Node ownership inconsistency!");

    double* lm  = cnode->MoData().lm();
    std::map<int,double>& augALinMap = cnode->CoData().GetAugALin();

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;

    for (int j=0; j<Dim();++j)
    {
      int rowId = cnode->Dofs()[j];
      double tmp = 0.0;
      // normal direction
      if (j==0)
        tmp = cn_inv*lm[j];
      // tangential direction
      else
        tmp = ct_inv*lm[j];

      for (CI p=augALinMap.begin();p!=augALinMap.end();++p)
      {
        int colId = p->first;
        double val = tmp*(p->second);

        if (abs(val)>1.0e-12) augInactiveLinMatrix.Assemble(val,rowId,colId);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate the nodal weighted gap                      hiermeier 07/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::WGap()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid  =snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGGrhs: Cannot find slave node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleGGrhs: Node ownership inconsistency!");

    GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->CoData().GetVarWGapSl();

    if (varWGapSlMap.size()!=0)
    {
      typedef std::map<int,std::pair<int,double> >::const_iterator CI;
      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CII;
      // *** Slave fraction ****************************************************
      for (CII p=varWGapSlMap.begin();p!=varWGapSlMap.end();++p)
      {
        int sGid = (p->second).first;

        CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(sGid));
        if (!snode) dserror("ERROR: Cannot find slave node with gid %",sGid);
        // Check if Slave node
        if (!snode->IsSlave()) dserror("ERROR: This has to be a slave node!");

        // Switch from the global DofId to the local DofId
        int sDof  = ((p->first)-snode->Dofs()[0])%Dim();
        // get the spatial nodal coordinate
        double xs = snode->xspatial()[sDof];

        double val = -(p->second).second*xs;
        cnode->AddWGapValue(val);
      }
      // *** Master fraction ***************************************************
      std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->CoData().GetVarWGapMa();
      if (varWGapMaMap.size()==0) dserror("ERROR: AugmentedInterface::WGap(): "
          "Variation of the slave side is not empty, but the master side is empty!");
      for (CI p=varWGapMaMap.begin();p!=varWGapMaMap.end();++p)
      {
        int mGid = (p->second).first;

        CoNode* mnode = dynamic_cast<CoNode*>(idiscret_->gNode(mGid));
        if (!mnode) dserror("ERROR: Cannot find master node with gid %",mGid);
        // Check if Slave node
        if (mnode->IsSlave()) dserror("ERROR: This has to be a master node!");

        // Switch from the global DofId to the local DofId
        int mDof  = ((p->first)-mnode->Dofs()[0])%Dim();
        // get the spatial nodal coordinate
        double xm = mnode->xspatial()[mDof];

        double val = (p->second).second*xm;
        cnode->AddWGapValue(val);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate the nodal averaged weighted gap             hiermeier 04/14|
 | linearization                                                        |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AWGapLin()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGGrhs: Cannot find slave node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleGGrhs: Node ownership inconsistency!");

    // get some pre-calculated results for the current slave node:
    double kappainv = 1.0/cnode->CoData().GetKappa();
    if (cnode->CoData().GetKappa()!=1.0e12)
    {
      std::map<int,double>& kappaLinMap  = cnode->CoData().GetKappaLin();

      // Get std::map to save the linearization of the weighted gap
      std::map<int,double>& wGapLinMap = cnode->CoData().GetWGapLin();
      // Get std::map to save the linearization of the averaged weighted gap
      std::map<int,double>& aWGapLinMap = cnode->CoData().GetAWGapLin();

      // typedef of the constant map iterators
      typedef std::map<int,double> ::const_iterator CI;
      typedef std::map<int,std::pair<int,double> >::const_iterator CII;
      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CIII;
      // *** Slave side **************************************************
      GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->CoData().GetVarWGapSl();
      // loop over ALL Dofs of the slave side
      for (CIII p=varWGapSlMap.begin();p!=varWGapSlMap.end();++p)
      {
        int sGid = (p->second).first;

        CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(sGid));
        if (!snode) dserror("ERROR: Cannot find slave node with gid %",sGid);
        // Check if Slave node
        if (!snode->IsSlave()) dserror("ERROR: This has to be a slave node!");

        // Switch from the global DofId to the local DofId
        int sDof  = ((p->first)-snode->Dofs()[0])%Dim();

        double xs = snode->xspatial()[sDof];

        // Lin(kappa)
        double val = (p->second).second*xs*kappainv*kappainv;
        for (CI pp=kappaLinMap.begin();pp!=kappaLinMap.end();++pp)
          aWGapLinMap[pp->first] += (pp->second)*val;

        // Lin(varWGapLinSl)
        std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[p->first];
        for (CI pp=varWGapLinSlMap.begin();pp!=varWGapLinSlMap.end();++pp)
        {
          val = (pp->second)*xs;
          wGapLinMap[pp->first]  -= val;
          aWGapLinMap[pp->first] -= val*kappainv;
        }

        // Lin(xs)
        wGapLinMap[p->first]  -= (p->second).second;
        aWGapLinMap[p->first] -= (p->second).second*kappainv;
      }

      // *** Master side *************************************************
      std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->CoData().GetVarWGapMa();
      // loop over ALL Dofs of the master side
      for (CII p=varWGapMaMap.begin();p!=varWGapMaMap.end();++p)
      {
        int mGid = (p->second).first;

        CoNode* mnode = dynamic_cast<CoNode*>(idiscret_->gNode(mGid));
        if (!mnode) dserror("ERROR: Cannot find master node with gid %",mGid);
        // Check if master node
        if (mnode->IsSlave()) dserror("ERROR: This has to be a master node!");

        // Switch from the global DofId to the local DofId
        int mDof = ((p->first)-mnode->Dofs()[0])%Dim();

        double xm = mnode->xspatial()[mDof];

        // Lin(kappa)
        double val = (p->second).second*xm*kappainv*kappainv;
        for (CI pp=kappaLinMap.begin();pp!=kappaLinMap.end();++pp)
          aWGapLinMap[pp->first] -= (pp->second)*val;

        // Lin(varWGapLinMa)
        std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[p->first];
        for (CI pp=varWGapLinMaMap.begin();pp!=varWGapLinMaMap.end();++pp)
        {
          val = (pp->second)*xm;
          wGapLinMap[pp->first]  += val;
          aWGapLinMap[pp->first] += val*kappainv;
        }

        // Lin(xm)
        wGapLinMap[p->first]  += (p->second).second;
        aWGapLinMap[p->first] += (p->second).second*kappainv;
      }
    }
    else
      dserror("ERROR: CONTACT::AugmentedInterface::AWGapLin(): kappa shouldn't be "
          "equal 1.0e12 for active nodes!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Build augmented active set (nodes / dofs)            hiermeier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::BuildAugActiveSet(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mydofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERRO: BuildAugActiveSet: Cannot find node with gid %",gid);

    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // flag for initialization of init active nodes with nodal gaps
      bool initcontactbygap = DRT::INPUT::IntegralValue<int>(IParams(),"INITCONTACTBYGAP");
      // value
      double initcontactval = IParams().get<double>("INITCONTACTGAPVALUE");

      // Either init contact by definition or by gap
      if(cnode->IsInitActive() and initcontactbygap)
        dserror("Init contact either by definition in condition or by gap!");

      // check if node is initially active or, if initialization with nodal, gap,
      // the gap is smaller than the prescribed value
      if (cnode->IsInitActive() or (initcontactbygap and cnode->CoData().GetWGap() < initcontactval))
      {
        cnode->AugActive()=true;
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }
    }
    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->AugActive())
      {
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }
    }
  } // end loop over all slave nodes

  // create interface local augmented active node map and augmented active dof map
  augActiveSlaveNodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)mynodegids.size(),&mynodegids[0],0,Comm()));
  augActiveSlaveDofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)mydofgids.size(),&mydofgids[0],0,Comm()));

  SplitAugActiveDofs();

  return;
}

/*----------------------------------------------------------------------*
 |  Split augmented active dof set in normal and         hiermeier 05/14|
 |  tangential direction                                                |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::SplitAugActiveDofs()
{
  // get out of here if augmented active set is empty
  if (augActiveSlaveNodes_==Teuchos::null or
      augActiveSlaveNodes_->NumGlobalElements() == 0)
  {
    augActiveSlaveNDofs_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    augActiveSlaveTDofs_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(augActiveSlaveNodes_->NumMyElements());
  std::vector<int> myTGids((Dim()-1)*augActiveSlaveNodes_->NumMyElements());

  // dimension check
  double dimcheck =(augActiveSlaveDofs_->NumGlobalElements())/(augActiveSlaveNodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitAugActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all augmented active row nodes
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    // add first dof to nMap
    myNGids[countN] = cnode->Dofs()[0];
    ++countN;

    //add reamining dofs to tMap
    for (int j=1;j<cnode->NumDof();++j)
    {
      myTGids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNGids.resize(countN);
  myTGids.resize(countT);

  // communicate countN and countT among procs
  int gCountN, gCountT;
  Comm().SumAll(&countN,&gCountN,1);
  Comm().SumAll(&countT,&gCountT,1);

  // check global dimensions
  if ((gCountN+gCountT)!=augActiveSlaveDofs_->NumGlobalElements())
    dserror("ERROR: SplitAugActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  augActiveSlaveNDofs_ = Teuchos::rcp(new Epetra_Map(gCountN,countN,&myNGids[0],0,Comm()));
  augActiveSlaveTDofs_ = Teuchos::rcp(new Epetra_Map(gCountT,countT,&myTGids[0],0,Comm()));

  return;
}

/*----------------------------------------------------------------------*
 | Update the augmented active set                      hiermeier 04/14 |
 *----------------------------------------------------------------------*/
bool CONTACT::AugmentedInterface::UpdateAugActiveSetSemiSmooth()
{
  // In the augmented Lagrange formulation we do the active set decision at
  // a different time compared to the standard lagrange formulation. Due to this
  // and because of the slightly different AVERAGED weighted gap definition we
  // need a new active set, called augmented active set.

  // get out of here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = DRT::INPUT::IntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<snoderowmap_->NumMyElements();++j)
    {
      int gid = snoderowmap_->GID(j);
      DRT::Node* node = idiscret_->gNode(gid);
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
      if (cnode->Active() && !cnode->HasSegment())
        dserror("ERROR: Active node %i without any segment/cell attached",cnode->Id());
    }
    return (true);
  }

  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weight the two violations via cn!
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // assume that active set has converged and check for opposite
  bool activesetconv=true;

  // loop over all slave nodes of the current interface
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AugmentedInterface::UpdateAugActiveSetSemiSmooth: "
        "Cannot find node with gid %",gid);

    // compute averaged weighted gap
    double kappa = cnode->CoData().GetKappa();
    double awgap = cnode->CoData().GetWGap();
    if (kappa != 1.0e12)
      awgap/= kappa;

    // get normal part of the Lagrange multiplier
    double nz = cnode->MoData().lm()[0];

    // check nodes of inactive set *************************************
    if (cnode->AugActive()==false)
    {
      // check for fulfillment of contact condition
      if (nz - cn*awgap > 0.0)
      {
        cnode->AugActive() = true;
        activesetconv = false;
      }
    }
    // check nodes of active set ***************************************
    else
    {
      if (nz-cn*awgap<=0.0)
      {
        cnode->AugActive() = false;
        activesetconv = false;
      }
    }
  }

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv;
  Comm().SumAll(&localcheck,&convcheck,1);

  // active set is only converged, if it is converged on all procs
  // if no, increase number of active set steps
  if (convcheck != Comm().NumProc())
    activesetconv = false;


  // update active set Epetra_Maps
  BuildAugActiveSet();

  return (activesetconv);
}


