/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_interface.cpp

\brief Augmented contact interface.

\level 2

\maintainer Michael Hiermeier

\date Apr 16, 2014

*/
/*---------------------------------------------------------------------*/
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
 *----------------------------------------------------------------------*/
CONTACT::AugmentedInterface::AugmentedInterface(const int id,
    const Epetra_Comm& comm,
    const int dim,
    const Teuchos::ParameterList& icontact,
    bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : CONTACT::CoInterface(id,comm,dim,icontact,selfcontact,redundant),
      penBound_(-1.0),
      sndofrowmap_(Teuchos::null),
      stdofrowmap_(Teuchos::null)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FillComplete(
    int maxdof,
    bool newghosting)
{
  // call the standard mortar routine first
  MortarInterface::FillComplete(maxdof,newghosting);

  /* check if a user defined penetration bound is provided, if not we
   * use an automatic routine. */
  if (penBound_ < 0.0)
  {
    // find smallest interface element edge length
    double myMinEdgeLength = 1.0e12;
    // loop over all slave elements
    for (int i=0; i<selerowmap_->NumMyElements(); ++i)
    {
      int gid = selerowmap_->GID(i);
      DRT::Element* ele = idiscret_->gElement(gid);
      if (!ele) dserror("ERROR: Cannot find slave element with gid %i",gid);

      MORTAR::MortarElement* sele = dynamic_cast<MORTAR::MortarElement*>(ele);
      if (myMinEdgeLength > sele->MinEdgeSize())
        myMinEdgeLength = sele->MinEdgeSize();
    }

    // loop over all master elements
    for (int i=0; i<melerowmap_->NumMyElements(); ++i)
    {
      int gid = melerowmap_->GID(i);
      DRT::Element* ele = idiscret_->gElement(gid);
      if (!ele) dserror("ERROR: Cannot find master element with gid %i",gid);

      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);
      if (myMinEdgeLength > mele->MinEdgeSize())
        myMinEdgeLength = mele->MinEdgeSize();
    }
    // communicate the minimal edge length over all procs
    double gMinEdgeLength = 1.0e12;
    Comm().MinAll(&myMinEdgeLength,&gMinEdgeLength,1);

    if (gMinEdgeLength == 1.0e12 or gMinEdgeLength < 0.0)
      dserror("ERROR: Global minimal interface edge length couldn't"
          " be calculated! (value: %d)",gMinEdgeLength);

    // set penetration bound to 50% of the minimal interface element edge length,
    // if no user defined value is provided.
    penBound_ = 0.5*gMinEdgeLength;
  }

  return;
}

/*----------------------------------------------------------------------*
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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::UpdateMasterSlaveSets()
{
  MORTAR::MortarInterface::UpdateMasterSlaveSets();
  SplitSlaveDofs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::RedEvaluate(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
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
    // create a Augmented Lagrangian integrator instance with correct NumGP and Dim
    CONTACT::AugmentedIntegrator augIntegrator(IParams(),selement->Shape(),Comm());
    if (Dim()==2)
      augIntegrator.IntegrateDerivSlEle2D((*selement),Comm(),mparams_ptr);
    else if (Dim()==3)
      augIntegrator.IntegrateDerivSlEle3D((*selement),Comm(),mparams_ptr);
    else
      dserror("ERROR: RedEvaluate: Dim value has to be 2 or 3!");
    /**************************************************************************
     *                       !!! end integration !!!                          *
     **************************************************************************/
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugDnMnMatrix(
    LINALG::SparseMatrix& augDnMatrix,
    LINALG::SparseMatrix& augMnMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);

    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleAugDnMn: Node ownership inconsistency!");

    // We use the wrong map here and replace it later! This has the advantage, that
    // no extra sndofrowmap_ is needed.
    int rowId = gid;

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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDGLmLinMatrix(
    LINALG::SparseMatrix& dGLmSlLinMatrix,
    LINALG::SparseMatrix& dGLmMaLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid  = activenodes_->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDGGLinMatrix(
    LINALG::SparseMatrix& dGGSlLinMatrix,
    LINALG::SparseMatrix& dGGMaLinMatrix,
    const Epetra_Vector& cnVec,
    const bool& completeLin)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // get cn
//  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;
    typedef std::map<int,std::pair<int,double> >::const_iterator CII;
    typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CIII;

    // get nodal values
    double cn     = cnVec[cnVec.Map().LID(gid)];
    double aWGap  = cnode->CoData().GetWGap();
    double ckappa = cnode->CoData().GetKappa();
    double ckappainv = 1/ckappa;
    if (aWGap!=1.0e12 and ckappa!=1.0e12)
      aWGap *= ckappainv;
    else
      dserror("ERROR: GID % seems to be no active slave node!",gid);
    std::map<int,double>& aWGapLinMap = cnode->CoData().GetAWGapLin();

    GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->CoData().GetVarWGapSl();
    std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->CoData().GetVarWGapMa();
    /* --------------------------- SLAVE SIDE ----------------------------------- */
    // iteration over ALL slave Dof Ids
    for (CIII p=varWGapSlMap.begin();p!=varWGapSlMap.end();++p)
    {
      int sRow = p->first;
      int col = 0;
      // *** asymptotic convergence phase: consistent linearization ***
      // This term is only considered, when the current node shows no large penetration.
      if (completeLin)
      {
        // *** linearization of varWGap w.r.t. displacements ***
        std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[sRow];
        for (CI pp=varWGapLinSlMap.begin();pp!=varWGapLinSlMap.end();++pp)
        {
          double val = cn*(pp->second)*aWGap;
          col = pp->first;
          if (abs(val)>1.0e-12) dGGSlLinMatrix.FEAssemble(val,sRow,col);
        }
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
    // iteration over ALL master Dof Ids
    for (CII p=varWGapMaMap.begin();p!=varWGapMaMap.end();++p)
    {
      int mRow = p->first;
      int col = 0;
      // *** asymptotic convergence phase: consistent linearization ***
      // This term is only considered, when the current node shows no large penetration.
      if (completeLin)
      {
        // *** linearization of varWGap w.r.t. displacements ***
        std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[mRow];
        for (CI pp=varWGapLinMaMap.begin();pp!=varWGapLinMaMap.end();++pp)
        {
          double val = -cn*(pp->second)*aWGap;
          col = pp->first;
          if (abs(val)>1.0e-12) dGGMaLinMatrix.FEAssemble(val,mRow,col);
        }
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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleResidualVectors(
    Epetra_Vector& lmNVec,
    Epetra_Vector& aWGapVec,
    Epetra_Vector& wGapVec)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get cn
//  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGLmrhs: Cannot find slave"
        " node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    double lmn  = cnode->MoData().lm()[0];

    // calculate averaged weighted gap
    double wGap = cnode->CoData().GetWGap();
    double aWGap = wGap;
    double kappa = cnode->CoData().GetKappa();
    if (kappa!=1.0e12 && aWGap!=1.0e12)
    {
      aWGap /= kappa;

      Epetra_SerialDenseVector lmNV(1);
      Epetra_SerialDenseVector aWGapV(1);
      Epetra_SerialDenseVector wGapV(1);
      lmNV[0]   = lmn;
      aWGapV[0] = aWGap;
      wGapV[0]  = wGap;

      std::vector<int> rgid(1,activen_->GID(i));
      std::vector<int> rowner(1,cnode->Owner());

      LINALG::Assemble(lmNVec,lmNV,rgid,rowner);
      LINALG::Assemble(aWGapVec,aWGapV,rgid,rowner);
      LINALG::Assemble(wGapVec,wGapV,rgid,rowner);
    }
    else
      dserror("ERROR: Kappa and/or the weighted gap should "
          "not be equal 1.0e12 for active nodes!");
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmNWGapLinMatrix(
    LINALG::SparseMatrix& dLmNWGapLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all active augmented slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find node with gid %",gid);

    // typedef of the map iterator
    typedef std::map<int,double>::const_iterator CI;

    double wGap = cnode->CoData().GetWGap();

    if (wGap==1.0e12)
      dserror("ERROR: The weighted gap should not be equal 1.0e12 "
          "for active nodes!");

    std::map<int,double>& wGapLinMap = cnode->CoData().GetWGapLin();
    double val = 0.0;

    // linearization of the weighted gap
    for (CI p=wGapLinMap.begin();p!=wGapLinMap.end();++p)
    {
      int rowId = activen_->GID(i);
      int colId = p->first;
      val = p->second;

      if (abs(val)>1.0e-12)
        dLmNWGapLinMatrix.Assemble(val,rowId,colId);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveRhs(
    Epetra_Vector& augInactiveRhs,
    Epetra_Vector& cnVec)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleInactiverhs: "
          "Node ownership inconsistency!");

    double cn_inv   = 1/cnVec[cnVec.Map().LID(gid)];
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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTRhs(
    Epetra_Vector& dLmTLmTRhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for(int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTRhs: Cannot find active"
        " slave node with gid %",gid);

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
     |        0             ==>           0                             |
     |        1             ==>           1                             |
     |                       :                                          |
     |                       :                                          |
     |                                                                  |
     |* 3-D case *******************************************************|
     | (Dim()-1) = 2 and j=0,1                                          |
     |                                                                  |
     |        i                      (Dim()-1)*i+j                      |
     |==================================================================|
     |        0             ==>          0,1                            |
     |        1             ==>          2,3                            |
     |                       :                                          |
     |                       :                                          |
     *------------------------------------------------------------------*/
    for (int j=0;j<(Dim()-1);++j)
    {
      rGid[j] = activet_->GID((Dim()-1)*i+j);
      rVal[j] = ct_inv*lm[j+1]*augA;
    }

    // Assemble
    LINALG::Assemble(dLmTLmTRhs,rVal,rGid,rOwner);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTMatrix(
    LINALG::SparseMatrix& dLmTLmTMatrix)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for(int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTrhs: Cannot find slave"
        " node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // get weighted element area
    double augA = cnode->CoData().GetAugA();

    for (int j=0;j<(Dim()-1);++j)
    {
      int rowId = activet_->GID((Dim()-1)*i+j);
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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleDLmTLmTLinMatrix(
    LINALG::SparseMatrix& dLmTLmTLinMatrix)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
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
      int rowId = activet_->GID((Dim()-1)*i+j);

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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveMatrix(
    LINALG::SparseMatrix& augInactiveMatrix,
    const Epetra_Vector& cnVec)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find"
        " inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
    double augA   = cnode->CoData().GetAugA();

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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugInactiveLinMatrix(
    LINALG::SparseMatrix& augInactiveLinMatrix,
    const Epetra_Vector& cnVec)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find"
        " inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
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
    if (!node) dserror("ERROR: AssembleDGGrhs: Cannot find slave"
        " node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap =
        cnode->CoData().GetVarWGapSl();

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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AWGapLin(bool completeLin)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
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

        double val = 0.0;
        // Lin(kappa)
        if (completeLin)
        {
          val = (p->second).second*xs*kappainv*kappainv;
          for (CI pp=kappaLinMap.begin();pp!=kappaLinMap.end();++pp)
            aWGapLinMap[pp->first] += (pp->second)*val;
        }

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

        double val = 0.0;
        // Lin(kappa)
        if (completeLin)
        {
          val = (p->second).second*xm*kappainv*kappainv;
          for (CI pp=kappaLinMap.begin();pp!=kappaLinMap.end();++pp)
            aWGapLinMap[pp->first] -= (pp->second)*val;
        }

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
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugContactPotential(
    double& conPot,
    double& augConPot,
    Epetra_Vector& cnVec)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get cn
//  double cn     = IParams().get<double>("SEMI_SMOOTH_CN");
//  double cn_inv = 1.0/cn;
  double ct_inv = 1.0/IParams().get<double>("SEMI_SMOOTH_CT");

  // *** Active part *************************************************
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    double cn = cnVec[cnVec.Map().LID(gid)];

    double wgap = cnode->CoData().GetWGap();
    if (wgap==1.0e12) dserror("ERROR: WGap is equal 1.e12 for a active node! (node-id: %d)",gid);
    double kappa = cnode->CoData().Kappa();
    if (kappa==1.0e12) dserror("ERROR: Kappa is equal 1.e12 for a active node! (node-id: %d)",gid);
    double augA = cnode->CoData().GetAugA();

    // get the Lagrange multiplier
    double* lm  = cnode->MoData().lm();

    // *** ACTIVE - NORMAL DIRECTION ***
    if (wgap != 0.0)
    {
      // ** - (zn_i * awgap_i * A_i) **
      conPot -= wgap*lm[0];
      // ** - (zn_i * awgap_i * A_i - cn/2 * awgap * awgap * A_i) **
      augConPot -= wgap*lm[0];
      augConPot += 0.5*cn*wgap*wgap/kappa;
    }
    // *** ACTIVE - TANGENTIAL DIRECTION ***
    for (int d=1;d<Dim();++d)
    {
      // ** - 1/(ct) * zt_i^T*zt_i * A_i **
      conPot -= ct_inv*lm[d]*lm[d]*augA;
      // ** - 1/(2*ct) * zt_i^T*zt_i * A_i **
      augConPot -= 0.5*ct_inv*lm[d]*lm[d]*augA;
    }
  }

  // *** Inactive part *************************************************
  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  =
      LINALG::SplitMap(*snoderowmap_, *activenodes_);
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<augInactiveSlaveNodes->NumMyElements();++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
    double augA = cnode->CoData().GetAugA();

    // get the lagrange multiplier
    double* lm  = cnode->MoData().lm();

    // *** INACTIVE - NORMAL DIRECTION ***
    // ** - 1/(cn) * zn_i * zn_i * A_i **
    augConPot -= cn_inv*lm[0]*lm[0]*augA;
    // ** - 1/(2*cn) * zn_i * zn_i * A_i **
    augConPot -= 0.5*cn_inv*lm[0]*lm[0]*augA;
    // *** INACTIVE - TANGENTIAL DIRECTION ***
    // ** - 1/(2*ct) * zt_i^T*zt_i * A_i
    for (int d=1;d<Dim();++d)
    {
      conPot -= ct_inv*lm[d]*lm[d]*augA;
      augConPot -= 0.5*ct_inv*lm[d]*lm[d]*augA;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AugmentedInterface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mydofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find node with gid %i",gid);

    const int numdof = cnode->NumDof();

    // -------------------------------------------------------------------
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // -------------------------------------------------------------------
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // -------------------------------------------------------------------
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
        cnode->Active()=true;
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }
    }
    // -------------------------------------------------------------------
    // RE-BUILDING OF THE ACTIVE SET
    // -------------------------------------------------------------------
    else
    {
      // check if node is active
      if (cnode->Active())
      {
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }
    }
  } // end loop over all slave nodes

  // create interface local augmented active node map and augmented active dof map
  activenodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)mynodegids.size(),&mynodegids[0],0,Comm()));
  activedofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)mydofgids.size(),&mydofgids[0],0,Comm()));

  SplitAugActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::SplitAugActiveDofs()
{
  // get out of here if augmented active set is empty
  if (activenodes_==Teuchos::null or
      activenodes_->NumGlobalElements() == 0)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(activenodes_->NumMyElements());
  std::vector<int> myTGids((Dim()-1)*activenodes_->NumMyElements());

  // dimension check
  double dimcheck =(activedofs_->NumGlobalElements())/(activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitAugActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all augmented active row nodes
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
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
  if ((gCountN+gCountT)!=activedofs_->NumGlobalElements())
    dserror("ERROR: SplitAugActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = Teuchos::rcp(new Epetra_Map(gCountN,countN,&myNGids[0],0,Comm()));
  activet_ = Teuchos::rcp(new Epetra_Map(gCountT,countT,&myTGids[0],0,Comm()));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::SplitSlaveDofs()
{
  // get out of here if augmented active set is empty
  if (snoderowmap_==Teuchos::null or
      snoderowmap_->NumGlobalElements() == 0)
  {
    sndofrowmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    stdofrowmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(snoderowmap_->NumMyElements());
  std::vector<int> myTGids((Dim()-1)*snoderowmap_->NumMyElements());

  // dimension check
  double dimcheck =(sdofrowmap_->NumGlobalElements())/(snoderowmap_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitSlaveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all augmented active row nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
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
  if ((gCountN+gCountT)!=sdofrowmap_->NumGlobalElements())
    dserror("ERROR: SplitSlaveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  sndofrowmap_ = Teuchos::rcp(new Epetra_Map(gCountN,countN,&myNGids[0],0,Comm()));
  stdofrowmap_ = Teuchos::rcp(new Epetra_Map(gCountT,countT,&myTGids[0],0,Comm()));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::AssembleAugAVector(
    Epetra_Vector& augAVec,
    Epetra_Vector& kappaVec)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    std::vector<int> rgid(1,gid);
    std::vector<int> rowner(1,cnode->Owner());

    // *** augmented Area ***
    double augA  = cnode->CoData().GetAugA();

    if (augA > 0.0)
    {
      Epetra_SerialDenseVector augAV(1);
      augAV[0] = augA;

      LINALG::Assemble(augAVec,augAV,rgid,rowner);
    }
    else
      dserror("ERROR: The augmented nodal area shouldn't be equal/lower than zero! "
          "(value= %.2e)",augA);

    // *** kappa ***
    if (cnode->Active())
    {
      double kappa  = cnode->CoData().GetKappa();

      if (kappa > 0.0)
      {
        Epetra_SerialDenseVector kappaV(1);
        kappaV[0] = kappa;

        LINALG::Assemble(kappaVec,kappaV,rgid,rowner);
      }
      else
        dserror("ERROR: The weighted area kappa shouldn't be equal/lower than zero! "
            "(value= %.2e)",augA);
    }
  }

  return;
}
