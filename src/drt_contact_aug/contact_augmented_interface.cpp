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
#include "../linalg/linalg_serialdensevector.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AUG::Interface::Interface(
    int id,
    const Epetra_Comm& comm,
    int dim,
    const Teuchos::ParameterList& icontact,
    bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : CONTACT::CoInterface(id,comm,dim,icontact,selfcontact,redundant),
      penBound_( -1.0 ),
      ct_( IParams().get<double>("SEMI_SMOOTH_CT") ),
      maxNumMasterElements_( 0 ),
      issetup_( false ),
      sndofrowmap_( Teuchos::null ),
      stdofrowmap_( Teuchos::null )
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FillComplete(
    int maxdof,
    bool newghosting)
{
  // call the standard mortar routine first
  MortarInterface::FillComplete(maxdof,newghosting);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Setup()
{
  if ( issetup_ )
    return;

  // find smallest interface element edge length
  double myMinEdgeLength = 1.0e12;

  // find largest slave element area
  double myMaxAreaSl = 0.0;

  // find smallest master element area
  double myMinAreaMa = 1.0e12;

  int myTriangleOnMaster = 0;

  // ------------------------------------------------------------------------
  // loop over all slave elements
  for (int i=0; i<selerowmap_->NumMyElements(); ++i)
  {
    int gid = selerowmap_->GID(i);
    DRT::Element* ele = idiscret_->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find slave element with gid %i",gid);

    MORTAR::MortarElement* sele = static_cast<MORTAR::MortarElement*>(ele);
    if (myMinEdgeLength > sele->MinEdgeSize())
      myMinEdgeLength = sele->MinEdgeSize();

    if ( myMaxAreaSl < sele->MoData().Area() )
      myMaxAreaSl = sele->MoData().Area();
  }

  // ------------------------------------------------------------------------
  // loop over all master elements
  for (int i=0; i<melerowmap_->NumMyElements(); ++i)
  {
    int gid = melerowmap_->GID(i);
    DRT::Element* ele = idiscret_->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find master element with gid %i",gid);

    MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);
    if (myMinEdgeLength > mele->MinEdgeSize())
      myMinEdgeLength = mele->MinEdgeSize();

    // no Mortar Data container on the master side
    double myarea = mele->ComputeArea();
    if ( myMinAreaMa > myarea )
      myMinAreaMa = myarea;

    switch ( ele->Shape() )
    {
      case DRT::Element::tri3:
      case DRT::Element::tri6:
        myTriangleOnMaster = 1;
        break;
      default:
        // do nothing
        break;
    }
  }

  // ------------------------------------------------------------------------
  // communicate the minimal edge length and element area over all procs
  double lMins[2] = { myMinEdgeLength, myMinAreaMa };
  double gMins[2] = { 1.0e12, 1.0e12 };
  Comm().MinAll(&lMins[0],&gMins[0],2);

  if (gMins[0] == 1.0e12 or gMins[0] < 0.0)
    dserror("ERROR: Global minimal interface edge length couldn't"
        " be calculated! (value: %d)",gMins[0]);

  // set penetration bound to 50% of the minimal interface element edge length,
  // if no user defined value is provided.
  if (penBound_ < 0.0)
    penBound_ = 0.5*gMins[0];

  if (gMins[1] == 1.0e12 or gMins[1] < 0.0)
    dserror("ERROR: Global minimal master element area couldn't"
        " be calculated! (value: %d)",gMins[1]);

  const double gMinAreaMa = gMins[1];

  // ------------------------------------------------------------------------
  // communicate the maximal slave element area over all procs
  double lMaxs[2] = { myMaxAreaSl, static_cast<double>( myTriangleOnMaster ) };
  double gMaxs[2] = { 0.0, 0.0 };
  Comm().MaxAll( &lMaxs[0], &gMaxs[0], 2 );

  const double gMaxAreaSl = gMaxs[0];
  const bool isTriangleOnMaster = static_cast<bool>( gMaxs[1] );

  maxNumMasterElements_ = std::ceil( gMaxAreaSl / gMinAreaMa );

  if ( isTriangleOnMaster )
  {
    // approximated number of surrounding elements
    const int maxNumSurrounding = 4 * std::ceil(
        std::sqrt( static_cast<double>( maxNumMasterElements_ ) / 2.0 ) );

    // two elements in the corner
    const int numCornerEles = 2;

    // complete number of approximated master elements
    maxNumMasterElements_ += maxNumSurrounding + numCornerEles;
  }
  else
  {
    // approximated number of surrounding elements
    const int maxNumSurrounding = 2 * std::ceil(
        std::sqrt( static_cast<double>( maxNumMasterElements_ ) ) );

    // one element in the corner
    const int numCornerEles = 1;

    maxNumMasterElements_ += maxNumSurrounding + numCornerEles;
  }

  std::cout << "maxNumMasterElements_ = " << maxNumMasterElements_ << std::endl;

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // call initialization routine of the contact_interface
  CoInterface::Initialize();

  // setup member variables (has to be done only once)
  Setup();

  // *** Reset all augmented quantities *********************************
  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  const int nummynodes = SlaveColNodesBound()->NumMyElements();
  int * mynodegids = SlaveColNodesBound()->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    cnode->InitializeAugDataContainer( maxNumMasterElements_ );

    // reset nodal weighted gap
    cnode->AugData().GetWGap() = 1.0e12;

    // reset nodal scaling factor
    cnode->AugData().GetKappa() = 1.0e12;
    cnode->AugData().GetAugA()  = 0.0;

    // reset variation of the weighted gap
    cnode->AugData().GetVarWGapSl().clear();
    cnode->AugData().GetVarWGapMa().clear();

    // reset kappaLin
    cnode->AugData().GetKappaLin().clear();
    cnode->AugData().GetAugALin().clear();

    // reset varWGapLin
    cnode->AugData().GetVarWGapLinSl().clear();
    cnode->AugData().GetVarWGapLinMa().clear();

    // reset linearization of the averaged weighted gap
    cnode->AugData().GetAWGapLin().clear();

    // reset linearization of the averaged weighted gap
    cnode->AugData().GetWGapLin().clear();
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::UpdateMasterSlaveSets()
{
  MORTAR::MortarInterface::UpdateMasterSlaveSets();
  SplitSlaveDofs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::RedEvaluate(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
// interface needs to be complete
  if (!Filled() && Comm().MyPID()==0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  const int nummyelements = selecolmap_->NumMyElements();
  int * myelementgids = selecolmap_->MyGlobalElements();

  for ( int i=0; i<nummyelements; ++i )
  {
    const int gid1 = myelementgids[i];

    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(ele1);

    if (selement->ZeroSized())
      continue;

    /**************************************************************************
     *    Integrate all remaining quantities only over the slave interface    *
     **************************************************************************/
    // create a Augmented Lagrangian integrator instance with correct NumGP and Dim
    CONTACT::AUG::IntegrationWrapper augIntegrationWrapper(IParams(),selement->Shape(),Comm());
    switch ( Dim() )
    {
      case 2:
      case 3:
        augIntegrationWrapper.IntegrateDerivSlaveElement((*selement),Comm(),mparams_ptr);
        break;
      default:
        dserror("ERROR: RedEvaluate: Dim value has to be 2 or 3!");
        exit( EXIT_FAILURE );
    }
    /**************************************************************************
     *                       !!! end integration !!!                          *
     **************************************************************************/
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugDnMnMatrix(
    LINALG::SparseMatrix& augDnMatrix,
    LINALG::SparseMatrix& augMnMatrix) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = snoderowmap_->NumMyElements();
  int * mynodegids = snoderowmap_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

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
    const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap =
        cnode->AugData().GetVarWGapSl();

    for(CII pp=varWGapSlMap.begin();pp!=varWGapSlMap.end();++pp)
    {
      int cDofId  = pp->first;
      double cval = -(pp->second).second;
      augDnMatrix.Assemble(cval,rowId,cDofId);
    }
    /************************************* augMn-Matrix ******/
    const std::map<int,std::pair<int,double> >& varWGapMaMap =
        cnode->AugData().GetVarWGapMa();

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
void CONTACT::AUG::Interface::AssembleDGLmLinMatrix(
    LINALG::SparseMatrix& dGLmSlLinMatrix,
    LINALG::SparseMatrix& dGLmMaLinMatrix) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid  = mynodegids[i];

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    // get Lagrange multiplier in normal direction
    double lm_n = cnode->MoData().lm()[0];

    // typedef of the map iterator
    typedef std::map<int,std::map<int,double> >::const_iterator CII;

    /* --------------------------- SLAVE SIDE ----------------------------------- */
    const std::map<int,std::map<int,double> >& varWGapLinSlMap =
        cnode->AugData().GetVarWGapLinSl();

    // iteration over ALL slave Dof Ids
    for (CII p=varWGapLinSlMap.begin();p!=varWGapLinSlMap.end();++p)
    {
      const int sRow = p->first;

      // *** linearization of varWGap w.r.t. displacements ***
      AssembleMapIntoMatrix( sRow, -lm_n, varWGapLinSlMap.at(sRow), dGLmSlLinMatrix );
    }
    /* --------------------------- MASTER SIDE ---------------------------------- */
    const std::map<int,std::map<int,double> >& varWGapLinMaMap = cnode->AugData().GetVarWGapLinMa();

    // iteration over ALL master Dof Ids
    for (CII p=varWGapLinMaMap.begin();p!=varWGapLinMaMap.end();++p)
    {
      const int mRow = p->first;

      // *** linearization of varWGap w.r.t. displacements ***
      AssembleMapIntoMatrix( mRow, lm_n, varWGapLinMaMap.at(mRow), dGLmMaLinMatrix );
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDGGLinMatrix(
    LINALG::SparseMatrix& dGGSlLinMatrix,
    LINALG::SparseMatrix& dGGMaLinMatrix,
    const Epetra_Vector& cnVec ) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %",gid);

    // get nodal values
    double cn     = cnVec[cnVec.Map().LID(gid)];
    double aWGap  = cnode->AugData().GetWGap();
    double ckappa = cnode->AugData().GetKappa();
    double ckappainv = 1/ckappa;

    if ( aWGap != 1.0e12 and ckappa != 1.0e12 )
      aWGap *= ckappainv;
    else
      dserror("ERROR: GID % seems to be no active slave node!",gid);

    const std::map<int,double>& aWGapLinMap = cnode->AugData().GetAWGapLin();

    const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap =
        cnode->AugData().GetVarWGapSl();
    const std::map<int,std::pair<int,double> >& varWGapMaMap =
        cnode->AugData().GetVarWGapMa();

    /* --------------------------- SLAVE SIDE --------------------------------*/
    AssembleDGGLinMatrixOnSlaveSide( *cnode, varWGapSlMap, aWGapLinMap, cn,
        aWGap, dGGSlLinMatrix );

    /* --------------------------- MASTER SIDE -------------------------------*/
    AssembleDGGLinMatrixOnMasterSide( *cnode, varWGapMaMap, aWGapLinMap, cn,
        aWGap, dGGMaLinMatrix );
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDGGLinMatrixOnSlaveSide(
    const CoNode& cnode,
    const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap,
    const std::map<int,double>& aWGapLinMap,
    double cn,
    double aWGap,
    LINALG::SparseMatrix& dGGSlLinMatrix ) const
{
  typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

  // iteration over ALL slave Dof Ids
  for ( CI p=varWGapSlMap.begin(); p!=varWGapSlMap.end(); ++p )
  {
    const int sRow = p->first;

    // *** linearization of varWGap w.r.t. displacements ***
    const std::map<int,double>& varWGapLinSlMap = cnode.AugData().GetVarWGapLinSl().at(sRow);
    AssembleMapIntoMatrix( sRow, cn*aWGap, varWGapLinSlMap, dGGSlLinMatrix );

    // *** linearization of the averaged weighted gap w.r.t. displacements ***
    AssembleMapIntoMatrix( sRow, cn*(p->second).second, aWGapLinMap, dGGSlLinMatrix );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDGGLinMatrixOnMasterSide(
    const CoNode& cnode,
    const std::map<int,std::pair<int,double> >& varWGapMaMap,
    const std::map<int,double>& aWGapLinMap,
    double cn,
    double aWGap,
    LINALG::SparseMatrix& dGGMaLinMatrix ) const
{
  typedef std::map<int,std::pair<int,double> >::const_iterator CI;

  // iteration over ALL master Dof Ids
  for ( CI p=varWGapMaMap.begin(); p!=varWGapMaMap.end(); ++p )
  {
    const int mRow = p->first;

    // *** linearization of varWGap w.r.t. displacements ***
    const std::map<int,double>& varWGapLinMaMap = cnode.AugData().GetVarWGapLinMa().at(mRow);
    AssembleMapIntoMatrix( mRow, -cn*aWGap, varWGapLinMaMap, dGGMaLinMatrix );

    // *** linearization of the averaged weighted gap w.r.t. displacements ***
    AssembleMapIntoMatrix( mRow, -cn*(p->second).second, aWGapLinMap, dGGMaLinMatrix );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleResidualVectors(
    Epetra_Vector& lmNVec,
    Epetra_Vector& aWGapVec,
    Epetra_Vector& wGapVec) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get cn
//  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int* mynodegids = activenodes_->MyGlobalElements();

  int* myactiven_gids = activen_->MyGlobalElements();
  if ( nummynodes != activen_->NumMyElements() )
    dserror( "Dimension mismatch!" );

  double* lmNVec_vals = lmNVec.Values();
  double* aWGapVec_vals = aWGapVec.Values();
  double* wGapVec_vals = wGapVec.Values();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGLmrhs: Cannot find slave"
        " node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    double lmn  = cnode->MoData().lm()[0];

    // calculate averaged weighted gap
    const double wGap = cnode->AugData().GetWGap();
    double aWGap = wGap;
    double kappa = cnode->AugData().GetKappa();
    if ( kappa != 1.0e12 and aWGap != 1.0e12 )
    {
      aWGap /= kappa;

      const int rgid = myactiven_gids[ i ];

      // --- normal lagrange multiplier vector
      int rlid = lmNVec.Map().LID( rgid );
      if ( rlid == -1 )
        dserror("Sparse vector lmNVec does not have global row %d", rgid);

      lmNVec_vals[ rlid ] += lmn;

      // --- averaged weighted gap vector
      rlid = aWGapVec.Map().LID( rgid );
      if ( rlid == -1 )
        dserror("Sparse vector aWGapVec does not have global row %d", rgid);

      aWGapVec_vals[ rlid ] += aWGap;

      // --- weighted gap vector
      rlid = wGapVec.Map().LID( rgid );
      if ( rlid == -1 )
        dserror("Sparse vector wGapVec does not have global row %d", rgid);

      wGapVec_vals[ rlid ] += wGap;
    }
    else
      dserror("ERROR: Kappa and/or the weighted gap should "
          "not be equal 1.0e12 for active nodes!");
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDLmNWGapLinMatrix(
    LINALG::SparseMatrix& dLmNWGapLinMatrix) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all active augmented slave nodes of the interface
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find node with gid %",gid);

    double wGap = cnode->AugData().GetWGap();

    if (wGap==1.0e12)
      dserror("ERROR: The weighted gap should not be equal 1.0e12 "
          "for active nodes!");

    const int rowId = activen_->GID(i);

    // linearization of the weighted gap
    const std::map<int,double>& wGapLinMap = cnode->AugData().GetWGapLin();
    AssembleMapIntoMatrix( rowId, 1.0, wGapLinMap, dLmNWGapLinMatrix );
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveRhs(
    Epetra_Vector& augInactiveRhs,
    Epetra_Vector& cnVec) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  const int nummynodes = augInactiveSlaveNodes->NumMyElements();
  int * mynodegids = augInactiveSlaveNodes->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleInactiverhs: "
          "Node ownership inconsistency!");

    double cn_inv   = 1/cnVec[cnVec.Map().LID(gid)];
    double* lm  = cnode->MoData().lm();
    double augA = cnode->AugData().GetAugA();

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
void CONTACT::AUG::Interface::AssembleDLmTLmTRhs(
    Epetra_Vector& dLmTLmTRhs) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTRhs: Cannot find active"
        " slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // Get the Lagrange multiplier and txi of the current node
    double* lm  = cnode->MoData().lm();
    // Get weighted element area
    double augA  = cnode->AugData().GetAugA();

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
void CONTACT::AUG::Interface::AssembleDLmTLmTMatrix(
    LINALG::SparseMatrix& dLmTLmTMatrix) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleDLmTLmTrhs: Cannot find slave"
        " node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // get weighted element area
    double augA = cnode->AugData().GetAugA();

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
void CONTACT::AUG::Interface::AssembleDLmTLmTLinMatrix(
    LINALG::SparseMatrix& dLmTLmTLinMatrix) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDGLmLin: Node ownership inconsistency!");

    double* lm  = cnode->MoData().lm();
    GEN::pairedvector<int,double>& augALinMap = cnode->AugData().GetAugALin();

    /*---------------------------------------------------------*
     | Attention:                                              |
     | D_{d}[-1/ct*(dlm_t * lm_t)] = D_{d}[-1/ct*[dlm_t*lm_t]] |
     *---------------------------------------------------------*/
    for (int j=0;j<(Dim()-1);++j)
    {
      const int rowId = activet_->GID((Dim()-1)*i+j);

      double tmp = ct_inv*lm[j+1];

      AssembleMapIntoMatrix( rowId, tmp, augALinMap, dLmTLmTLinMatrix );
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveDiagMatrix(
    Epetra_Vector& augInactiveDiagMatrix,
    const Epetra_Vector& cnVec) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  const int nummynodes = augInactiveSlaveNodes->NumMyElements();
  int * mynodegids = augInactiveSlaveNodes->MyGlobalElements();

  LINALG::SerialDenseVector vals;
  std::vector<int> rowIds(0);
  std::vector<int> rowner(0);

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find"
        " inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
    double augA   = cnode->AugData().GetAugA();

    const int numdof = cnode->NumDof();

    if ( vals.Length() != numdof )
    {
      vals.Resize( numdof );
      rowIds.resize( numdof, -1 );
      rowner.resize( numdof, -1 );
    }

    // normal direction
    vals( 0 ) = cn_inv*augA;

    // tangential directions
    std::fill( vals.A() + 1, vals.A() + numdof, ct_inv*augA );

    // copy dof ids
    std::copy( cnode->Dofs(), cnode->Dofs() + numdof, &rowIds[0] );

    // insert owner
    std::fill( &rowner[0], &rowner[0] + numdof, cnode->Owner() );

    LINALG::Assemble( augInactiveDiagMatrix, vals, rowIds, rowner );
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveLinMatrix(
    LINALG::SparseMatrix& augInactiveLinMatrix,
    const Epetra_Vector& cnVec) const
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes  = LINALG::SplitMap(
      *snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = IParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1/ct_inv;

  const int nummynodes = augInactiveSlaveNodes->NumMyElements();
  int * mynodegids = augInactiveSlaveNodes->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: AssembleAugInactiveMatrix: Cannot find"
        " inactive slave node with gid %",gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
    double* lm  = cnode->MoData().lm();
    GEN::pairedvector<int,double>& augALinMap = cnode->AugData().GetAugALin();

    // typedef of the map iterator
    typedef GEN::pairedvector<int,double>::const_iterator CI;

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

        if (abs(val)>1.0e-12)
          augInactiveLinMatrix.Assemble(val,rowId,colId);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::WGap() const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = snoderowmap_->NumMyElements();
  int * mynodegids = snoderowmap_->MyGlobalElements();

  for ( int i=0;i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGGrhs: Cannot find slave"
        " node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap =
        cnode->AugData().GetVarWGapSl();

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
      const std::map<int,std::pair<int,double> >& varWGapMaMap =
          cnode->AugData().GetVarWGapMa();
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
void CONTACT::AUG::Interface::AWGapLin()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: AssembleDGGrhs: Cannot find slave node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleGGrhs: Node ownership inconsistency!");

    // get some pre-calculated results for the current slave node:
    double kappainv = 1.0/cnode->AugData().GetKappa();

    if (cnode->AugData().GetKappa() == 1.0e12)
      dserror("ERROR: CONTACT::AugmentedInterface::AWGapLin(): kappa shouldn't be "
        "equal 1.0e12 for active nodes!");

    const std::map<int,double>& kappaLinMap  = cnode->AugData().GetKappaLin();

    // Get std::map to save the linearization of the weighted gap
    std::map<int,double>& wGapLinMap = cnode->AugData().GetWGapLin();
    // Get std::map to save the linearization of the averaged weighted gap
    std::map<int,double>& aWGapLinMap = cnode->AugData().GetAWGapLin();

    // typedef of the constant map iterators
    typedef std::map<int,double> ::const_iterator CI;
    typedef std::map<int,std::pair<int,double> >::const_iterator CII;
    typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CIII;

    // *** Slave side **************************************************
    GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap = cnode->AugData().GetVarWGapSl();
    // loop over ALL Dofs of the slave side
    for (CIII p=varWGapSlMap.begin();p!=varWGapSlMap.end();++p)
    {
      int sGid = (p->second).first;

      CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(sGid));
      if (!snode) dserror("ERROR: Cannot find slave node with gid %",sGid);
      // Check if Slave node
      if (!snode->IsSlave()) dserror("ERROR: This has to be a slave node!");

      // Switch from the global DofId to the local DofId
      int sDof  = ((p->first)-snode->Dofs()[0])%Dim();

      const double xs = snode->xspatial()[sDof];

      // Lin(kappa)
      AddKappaLinToGapLinearization( kappaLinMap, xs, kappainv,
          (p->second).second, 1.0, aWGapLinMap );

      // Lin(varWGapLinSl)
      std::map<int,double>& varWGapLinSlMap = cnode->AugData().GetVarWGapLinSl()[p->first];
      for (CI pp=varWGapLinSlMap.begin();pp!=varWGapLinSlMap.end();++pp)
      {
        const double val = (pp->second)*xs;
        wGapLinMap[pp->first]  -= val;
        aWGapLinMap[pp->first] -= val*kappainv;
      }

      // Lin(xs)
      wGapLinMap[p->first]  -= (p->second).second;
      aWGapLinMap[p->first] -= (p->second).second*kappainv;
    }

    // *** Master side *************************************************
    std::map<int,std::pair<int,double> >& varWGapMaMap = cnode->AugData().GetVarWGapMa();
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

      const double xm = mnode->xspatial()[mDof];

      // Lin(kappa)
      AddKappaLinToGapLinearization( kappaLinMap, xm, kappainv,
          (p->second).second, -1.0, aWGapLinMap );


      // Lin(varWGapLinMa)
      std::map<int,double>& varWGapLinMaMap = cnode->AugData().GetVarWGapLinMa()[p->first];
      for (CI pp=varWGapLinMaMap.begin();pp!=varWGapLinMaMap.end();++pp)
      {
        const double val = (pp->second)*xm;
        wGapLinMap[pp->first]  += val;
        aWGapLinMap[pp->first] += val*kappainv;
      }

      // Lin(xm)
      wGapLinMap[p->first]  += (p->second).second;
      aWGapLinMap[p->first] += (p->second).second*kappainv;
    }

  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AddKappaLinToGapLinearization(
    const std::map<int,double>& kappaLinMap,
    double x,
    double kappainv,
    double varWGap,
    double scale,
    std::map<int,double>& aWGapLinMap ) const
{
  const double val = scale * varWGap * x * kappainv * kappainv;
  for ( std::map<int,double> ::const_iterator pp=kappaLinMap.begin();
        pp!=kappaLinMap.end(); ++pp )
    aWGapLinMap[pp->first] += (pp->second)*val;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleContactPotentialTerms(
    const Epetra_Vector& cnVec,
    double& zn_gn,
    double& gn_gn,
    double& zn_zn,
    double& zt_zt ) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  double ct_inv = 1.0 / ct_;

  // *** Active part *************************************************
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  int * mynodegids = activenodes_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    const double cn = cnVec[cnVec.Map().LID(gid)];

    const double wgap = cnode->AugData().GetWGap();
    dsassert( wgap!=1.0e12, "ERROR: WGap is equal 1.e12 for a active node!" );

    const double kappa = cnode->AugData().GetKappa();
    dsassert( kappa!=1.0e12, "ERROR: Kappa is equal 1.e12 for a active node! (node-id: %d)" );

    const double augA = cnode->AugData().GetAugA();

    // get the Lagrange multiplier
    const double* lm  = cnode->MoData().lm();

    // *** ACTIVE - NORMAL DIRECTION ***
    if (wgap != 0.0)
    {
      // ** zn_i * awgap_i * A_i **
      zn_gn += wgap*lm[0];
      // ** cn/2 * awgap * awgap * A_i **
      gn_gn += 0.5 * cn * wgap * wgap / kappa;
    }

    // *** ACTIVE - TANGENTIAL DIRECTION ***
    for (int d=1;d<Dim();++d)
    {
      // ** 1/(ct) * zt_i^T*zt_i * A_i **
      zt_zt += ct_inv * lm[d] * lm[d] * augA;
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

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    const double cn_inv = 1/cnVec[cnVec.Map().LID(gid)];
    const double augA = cnode->AugData().GetAugA();

    // get the lagrange multiplier
    const double* lm  = cnode->MoData().lm();

    // *** INACTIVE - NORMAL DIRECTION ***
    // ** 1/(cn) * zn_i * zn_i * A_i **
    zn_zn += cn_inv * lm[0] * lm[0] * augA;

    // *** INACTIVE - TANGENTIAL DIRECTION ***
    for (int d=1;d<Dim();++d)
    {
      // ** 1/(ct) * zt_i^T*zt_i * A_i **
      zt_zt += ct_inv * lm[d] * lm[d] * augA;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Interface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> myactivenodegids(0);
  std::vector<int> myactivedofgids(0);

  // loop over all slave nodes
  const int nummynodes = snoderowmap_->NumMyElements();
  int * mynodegids = snoderowmap_->MyGlobalElements();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

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
      if (cnode->IsInitActive() or (initcontactbygap and cnode->AugData().GetWGap() < initcontactval))
      {
        cnode->Active()=true;
        myactivenodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          myactivedofgids.push_back(cnode->Dofs()[j]);
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
        myactivenodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          myactivedofgids.push_back(cnode->Dofs()[j]);
      }
    }
  } // end loop over all slave nodes

  // create interface local augmented active node map and augmented active dof map
  activenodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)myactivenodegids.size(),
      &myactivenodegids[0],0,Comm()));
  activedofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)myactivedofgids.size(),
      &myactivedofgids[0],0,Comm()));

  SplitAugActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SplitAugActiveDofs()
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

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    CoNode* cnode = static_cast<CoNode*>( node );

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
void CONTACT::AUG::Interface::SplitSlaveDofs()
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
void CONTACT::AUG::Interface::AssembleAugAVector(
    Epetra_Vector& augAVec,
    Epetra_Vector& kappaVec) const
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = snoderowmap_->NumMyElements();
  int * mynodegids = snoderowmap_->MyGlobalElements();

  double* augA_values = augAVec.Values();
  double* kappa_values = kappaVec.Values();

  for ( int i=0; i<nummynodes; ++i )
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    // *** augmented Area ***
    const double augA  = cnode->AugData().GetAugA();

    if (augA > 0.0)
    {
      const int lid = augAVec.Map().LID( gid );
      if ( lid == -1 )
        dserror("Sparse vector augAVec does not have global row %d", gid);

      augA_values[ lid ] += augA;
    }
    else
      dserror("ERROR: The augmented nodal area shouldn't be equal/lower than zero! "
          "(value= %.2e)",augA);

    // *** kappa ***
    if (cnode->Active())
    {
      const double kappa  = cnode->AugData().GetKappa();

      if (kappa > 0.0)
      {
        const int lid = kappaVec.Map().LID( gid );
        if ( lid == -1 )
          dserror("Sparse vector kappaVec does not have global row %d", gid);

        kappa_values[ lid ] += kappa;
      }
      else
        dserror("ERROR: The weighted area kappa shouldn't be equal/lower than zero! "
            "(value= %.2e)",augA);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < class T >
void CONTACT::AUG::Interface::AssembleMapIntoMatrix(
    int row,
    double scal,
    const T& values,
    LINALG::SparseMatrix& mat,
    double threshold ) const
{
  dsassert( threshold >= 0.0, "The threshold value has to be positive!" );

  for ( typename T::const_iterator pp = values.begin();
        pp != values.end(); ++pp )
  {
    double val = scal * ( pp->second );
    const int col = pp->first;

    if ( std::abs( val ) < threshold )
      continue;

    switch ( mat.GetMatrixtype() )
    {
      case LINALG::SparseMatrix::FE_MATRIX:
        mat.FEAssemble( val, row, col );
        break;
      case LINALG::SparseMatrix::CRS_MATRIX:
        mat.Assemble( val, row, col );
        break;
    }
  }
}

template void CONTACT::AUG::Interface::AssembleMapIntoMatrix< GEN::pairedvector<int,double> >(
    int row,
    double scal,
    const GEN::pairedvector<int,double>& values,
    LINALG::SparseMatrix& mat,
    double threshold ) const;
template void CONTACT::AUG::Interface::AssembleMapIntoMatrix< std::map<int,double> >(
    int row,
    double scal,
    const std::map<int,double>& values,
    LINALG::SparseMatrix& mat,
    double threshold ) const;
