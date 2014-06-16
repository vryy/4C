/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsi_crack.cpp

\brief Adapter Layer for FSI with cracking structure

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15267
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsi_crack.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_crack/InsertCohesiveElements.H"
#include "../drt_crack/dcohesive.H"
#include "../drt_crack/crackUtils.H"

#include "../linalg/linalg_utils.H"

/*======================================================================*/
/* constructor */
ADAPTER::FSICrackingStructure::FSICrackingStructure
(
  Teuchos::RCP<Structure> stru
)
:FSIStructureWrapper(stru)
{
  firstTime_ = true;

  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create the underlying structural adapter");

  structdis_ = DRT::Problem::Instance()->GetDis("structure");
  if( structdis_ == Teuchos::null )
    dserror("Structure discretization not found\n");

  //TODO: Check whether fluid discretization is required here..
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  if( fluiddis_ == Teuchos::null )
    dserror("Fluid discretization not found\n");

  //TODO: should be read from input file
  minOpeningDist_ = 0.001;
}


/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 03/14
 *------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::addNodesToConditions( DRT::Condition * cond,
                                                          std::map<int,int> oldnew )
{
  std::vector<int> add;

  for( std::map<int,int>::iterator itm = oldnew.begin(); itm != oldnew.end(); itm++ )
  {
    add.push_back(itm->first);
    add.push_back(itm->second);
  }
  this->addNodesToConditions( cond, add );
}

/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 03/14
 *------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::addNodesToConditions( DRT::Condition * cond,
                                                          std::vector<int> add )
{
  if( add.size() == 0 )
    dserror("No nodes passed for adding to condition\n");

  const std::vector<int>* conNodes = cond->Nodes();
  for( unsigned ii = 0; ii < conNodes->size(); ii++ )
      add.push_back( (*conNodes)[ii] );

  Teuchos::RCP<std::vector<int> > storage = Teuchos::rcp( new std::vector<int>(add) );

  std::sort( storage->begin(), storage->end() );
  cond->Add( "Node Ids", storage );
}

/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::addNodesToConditions( DRT::Condition * cond,
                                                          std::set<int> add )
{
  std::vector<int> vec (add.size());

  int count = 0;
  for( std::set<int>::iterator it = add.begin(); it != add.end(); it++ )
  {
    vec[count] = *it;
    count++;
  }

  this->addNodesToConditions( cond, vec );
}

/*------------------------------------------------------------------------------------*
 * Add newly formed crack surfaces into cut discretization                    sudhakar 03/14
 * This is the only function that communicates with FSI implementation
 *------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                                std::map<int, LINALG::Matrix<3,1> >& tip_nodes )
{
  //RebuildInterfaceWithConditionCheck( boundary_dis, tip_nodes );

  RebuildInterfaceWithoutConditionCheck( boundary_dis, tip_nodes );
}

/*------------------------------------------------------------------------------------*
 * Rebuild FSI interface only when the crack mouth opening displacement        sudhakar 03/14
 * reaches a predefined given value.
 *------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::RebuildInterfaceWithConditionCheck( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                                        std::map<int, LINALG::Matrix<3,1> >& tip_nodes )
{
  std::map<int,int> tips = getOldNewCrackNodes();
  std::vector<int> crtip = GetCrackTipNodes();

  if( tips.size() == 0 )
    dserror("atleast one new node should have been added during the crack propagation\n");

  tipHistory_.push_back( tips );

  if( not isCrackMouthOpeningConditionSatisfied( tipHistory_[0] ) )
    return;

  // We do not explicitly add new elements to the discretization
  // Instead, we add the new nodes to FSI condition in structual discretization
  // And we build a new discretization based on this new condition
  DRT::Condition* cond_fsi = structdis_->GetCondition("FSICoupling");
  DRT::Condition* cond_xfem = structdis_->GetCondition("XFEMCoupling");

  if( cond_fsi == NULL or cond_xfem == NULL )
    dserror( "XFEM or FSI coupling conditions undefined in XFSI problem?\n" );

  bool addtip = false;

  //if( not firstTime_ )
  {
    addNodesToConditions( cond_fsi, tipHistory_[0] );
    addNodesToConditions( cond_xfem, tipHistory_[0] );
    tipHistory_.erase( tipHistory_.begin() );
    firstTime_ = false;
  }

  if( tipHistory_.size() == 0 )
  {
    addNodesToConditions( cond_fsi, crtip );
    addNodesToConditions( cond_xfem, crtip );
    addtip = true;
  }


  // Fillcomplete() is necessary here because we introduced new nodes into conditions
  // In order to build discretization based on this condition, it is necessary to rebuild
  // the geometry of the conditions
  structdis_->FillComplete();

  // build boundary discretization based on FSI condition
  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundary_dis = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentIndependentDofSet(structdis_,true,Teuchos::null));
  boundary_dis->ReplaceDofSet(newdofset);

  #if 0 // try to add element
    tip_nodes.clear();

    int finalids[4] = { 132, 133, 39, 38 };


    int neweleid = boundary_dis->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3_3","quad4", neweleid, 0/*boundary_dis2_->gNode(nodeids[0])->Owner()*/ );
    spr->SetNodeIds( 4, finalids );
    boundary_dis->AddElement( spr );

  #endif

    boundary_dis->FillComplete();

    tip_nodes.clear();
  if( not addtip )
  {
    //TODO: parallelize this process
    Teuchos::RCP<const Epetra_Vector> disp = Dispnp();
    for( std::map<int,int>::iterator itmap = tipHistory_[0].begin(); itmap != tipHistory_[0].end(); itmap++ )
    {
      int nid1 = itmap->first;
      int nid2 = itmap->second;

      DRT::Node* node1 = structdis_->gNode( nid1 );
      DRT::Node* node2 = structdis_->gNode( nid2 );

      std::vector<int> lm1,lm2;
      std::vector<double> myval1, myval2;
      structdis_->Dof( node1, lm1 );
      structdis_->Dof( node2, lm2 );
      DRT::UTILS::ExtractMyValues( *disp, myval1, lm1 );
      DRT::UTILS::ExtractMyValues( *disp, myval2, lm2 );

      LINALG::Matrix<3,1> displ(true);
      for( unsigned i=0; i<3; i++ )
        displ(i,0) = 0.5 * (myval1[i]+myval2[i]);

      std::cout<<"disp 1 = "<<myval1[0]<<"\t"<<myval1[1]<<"\t"<<myval1[2]<<"\n";
      std::cout<<"disp 2 = "<<myval2[0]<<"\t"<<myval2[1]<<"\t"<<myval2[2]<<"\n";

      tip_nodes[nid1] = displ;
      tip_nodes[nid2] = displ;
    }
  }

  // TODO: This has not been checked in parallel
  // Now it is mandatory to distribute all the nodes and elements to all processors
  const Epetra_Map noderowmap = *boundary_dis->NodeRowMap();
  const Epetra_Map elemrowmap = *boundary_dis->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  boundary_dis->ExportColumnNodes(nodecolmap);
  boundary_dis->ExportColumnElements(elemcolmap);

  boundary_dis->FillComplete();
}

/*------------------------------------------------------------------------------------*
 * Check whether crack mouth opening displacement reached the predefined value    sudhakar 03/14
 * Only when this condition is satisfied, new crack surfaces are added to the
 * cut discretization
 *------------------------------------------------------------------------------------*/
bool ADAPTER::FSICrackingStructure::isCrackMouthOpeningConditionSatisfied( const std::map<int,int>& tips )
{
  if( tipHistory_.size() == 0 )
    dserror( "No tip nodes stored until now\n" );

  Teuchos::RCP<const Epetra_Vector> disp = Dispnp();
  int nid1 = tips.begin()->first;
  int nid2 = tips.begin()->second;

  DRT::Node* node1 = structdis_->gNode( nid1 );
  DRT::Node* node2 = structdis_->gNode( nid2 );

  std::vector<int> lm1,lm2;
  std::vector<double> myval1, myval2;
  structdis_->Dof( node1, lm1 );
  structdis_->Dof( node2, lm2 );
  DRT::UTILS::ExtractMyValues( *disp, myval1, lm1 );
  DRT::UTILS::ExtractMyValues( *disp, myval2, lm2 );

  double dist = 0.0;
  for( unsigned dim = 0; dim < 3; dim++ )
    dist += pow(( myval1[dim] - myval2[dim] ),2);

  dist = sqrt(dist);

  if( dist >= minOpeningDist_ )
    return true;
  return false;
}

/*------------------------------------------------------------------------------------------------------*
 * Rebuild FSI interface without checking any conditions                                          sudhakar 04/14
 * This means even if the crack mouth opening is infinitely small, the crack surfaces
 * are added to cut discretization; so fluid will enter into this which may
 * create some problems
 *------------------------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::RebuildInterfaceWithoutConditionCheck( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                                            std::map<int, LINALG::Matrix<3,1> >& tip_nodes )
{
  std::map<int,int> tips = getOldNewCrackNodes();
  std::vector<int> crtip = GetCrackTipNodes();

  if( tips.size() == 0 )
    dserror("atleast one new node should have been added during the crack propagation\n");

  // We do not explicitly add new elements to the discretization
  // Instead, we add the new nodes to FSI condition in structual discretization
  // And we build a new discretization based on this new condition
  DRT::Condition* cond_fsi = structdis_->GetCondition("FSICoupling");
  DRT::Condition* cond_xfem = structdis_->GetCondition("XFEMCoupling");

  if( cond_fsi == NULL or cond_xfem == NULL )
    dserror( "XFEM or FSI coupling conditions undefined in XFSI problem?\n" );

  {
    if( not firstTime_ )
    {
    addNodesToConditions( cond_fsi, tips );
    addNodesToConditions( cond_xfem, tips );
    }

    // activate this if you are adding crack nodes irrespective of the opening displacement
    addNodesToConditions( cond_fsi, crtip );
    addNodesToConditions( cond_xfem, crtip );
    firstTime_ = false;
  }

  // Fillcomplete() is necessary here because we introduced new nodes into conditions
  // In order to build discretization based on this condition, it is necessary to rebuild
  // the geometry of the conditions
  structdis_->FillComplete();

  Teuchos::RCP<DRT::Discretization> temp_dis = Teuchos::null;

  // build a temporary discretization based on FSI condition
  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  temp_dis = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentIndependentDofSet(structdis_,true,Teuchos::null));
  temp_dis->ReplaceDofSet(newdofset);

  temp_dis->FillComplete();

  tip_nodes.clear();


#if 0 //seems to be easier implementation. but with some strange errors

  //distributeDisToAllProcs( temp_dis );

  std::map<int,int> changeEleIds;
  std::vector<int> notFound;
  const int numcolele = temp_dis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = temp_dis->lColElement(i);
    int ele_other_dis_id=0; //id of the element which has the same node ids in other discretization
    if( checkElementExist( boundary_dis, actele, ele_other_dis_id ) )
      changeEleIds[actele->Id()] = ele_other_dis_id;
    else
      notFound.push_back(actele->Id());
  }

  if( notFound.size() > 3*(crtip.size()-1) )
  {
    std::cout<<"the not found elements are\t";
    for(unsigned i=0;i<notFound.size();i++)
      std::cout<<notFound[i]<<"\t";
    std::cout<<"\n";
    dserror( "more elements not found on original discretization %d",notFound.size() );
  }

  // std::map<key,val> tips = std::map<val,key> rev_tips;
  std::map<int,int> rev_tips;
  for( std::map<int,int>::iterator it = tips.begin(); it != tips.end(); it++ )
    rev_tips[it->second] = it->first;

  std::set<int> tipeleIds;
  for( unsigned inot=0; inot<notFound.size(); inot++ )
  {
    if( temp_dis->HaveGlobalElement( notFound[inot] ) )
    {
      DRT::Element * ele = temp_dis->gElement( notFound[inot] );
      const int * elenodes = ele->NodeIds();
      int modinodes[4];

      bool tipele = false;
      for( int ino = 0; ino < ele->NumNode(); ino++ )
      {
        int elenodeid = elenodes[ino];
        if( std::find( crtip.begin(), crtip.end(), elenodeid ) != crtip.end() )
        {
          tipele = true;
          tipeleIds.insert( notFound[inot] );
          break;
        }
        else
        {
          std::map<int, int>::iterator delnod = rev_tips.find( elenodeid );
          if( delnod == rev_tips.end() )
            modinodes[ino] = elenodeid;
          else
            modinodes[ino] = delnod->second;
        }
      }
      if( not tipele )
      {
        int ele_other_dis_id;
        if( checkElementExist( boundary_dis, modinodes, ele_other_dis_id ) )
          changeEleIds[notFound[inot]] = ele_other_dis_id;
        else
          dserror("if an element is not found, it should be a tip element");
      }
    }
  }

  for( std::map<int,int>::iterator it = changeEleIds.begin(); it != changeEleIds.end(); it++ )
  {
    int old_eleid = it->first;
    DRT::Element* ele = temp_dis->gElement( old_eleid );
    ele->SetId( it->second );
  }

  int count_rev = 1;
  for(std::set<int>::iterator it = tipeleIds.begin(); it!=tipeleIds.end(); it++)
  {
    if( temp_dis->HaveGlobalElement( *it ) )
    {
      DRT::Element* ele = temp_dis->gElement( *it );
      ele->SetId( temp_dis->NumGlobalElements() - count_rev );
    }
    count_rev = count_rev + 1;
  }

  boundary_dis = Teuchos::null;
  boundary_dis = temp_dis;

  boundary_dis->FillComplete();

#else // checked already. working fine
  // --------------
  // Now we have two discretizations here
  // boundary_dis --> old discretization
  // temp_dis --> new one with updated elements - Fine to work with further, except for xfem time integration
  // In order to have xfem time integration working, we need to have same side ids in both old and new
  // discretizations. To have this, now we copy all newly created sides in new discretization to
  // the old one.
  // --------------
  //
  //       \       /               \                 /
  //        \     /                 \               /
  //         \   /                   \             /
  //          \ /                     \           /
  //           v                       *         *
  //                                    \       /
  //     BOUNDARY_DIS                    \     /
  //                                      \   /
  //                                       \ /
  //                                        v
  //
  //                                     TEMP_DIS
  //
  // Here, we modify boundary_dis to become same as that of temp_dis
  //
  //

  int total_num = boundary_dis->NumGlobalElements();

  //--------------------------------------------------------------------------------------------------
  // STEP 1: Add newly created nodes into boundary discretization
  //
  //       \       /                    \       /
  //        \     /                      \     /              * added nodes
  //         \   /        =======>        \   /
  //          \ /                          \ /
  //           v                     *      v
  //
  //
  //
  //                                     *
  //
  //--------------------------------------------------------------------------------------------------
  std::vector<int> addnodes;
  for( unsigned i=0; i<crtip.size(); i++ )
    addnodes.push_back(crtip[i]);
  for( std::map<int,int>::iterator it = tips.begin(); it != tips.end(); it++)
    addnodes.push_back( it->second );

  for( unsigned i=0; i< addnodes.size(); i++ )
  {
    int lmaster = 0, gmaster = 0;
    double xx[3]={0.0, 0.0, 0.0};
    if( temp_dis->HaveGlobalNode( addnodes[i] ) )
    {
      DRT::Node * nod = temp_dis->gNode( addnodes[i] );
      if( nod->Owner() == temp_dis->Comm().MyPID() )
      {
        lmaster = temp_dis->Comm().MyPID();
        for( unsigned j=0; j<3; j++ )
          xx[j] = nod->X()[j];
      }
    }
    temp_dis->Comm().SumAll( &lmaster, &gmaster, 1 );
    temp_dis->Comm().Broadcast( &xx[0], 3, gmaster );
    boundary_dis->AddNode( Teuchos::rcp(new DRT::Node(addnodes[i], xx, gmaster )));
  }

  //--------------------------------------------------------------------------------------------------
  // STEP 2: Add newly created elements into boundary discretization
  //
  //          \       /                                    \       /
  //           \     /                                      \     /
  //            \   /                                        \   /
  //             \ /            =======>                      \ /
  //       *      v                                    *       v
  //                                                    \     /
  //                                                     \   /
  //                                                      \ /
  //          *                                            V
  //--------------------------------------------------------------------------------------------------
  std::set<int> neweles; //new elements that are to be added to old discret.

  const int numrowele = temp_dis->NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = temp_dis->lRowElement(i);

    for( unsigned crid = 0; crid < crtip.size(); crid++ )
    {
      if( DRT::CRACK::UTILS::ElementHasThisNodeId( actele, crtip[crid] ) )
        neweles.insert( actele->Id() );
    }
  }

  LINALG::GatherAll( neweles, temp_dis->Comm() );

  for( std::set<int>::iterator it = neweles.begin(); it != neweles.end(); it++ )
  {
    int incr_id = 0; // this is necessary to make sure that a particular element has same id on all processors
    int lmaster = 0, gmaster = 0, numnodes = 0;
    int nodeids[4] = {0, 0, 0, 0};
    if( temp_dis->HaveGlobalElement( *it ) )
    {
      DRT::Element* elem = temp_dis->gElement( *it );
      //Teuchos::RCP<DRT::Element> adelem = DRT::UTILS::Factory("BELE3_3","quad4", total_num+incr_id, elem->Owner() );
      if( elem->Owner() == temp_dis->Comm().MyPID() )
      {
        lmaster = temp_dis->Comm().MyPID();
        numnodes = elem->NumNode();
        for( unsigned j=0; j<4; j++ )
          nodeids[j] = elem->NodeIds()[j];
      }
    }
    temp_dis->Comm().SumAll( &lmaster, &gmaster, 1 );
    temp_dis->Comm().Broadcast( &numnodes, 1, gmaster );
    temp_dis->Comm().Broadcast( &nodeids[0], 4, gmaster );
    Teuchos::RCP<DRT::Element> adelem = DRT::UTILS::Factory("BELE3_3","quad4", total_num++, gmaster );
    adelem->SetNodeIds( numnodes, nodeids );
    boundary_dis->AddElement( adelem );
    incr_id++;
  }

  //--------------------------------------------------------------------------------------------------
  // STEP 3: Modify the connectivity of existing elements
  //
  //            \       /                    \                 /
  //             \     /                      \               /
  //              \   /                        \             /
  //               \ /                          \           /
  //        *       v        =======>            *         *
  //         \     /                              \       /
  //          \   /                                \     /
  //           \ /                                  \   /
  //            v                                    \ /
  //                                                  v
  //
  //--------------------------------------------------------------------------------------------------
  distributeDisToAllProcs( temp_dis );

  std::set<int> modifyEle;
  for( std::map<int,int>::iterator it = tips.begin(); it != tips.end(); it++ )
  {
    int nodid = it->first;
    if( boundary_dis->HaveGlobalNode( nodid ) )
    {
      DRT::Node* nod = boundary_dis->gNode( nodid );
      if( nod->Owner() == boundary_dis->Comm().MyPID() )
      {
        DRT::Element** attele = nod->Elements();
        int numattele = nod->NumElement();

        for( int elno = 0; elno < numattele; elno++ )
        {
          DRT::Element * thisele = attele[elno];
          int thiseleid = thisele->Id();

          if( std::find( modifyEle.begin(), modifyEle.end(), thiseleid ) != modifyEle.end() )
            continue;

          int unused_id;
          if( not checkElementExist( temp_dis, thisele, unused_id ) )
          {
            modifyEle.insert( thiseleid );
          }
        }

      }
    }
  }

  LINALG::GatherAll( modifyEle, boundary_dis->Comm() );

  for( std::set<int>::iterator it = modifyEle.begin(); it != modifyEle.end(); it++ )
  {
    int eleid = *it;

    if( boundary_dis->HaveGlobalElement( eleid ) )
    {
      bool del = false;
      DRT::Element* ele = boundary_dis->gElement( eleid );
      const int * oldnodes = ele->NodeIds();
      std::vector<int> newnodes( ele->NumNode() );

      for( int i = 0; i < ele->NumNode(); i++ )
      {
        std::map<int, int>::iterator delnod = tips.find( oldnodes[i] );
        if( delnod != tips.end() )
        {
          del = true;
          newnodes[i] = delnod->second;
        }
        else
          newnodes[i] = oldnodes[i];
      }

      if( not del )
        dserror("This element should have atleast one replaceable node\n");

      if( newnodes.size() != static_cast<std::size_t>(ele->NumNode()) )
        dserror("Check the number of new nodes\n");

      if( del )
        ele->SetNodeIds(newnodes.size(), &newnodes[0]);
    }
  }

  //--------------------------------------------------------------------------------------------------
  // STEP 4: Set XFEM and FSI conditions to new nodes
  //--------------------------------------------------------------------------------------------------
  DRT::Condition* cond_fsi_boun = boundary_dis->GetCondition("FSICoupling");
  DRT::Condition* cond_xfem_boun = boundary_dis->GetCondition("XFEMCoupling");

  if( cond_fsi_boun == NULL or cond_xfem_boun == NULL )
    dserror( "XFEM or FSI coupling conditions undefined in old boundary discretization?\n" );

  addNodesToConditions( cond_fsi_boun, addnodes );
  addNodesToConditions( cond_xfem_boun, addnodes );

  boundary_dis->FillComplete();
#endif

  //--------------------------------------------------------------------------------------------------
  // STEP 5: Distribute boundary discretization to all processors
  //--------------------------------------------------------------------------------------------------
  distributeDisToAllProcs( boundary_dis );
}

/*---------------------------------------------------------------------------------------------------*
 * Distribute the discretization to all processors within its communicator                  sudhakar 05/14
 *---------------------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::distributeDisToAllProcs( Teuchos::RCP<DRT::Discretization>& dis )
{
  const Epetra_Map noderowmap = *dis->NodeRowMap();
  const Epetra_Map elemrowmap = *dis->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  dis->ExportColumnNodes(nodecolmap);
  dis->ExportColumnElements(elemcolmap);

  dis->FillComplete();
}

/*-----------------------------------------------------------------------------------------------------------*
 * Add the given structural node to the given discretization                                      sudhakar 05/14
 * Remember : the given node should be available in structdis_
 *-----------------------------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::addThisStructNodeToDis( int nodeid, Teuchos::RCP<DRT::Discretization>& dis )
{
  if( structdis_->HaveGlobalNode( nodeid ) )
  {
    DRT::Node* node = structdis_->gNode( nodeid );

    //Teuchos::RCP<DRT::Node> temp = Teuchos::RCP(*node);

    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( new DRT::Node( nodeid, node->X(), node->Owner() ) );
    dis->AddNode( newnode );
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 * Check whether an element exists in the discretization eith all nodal ids of given element            sudhakar 05/14
 * Assumes that 1. both elements are ordered in same direction
 * 2. given discretization is redundant on all processors
 *-----------------------------------------------------------------------------------------------------------*/
bool ADAPTER::FSICrackingStructure::checkElementExist( Teuchos::RCP<DRT::Discretization>& dis,
                                                       DRT::Element* ele,
                                                       int& eleid )
{

  const int* elenod = ele->NodeIds();
  if( not ele->NumNode() == 4 )
    dserror("we can handle only quad4 in crack-fsi\n");

  return checkElementExist( dis, elenod, eleid );
}

/*-----------------------------------------------------------------------------------------------------------*
 * Check whether an element exists in the discretization eith all nodal ids                         sudhakar 05/14
 * Assumes that 1. both elements are ordered in same direction
 * 2. given discretization is redundant on all processors
 *-----------------------------------------------------------------------------------------------------------*/
bool ADAPTER::FSICrackingStructure::checkElementExist( Teuchos::RCP<DRT::Discretization>& dis,
                                                       const int* elenod,
                                                       int& eleid )
{
  const Epetra_Map* elecolmap = dis->ElementColMap();

  for (int i=0; i<elecolmap->NumMyElements(); ++i)
  {
    const int gid = elecolmap->GID(i);
    const DRT::Element* sourele = dis->gElement( gid );

    const int* sournod = sourele->NodeIds();

    if( not sourele->NumNode() == 4 )
      dserror("we can handle only quad4 in crack-fsi\n");

    eleid = sourele->Id();

    if( sournod[0] == elenod[0] and sournod[1] == elenod[1] and sournod[2] == elenod[2] and sournod[3] == elenod[3] )
      return true;
    if( sournod[0] == elenod[1] and sournod[1] == elenod[2] and sournod[2] == elenod[3] and sournod[3] == elenod[0] )
      return true;
    if( sournod[0] == elenod[2] and sournod[1] == elenod[3] and sournod[2] == elenod[0] and sournod[3] == elenod[1] )
      return true;
    if( sournod[0] == elenod[3] and sournod[1] == elenod[0] and sournod[2] == elenod[1] and sournod[3] == elenod[2] )
      return true;
  }
  return false;
}
