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
  boundary_dis = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "FSICoupling", "boundary", "BELE3", conditions_to_copy);

  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentIndependentDofSet(structdis_,true,Teuchos::null));
  boundary_dis->ReplaceDofSet(newdofset);

  #if 0 // try to add element
    tip_nodes.clear();

    int finalids[4] = { 132, 133, 39, 38 };


    int neweleid = boundary_dis->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3","quad4", neweleid, 0/*boundary_dis2_->gNode(nodeids[0])->Owner()*/ );
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

  // build boundary discretization based on FSI condition
  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundary_dis = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "FSICoupling", "boundary", "BELE3", conditions_to_copy);

  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentIndependentDofSet(structdis_,true,Teuchos::null));
  boundary_dis->ReplaceDofSet(newdofset);

  boundary_dis->FillComplete();

  tip_nodes.clear();

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

