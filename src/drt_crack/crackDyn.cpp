/*----------------------------------------------------------------------*/
/*!
\file PropagateCrack.cpp

\brief After each time step, check the structure field and propagate
crack in the structure if necessary.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "crack_tolerance.H"
#include "crackUtils.H"
#include "aleCrack.H"
#include "SplitHexIntoTwoWedges.H"
#include "initiateCrack.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_mat/elasthyper.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_matelast/elast_coupneohooke.H"

#include "../linalg/linalg_utils.H"


/*------------------------------------------------------------------------------------*
 * Constructor                                                                     sudhakar 12/13
 * Read crack tip nodes and nodes falling on the crack surface from input data
 *------------------------------------------------------------------------------------*/
DRT::CRACK::CrackDyn::CrackDyn( Teuchos::RCP<DRT::Discretization>& discret )
:discret_(discret),
 comm_( discret_->Comm() ),
 clearCondns_( false ),
 myrank_( comm_.MyPID() )
{
  strIsSplit_ = false;
  crackInitiated_ = false;
  isRestart_ = false;

  //---
  // Read crack tip nodes from the input file
  // If these nodes exist, then crack is already initiated in the structure
  //---
  std::vector<std::vector<int> > all_crack_tips;
  std::vector<Teuchos::RCP<Condition> > crackpts;
  discret_->GetCondition( "CrackInitiationPoints", crackpts );
  for( std::vector<Teuchos::RCP<Condition> >::iterator conit = crackpts.begin(); conit != crackpts.end(); conit++ )
  {
    std::vector<int> this_tip;
    Teuchos::RCP<Condition> tip_pts = *conit;
    if( tip_pts != Teuchos::null )
    {
      const std::vector<int>* tipnodes = const_cast<std::vector<int>* >(tip_pts->Nodes());

      for( unsigned i=0; i< tipnodes->size(); i++ )
      {
        int val = (*tipnodes)[i];
        this_tip.push_back(val);
      }
      all_crack_tips.push_back(this_tip);
      crackInitiated_ = true;
    }
  }

  //---
  // If crack is initiated master and slave crack surface nodes exist in the input fil
  // Read them and store for further calculations
  //---
  if( crackInitiated_ )
  {
    DRT::Condition* maspts = discret_->GetCondition( "masterCrackSurface" );
    DRT::Condition* slapts = discret_->GetCondition( "slaveCrackSurface" );

    const std::vector<int>* masternodes = const_cast<std::vector<int>* >(maspts->Nodes());
    const std::vector<int>* slavenodes = const_cast<std::vector<int>* >(slapts->Nodes());

    //if( masternodes->size() != slavenodes->size() )
    //  dserror("There should be equal number of master and slave nodes\n");

    if( masternodes->size() == 0 )
      dserror("No master nodes defined. Are you dreaming of simulating crack initiation?\n");

    for( unsigned i=0; i<masternodes->size(); i++ )
    {
      cracknodes_.insert( (*masternodes)[i] );
      cracknodes_.insert( (*slavenodes)[i] );
    }
  }

  //---
  // construct propagatetip vector
  // Each handles propagation of the associated crack tip segment
  //---
  for( unsigned tipno = 0; tipno < all_crack_tips.size(); tipno++ )
  {
    //PropagateTip pt( discret, all_crack_tips[tipno], cracknodes_, tipno );
    //tip_segments_.push_back( &pt );

    Teuchos::RCP<PropagateTip> pt = Teuchos::rcp(new PropagateTip( discret, all_crack_tips[tipno], cracknodes_, tipno ));
    tip_segments_.push_back( pt );
  }

  //---
  // Read necessary parameters
  //---
  const Teuchos::ParameterList& crackparam = DRT::Problem::Instance()->CrackParams();
  startNewNodeId_ = crackparam.get<int>("START_NEW_NODE_ID");
  startNewEleId_ = crackparam.get<int>("START_NEW_ELE_ID");

  if( startNewNodeId_ == 0 )
    dserror("Set correct value of START_NEW_NODE_ID in the input file\n");

  if( startNewEleId_ == 0 )
    dserror("Set correct value of START_NEW_ELE_ID in the input file\n");

  //---
  // Store ALE line condition nodes
  //---
  std::vector<DRT::Condition*> linecond;
  discret_->GetCondition( "ALEDirichlet", linecond );

  if( linecond.size() == 0 )
    dserror("No Dirichlet condition found for this discretization\n");

  for(unsigned i=0;i<linecond.size();i++)
  {
    DRT::Condition* con = linecond[i];
    if( con->Type() == DRT::Condition::LineDirichlet )
    {
      const std::vector<int>* nodes = con->Nodes();
      for( unsigned nodno = 0; nodno < nodes->size(); nodno++ )
        ALE_line_nodes_.insert( (*nodes)[nodno] );
    }
  }

  //---
  //Store inner and outer layer elements for crack propagation
  //---
  if( not crackInitiated_ )
  {
    DRT::Condition* innpts = discret_->GetCondition( "CrackInnerLayerPoints" );
    DRT::Condition* outpts = discret_->GetCondition( "CrackOuterLayerPoints" );

    const std::vector<int>* inner_nodes = const_cast<std::vector<int>* >(innpts->Nodes());
    const std::vector<int>* outer_nodes = const_cast<std::vector<int>* >(outpts->Nodes());

    if( inner_nodes->size() == 0 )
      dserror( "inner layer of nodes should be specified for crack initiation\n" );
    if( outer_nodes->size() == 0 )
      dserror( "outer layer of nodes should be specified for crack initiation\n" );

    for( unsigned i=0; i<inner_nodes->size(); i++ )
      inner_layer_nodes_.push_back( (*inner_nodes)[i] );

    for( unsigned i=0; i<outer_nodes->size(); i++ )
      outer_layer_nodes_.push_back( (*outer_nodes)[i] );
  }

  crackInitiatedThisStep_ = false;

  //this->WriteNodesDebug();
}


/*------------------------------------------------------------------------------------*
 * Perform all the operations related to crack propagation                  sudhakar 11/13
 * This calculates stress intensity factor (K), and if it is  higher than
 * the critical value, crack is introduced into discretization
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::propagateOperations( Teuchos::RCP<const Epetra_Vector>& displace,
                                                Teuchos::RCP<std::vector<char> >& stressdata,
                                                Teuchos::RCP<std::vector<char> >& straindata )
{
  // general geometric paramters
  oldnew_.clear();
  tipnodes_.clear();

  // parameters related to clearing all conditions
  //TODO:: if( not restart at this step )
  justClearedCondns_ = false;
  if( clearCondns_ )
  {
    //DeleteConditions();
    justClearedCondns_ = true;
    return;
  }

  // Crack has already propagated that the structure is completely split into two
  if( strIsSplit_ )
    return;

  if( crackInitiated_ )
  {
    crackInitiatedThisStep_ = false;
  }

  //-------------
  // STEP 1 : Initiate crack if not already existing
  //-------------
  InitiateCrack( stressdata, straindata );

  // crack has not yet initiated. nothing more to do here
  if( not crackInitiated_ )
    return;

  // crack just initated. not propagate at this step
  if( crackInitiatedThisStep_ )
    return;

  //------------
  // STEP 2: Propagate each tip segments
  //------------
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    tip_segments_[segno]->set_stress_values( stressdata );
    tip_segments_[segno]->set_strain_values( straindata );
    tip_segments_[segno]->propagateThisTip( displace );
  }

  //--------------
  // STEP 3: Perform ALE step of moving nodes
  //--------------
  // it is performed here since each tip segments know only about their data
  // Here we collect all the  BC from each segments and perform the ALE operations
  std::map<int, std::vector<double> > tipbc_disp;
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    if( not tip_segments_[segno]->isCrackPropagatedThisStep() )
      return;
    std::map<int, std::vector<double> > temp_bc_disp;
    temp_bc_disp = tip_segments_[segno]->getALEtipBC();
    tipbc_disp.insert(temp_bc_disp.begin(), temp_bc_disp.end());
  }
  Perform_ALE_Step( displace, tipbc_disp );

  std::cout<<"---------combined ALE boundary conditions----------------\n";
  for(std::map<int,std::vector<double> >::iterator it = tipbc_disp.begin(); it != tipbc_disp.end(); it++)
    std::cout<<"node id = "<<it->first<<" ALE disp = "<<(it->second)[0]<<" "<<(it->second)[1]<<" "<<(it->second)[2]<<"\n";
  std::cout<<"---------------------------------------------------------------------------\n";

  //-----------------------
  // STEP 4: Perform nodal releasing process
  //-----------------------
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    if( not tip_segments_[segno]->isCrackPropagatedThisStep() )
      return;
    tip_segments_[segno]->updateCrack( startNewEleId_, startNewNodeId_ );
  }

  //-----------------------
  // STEP 5: Perform complete splitting of structure if needed
  //-----------------------
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    if( not tip_segments_[segno]->isCrackPropagatedThisStep() )
      return;
    tip_segments_[segno]->CheckCompleteSplit( startNewNodeId_ );
  }

  //-----------------------
  //STEP 6: Update crack tip nodes
  //-----------------------
  update_crack_tip_nodes();

  //-----------------------
  // STEP 7: Update ALE boundarycondition nodes
  //-----------------------
  updateALE_BC_nodes();

  //-----------------------
  // STEP 8: Update boundary nodes
  //-----------------------
  update_boundary_nodes();

  //-------------
  // STEP 9: Analyze the discretization, and delete unnecessary nodes in each processor
  //-------------
  analyzeAndCleanDiscretization();

  //---------------
  // STEP 10: Check if all tip segments completely break the body into two
  //---------------
  strIsSplit_ = true;
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    if ( not tip_segments_[segno]->isSplit() )
    {
      strIsSplit_ = false;
      break;
    }
  }

  step_++;

  isRestart_ = false;
}

/*------------------------------------------------------------------------------------*
 * If not already present, initiate crack into the structure                sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::InitiateCrack( Teuchos::RCP<std::vector<char> >& stressdata,
                                          Teuchos::RCP<std::vector<char> >& straindata )
{

  if( crackInitiated_ )
    return;

  if( stressdata == Teuchos::null )
  {
    dserror( "stress data is empty!" );
  }

  if( straindata == Teuchos::null )
  {
    dserror( "strain data is empty!" );
  }

  Teuchos::RCP<DRT::CRACK::InitiateCrack> ini = Teuchos::null;

  if( 1 ) // strain cases
  {
    ini = Teuchos::rcp(new DRT::CRACK::InitiateCrack(discret_, straindata, "strain" ));
  }
  else // stress cases
  {
    ini = Teuchos::rcp(new DRT::CRACK::InitiateCrack(discret_, stressdata, "stress" ));
  }

  //TODO: read failure value from input file

  ini->set_failure_value( 0.045/*100*/ );
  ini->set_inner_layer_nodes( inner_layer_nodes_ );
  ini->set_outer_layer_nodes( outer_layer_nodes_ );

  if( ini->initiateOperations( startNewNodeId_ ) )
  {
    crackInitiated_ = true;
    crackInitiatedThisStep_ = true;
    analyzeAndCleanDiscretization();
    oldnew_ = ini->getOldNewNodes();

    Teuchos::RCP<PropagateTip> ptpos = Teuchos::rcp(new PropagateTip( discret_, ini->getPositiveCrackTipNodes(), ini->getCrackSurfaceNodes(), 0 ));
    Teuchos::RCP<PropagateTip> ptneg = Teuchos::rcp(new PropagateTip( discret_, ini->getNegativeCrackTipNodes(), ini->getCrackSurfaceNodes(), 1 ));

    tip_segments_.push_back( ptpos );
    tip_segments_.push_back( ptneg );
  }
}

/*------------------------------------------------------------------------------------*
 * Perform all operations related to ALE step                                 sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::Perform_ALE_Step(  Teuchos::RCP<const Epetra_Vector>& displace,
                                              std::map<int, std::vector<double> >& tipbc_disp  )
{
  if( tipbc_disp.size() == 0 )
    return;

  if( 0 )
  {
    for( std::map<int,std::vector<double> >::iterator it = tipbc_disp.begin(); it != tipbc_disp.end(); it++ )
    {
      int tipid = it->first;
      if( discret_->HaveGlobalNode( tipid ) )
      {
        DRT::Node * tipnode = discret_->gNode( tipid );
        tipnode->ChangePos( it->second );
      }
    }
  }
  else
  {
    aleCrack ale( discret_ );
    ale.ALE_step( tipbc_disp, new_ale_bc_nodes_ );


    ALE_line_nodes_ = ale.getLineDirichNodes();
    if( ( step_ == 1 ) and ALE_line_nodes_.size() == 0 )
      dserror( "No Line Dirichlet found on ALE discretization\n" );

    ale.clearALE_discret();
  }
}

/*---------------------------------------------------------------------------------------------*
 * Update nodes those holds ALE boundary condition after crack propagation             sudhakar 09/14
 *---------------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::updateALE_BC_nodes()
{
  // get oldtip nodes from all segments
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    std::map<int,int> o_n = tip_segments_[segno]->getOldNewTipNodes();
    oldnew_.insert( o_n.begin(), o_n.end() );
  }

  // update all data structure
  for( std::map<int,int>::iterator it = oldnew_.begin(); it != oldnew_.end(); it++ )
  {
    cracknodes_.insert( it->first );
    cracknodes_.insert( it->second );
    new_ale_bc_nodes_.insert( it->second );
    ALE_line_nodes_.insert( it->second );
  }
  for( std::vector<int>::iterator it = tipnodes_.begin(); it != tipnodes_.end(); it++ )
  {
    new_ale_bc_nodes_.insert( *it );
    ALE_line_nodes_.insert( *it );
  }

  // after updating all crack nodes, let all crack segments know about the updated crack nodes
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    tip_segments_[segno]->setCrackNodes( cracknodes_ );
  }

}

/*-----------------------------------------------------------------------*
 * Update boundary nodes after crack propagation                 sudhakar 09/14
 *-----------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::update_boundary_nodes()
{
  // update boundary nodes of all crack segments
  // This is mandatory because each crack segment does not know about other
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    for( std::map<int,int>::iterator it = oldnew_.begin(); it != oldnew_.end(); it++ )
    {
      tip_segments_[segno]->addBoundaryNode( it->first );
      tip_segments_[segno]->addBoundaryNode( it->second );
    }

    for( std::vector<int>::iterator it = tipnodes_.begin(); it != tipnodes_.end(); it++ )
      tip_segments_[segno]->addBoundaryNode( *it );
  }
}

/*--------------------------------------------------------------------------------------------------*
 * After propagating the crack tips, store all new tip nodes                                sudhakar 09/14
 *--------------------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::update_crack_tip_nodes()
{
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    std::vector<int>segtip = tip_segments_[segno]->getSegmentTipNodes();
    tipnodes_.insert(tipnodes_.end(), segtip.begin(), segtip.end());
  }
}

/*------------------------------------------------------------------------------------*
 * Analyze --> Make sure all nodes are available on each proc for considered
 *             element distribution                                             sudhakar 06/14
 * Clean   --> Delete all extra nodes on each proc that are not needed for
 *             given element distribution
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::analyzeAndCleanDiscretization()
{
  std::set<int> colnodesReqd;  // all nodes that are required in this processor for available element distribution
  std::set<int> colnodesAvail; // all nodes that are availble in this processor for available element distribution

  // get all required nodes on this proc
  int numelements = discret_->NumMyColElements();
  for (int i=0; i<numelements; i++ )
  {
    DRT::Element* actele = discret_->lColElement(i);
    const int* nodes = actele->NodeIds();
    for( int nno = 0; nno < actele->NumNode(); nno++ )
      colnodesReqd.insert( nodes[nno] );
  }

  // get all available nodes on this proc
  int numnodes = discret_->NumMyColNodes();
  const Epetra_Map * colnodemap = discret_->NodeColMap();
  for( int i=0; i<numnodes; i++ )
  {
    int nid = colnodemap->GID( i );
    colnodesAvail.insert( nid );
  }

  std::vector<int> delNodeIds;
  if( colnodesReqd == colnodesAvail )
  {
  }
  else
  {
    std::set<int>::iterator itavail, itreqd;

    // check whether all node ids are available on this proc
    std::vector<int> notavail;
    for( itreqd = colnodesReqd.begin(); itreqd != colnodesReqd.end(); itreqd++ )
    {
      int nid = *itreqd;
      itavail = colnodesAvail.find( nid );
      if( itavail == colnodesAvail.end() )
        notavail.push_back( nid );
    }

    if( not (notavail.size() == 0) )
    {
      std::cout<<"The following nodeids are not available in proc "<<myrank_<<" : ";
      for( unsigned i=0; i< notavail.size(); i++ )
        std::cout<<notavail[i]<<" ";
      std::cout<<"\n";
      dserror("failed!\n");
    }

    // node ids that is not required but available in this proc
    // Hence these will  be erased

    for( itavail = colnodesAvail.begin(); itavail != colnodesAvail.end(); itavail++ )
    {
      int nid = *itavail;
      itreqd = colnodesReqd.find( nid );
      if( itreqd == colnodesReqd.end() )
        delNodeIds.push_back( nid );
    }
  }

  for( std::vector<int>::iterator it = delNodeIds.begin(); it != delNodeIds.end(); it++ )
  {
    discret_->DeleteNode( *it );
  }

  discret_->FillComplete();
}

/*-------------------------------------------------------------------------------------------------*
 * Write discretization with nodal Ids for debugging reasons                               sudhakar 06/14
 *-------------------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::WriteNodesDebug()
{
  std::ofstream f;
  std::string nodfil = "nodeids";

  f.open( nodfil.c_str());
  f<<"View \"Element geo\" {\n";

  for( int iele=0; iele<discret_->ElementRowMap()->NumMyElements(); ++iele )
  {
    DRT::Element * curele = discret_->lRowElement( iele );
    f<<"SH(";
    for(int j=0; j<curele->NumNode(); ++j)
    {
      DRT::Node* node = curele->Nodes()[j];
      const double * co = node->X();
      f<<co[0]<<","<<co[1]<<","<<co[2];
      if( j != curele->NumNode()-1 )
        f<<",";
    }
    f<<"){";
    for(int j=0; j<curele->NumNode(); ++j)
    {
      f<<curele->NodeIds()[j];
      if( j != curele->NumNode()-1 )
        f<<",";
    }
    f<<"};\n";
  }
  f<<"};\n";



  f<<"View \"Node Ids\" {\n";
  for( int inod = 0; inod < discret_->NodeRowMap()->NumMyElements(); inod++ )
  {
    DRT::Node * curnod = discret_->lRowNode( inod );
    f<<"SP(";
    const double * co=curnod->X();
    f<<co[0]<<","<<co[1]<<","<<co[2]<<"){"<<curnod->Id()<<"};\n";
  }
  f<<"};\n";

  f.close();

  dserror("done");
}

/*------------------------------------------------------------------------------------*
 * Print all conditions of the discretization (just for debugging)           sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::printConditions(std::multimap<std::string,Teuchos::RCP<Condition> > allcondn)
{
  std::cout<<"number of conditions = "<<allcondn.size()<<"\n";
  for( std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit = allcondn.begin();
                                                                       conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;
    std::cout<<"Id = "<<cond->Id()<<" condn type = "<<cond->Type()<<" geom disc = "<<cond->GeometryDescription()<<"\n";
    const std::vector<int>* nodeids = cond->Nodes();

    std::cout<<"condition nodes are = ";
    for( unsigned i=0;i<nodeids->size();i++ )
    {
      std::cout<<(*nodeids)[i]<<"  ";
    }
    std::cout<<"\n";
  }
}

/*---------------------------------------------------------------------------------------------------*
 * Write all necessary data for restarting the simulation                                    sudhakar 12/14
 *---------------------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::WriteRestartCrack( Teuchos::RCP<IO::DiscretizationWriter> outputWriter )
{
  if( myrank_ == 0 )
    std::cout<<"~~~ Writing restart data for crack ~~~\n";

  // all nodes that fall on the crack surface
  Teuchos::RCP<std::vector<int> > dummy_crack_nodes = Teuchos::rcp( new std::vector<int>() );
  for( std::set<int>::iterator it = cracknodes_.begin(); it != cracknodes_.end(); it++ )
    dummy_crack_nodes->push_back( *it );
  outputWriter->WriteRedundantIntVector( "cracksurfacenodes", dummy_crack_nodes );

  // all nodes that hold ALE line BC
  Teuchos::RCP<std::vector<int> > dummy_ale_line_nodes = Teuchos::rcp( new std::vector<int>() );
  for( std::set<int>::iterator it = ALE_line_nodes_.begin(); it != ALE_line_nodes_.end(); it++ )
    dummy_ale_line_nodes->push_back( *it );
  outputWriter->WriteRedundantIntVector( "alelinenodes", dummy_ale_line_nodes );

  // new Nodes that are added to ALE boundary conditions
  Teuchos::RCP<std::vector<int> > dummy_ale_new_nodes = Teuchos::rcp( new std::vector<int>() );
  for( std::set<int>::iterator it = new_ale_bc_nodes_.begin(); it != new_ale_bc_nodes_.end(); it++ )
    dummy_ale_new_nodes->push_back( *it );
  outputWriter->WriteRedundantIntVector( "newalebcnodes", dummy_ale_new_nodes );

  // current crack tip nodes
  Teuchos::RCP<std::vector<int> > dummy_tip_ids = Teuchos::rcp( new std::vector<int>() );
  int num_nodes_per_tip = 0;
  for( unsigned segno = 0 ; segno < tip_segments_.size(); segno++ )
  {
    std::vector<int>segtip = tip_segments_[segno]->getSegmentTipNodes();
    num_nodes_per_tip = segtip.size();
    dummy_tip_ids->insert(dummy_tip_ids->end(), segtip.begin(), segtip.end());
  }
  outputWriter->WriteInt( "numnodepertip", num_nodes_per_tip );
  outputWriter->WriteRedundantIntVector( "alltipsegmentnodes", dummy_tip_ids );

  // ids of new node and ele ids to be added to discret
  outputWriter->WriteInt( "startnewnodeid", startNewNodeId_ );
  outputWriter->WriteInt( "startneweleid", startNewEleId_ );

  // Write boundary nodes
  std::set<int> bounNodes = tip_segments_[0]->getAllBoundaryNodes();
  Teuchos::RCP<std::vector<int> > dummy_bounNodes = Teuchos::rcp( new std::vector<int>() );
  for( std::set<int>::iterator it = bounNodes.begin(); it != bounNodes.end(); it++ )
    dummy_bounNodes->push_back( *it );
  outputWriter->WriteRedundantIntVector( "boundaryNodes", dummy_bounNodes );

  // Other boolean variables
  int boolclearcondn = 0;
  if( clearCondns_ )
    boolclearcondn = 1;

  int boolstrsplit = 0;
  if( strIsSplit_ )
    boolstrsplit = 1;

  int justclearedcondn = 0;
  if( justClearedCondns_ )
    justclearedcondn = 1;

  outputWriter->WriteInt( "boolclearcondn", boolclearcondn );
  outputWriter->WriteInt( "boolstrsplit", boolstrsplit );
  outputWriter->WriteInt( "justclearedcondn", justclearedcondn );
}

/*---------------------------------------------------------------------------------------------------*
 * Read restart data and properly initialize all variables                                   sudhakar 12/14
 *---------------------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::ReadRestartCrack( IO::DiscretizationReader& reader )
{
  if( myrank_ == 0 )
    std::cout<<"~~~ Reading restart data for crack ~~~\n";

  //---
  // Read crack surface nodes
  //--
  Teuchos::RCP<std::vector<int> > dummy_crack_nodes = Teuchos::rcp( new std::vector<int>() );
  reader.ReadRedundantIntVector( dummy_crack_nodes, "cracksurfacenodes" );
  cracknodes_.clear();
  for( unsigned i=0; i< dummy_crack_nodes->size(); i++ )
    cracknodes_.insert( (*dummy_crack_nodes)[i] );

  //---
  // Read all ALE line nodes
  //---
  Teuchos::RCP<std::vector<int> > dummy_ale_line_nodes = Teuchos::rcp( new std::vector<int>() );
  reader.ReadRedundantIntVector( dummy_ale_line_nodes, "alelinenodes" );
  ALE_line_nodes_.clear();
  for( unsigned i=0; i<dummy_ale_line_nodes->size(); i++ )
    ALE_line_nodes_.insert( (*dummy_ale_line_nodes)[i] );

  //---
  // Read all new nodes that are added to discretization
  //---
  Teuchos::RCP<std::vector<int> > dummy_ale_new_nodes = Teuchos::rcp( new std::vector<int>() );
  reader.ReadRedundantIntVector( dummy_ale_new_nodes, "newalebcnodes" );
  new_ale_bc_nodes_.clear();
  for( unsigned i=0; i<dummy_ale_new_nodes->size(); i++ )
    new_ale_bc_nodes_.insert( (*dummy_ale_new_nodes)[i] );

  //---
  // Read boundary nodes
  //---
  Teuchos::RCP<std::vector<int> > dummy_bounNodes = Teuchos::rcp( new std::vector<int>() );
  reader.ReadRedundantIntVector( dummy_bounNodes, "boundaryNodes" );
  std::set<int> bounNodes( dummy_bounNodes->begin(), dummy_bounNodes->end() );

  //---
  // Reading starting node and element ids
  //---
  startNewNodeId_ = reader.ReadInt( "startnewnodeid" );
  startNewEleId_ = reader.ReadInt( "startneweleid" );

  //---
  // Set crack tip segments
  //---
  int num_nodes_per_tip = reader.ReadInt( "numnodepertip" );
  Teuchos::RCP<std::vector<int> > allcracktips = Teuchos::rcp( new std::vector<int>() );
  reader.ReadRedundantIntVector( allcracktips, "alltipsegmentnodes" );

  int totaltips = (int)allcracktips->size();
  if( totaltips % num_nodes_per_tip != 0 )
    dserror( "Different tip segments have different number of tip nodes? Not supported currently!" );
  if( totaltips != 0 )
    crackInitiated_ = true;

  int no_segments = totaltips/num_nodes_per_tip;
  tip_segments_.clear();
  for( int tipno=0; tipno<no_segments; tipno++ )
  {
    std::vector<int> segment_tip;
    int begin = tipno*num_nodes_per_tip;
    for( int jj=begin; jj< begin+num_nodes_per_tip; jj++ )
      segment_tip.push_back( (*allcracktips)[jj] );

    Teuchos::RCP<PropagateTip> pt = Teuchos::rcp(new PropagateTip( discret_, segment_tip, cracknodes_, tipno ));
    pt->InitializeBoundaryNodes( bounNodes );
    tip_segments_.push_back( pt );
  }

  //---
  // Read other boolean parameters
  //---
  int boolclearcondn = reader.ReadInt( "boolclearcondn" );
  if( boolclearcondn == 0 )
    clearCondns_ = false;
  else
    clearCondns_ = true;

  int boolstrsplit = reader.ReadInt( "boolstrsplit" );
  if( boolstrsplit == 0 )
    strIsSplit_ = false;
  else
    strIsSplit_ = true;

  int justclearedcondn = reader.ReadInt( "justclearedcondn" );
  if( justclearedcondn == 0 )
    justClearedCondns_ = false;
  else
    justClearedCondns_ = true;

  isRestart_ = true;

  std::vector<DRT::Condition*> allcond;
  discret_->GetCondition( "Dirichlet", allcond );

  if( allcond.size() == 0 )
    dserror("No Dirichlet condition found for this discretization\n");

  for(unsigned i=0;i<allcond.size();i++)
  {
    DRT::Condition* con = allcond[i];
    if( con->Type() == DRT::Condition::VolumeDirichlet )
    {
      DRT::CRACK::UTILS::addNodesToConditions( con, cracknodes_ );
    }
  }

  if( DRT::Problem::Instance()->ProblemType() == prb_fsi_crack )
  {
    DRT::Condition* cond_fsi = discret_->GetCondition("FSICoupling");
    DRT::Condition* cond_xfem = discret_->GetCondition("XFEMCoupling");

    if( not cond_fsi->Nodes()->size() == dummy_ale_new_nodes->size() )
    {
      if( cond_fsi == NULL or cond_xfem == NULL )
        dserror( "XFEM or FSI coupling conditions undefined in XFSI problem?\n" );

      DRT::CRACK::UTILS::addNodesToConditions( cond_fsi, *dummy_ale_new_nodes );
      DRT::CRACK::UTILS::addNodesToConditions( cond_xfem, *dummy_ale_new_nodes );
    }
  }

  discret_->FillComplete();
}

#if 0

/*------------------------------------------------------------------------------------*
 * For the given element, get the surface that lies in the z-plane of "coord      sudhakar 12/13
 *------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::CRACK::CrackDyn::getSurfaceSameZplane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
                                                                             const double * coord )
{
  Teuchos::RCP<DRT::Element> surf = Teuchos::null;
  bool foundsurf = true;

  // which surfaces has all the nodes in considered z-plane
  for (unsigned int surfele = 0; surfele < Surfaces.size(); surfele++)
  {
    foundsurf = true;
    surf = Surfaces[surfele];

    DRT::Node ** searchnodes = surf->Nodes();
#ifdef DEBUG
    if (searchnodes==NULL) dserror("No Nodes in Surface Element");
#endif
    for (int num = 0; num < surf->NumNode(); num++)
    {
      DRT::Node * nod = searchnodes[num];
      const double * atcord = nod->X();

      if( fabs( atcord[2] - coord[2] ) > 1e-12 )
      {
        foundsurf = false;
        break;
      }
    }

    if( foundsurf )
      break;
  }

  if( not foundsurf )
    dserror( "the required surface is not found\n" );

  if( surf == Teuchos::null )
    dserror( "the required surface is not found\n" );

  return surf;
}

/*------------------------------------------------------------------------------------*
 * Delete the Dirichlet conditions existing at the previous crack tip nodes      sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::CrackDyn::DeleteConditions()
{
  if( not clearCondns_ )
    return;

  DRT::CRACK::UTILS::deleteConditions( discret_ );

  /*std::vector<std::multimap<std::string,Teuchos::RCP<Condition> >::iterator> del;

  std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit;
  std::multimap<std::string,Teuchos::RCP<Condition> >& allcondn = discret_->GetAllConditions();
  for( conit = allcondn.begin(); conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;
    if( cond->Nodes()->size() == 1 )
      del.push_back( conit );
  }

  for( unsigned i=0; i< del.size(); i++ )
  {
    conit = del[i];
    allcondn.erase( conit );
  }*/

  // already we have cleared it now. unless new condns are set, no need to delete anyting
  clearCondns_ = false;
}
#endif
