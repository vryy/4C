/*----------------------------------------------------------------------*/
/*!
\file initiateCrack.cpp

\brief Crack initiation procedure in a homogeneous flawless
structure.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#include "initiateCrack.H"
#include "crackUtils.H"

#include "Epetra_Map.h"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils.H"

//------------------------------------------------------------------------------------------------
//
//
//                   ----------------*----------------
//                  /|              /|               /|
//                 / |             / |              / |
//                /  |            /  |             /  |
//               /   |           /   |            /   |
//              /    |          /    |           /    |
//             |---------------|----------------|     |                  * split_nodes
//             |     |         |     |          |     |                  o duplicate nodes
//             |     | ........|.....o..........|...  |                  + cracktip 1
//             |     |.        |                |   . |                  # cracktip2
//             |     +---------------*----------------#
//             |    /|         |    /|          |    /|
//             |   / |         |   / |          |   / |
//             |  /  |         |  /  |          |  /  |
//             | / ..|.........o./...|......... | /   |
//             |/ .  |          /    |         .|/    |
//             +---------------*----------------#     |
//             |     |---------|-----|----------|-----|
//             |     /         |     /          |    /
//             |    /          |    /           |   /
//             |   /           |   /            |  /
//             |  /            |  /             | /
//             | /             | /              |/
//             |---------------|----------------|
//
//  The crack initiation operation has the following steps
//  1. Find maximum stress/strain node --> this is one of the split_nodes
//  2. If the max stress/strain is less than failure value defined, nothing to do --> leave. If
//     it is more then proceed with the following steps
//  3. Find another split_node --> At these nodes, the connectivity will be modified later
//  4. Find crack tip nodes --> At the moment, we look for nodes that are away in x-direction
//  5. Apply nodal releasing --> create duplicate nodes in place of split_nodes, and
//     modify the connectivity locally
//  6. Add split_nodes, duplicate nodes and cracktip nodes as crack nodes
//------------------------------------------------------------------------------------------------
bool DRT::CRACK::InitiateCrack::initiateOperations( int & start_new_node_id )
{
  //TODO: Read from input file
  dist_index_ = 1;

  //-------------
  // STEP 1: Determine the node at which von Mises stress is maximum
  //-------------
  int max_von_node = 0;
  double max_von_stress = 0.0;
  Determine_node_max_von_mises_stress( max_von_node, max_von_stress );

  //-------------
  // STEP 2: Check whether failure condition satisfied
  //-------------
  if( max_von_stress < failure_value_ )
    return false;

  //-------------
  // STEP 3: Determine the split nodes
  //-------------
  Determine_split_nodes( max_von_node );

  //-------------
  // STEP 4: Determine Crack tip nodes
  //-------------
  Determine_crack_tip_nodes();

  //-------------
  // STEP 5: Determine Crack tip nodes
  //-------------
  nodal_release( start_new_node_id );

  //-------------
  // STEP 6: Build crack tip and crack surface nodes
  //-------------
  build_tip_surface_nodes();

  return true;
}

/*------------------------------------------------------------------------------------------------------------------------*
 * Find the node at which max stress/strain criterion is met.                                                     sudhakar 09/14
 * This node is always located either in inner or outer layer given in the input file
 *------------------------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::Determine_node_max_von_mises_stress( int & max_von_node, double & max_von_stress )
{
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> nodal_str = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,6,true));
  DRT::CRACK::UTILS::get_nodal_values_from_gauss_point_val( discret_, gausspts_str_, nodal_str );

  double max_str = 0.0;
  std::map<int,double> node_von_mis_stress;


  const int numnodes = discret_->NumMyRowNodes();

  for (int i=0;i<numnodes;++i)
  {

    int this_node = noderowmap->GID(i);
    //if( std::find( boun_nodes_.begin(), boun_nodes_.end(), this_node ) == boun_nodes_.end() )
    //  continue;

    double sxx = (*((*nodal_str)(0)))[i];
    double syy = (*((*nodal_str)(1)))[i];
    double szz = (*((*nodal_str)(2)))[i];
    double sxy = (*((*nodal_str)(3)))[i];
    double syz = (*((*nodal_str)(4)))[i];
    double sxz = (*((*nodal_str)(5)))[i];

    //calculate von-mises stress
    double von_mis = 0.0;

    if( criterion_ == "stress" )
    {
      von_mis = pow( (sxx-syy), 2 ) + pow( (syy-szz), 2 ) + pow( (szz-sxx), 2 ) +
                       6.0 * ( sxy*sxy + syz*syz + sxz*sxz );
      von_mis = sqrt( von_mis * 0.5 );
    }
    else
    {
      von_mis = fabs( syy );
    }

    //double von_mis = fabs( syy );

    if( von_mis > max_str )
    {
      max_str = von_mis;
      node_von_mis_stress.clear();
      node_von_mis_stress[this_node] = von_mis;
    }
    //else if( fabs(von_mis - max_str) < 1e-12 )
      //node_no.push_back( this_node );
  }

  // gather max stress nodes from all processors
  LINALG::GatherAll( node_von_mis_stress, comm_ );

  max_von_node = -1; // initialize with some negative value
  max_von_stress = 0.0;

  for( std::map<int,double>::iterator it = node_von_mis_stress.begin(); it != node_von_mis_stress.end(); it++ )
  {
    if( it->second > max_von_stress )
    {
      max_von_stress = it->second;
      max_von_node = it->first;
    }
  }

  std::cout<<"maximum criterion node = "<<max_von_node<<" value = "<<max_von_stress<<"\n";

  if( max_von_stress == -1 )
    dserror( "von Mises stress throughout the domain is zero?" );
}

/*------------------------------------------------------------------------------------------------*
 * From inner and outer layer of nodes, and using the node at which
 * stress/strain criterion is maximum, determine the splitnodes; nodes                    sudhakar 09/14
 * at which structure is initially split to introduce crack
 *------------------------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::Determine_split_nodes( int& max_von_node )
{
  //---
  // max_von_node is a split node; decide whether it lies in inner or outer nodes
  //---
  if( std::find( inner_nodes_.begin(), inner_nodes_.end(), max_von_node ) != inner_nodes_.end() )
  {
    split_nodes_[ max_von_node ] = "inner";
  }
  else if( std::find( outer_nodes_.begin(), outer_nodes_.end(), max_von_node ) != outer_nodes_.end() )
  {
    split_nodes_[ max_von_node ] = "outer";
  }
  else
    dserror("nodes of maximum stress/strain not found either in inner or outer layer nodes\n");

  //---
  // After deciding the inner/outer location of max_von_node, the next step is to find another
  // split node: node on the outer/inner layer that is directly connected to max_von_node.
  // Idea : get all elements that are connected to max_von_node; the other split node must have all these
  // elements connected to it
  //---

  // stores id of all elements connected to this node
  std::set<int> connected_ele;
  // stores all node ids of one of the elements connected to this node
  std::set<int> ele_nodes;
  if( discret_->HaveGlobalNode( max_von_node ) )
  {
    DRT::Node* max_node = discret_->gNode( max_von_node );
    if( max_node->Owner() == myrank_ )
    {
      const DRT::Element* const * ElementsPtr = max_node->Elements();
      int NumElement = max_node->NumElement();

      for (int jEle = 0; jEle < NumElement; jEle++)
        connected_ele.insert(ElementsPtr[jEle]->Id());

      // since another node also contains all the elements connected to this node,
      // we store nodes of one element, and another node should be one among these nodes
      const DRT::Element* Element = ElementsPtr[0];
      const int* nodeids = Element->NodeIds();
      int numnode = Element->NumNode();
      for(int i=0; i<numnode; ++i)
        ele_nodes.insert( nodeids[i] );
    }
  }

  LINALG::GatherAll( connected_ele, comm_ );
  LINALG::GatherAll( ele_nodes, comm_ );

  for( std::set<int>::iterator it = ele_nodes.begin(); it != ele_nodes.end(); it++ )
  {
    int eleid = *it;

    if( eleid == max_von_node )
      continue;

    if( discret_->HaveGlobalNode( eleid ) )
    {
      DRT::Node* node = discret_->gNode( eleid );
      if( node->Owner() == myrank_ )
      {
        std::set<int> thisnodeele;
        const DRT::Element* const * ElementsPtr = node->Elements();
        int NumElement = node->NumElement();

        for (int jEle = 0; jEle < NumElement; jEle++)
          thisnodeele.insert(ElementsPtr[jEle]->Id());

        if( connected_ele == thisnodeele )
        {
          if( split_nodes_[ max_von_node ] == "inner" )
          {
            split_nodes_[eleid] = "outer";
            break;
          }
          else if( split_nodes_[ max_von_node ] == "outer" )
          {
            split_nodes_[eleid] = "inner";
            break;
          }
          else
            dserror( "maximum von MIses stress node not found?\n" );
        }
      }
    }
  }

  LINALG::GatherAll( split_nodes_, comm_ );

  if( split_nodes_.size() != 2 )
    dserror("There should be 2 split nodes in any configuration\n");
}

/*----------------------------------------------------------------------------------*
 * Determine nodes that will be assigned as tip nodes                       sudhakar 09/14
 *----------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::Determine_crack_tip_nodes()
{
  //---
  // Go through each surface of the elements
  // Find an attached node that has max x-distance from the split nodes
  // in either direction
  // Two nodes on +ve x-direction --> cracktip1
  // Two nodes on -ve x-direction --> cracktip2
  //---
  for( std::map<int,std::string>::iterator it = split_nodes_.begin(); it != split_nodes_.end(); it++ )
  {
    int splnode_id = it->first;
    int posnode = 0, negnode = 0;

    if( it->second == "inner" )
      get_positive_negative_xnode( splnode_id, inner_nodes_, posnode, negnode );
    else
      get_positive_negative_xnode( splnode_id, outer_nodes_, posnode, negnode );

    std::pair<int,int> posneg = std::make_pair( posnode, negnode );
    split_pos_neg_[ it->first ] = posneg;
  }
}

/*--------------------------------------------------------------------------------------------*
 * Get nodes that are connected to, and in max and min x-distance from                sudhakar 09/14
 * split nodes
 *--------------------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::get_positive_negative_xnode( int & splnodeid,
                                                             std::vector<int>& nodes,
                                                             int& pos_node,
                                                             int& neg_node )
{
  int lpos_node = 0, lneg_node = 0;

  if( discret_->HaveGlobalNode( splnodeid ) )
  {
    DRT::Node * splnode = discret_->gNode( splnodeid );
    if( splnode->Owner() == myrank_ )
    {
      double max_xdist = 0.0, min_xdist = 0.0;
      const double * splcord = splnode->X();

      DRT::Element** ElementsPtr = splnode->Elements();
      int NumElement = splnode->NumElement();

      for (int jEle = 0; jEle < NumElement; jEle++)
      {
        DRT::Element * ele = ElementsPtr[jEle];
        Teuchos::RCP<DRT::Element> surf = get_surface_all_in_nodes( ele, nodes );

        const int * surfnodeids = surf->NodeIds();
        int numnode = surf->NumNode();
        const DRT::Node* const * surfnodes = surf->Nodes();

        int spl_index = 0;
        for (int num = 0; num < numnode; num++)
        {
          if( splnodeid == surfnodeids[num] )
          {
            spl_index = num;
            break;
          }
        }
        const DRT::Node * conn1 = surfnodes[(spl_index+1)%numnode];
        int conn2_index = 0;
        if( spl_index == 0 )
          conn2_index = numnode-1;
        else
          conn2_index = spl_index-1;
        const DRT::Node * conn2 = surfnodes[conn2_index];

        const double * conn1_cord = conn1->X();
        const double * conn2_cord = conn2->X();

        double xdist = conn1_cord[dist_index_] - splcord[dist_index_];
        if( xdist > 0.0 and xdist > max_xdist )
        {
          max_xdist = xdist;
          lpos_node = conn1->Id();
        }
        else if( xdist < 0.0 and xdist < min_xdist )
        {
          min_xdist = xdist;
          lneg_node = conn1->Id();
        }

        xdist = conn2_cord[dist_index_] - splcord[dist_index_];
        if( xdist > 0.0 and xdist > max_xdist )
        {
          max_xdist = xdist;
          lpos_node = conn2->Id();
        }
        else if( xdist < 0.0 and xdist < min_xdist )
        {
          min_xdist = xdist;
          lneg_node = conn2->Id();
        }
      }
    }
  }

  comm_.SumAll( &lpos_node, &pos_node, 1 );
  comm_.SumAll( &lneg_node, &neg_node, 1 );
}

/*------------------------------------------------------------------------------------------------------------*
 * Get a surface of the given element that has all of its nodes in the given node set                 sudhakar 09/14
 *------------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::CRACK::InitiateCrack::get_surface_all_in_nodes( DRT::Element * ele, std::vector<int>& nodes )
{
  Teuchos::RCP<DRT::Element> surele = Teuchos::null;
  const std::vector<Teuchos::RCP<DRT::Element> > allsurf = ele->Surfaces();

  bool foundsurf = true;
  for( unsigned i=0; i<allsurf.size(); i++ )
  {
    foundsurf = true;
    Teuchos::RCP<DRT::Element> thissur = allsurf[i];

    const int * thissurnodes = thissur->NodeIds();
    int numnode = thissur->NumNode();
    for( int nn=0; nn < numnode; nn++ )
    {
      if( std::find( nodes.begin(), nodes.end(), thissurnodes[nn] ) == nodes.end() )
        foundsurf = false;
    }

    if( foundsurf )
    {
      surele = thissur;
      break;
    }
  }

  if( not foundsurf )
    dserror("Surface element not found\n");

  return surele;
}

/*---------------------------------------------------------------------------------------*
 * Perform nodal releasing step
 * At the points of maximum stress/strain, and its corresponding node            sudhakar 09/14
 * in inner/outer layer, a duplicate node is generated and the connectivity
 * is modified to include crack
 *---------------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::nodal_release( int & start_new_node_id )
{
  // plane over which to project 3D object
  // here [0] -- x, [1] -- y and [2] -- z
  int projPlane = getProjectionPlane();

  int ind1=0,ind2=0;
  if( projPlane == 0 )
  {
    ind1 = 1;
    ind2 = 2;
  }
  else if( projPlane == 1 )
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if( projPlane == 2 )
  {
    ind1 = 0;
    ind2 = 1;
  }

  std::vector<double> limitAngles = FindLimitingAngles( ind1, ind2 );

  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
  std::map<int, int> delEle;

  for( std::map<int, std::pair<int,int> >::iterator it = split_pos_neg_.begin(); it != split_pos_neg_.end(); it++ )
  {
    int splid = it->first;

    if( discret_->HaveGlobalNode( splid ) )
    {
      DRT::Node * splnode = discret_->gNode( splid );
      const double * splcord = splnode->X();

      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( start_new_node_id, splnode->X(), splnode->Owner() ) );
      oldnew_[splnode->Id()] = dupnode->Id();

      DRT::Element** ElementsPtr = splnode->Elements();
      int NumElement = splnode->NumElement();

      for (int jEle = 0; jEle < NumElement; jEle++)
      {
        DRT::Element* Element = ElementsPtr[jEle];
        std::vector<double> cencord = DRT::UTILS::ElementCenterRefeCoords(Element);

        double cenAng = atan2( (cencord[ind2]-splcord[ind2]), (cencord[ind1]-splcord[ind1]) );
        DRT::CRACK::UTILS::convertAngleTo_02PI_range( cenAng );

        if( cenAng > limitAngles[0] and cenAng < limitAngles[1] )
        {
          delEle[ Element->Id() ] = 0;
        }
      }
      discret_->AddNode( dupnode );
    }
    start_new_node_id++;
  }
  LINALG::GatherAll( oldnew_, comm_ );
  LINALG::GatherAll( delEle, comm_ );

  if( delEle.size() == 0 )
    dserror("initiated crack!!! but no elements found to attach new nodes. This leads to incorrect discretization\n");

  // modify necessary connectivities
  DRT::CRACK::UTILS::ModifyElementConnectivity( discret_, delEle, oldnew_ );

  // add newly generated nodes to appropriate conditions
  DRT::CRACK::UTILS::AddConditions( discret_, oldnew_ );

  discret_->FillComplete();
}

/*---------------------------------------------------------------------------------------*
 * Since real surfaces are 3D, to make easier calculation of angles
 * we project them into appropriate 2D plane This function decides to            sudhakar 09/14
 * which plane to project 3D surfaces
 *---------------------------------------------------------------------------------------*/
int DRT::CRACK::InitiateCrack::getProjectionPlane()
{
  std::map<int,std::string>::iterator it1 = split_nodes_.begin();
  int spl1 = it1->first;
  it1++;
  int spl2 = it1->first;

  // plane over which to project 3D object
  // here [0] -- x, [1] -- y and [2] -- z
  int lprojPlane = 0, gprojPlane = 0;

  // we consider the distance between the two split nodes in each direction
  // in which direction the dist is maximum, we consider that direction to
  // be the projection plane
  if( discret_->HaveGlobalNode( spl1 ) )
  {
    DRT::Node * splnode1 = discret_->gNode( spl1 );
    if( splnode1->Owner() == myrank_ )
    {
      bool planefound = false;

      const double * splcord1 = splnode1->X();

      if( not discret_->HaveGlobalNode( spl2 ) )
        dserror("owner of splitnode1 must contain splitnode2\n");

      const double * splcord2 = discret_->gNode( spl2 )->X();

      double dist[3];
      for( int dim = 0; dim < 3; dim++ )
        dist[dim] = fabs( splcord1[dim] - splcord2[dim] );

      if( dist[0] >= dist[1] and dist[0] >= dist[2] )
      {
        lprojPlane = 0;
        planefound = true;
      }
      else if( dist[1] >= dist[0] and dist[1] >= dist[2] )
      {
        lprojPlane = 1;
        planefound = true;
      }
      else if( dist[2] >= dist[1] and dist[2] >= dist[0] )
      {
        lprojPlane = 2;
        planefound = true;
      }

      if( not planefound )
        dserror( "projection plane not found\n" );
    }
  }

  comm_.SumAll( &lprojPlane, &gprojPlane, 1 );
  return gprojPlane;
}

/*------------------------------------------------------------------------------------------------------*
 * Build crack surface and crack tip segment nodes
 * There are always two tip segments, one for positive and another for negative x-nodes         sudhakar 09/14
 * Split nodes, their duplicates and the tip segment nodes fall on crack surfaces
 *------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::InitiateCrack::build_tip_surface_nodes()
{
  for( std::map<int,int>::iterator it = oldnew_.begin(); it != oldnew_.end(); it++ )
  {
    crackSurfNodes_.insert( it->first );
    crackSurfNodes_.insert( it->second );
  }

  for( std::map<int, std::pair<int,int> >::iterator it = split_pos_neg_.begin(); it != split_pos_neg_.end(); it++ )
  {
    std::pair<int,int> posneg = it->second;
    posCrackTip_.push_back( posneg.first );
    negCrackTip_.push_back( posneg.second );
  }
}

/*------------------------------------------------------------------------------------------------------*
 * Find limiting angles which will is used in nodal release                                     sudhakar 09/14
 *------------------------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::InitiateCrack::FindLimitingAngles( const int &ind1, const int & ind2 )
{
  std::vector<double> limits;

  std::map<int, std::pair<int,int> >::iterator it = split_pos_neg_.begin();

  int splid = it->first;
  std::pair<int,int> posnegid = it->second;
  int posid = posnegid.first;
  int negid = posnegid.second;

  // temporary storage of angles from different processors
  // "key" --> stores angles, and "value" --> dummy
  std::map<double,int> angles;

  if( discret_->HaveGlobalNode( splid ) )
  {
    DRT::Node * splnode = discret_->gNode( splid );
    const double * splcord = splnode->X();

    if( splnode->Owner() == myrank_ )
    {
      if( ( not discret_->HaveGlobalNode( posid ) ) or ( not discret_->HaveGlobalNode( negid  ) ) )
        dserror("positive or negative node not found on this processor\n");

      const double * poscord = discret_->gNode( posid )->X();
      double posAng = atan2( (poscord[ind2]-splcord[ind2]), (poscord[ind1]-splcord[ind1]) );
      DRT::CRACK::UTILS::convertAngleTo_02PI_range( posAng );

      const double * negcord = discret_->gNode( negid )->X();
      double negAng = atan2( (negcord[ind2]-splcord[ind2]), (negcord[ind1]-splcord[ind1]) );
      DRT::CRACK::UTILS::convertAngleTo_02PI_range( negAng );

      angles[posAng] = 0;
      angles[negAng] = 0;
    }
  }

  LINALG::GatherAll( angles, comm_ );

  for( std::map<double,int>::iterator it = angles.begin(); it != angles.end(); it++ )
  {
    double ang = it->first;
    limits.push_back( ang );
  }
  std::sort( limits.begin(), limits.end() );

  return limits;
}
