/*!----------------------------------------------------------------------
\file drt_discret_conditions.cpp
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "drt_parobject.H"

#include "drt_condition_utils.H"
#include "../linalg/linalg_utils.H"

#include <numeric>
#include <algorithm>


/*----------------------------------------------------------------------*
 |  Build boundary condition geometries (public)             mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BoundaryConditionsGeometry()
{
  // As a first step we delete ALL references to any conditions
  // in the discretization
  for (int i=0; i<NumMyColNodes(); ++i)
    lColNode(i)->ClearConditions();
  for (int i=0; i<NumMyColElements(); ++i)
    lColElement(i)->ClearConditions();

  // now we delete all old geometries that are attached to any conditions
  // and set a communicator to the condition
  std::multimap<std::string,RCP<DRT::Condition> >::iterator fool;
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    fool->second->ClearGeometry();
    fool->second->SetComm(comm_);
  }

  // for all conditions, we set a ptr in the nodes to the condition
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    const std::vector<int>* nodes = fool->second->Nodes();
    // There might be conditions that do not have a nodal cloud
    if (!nodes) continue;
    int nnode = nodes->size();
    for (int i=0; i<nnode; ++i)
    {
      if (!NodeColMap()->MyGID((*nodes)[i])) continue;
      DRT::Node* actnode = gNode((*nodes)[i]);
      if (!actnode) dserror("Cannot find global node");
      actnode->SetCondition(fool->first,fool->second);
    }
  }

  // create a map that holds the overall number of created elements
  // associated with a specific condition type
  std::map<std::string, int> numele;

  // Loop all conditions and build geometry description if desired
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    // do not build geometry description for this condition
    if (fool->second->GeometryDescription()==false)         continue;
    // do not build geometry description for this condition
    else if (fool->second->GType()==DRT::Condition::NoGeom) continue;
    // do not build anything for point wise conditions
    else if (fool->second->GType()==DRT::Condition::Point)  continue;
    // build a line element geometry description
    else if (fool->second->GType()==DRT::Condition::Line)
      BuildLinesinCondition(fool->first,fool->second);
    // build a surface element geometry description
    else if (fool->second->GType()==DRT::Condition::Surface)
      BuildSurfacesinCondition(fool->first,fool->second);
    // build a volume element geometry description
    else if (fool->second->GType()==DRT::Condition::Volume)
      BuildVolumesinCondition(fool->first,fool->second);

    if (fool->second->GType()!=DRT::Condition::Volume)
    {
      // determine the local number of created elements associated with
      // the active condition
      int localcount=0;
      for (std::map<int,RCP<DRT::Element> >::iterator iter=fool->second->Geometry().begin();
           iter!=fool->second->Geometry().end();
           ++iter)
      {
        // do not count ghosted elements
        if (iter->second->Owner()==Comm().MyPID())
        {
          localcount += 1;
        }
      }

      // determine the global number of created elements associated with
      // the active condition
      int count;
      Comm().SumAll(&localcount, &count, 1);

      if (numele.find(fool->first)==numele.end())
      {
        numele[fool->first] = 0;
      }

      // adjust the IDs of the elements associated with the active
      // condition in order to obtain unique IDs within one condition type
      fool->second->AdjustId(numele[fool->first]);

      // adjust the number of elements associated with the current
      // condition type
      numele[fool->first]+=count;
    }
  }
  return;
}


/*
 *  A helper function for BuildLinesinCondition and
 *  BuildSurfacesinCondition, below.
 *  Gets a map (vector_of_nodes)->Element that maps
 *
 *  (A map with globally unique ids.)
 *
 *  \param comm (i) communicator
 *  \param elementmap (i) map (vector_of_nodes_ids)->(element) that maps
 *  the nodes of an element to the element itself.
 *
 *  \param finalelements (o) map (global_id)->(element) that can be
 *  added to a condition.
 *
 *  h.kue 09/07
 */
static void AssignGlobalIDs( const Epetra_Comm& comm,
                             const std::map< std::vector<int>, RCP<DRT::Element> >& elementmap,
                             std::map< int, RCP<DRT::Element> >& finalelements )
{
#if 0
  // First, give own elements a local id and find out
  // which ids we need to get from other processes.

  std::vector< RCP<DRT::Element> > ownelements;

  // ghostelementnodes, the vector we are going to communicate.
  // Layout:
  //  [
  //    // nodes of elements to ask process 0, separated by -1:
  //    [ node011,node012,node013, ... , -1, node021, ... , -1, ... ],
  //    // nodes of elements to ask process 1, separated by -1:
  //    [ node111,node112,node113, ... , -1, node121, ... , -1, ... ],
  //    ... // etc.
  //  ]
  std::vector< std::vector<int> > ghostelementnodes( comm.NumProc() );
  // corresponding elements objects
  std::vector< std::vector< RCP<DRT::Element> > > ghostelements( comm.NumProc() );

  std::map< std::vector<int>, RCP<DRT::Element> >::const_iterator elemsiter;
  for( elemsiter = elementmap.begin(); elemsiter != elementmap.end(); ++elemsiter )
  {
      const RCP<DRT::Element> element = elemsiter->second;
      if ( element->Owner() == comm.MyPID() )
      {
          ownelements.push_back( element );
      }
      else // This is not our element, but we know it. We'll ask its owner for the id, later.
      {
          copy( elemsiter->first.begin(), elemsiter->first.end(),
                back_inserter( ghostelementnodes[element->Owner()] ) );
          ghostelementnodes[element->Owner()].push_back( -1 );

          ghostelements[element->Owner()].push_back( element );
      }
  }

  // Find out which ids our own elements are supposed to get.
  std::vector<int> snelements( comm.NumProc() );
  std::vector<int> rnelements( comm.NumProc() );
  fill( snelements.begin(), snelements.end(), 0 );
  snelements[ comm.MyPID() ] = ownelements.size();
  comm.SumAll( &snelements[0], &rnelements[0], comm.NumProc() );
  int sum = accumulate( &rnelements[0], &rnelements[comm.MyPID()], 0 );

  // Add own elements to finalelements (with right id).
  for ( unsigned i = 0; i < ownelements.size(); ++i )
  {
      ownelements[i]->SetId( i + sum );
      finalelements[i + sum] = ownelements[i];
  }
  ownelements.clear();

  // Last step: Get missing ids.
  std::vector< std::vector<int> > requests;
  LINALG::AllToAllCommunication( comm, ghostelementnodes, requests );

  std::vector< std::vector<int> > sendids( comm.NumProc() );

  std::vector<int>::iterator keybegin;
  for ( int proc = 0; proc < comm.NumProc(); ++proc )
  {
      keybegin = requests[proc].begin();
      for ( ;; ) {
          std::vector<int>::iterator keyend = find( keybegin, requests[proc].end(), -1 );
          if ( keyend == requests[proc].end() )
              break;
          std::vector<int> nodes = std::vector<int>( keybegin, keyend );
          elemsiter = elementmap.find( nodes );
          if ( elemsiter == elementmap.end() )
              dserror( "Got request for unknown element" );
          sendids[proc].push_back( elemsiter->second->Id() );

          ++keyend;
          keybegin = keyend;
      }
  }
  requests.clear();

#if 0 // Debug
  cout << "This is process " << comm.MyPID() << "." << endl;
  for ( int proc = 0; proc < comm.NumProc(); ++proc )
  {
      cout << "Send to process " << proc << ": ";
      for ( unsigned i = 0; i < sendids[proc].size(); ++i )
      {
          cout << sendids[proc][i] << ", ";
      }
      cout << endl;
  }
#endif // Debug

  LINALG::AllToAllCommunication( comm, sendids, requests );

#if 0 // Debug
  cout << "This is process " << comm.MyPID() << "." << endl;
  for ( int proc = 0; proc < comm.NumProc(); ++proc )
  {
      cout << "Got from process " << proc << ": ";
      for ( unsigned i = 0; i < requests[proc].size(); ++i )
      {
          cout << requests[proc][i] << ", ";
      }
      cout << endl;
  }
#endif // Debug

  for ( int proc = 0; proc < comm.NumProc(); ++proc )
  {
      if ( requests[proc].size() != ghostelements[proc].size() )
          dserror( "Wrong number of element ids from proc %i: expected %i, got %i",
                   proc, ghostelements[proc].size(), requests[proc].size() );

      for ( unsigned i = 0; i < ghostelements[proc].size(); ++i )
      {
          if ( finalelements.find( requests[proc][i] ) != finalelements.end() )
              dserror( "Received already known id %i", requests[proc][i] );

          ghostelements[proc][i]->SetId( requests[proc][i] );
          finalelements[ requests[proc][i] ] = ghostelements[proc][i];
      }
  }

#else

  // The point here is to make sure the element gid are the same on any
  // parallel distribution of the elements. Thus we allreduce thing to
  // processor 0 and sort the element descriptions (vectors of nodal ids)
  // there.
  //
  // This routine has not been optimized for efficiency. I don't think that is
  // needed.
  //
  // pack elements on all processors

  int size = 0;
  std::map<std::vector<int>, Teuchos::RCP<DRT::Element> >::const_iterator elemsiter;
  for (elemsiter=elementmap.begin();
       elemsiter!=elementmap.end();
       ++elemsiter)
  {
    size += elemsiter->first.size()+1;
  }
  std::vector<int> sendblock;
  sendblock.reserve(size);
  for (elemsiter=elementmap.begin();
       elemsiter!=elementmap.end();
       ++elemsiter)
  {
    sendblock.push_back(elemsiter->first.size());
    std::copy(elemsiter->first.begin(), elemsiter->first.end(), std::back_inserter(sendblock));
  }

  // communicate elements to processor 0

  int mysize = sendblock.size();
  comm.SumAll(&mysize,&size,1);
  int mypos = LINALG::FindMyPos(sendblock.size(),comm);

  std::vector<int> send(size);
  std::fill(send.begin(),send.end(),0);
  std::copy(sendblock.begin(),sendblock.end(),&send[mypos]);
  sendblock.clear();
  std::vector<int> recv(size);
  comm.SumAll(&send[0],&recv[0],size);

  send.clear();

  // unpack, unify and sort elements on processor 0

  if (comm.MyPID()==0)
  {
    std::set<std::vector<int> > elements;
    int index = 0;
    while (index < static_cast<int>(recv.size()))
    {
      int esize = recv[index];
      index += 1;
      std::vector<int> element;
      element.reserve(esize);
      std::copy(&recv[index], &recv[index+esize], std::back_inserter(element));
      index += esize;
      elements.insert(element);
    }
    recv.clear();

    // pack again to distribute pack to all processors

    send.reserve(index);
    for (std::set<std::vector<int> >::iterator i=elements.begin();
         i!=elements.end();
         ++i)
    {
      send.push_back(i->size());
      std::copy(i->begin(), i->end(), std::back_inserter(send));
    }
    size = send.size();
  }
  else
  {
    recv.clear();
  }

  // broadcast sorted elements to all processors

  comm.Broadcast(&size,1,0);
  send.resize(size);
  comm.Broadcast(&send[0],send.size(),0);

  // Unpack sorted elements. Take element position for gid.

  int index = 0;
  int gid = 0;
  while (index < static_cast<int>(send.size()))
  {
    int esize = send[index];
    index += 1;
    std::vector<int> element;
    element.reserve(esize);
    std::copy(&send[index], &send[index+esize], std::back_inserter(element));
    index += esize;

    // set gid to my elements
    std::map<std::vector<int>, RCP<DRT::Element> >::const_iterator iter = elementmap.find(element);
    if (iter!=elementmap.end())
    {
      iter->second->SetId(gid);
      finalelements[gid] = iter->second;
    }

    gid += 1;
  }

#endif
} // AssignGlobalIDs


/*----------------------------------------------------------------------*
 |  Build line geometry in a condition (public)              mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
void DRT::Discretization::BuildLinesinCondition( const std::string name,
                                                 RCP<DRT::Condition> cond )
{
  /* First: Create the line objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->Nodes();
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  std::map<int,DRT::Node*> rownodes;
  std::map<int,DRT::Node*> colnodes;

  for( int i=0; i<ngnode; ++i )
  {
    if (NodeColMap()->MyGID((*nodeids)[i]))
    {
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node");
      colnodes[actnode->Id()] = actnode;
    }
  }
  for (int i=0; i<ngnode; ++i)
  {
    if (NodeRowMap()->MyGID((*nodeids)[i]))
    {
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node");
      rownodes[actnode->Id()] = actnode;
    }
  }

  // map of lines in our cloud: (node_ids) -> line
  std::map< std::vector<int>, RCP<DRT::Element> > linemap;
  // loop these nodes and build all lines attached to them
  std::map<int,DRT::Node*>::iterator fool;
  for( fool = rownodes.begin(); fool != rownodes.end(); ++fool )
  {
    // currently looking at actnode
    DRT::Node*     actnode  = fool->second;
    // loop all elements attached to actnode
    DRT::Element** elements = actnode->Elements();
    for( int i = 0; i < actnode->NumElement(); ++i )
    {
      // loop all lines of all elements attached to actnode
      const int numlines = elements[i]->NumLine();
      if( !numlines ) continue;
      std::vector<Teuchos::RCP<DRT::Element> >  lines = elements[i]->Lines();
      if(lines.size()==0) dserror("Element returned no lines");
      for( int j = 0; j < numlines; ++j )
      {
        RCP<DRT::Element> actline = lines[j];
        // find lines that are attached to actnode
        const int nnodeperline   = actline->NumNode();
        DRT::Node** nodesperline = actline->Nodes();
        if( !nodesperline ) dserror("Line returned no nodes");
        for( int k = 0; k < nnodeperline; ++k )
          if( nodesperline[k] == actnode )
          {
            // line is attached to actnode
            // see whether all nodes on the line are in our nodal cloud
            bool allin = true;
            for( int l=0; l < nnodeperline; ++l )
            {
              std::map<int,DRT::Node*>::iterator test = colnodes.find(nodesperline[l]->Id());
              if( test==colnodes.end() )
              {
                allin = false;
                break;
              }
            } // for (int l=0; l<nnodeperline; ++l)
            // if all nodes on line are in our cloud, add line
            if( allin )
            {
              std::vector<int> nodes( actline->NumNode() );
              transform( actline->Nodes(), actline->Nodes() + actline->NumNode(),
                         nodes.begin(), std::mem_fun( &DRT::Node::Id ) );
              sort( nodes.begin(), nodes.end() );

              if ( linemap.find( nodes ) == linemap.end() )
              {
                  RCP<DRT::Element> line = Teuchos::rcp( actline->Clone() );
                  // Set owning process of line to node with smallest gid.
                  line->SetOwner( gNode( nodes[0] )->Owner() );
                  linemap[nodes] = line;
              }
            }
            break;
          } // if (nodesperline[k] == actnode)
      } // for (int j=0; j<numlines; ++j)
    } // for (int i=0; i<actnode->NumElement(); ++i)
  } // for (fool=nodes.begin(); fool != nodes.end(); ++fool)


  // Lines be added to the condition: (line_id) -> (line).
  std::map< int, RCP<DRT::Element> > finallines;

  AssignGlobalIDs( Comm(), linemap, finallines );
  cond->AddGeometry( finallines );
} // DRT::Discretization::BuildLinesinCondition


/*----------------------------------------------------------------------*
 |  Build surface geometry in a condition (public)           mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
void DRT::Discretization::BuildSurfacesinCondition(
                                        const std::string name,
                                        RCP<DRT::Condition> cond)
{
  // these conditions are special since associated volume conditions also need
  // to be considered
  std::set<int> VolEleIDs;
  if (cond->Type() == DRT::Condition::StructFluidSurfCoupling)
  {
    if (*(cond->Get<std::string>("field")) == "structure")
    {
      FindAssociatedEleIDs(cond, VolEleIDs, "StructFluidVolCoupling");
    }
  }
  else if (cond->Type() == DRT::Condition::RedAirwayTissue)
    FindAssociatedEleIDs(cond, VolEleIDs, "StructFluidVolCoupling");

  /* First: Create the surface objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->Nodes();
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  std::map<int,DRT::Node*> rownodes;
  std::map<int,DRT::Node*> colnodes;
  for (int i=0; i<ngnode; ++i)
  {
    if (NodeColMap()->MyGID((*nodeids)[i]))
    {
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node");
      colnodes[actnode->Id()] = actnode;
    }
    if (NodeRowMap()->MyGID((*nodeids)[i]))
    {
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node");
      rownodes[actnode->Id()] = actnode;
    }
  }

  // map of surfaces in this cloud: (node_ids) -> (surface)
  std::map< std::vector<int>, RCP<DRT::Element> > surfmap;

  // loop these row nodes and build all surfs attached to them
  std::map<int,DRT::Node*>::iterator fool;
  for (fool=rownodes.begin(); fool != rownodes.end(); ++fool)
  {
    // currently looking at actnode
    DRT::Node*     actnode  = fool->second;
    // loop all elements attached to actnode
    DRT::Element** elements = actnode->Elements();
    bool foundvolele = false;
    for (int i=0; i<actnode->NumElement(); ++i)
    {
      if (VolEleIDs.size())
        if (VolEleIDs.find(elements[i]->Id()) == VolEleIDs.end())
          continue;

      // loop all surfaces of all elements attached to actnode
      const int numsurfs = elements[i]->NumSurface();
      if (!numsurfs) continue;
      std::vector<Teuchos::RCP<DRT::Element> >  surfs = elements[i]->Surfaces();
      if (surfs.size()==0) dserror("Element does not return any surfaces");
      for (int j=0; j<numsurfs; ++j)
      {
        RCP<DRT::Element> actsurf = surfs[j];
        // find surfs attached to actnode
        const int nnodepersurf = actsurf->NumNode();
        DRT::Node** nodespersurf = actsurf->Nodes();
        if (!nodespersurf) dserror("Surface returned no nodes");
        for (int k=0; k<nnodepersurf; ++k)
        {
          if (nodespersurf[k]==actnode)
          {
            // surface is attached to actnode
            // see whether all  nodes on the surface are in our cloud
            bool allin = true;
            for (int l=0; l<nnodepersurf; ++l)
            {
              std::map<int,DRT::Node*>::iterator test = colnodes.find(nodespersurf[l]->Id());
              if (test==colnodes.end())
              {
                allin = false;
                break;
              }
            }
            // if all nodes are in our cloud, add surface
            if (allin)
            {
              //If condition type is c volume coupling condition, remove internal surfaces
              if ((cond->Type() == DRT::Condition::StructFluidSurfCoupling) or (cond->Type() == DRT::Condition::RedAirwayTissue))
              //if (true)
              {
                std::cout << "Volume Coupling Condition ID: " << cond->Id()+1 << std::endl;
                //cout << DRT::Problem::Instance()->ProblemType() << endl;

                // remove internal surfaces that are connected to two volume elements
                DRT::Node** actsurfnodes = actsurf->Nodes();
                const int actsurfnumnode = actsurf->NumNode();

                // to be tested
                std::vector<int> nodes( actsurf->NumNode() );
                transform( actsurf->Nodes(), actsurf->Nodes() + actsurf->NumNode(),
                           nodes.begin(), std::mem_fun( &DRT::Node::Id ) );
                sort( nodes.begin(), nodes.end() );

                // get all volume elements connected to the nodes of this surface element
                std::set<DRT::Element*> adjacentvoleles;
                for(int n=0; n<actsurfnumnode; ++n)
                {
                  DRT::Node* actsurfnode = actsurfnodes[n];

                  //std::cout << "Current active Surface: " << actsurfnodes[0]->Id() <<" " << actsurfnodes[1]->Id() << " " << actsurfnodes[2]->Id() << std::endl;

                  DRT::Element** eles = actsurfnode->Elements();
                  int numeles = actsurfnode->NumElement();
                  for(int e=0; e<numeles; ++e)
                  {
                    // do not consider volume elements that do not belong to VolEleIDs
                    if (VolEleIDs.size())
                      if (VolEleIDs.find(eles[e]->Id()) == VolEleIDs.end())
                        continue;
                    adjacentvoleles.insert(eles[e]);
                  }

                  //std::cout << "Number of adjacentvoleles: " << adjacentvoleles.size() << std::endl;
                  //std::cout << "IDs: ";
                  //for(std::set<DRT::Element*>::const_iterator iter1=adjacentvoleles.begin(); iter1!=adjacentvoleles.end(); ++iter1)
                  //  std::cout << "   " << (*iter1)->Id() << std::endl;

                }

                int identical = 0;
                // get surfaces of all adjacent vol eles and check how often actsurf is included via comparison of node ids
                for(std::set<DRT::Element*>::const_iterator iter=adjacentvoleles.begin(); iter!=adjacentvoleles.end(); ++iter)
                {
                  std::vector<Teuchos::RCP<DRT::Element> >  adjacentvolelesurfs = (*iter)->Surfaces();
                  const int adjacentvolelenumsurfs = (*iter)->NumSurface();
                  for(int n=0; n<adjacentvolelenumsurfs; ++n)
                  {
                    // current surf of adjacent vol ele
                    std::vector<int> nodesadj( adjacentvolelesurfs[n]->NumNode() );
                    transform( adjacentvolelesurfs[n]->Nodes(), adjacentvolelesurfs[n]->Nodes() + adjacentvolelesurfs[n]->NumNode(),
                               nodesadj.begin(), std::mem_fun( &DRT::Node::Id ) );
                    sort( nodesadj.begin(), nodesadj.end() );

                    if (nodes.size() == nodesadj.size())
                    {
                      if ( std::equal (nodes.begin(), nodes.end(), nodesadj.begin()) )
                        identical++;
                    }
                  }
                }
                if(identical == 0)
                  dserror("surface found with missing underlying volume element");
                else if(identical > 1)
                {
                  /*const int* actsurfnodeids = actsurf->NodeIds();
                  std::cout << "A surface element which was included unintentionally was removed (node ids: ";
                  for(int n=0; n<actsurfnumnode; ++n)
                    std::cout << actsurfnodeids[n] << " ";
                  std::cout << " )" << std::endl;*/
                }
                else
                {
                  // now we can add the surface
                  if ( surfmap.find( nodes ) == surfmap.end() )
                  {
                    Teuchos::RCP<DRT::Element> surf = Teuchos::rcp( actsurf->Clone() );
                    // Set owning process of surface to node with smallest gid.
                    surf->SetOwner( gNode( nodes[0] )->Owner() );
                    surfmap[nodes] = surf;
                  }
                  foundvolele = true;
                }

              } //end of special volume condition checking
              else
              {
                std::cout << "No Volume Coupling Condition: No internal boundary Faces check ID: " << cond->Id()+1 << std::endl;

                std::vector<int> nodes( actsurf->NumNode() );
                transform( actsurf->Nodes(), actsurf->Nodes() + actsurf->NumNode(),
                           nodes.begin(), std::mem_fun( &DRT::Node::Id ) );
                sort( nodes.begin(), nodes.end() );

                if ( surfmap.find( nodes ) == surfmap.end() )
                {
                  RCP<DRT::Element> surf = Teuchos::rcp( actsurf->Clone() );
                  // Set owning process of surface to node with smallest gid.
                  surf->SetOwner( gNode( nodes[0] )->Owner() );
                  surfmap[nodes] = surf;
                }

                foundvolele = true;
              }
            }
            break;
          }
        }
      }
    }
    if (VolEleIDs.size() and foundvolele == false)
      dserror("special surface condition: missing associated volume element");
  }

  //Write output for Gmsh format for debugging
  //if(true)
  if ((cond->Type() == DRT::Condition::StructFluidSurfCoupling) or (cond->Type() == DRT::Condition::RedAirwayTissue))
  {
    std::cout << "Surfaces in Surfmap for Coupling Condition No. " << cond->Id()+1 << ":" << std::endl;
    std::cout << "Format: [Node1, Node2, Node3, CondID] " << std::endl;
    for(std::map< std::vector<int>, RCP<DRT::Element> >::const_iterator iterat=surfmap.begin(); iterat!=surfmap.end(); ++iterat)
    {
      cout << iterat->first[0] << " " << iterat->first[1] << " " << iterat->first[2] << " " << cond->Id()+1 << endl;
    }
    std::cout << "End " << std::endl;
  }

  // Surfaces be added to the condition: (line_id) -> (surface).
  std::map< int, RCP<DRT::Element> > finalsurfs;

  AssignGlobalIDs( Comm(), surfmap, finalsurfs );
  cond->AddGeometry( finalsurfs );
} // DRT::Discretization::BuildSurfacesinCondition


/*----------------------------------------------------------------------*
 |  Build volume geometry in a condition (public)            mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildVolumesinCondition(
                                        const std::string name,
                                        RCP<DRT::Condition> cond)
{
  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->Nodes();
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // extract colnodes on this proc from condition
  const Epetra_Map* colmap = NodeColMap();
  std::set<int> mynodes;

  std::remove_copy_if(nodeids->begin(), nodeids->end(),
                      std::inserter(mynodes, mynodes.begin()),
                      std::not1(DRT::UTILS::MyGID(colmap)));

  // this is the map we want to construct
  std::map<int,RCP<DRT::Element> > geom;

  for (std::map<int,RCP<DRT::Element> >::iterator actele=element_.begin();
       actele!=element_.end();
       ++actele)
  {
    std::vector<int> myelenodes(actele->second->NodeIds(),actele->second->NodeIds()+actele->second->NumNode());

    // check whether all node ids of the element are nodes belonging
    // to the condition and stored on this proc
    bool allin=true;
    for(std::vector<int>::iterator myid=myelenodes.begin();myid!=myelenodes.end();++myid)
    {
      if(mynodes.find(*myid)==mynodes.end())
      {
        // myid is not in the condition
        allin=false;
        break;
      }
    }

    if(allin)
    {
      geom[actele->first] = actele->second;
    }
  }

  cond->AddGeometry(geom);

  return;
} // DRT::Discretization::BuildVolumesinCondition



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::FindAssociatedEleIDs(RCP<DRT::Condition> cond, std::set<int>& VolEleIDs, const std::string name)
{
  // determine constraint number
  int condID = cond->GetInt("coupling id");

  std::vector<DRT::Condition*> volcond;
  GetCondition(name,volcond);

  for (unsigned int i = 0; i < volcond.size(); ++i)
  {
    DRT::Condition& actvolcond = *(volcond[i]);

    if (actvolcond.GetInt("coupling id") == condID)
    {
      // get ptrs to all node ids that have this condition
      const std::vector<int>* nodeids = actvolcond.Nodes();
      if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

      // extract colnodes on this proc from condition
      const Epetra_Map* colmap = NodeColMap();
      std::set<int> mynodes;

      std::remove_copy_if(nodeids->begin(), nodeids->end(),
                          std::inserter(mynodes, mynodes.begin()),
                          std::not1(DRT::UTILS::MyGID(colmap)));

      for (std::map<int,RCP<DRT::Element> >::iterator actele=element_.begin();
           actele!=element_.end();
           ++actele)
      {
        std::vector<int> myelenodes(actele->second->NodeIds(),actele->second->NodeIds()+actele->second->NumNode());

        // check whether all node ids of the element are nodes belonging
        // to the condition and stored on this proc
        bool allin=true;
        for(std::vector<int>::iterator myid=myelenodes.begin();myid!=myelenodes.end();++myid)
        {
          if(mynodes.find(*myid)==mynodes.end())
          {
            // myid is not in the condition
            allin=false;
            break;
          }
        }

        if(allin)
        {
          VolEleIDs.insert(actele->second->Id());
        }
      }
    }
  }
}

