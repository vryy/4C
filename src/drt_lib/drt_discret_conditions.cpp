/*!----------------------------------------------------------------------
\file drt_discret_conditions.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"

#include "drt_utils.H"

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
  multimap<string,RefCountPtr<DRT::Condition> >::iterator fool;
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    fool->second->ClearGeometry();
    fool->second->SetComm(comm_);
  }

  // for all conditions, we set a ptr in the nodes to the condition
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    const vector<int>* nodes = fool->second->Get<vector<int> >("Node Ids");
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
  }

  return;
}

#if 0 // old version
/*----------------------------------------------------------------------*
 |  Build line geometry in a condition (public)              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildLinesinCondition(
                                        const string name,
                                        RefCountPtr<DRT::Condition> cond)
{

  // get ptrs to all node ids that have this condition
  const vector<int>* nodeids = cond->Get<vector<int> >("Node Ids");
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  map<int,DRT::Node*> rownodes;
  map<int,DRT::Node*> colnodes;

  for (int i=0; i<ngnode; ++i)
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

  // multimap of lines in our cloud
  multimap<int,RefCountPtr<DRT::Element> > linemap;
  // loop these nodes and build all lines attached to them
  map<int,DRT::Node*>::iterator fool;
  for (fool=rownodes.begin(); fool != rownodes.end(); ++fool)
  {
    // currently looking at actnode
    DRT::Node*     actnode  = fool->second;
    // loop all elements attached to actnode
    DRT::Element** elements = actnode->Elements();
    for (int i=0; i<actnode->NumElement(); ++i)
    {
      // loop all lines of all elements attached to actnode
      const int numlines = elements[i]->NumLine();
      if (!numlines) continue;
      DRT::Element** lines = elements[i]->Lines();
      if (!lines) dserror("Element returned no lines");
      for (int j=0; j<numlines; ++j)
      {
        DRT::Element* actline = lines[j];
        // find lines that are attached to actnode
        const int nnodeperline   = actline->NumNode();
        DRT::Node** nodesperline = actline->Nodes();
        if (!nodesperline) dserror("Line returned no nodes");
        for (int k=0; k<nnodeperline; ++k)
          if (nodesperline[k] == actnode)
          {
            // line is attached to actnode
            // see whether all nodes on the line are in our nodal cloud
            bool allin = true;
            for (int l=0; l<nnodeperline; ++l)
            {
              map<int,DRT::Node*>::iterator test = colnodes.find(nodesperline[l]->Id());
              if (test==colnodes.end())
              {
                allin = false;
                break;
              }
            } // for (int l=0; l<nnodeperline; ++l)
            // if all nodes on line are in our cloud, add line
            if (allin)
            {
              RefCountPtr<DRT::Element> tmp = rcp(actline->Clone());
              linemap.insert(pair<int,RefCountPtr<DRT::Element> >(actnode->Id(),tmp));
            }
            break;
          } // if (nodesperline[k] == actnode)
      } // for (int j=0; j<numlines; ++j)
    } // for (int i=0; i<actnode->NumElement(); ++i)
  } // for (fool=nodes.begin(); fool != nodes.end(); ++fool)

  // linemap contains all lines in our cloud, but it also contains a lot
  // of duplicates for now which need to be detected and deleted
  multimap<int,RefCountPtr<DRT::Element> >::iterator linecurr;
  for (linecurr=linemap.begin(); linecurr!=linemap.end(); ++linecurr)
  {
    // this lines was already deleted
    if (linecurr->second == null) continue;

    // get the lines
    RefCountPtr<DRT::Element> actline = linecurr->second;

    // get all nodal ids on this lineace
    const int  nnode   = actline->NumNode();
    const int* nodeids = actline->NodeIds();

    // loop all lines associated with entries of nodeids
    for(int nid=0;nid<nnode;nid++)
    {
      multimap<int,RefCountPtr<DRT::Element> >::iterator startit =
        linemap.lower_bound(nodeids[nid]);
      multimap<int,RefCountPtr<DRT::Element> >::iterator endit   =
        linemap.upper_bound(nodeids[nid]);

      multimap<int,RefCountPtr<DRT::Element> >::iterator curr;
      for (curr=startit; curr!=endit; ++curr)
      {
        if(curr->second == null   ) continue;
        if(curr         == linecurr) continue;

        const int nn    = curr->second->NumNode();
        if (nn != nnode) continue;
        const int* nids = curr->second->NodeIds();

        // nids must contain same ids as nodeids,
        // where ordering is arbitrary
        bool ident = true;
        for (int i=0; i<nnode; ++i)
        {
          bool foundit = false;
          for (int j=0; j<nnode; ++j)
            if (nodeids[i]==nids[j])
            {
              foundit = true;
              break;
            }
          if (!foundit)
          {
            ident = false;
            break;
          }
        }
        if (ident)
          curr->second = null;
        else
          continue;
      }
    }
  }

  // Build a clean map of the remaining now unique lines
  // and add it to the condition
  int count=0;
  map<int,RefCountPtr<DRT::Element> > finallines;
  for (linecurr=linemap.begin(); linecurr!=linemap.end(); ++linecurr)
  {
    if (linecurr->second == null) continue;
    linecurr->second->SetId(count);
    finallines[count] = linecurr->second;
    ++count;
  }

  // Build a global numbering for these elements
  // the elements are in a column map state but the numbering is unique anyway
  // and does NOT reflect the overlap!
  // This is somehow dirty but works for the moment
  vector<int> snelements(Comm().NumProc());
  vector<int> rnelements(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) snelements[i] = 0;
  snelements[Comm().MyPID()] = finallines.size();
  Comm().SumAll(&snelements[0],&rnelements[0],Comm().NumProc());
  int sum=0;
  for (int i=0; i<Comm().MyPID(); ++i) sum += rnelements[i];
  map<int,RefCountPtr<DRT::Element> > finalfinallines;
  count=0;
  for (linecurr=finallines.begin(); linecurr!=finallines.end(); ++linecurr)
  {
    linecurr->second->SetId(count+sum);
    finalfinallines[count+sum] = linecurr->second;
    ++count;
  }
  finallines.clear();

  cond->AddGeometry(finalfinallines);


  return;
} // DRT::Discretization::BuildLinesinCondition
#endif // DEBUG

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
                             const map< vector<int>, RefCountPtr<DRT::Element> >& elementmap,
                             map< int, RefCountPtr<DRT::Element> >& finalelements )
{
  // First, give own elements a local id and find out
  // which ids we need to get from other processes.

  vector< RefCountPtr<DRT::Element> > ownelements;

  // ghostelementnodes, the vector we are going to communicate.
  // Layout:
  //  [
  //    // nodes of elements to ask process 0, separated by -1:
  //    [ node011,node012,node013, ... , -1, node021, ... , -1, ... ],
  //    // nodes of elements to ask process 1, separated by -1:
  //    [ node111,node112,node113, ... , -1, node121, ... , -1, ... ],
  //    ... // etc.
  //  ]
  vector< vector<int> > ghostelementnodes( comm.NumProc() );
  // corresponding elements objects
  vector< vector< RefCountPtr<DRT::Element> > > ghostelements( comm.NumProc() );

  map< vector<int>, RefCountPtr<DRT::Element> >::const_iterator elemsiter;
  for( elemsiter = elementmap.begin(); elemsiter != elementmap.end(); ++elemsiter )
  {
      const RefCountPtr<DRT::Element> element = elemsiter->second;
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
  vector<int> snelements( comm.NumProc() );
  vector<int> rnelements( comm.NumProc() );
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
  vector< vector<int> > requests;
  DRT::Utils::AllToAllCommunication( comm, ghostelementnodes, requests );

  vector< vector<int> > sendids( comm.NumProc() );

  vector<int>::iterator keybegin;
  for ( int proc = 0; proc < comm.NumProc(); ++proc )
  {
      keybegin = requests[proc].begin();
      for ( ;; ) {
          vector<int>::iterator keyend = find( keybegin, requests[proc].end(), -1 );
          if ( keyend == requests[proc].end() )
              break;
          vector<int> nodes = vector<int>( keybegin, keyend );
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

  DRT::Utils::AllToAllCommunication( comm, sendids, requests );

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
} // AssignGlobalIDs

#if 1 // new version
/*----------------------------------------------------------------------*
 |  Build line geometry in a condition (public)              mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
void DRT::Discretization::BuildLinesinCondition( const string name,
                                                 RefCountPtr<DRT::Condition> cond )
{
  /* First: Create the line objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const vector<int>* nodeids = cond->Get< vector<int> >("Node Ids");
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  map<int,DRT::Node*> rownodes;
  map<int,DRT::Node*> colnodes;

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
  map< vector<int>, RefCountPtr<DRT::Element> > linemap;
  // loop these nodes and build all lines attached to them
  map<int,DRT::Node*>::iterator fool;
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
      DRT::Element** lines = elements[i]->Lines();
      if( !lines ) dserror("Element returned no lines");
      for( int j = 0; j < numlines; ++j )
      {
        DRT::Element* actline = lines[j];
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
              map<int,DRT::Node*>::iterator test = colnodes.find(nodesperline[l]->Id());
              if( test==colnodes.end() )
              {
                allin = false;
                break;
              }
            } // for (int l=0; l<nnodeperline; ++l)
            // if all nodes on line are in our cloud, add line
            if( allin )
            {
              vector<int> nodes( actline->NumNode() );
              transform( actline->Nodes(), actline->Nodes() + actline->NumNode(),
                         nodes.begin(), mem_fun( &DRT::Node::Id ) );
              sort( nodes.begin(), nodes.end() );

              if ( linemap.find( nodes ) == linemap.end() )
              {
                  RefCountPtr<DRT::Element> line = rcp( actline->Clone() );
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
  map< int, RefCountPtr<DRT::Element> > finallines;

  AssignGlobalIDs( Comm(), linemap, finallines );
  cond->AddGeometry( finallines );
} // DRT::Discretization::BuildLinesinCondition
#endif // new version of this method

/*----------------------------------------------------------------------*
 |  Build surface geometry in a condition (public)           mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
void DRT::Discretization::BuildSurfacesinCondition(
                                        const string name,
                                        RefCountPtr<DRT::Condition> cond)
{
  /* First: Create the surface objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const vector<int>* nodeids = cond->Get<vector<int> >("Node Ids");
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  map<int,DRT::Node*> rownodes;
  map<int,DRT::Node*> colnodes;
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
  map< vector<int>, RefCountPtr<DRT::Element> > surfmap;

  // loop these row nodes and build all surfs attached to them
  map<int,DRT::Node*>::iterator fool;
  for (fool=rownodes.begin(); fool != rownodes.end(); ++fool)
  {
    // currently looking at actnode
    DRT::Node*     actnode  = fool->second;
    // loop all elements attached to actnode
    DRT::Element** elements = actnode->Elements();
    for (int i=0; i<actnode->NumElement(); ++i)
    {
      // loop all surfaces of all elements attached to actnode
      const int numsurfs = elements[i]->NumSurface();
      if (!numsurfs) continue;
      DRT::Element** surfs = elements[i]->Surfaces();
      if (!surfs) dserror("Element does not return any surfaces");
      for (int j=0; j<numsurfs; ++j)
      {
        DRT::Element* actsurf = surfs[j];
        // find surfs attached to actnode
        const int nnodepersurf = actsurf->NumNode();
        DRT::Node** nodespersurf = actsurf->Nodes();
        if (!nodespersurf) dserror("Surface returned no nodes");
        for (int k=0; k<nnodepersurf; ++k)
          if (nodespersurf[k]==actnode)
          {
            // surface is attached to actnode
            // see whether all  nodes on the surface are in our cloud
            bool allin = true;
            for (int l=0; l<nnodepersurf; ++l)
            {
              map<int,DRT::Node*>::iterator test = colnodes.find(nodespersurf[l]->Id());
              if (test==colnodes.end())
              {
                allin = false;
                break;
              }
            }
            // if all nodes are in our cloud, add surface
            if (allin)
            {
                vector<int> nodes( actsurf->NumNode() );
                transform( actsurf->Nodes(), actsurf->Nodes() + actsurf->NumNode(),
                           nodes.begin(), mem_fun( &DRT::Node::Id ) );
                sort( nodes.begin(), nodes.end() );

                if ( surfmap.find( nodes ) == surfmap.end() )
                {
                    RefCountPtr<DRT::Element> surf = rcp( actsurf->Clone() );
                    // Set owning process of surface to node with smallest gid.
                    surf->SetOwner( gNode( nodes[0] )->Owner() );
                    surfmap[nodes] = surf;
                }
            }
            break;
          }
      }
    }
  }

  // Surfaces be added to the condition: (line_id) -> (surface).
  map< int, RefCountPtr<DRT::Element> > finalsurfs;

  AssignGlobalIDs( Comm(), surfmap, finalsurfs );
  cond->AddGeometry( finalsurfs );
} // DRT::Discretization::BuildSurfacesinCondition



/*----------------------------------------------------------------------*
 |  Build volume geometry in a condition (public)            mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildVolumesinCondition(
                                        const string name,
                                        RefCountPtr<DRT::Condition> cond)
{

  // get ptrs to all node ids that have this condition
  const vector<int>* nodeids = cond->Get<vector<int> >("Node Ids");
  if (!nodeids) dserror("Cannot find array 'Node Ids' in condition");

  // number of global nodes in this cloud
  const int ngnode = nodeids->size();

  // ptrs to my row/column nodes of those
  map<int,DRT::Node*> rownodes;
  map<int,DRT::Node*> colnodes;
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

  // multimap of volumes in our cloud
  multimap<int,RefCountPtr<DRT::Element> > volmap;

  // loop these nodes and build all volumes attached to them
  map<int,DRT::Node*>::iterator fool;
  for (fool=rownodes.begin(); fool != rownodes.end(); ++fool)
  {
    // currently looking at actnode
    DRT::Node*     actnode  = fool->second;
    // loop all elements attached to actnode
    DRT::Element** elements = actnode->Elements();
    for (int i=0; i<actnode->NumElement(); ++i)
    {
      // loop all volumes of all elements attached to actnode
      const int numvols = elements[i]->NumVolume();
      if (!numvols) continue;
      DRT::Element** volumes = elements[i]->Volumes();
      if (!volumes) dserror("Element returned no volumes");
      for (int j=0; j<numvols; ++j)
      {
        DRT::Element* actvol = volumes[j];
        // find volumes that are attached to actnode
        const int nnodepervol   = actvol->NumNode();
        DRT::Node** nodespervol = actvol->Nodes();
        if (!nodespervol) dserror("Volume returned no nodes");
        for (int k=0; k<nnodepervol; ++k)
          if (nodespervol[k] == actnode)
          {
            // volume is attached to actnode
            // see whether all nodes on the volume are in our nodal cloud
            bool allin = true;
            for (int l=0; l<nnodepervol; ++l)
            {
              map<int,DRT::Node*>::iterator test = colnodes.find(nodespervol[l]->Id());
              if (test==colnodes.end())
              {
                allin = false;
                break;
              }
            } // for (int l=0; l<nnodepervol; ++l)
            // if all nodes on volume are in our cloud, add volume
            if (allin)
            {
              RefCountPtr<DRT::Element> tmp = rcp(actvol->Clone());
              volmap.insert(pair<int,RefCountPtr<DRT::Element> >(actnode->Id(),tmp));
            }
            break;
          } // if (nodespervol[k] == actnode)
      } // for (int j=0; j<numvols; ++j)
    } // for (int i=0; i<actnode->NumElement(); ++i)
  } // for (fool=nodes.begin(); fool != nodes.end(); ++fool)

  // volmap contains a lot of duplicates which need to be detected and deleted
  multimap<int,RefCountPtr<DRT::Element> >::iterator volcurr;
  for (volcurr=volmap.begin(); volcurr!=volmap.end(); ++volcurr)
  {
    // this volume was already deleted
    if (volcurr->second == null) continue;

    // get the volume
    RefCountPtr<DRT::Element> actvol = volcurr->second;

    // get all nodal ids on this volume
    const int  nnode   = actvol->NumNode();
    const int* nodeids = actvol->NodeIds();

    // loop all volumes associated with entries of nodeids

    for(int nid=0;nid<nnode;nid++)
    {
      multimap<int,RefCountPtr<DRT::Element> >::iterator startit =
        volmap.lower_bound(nodeids[nid]);
      multimap<int,RefCountPtr<DRT::Element> >::iterator endit   =
        volmap.upper_bound(nodeids[nid]);

      multimap<int,RefCountPtr<DRT::Element> >::iterator curr;
      for (curr=startit; curr!=endit; ++curr)
      {
        if(curr->second == null   ) continue;
        if(curr         == volcurr) continue;

        const int nn    = curr->second->NumNode();
        if (nn != nnode) continue;
        const int* nids = curr->second->NodeIds();

        // nids must contain same ids as nodeids,
        // where ordering is arbitrary
        bool ident = true;
        for (int i=0; i<nnode; ++i)
        {
          bool foundit = false;
          for (int j=0; j<nnode; ++j)
            if (nodeids[i]==nids[j])
            {
              foundit = true;
              break;
            }
          if (!foundit)
          {
            ident = false;
            break;
          }
        }
        if (ident)
          curr->second = null;
        else
          continue;
      }
    }
  }

  // Build a clean map of the remaining now unique lines
  // and add it to the condition
  int count=0;
  map<int,RefCountPtr<DRT::Element> > finalvols;
  for (volcurr=volmap.begin(); volcurr!=volmap.end(); ++volcurr)
  {
    if (volcurr->second == null) continue;
    volcurr->second->SetId(count);
    finalvols[count] = volcurr->second;
    ++count;
  }

  // Build a global numbering for these elements
  // the elements are in a column map state but the numbering is unique anyway
  // and does NOT reflect the overlap!
  // This is somehow dirty but works for the moment (gee)
  vector<int> snelements(Comm().NumProc());
  vector<int> rnelements(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) snelements[i] = 0;
  snelements[Comm().MyPID()] = finalvols.size();
  Comm().SumAll(&snelements[0],&rnelements[0],Comm().NumProc());
  int sum=0;
  for (int i=0; i<Comm().MyPID(); ++i) sum += rnelements[i];
  map<int,RefCountPtr<DRT::Element> > finalfinalvols;
  count=0;
  for (volcurr=finalvols.begin(); volcurr!=finalvols.end(); ++volcurr)
  {
    volcurr->second->SetId(count+sum);
    finalfinalvols[count+sum] = volcurr->second;
    ++count;
  }
  finalvols.clear();

  cond->AddGeometry(finalfinalvols);

  // Normally, we would need to do a ton of cleanup here, but due to
  // RefCountPtr (and STL) we can forget about everything!
  // Great, isn't it?

  return;
} // DRT::Discretization::BuildVolumesinCondition











#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
