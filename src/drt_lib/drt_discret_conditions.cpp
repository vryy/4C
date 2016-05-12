/*!----------------------------------------------------------------------
\file drt_discret_conditions.cpp
\brief Implementation

<pre>
\brief Implementation
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "drt_parobject.H"

#include "drt_globalproblem.H"
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
  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator fool;
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
    // flag whether new elements have been created for this condition
    bool havenewelements = false;

    // do not build geometry description for this condition
    if (fool->second->GeometryDescription()==false)         continue;
    // do not build geometry description for this condition
    else if (fool->second->GType()==DRT::Condition::NoGeom) continue;
    // do not build anything for point wise conditions
    else if (fool->second->GType()==DRT::Condition::Point)  continue;
    // build element geometry description without creating new elements if the
    // condition is not a boundary condition; this would be:
    //  - line conditions in 1D
    //  - surface conditions in 2D
    //  - volume conditions in 3D
    else if ((int)(fool->second->GType())==DRT::Problem::Instance()->NDim())
      havenewelements = BuildVolumesinCondition(fool->first,fool->second);
    // dimension of condition must not larger than the one of the problem itself
    else if ((int)(fool->second->GType())>DRT::Problem::Instance()->NDim())
      dserror("Dimension of condition is larger than the problem dimension.");
    // build a line element geometry description
    else if (fool->second->GType()==DRT::Condition::Line)
      havenewelements = BuildLinesinCondition(fool->first,fool->second);
    // build a surface element geometry description
    else if (fool->second->GType()==DRT::Condition::Surface)
      havenewelements = BuildSurfacesinCondition(fool->first,fool->second);
    // this should be it. if not: dserror.
    else
      dserror("Somehow the condition geometry does not fit to the problem dimension.");

    if (havenewelements)
    {
      // determine the local number of created elements associated with
      // the active condition
      int localcount=0;
      for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator iter=fool->second->Geometry().begin();
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
void DRT::Discretization::AssignGlobalIDs(
  const Epetra_Comm& comm,
  const std::map< std::vector<int>, Teuchos::RCP<DRT::Element> >& elementmap,
  std::map< int, Teuchos::RCP<DRT::Element> >& finalelements )
{
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
    std::map<std::vector<int>, Teuchos::RCP<DRT::Element> >::const_iterator iter = elementmap.find(element);
    if (iter!=elementmap.end())
    {
      iter->second->SetId(gid);
      finalelements[gid] = iter->second;
    }

    gid += 1;
  }
} // AssignGlobalIDs


/*----------------------------------------------------------------------*
 |  Build line geometry in a condition (public)              mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
bool DRT::Discretization::BuildLinesinCondition( const std::string name,
                                                 Teuchos::RCP<DRT::Condition> cond )
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
  std::map< std::vector<int>, Teuchos::RCP<DRT::Element> > linemap;
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
        Teuchos::RCP<DRT::Element> actline = lines[j];
        // find lines that are attached to actnode
        const int nnodeperline   = actline->NumNode();
        DRT::Node** nodesperline = actline->Nodes();
        if( !nodesperline ) dserror("Line returned no nodes");
        for( int k = 0; k < nnodeperline; ++k )
          if( nodesperline[k]->Id() == actnode->Id())
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
                  Teuchos::RCP<DRT::Element> line = Teuchos::rcp( actline->Clone() );
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
  Teuchos::RCP<std::map< int, Teuchos::RCP<DRT::Element> > > finallines = Teuchos::rcp(new std::map<int,Teuchos::RCP<DRT::Element> >() );

  AssignGlobalIDs( Comm(), linemap, *finallines );

  cond->AddGeometry( finallines );

  // elements where created that need new unique ids
  bool havenewelements = true;
  // note: this seems useless since this BuildLinesinCondition always creates
  //       elements. However, this function is overloaded in
  //       MeshfreeDiscretization where it does not necessarily build elements.
  return havenewelements;

} // DRT::Discretization::BuildLinesinCondition


/*----------------------------------------------------------------------*
 |  Build surface geometry in a condition (public)           mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
bool DRT::Discretization::BuildSurfacesinCondition(
                                        const std::string name,
                                        Teuchos::RCP<DRT::Condition> cond)
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

  //Just for information - remove later
  if ((cond->Type() == DRT::Condition::StructFluidSurfCoupling) or (cond->Type() == DRT::Condition::RedAirwayTissue))
  {
    if (cond->comm_->MyPID()==0)
    {
      std::cout << "Volume Coupling Condition ID: " << cond->Id()+1 << " found and checked..."<< std::endl;
    }
  }

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
  std::map< std::vector<int>, Teuchos::RCP<DRT::Element> > surfmap;

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
        Teuchos::RCP<DRT::Element> actsurf = surfs[j];
        // find surfs attached to actnode
        const int nnodepersurf = actsurf->NumNode();
        DRT::Node** nodespersurf = actsurf->Nodes();
        if (!nodespersurf) dserror("Surface returned no nodes");
        for (int k=0; k<nnodepersurf; ++k)
        {
          if (nodespersurf[k]->Id()==actnode->Id())
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
              //If condition type is a volume coupling condition, remove internal surfaces
              if ((cond->Type() == DRT::Condition::StructFluidSurfCoupling) or (cond->Type() == DRT::Condition::RedAirwayTissue))
              {
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
                std::vector<int> nodes( actsurf->NumNode() );
                transform( actsurf->Nodes(), actsurf->Nodes() + actsurf->NumNode(),
                           nodes.begin(), std::mem_fun( &DRT::Node::Id ) );
                sort( nodes.begin(), nodes.end() );

                if ( surfmap.find( nodes ) == surfmap.end() )
                {
                  Teuchos::RCP<DRT::Element> surf = Teuchos::rcp( actsurf->Clone() );
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
    {
      //std::cout << " Warning: special surface condition: missing associated volume element" << std::endl;
    }
  }

  //Write output for Gmsh format for debugging of StructFluidSurfCoupling surface correction,
  //note that this can only be done on one proc
  if ((cond->Type() == DRT::Condition::StructFluidSurfCoupling) or (cond->Type() == DRT::Condition::RedAirwayTissue))
  {
    DRT::UTILS::WriteBoundarySurfacesVolumeCoupling(surfmap,cond->Id(),cond->comm_->NumProc(),cond->comm_->MyPID());
  }

  // Surfaces be added to the condition: (line_id) -> (surface).
  Teuchos::RCP<std::map< int, Teuchos::RCP<DRT::Element> > > finalsurfs = Teuchos::rcp(new std::map<int,Teuchos::RCP<DRT::Element> >() );

  AssignGlobalIDs( Comm(), surfmap, *finalsurfs );
  cond->AddGeometry(finalsurfs);

  // elements where created that need new unique ids
  bool havenewelements = true;
  // note: this seems useless since this BuildSurfacesinCondition always
  //       creates elements. However, this function is overloaded in
  //       MeshfreeDiscretization where it does not necessarily build elements.
  return havenewelements;

} // DRT::Discretization::BuildSurfacesinCondition


/*----------------------------------------------------------------------*
 |  Build volume geometry in a condition (public)            mwgee 01/07|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::BuildVolumesinCondition(
                                        const std::string name,
                                        Teuchos::RCP<DRT::Condition> cond)
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
  Teuchos::RCP<std::map<int,Teuchos::RCP<DRT::Element> > >  geom = Teuchos::rcp(new std::map<int,Teuchos::RCP<DRT::Element> >() );

  for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator actele=element_.begin();
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
      (*geom)[actele->first] = actele->second;
    }
  }

  cond->AddGeometry(geom);

  // no elements where created to assign new unique ids to
  return false;
} // DRT::Discretization::BuildVolumesinCondition



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::FindAssociatedEleIDs(Teuchos::RCP<DRT::Condition> cond, std::set<int>& VolEleIDs, const std::string name)
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

      for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator actele=element_.begin();
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

/*--------------------------------------------------------------------------*
|  Set new nodes to a condition & build geometry new            (public)    |
|                                                                ager 11/14 |
 *-------------------------------------------------------------------------*/
//Do not use this function to manipulate conditions at the beginning of a calculation!
//Should just be used for dynamically changing conditions!
//Will always act on the first condition with this condname (if there is just one, this is ok!)
void DRT::Discretization::ModifyCondition(const std::string& condname, std::vector<int>& nodegid, bool buildgeom)
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  DRT::Condition* cond = GetCondition(condname);
  cond->Delete("Node Ids");
  cond->Add("Node Ids", nodegid);
  if (buildgeom) BoundaryConditionsGeometry();
}

