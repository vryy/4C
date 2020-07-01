/*----------------------------------------------------------------------*/
/*! \file

\brief class that provides all routines to handle cutsides which cut each other

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "cut_kernel.H"
#include "cut_point_impl.H"
#include "cut_pointgraph.H"
#include "cut_mesh.H"
#include "cut_meshhandle.H"
#include "cut_triangulateFacet.H"
#include "cut_side.H"
#include "cut_pointpool.H"
#include "cut_options.H"
#include "cut_output.H"

#include "../drt_geometry/searchtree.H"

#include "cut_selfcut.H"

#include <Teuchos_TimeMonitor.hpp>

//#define DEBUG_SELFCUT

/*-------------------------------------------------------------------------------------*
 * constructor                                                              wirtz 05/13
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SelfCut::SelfCut(MeshHandle& cut_mesh_handle, int myrank)
    : myrank_(myrank), mesh_(cut_mesh_handle.LinearMesh()), meshhandle_(cut_mesh_handle)
{
  meshsizeparam_ = 1e200;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator sit = mesh_.Sides().begin();
       sit != mesh_.Sides().end(); ++sit)
  {
    Teuchos::RCP<BoundingBox> sbb = Teuchos::rcp(BoundingBox::Create(*(sit->second)));
    meshsizeparam_ = std::min(meshsizeparam_, sbb->Diagonal());
  }
}

void GEO::CUT::SelfCut::PerformSelfCut()
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_Mesh --- SELFCUT");
  if (CollisionDetection())
  {
    MeshIntersection();
    ElementSelection();
    if (mesh_.GetOptions().SelfCut_Do_MeshCorrection()) MeshCorrection();
  }
}

/*-------------------------------------------------------------------------------------*
 * detects cutsides which cut each other by finding the respective cutpoints
 *                                                                          wirtz 05/13
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::SelfCut::CollisionDetection()
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_Mesh --- SELFCUT --- CollisionDetection");

  FindCuttingSides();
  FindSelfCutPoints();
  GetSelfCutObjects();

  return (selfcut_sides_.size() > 0);
}

/*-------------------------------------------------------------------------------------*
 * replaces cutted sides by creating new nodes, edges and sides             wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::MeshIntersection()
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_Mesh --- SELFCUT --- MeshIntersection");

  CreateSelfCutNodes();
  CreateSelfCutEdges();
  FindSelfCutTriangulation();
  CreateSelfCutSides();
  EraseCuttedSides();
  EraseCuttedEdges();
}

/*-------------------------------------------------------------------------------------*
 * erases nodes, edges and sides which lies inside a structure body by
 * locating there position                                                  wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::ElementSelection()
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_Mesh --- SELFCUT --- ElementSelection");

  DetermineSelfCutPosition();
  PropagateSelfCutPosition();
  EraseInsideSides();
  EraseInsideEdges();
  EraseInsideNodes();
}

/*-------------------------------------------------------------------------------------*
 * repair the mesh from potential islands                                   wirtz 07/16
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::MeshCorrection()
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_Mesh --- SELFCUT --- MeshCorrection");

  ConstructConnectivity();
  FindIslands();
  EraseInsideSides();
  EraseInsideEdges();
  EraseInsideNodes();
}

/*-------------------------------------------------------------------------------------*
 * represents the result by using console or gmsh                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::StatusInspection()
{
  CuttedSideStatusText();
  CutMeshStatusText();
  CuttedSideStatusGmsh("selfcut.pos");
  SCmgmGmsh("SCmgmGmsh.pos");
}

/*-------------------------------------------------------------------------------------*
 * finds cutsides which possibly cuts the considered side                   wirtz 05/13
 * we construct bounding box over each side and collect all cut sides that intersect
 * this bounding box. Hence possible self cut sides that we collect here may not
 * actually produce self-cut --> we check this in subsequent procedure
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::FindCuttingSides()
{
  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();
  const std::map<int, LINALG::Matrix<3, 2>>& selfcutbvs = mesh_.SelfCutBvs();
  const Teuchos::RCP<GEO::SearchTree>& selfcuttree = mesh_.SelfCutTree();
  const std::map<int, Side*>& shadowsides = mesh_.ShadowSides();

  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    LINALG::Matrix<3, 2> cutsideBV = cutside->GetBoundingVolume().GetBoundingVolume();
    if (selfcutbvs.size() != 0)  // ******************************************************* possible
                                 // in case of parallel computing and using only relevant elements
    {
      std::set<int> collisions;
      // search collision between cutside and elements in the static search tree
      selfcuttree->searchCollisions(selfcutbvs, cutsideBV, 0, collisions);

      // Quick check if another side consists of the same points
      bool skip_cutside = false;
      for (std::set<int>::iterator ic = collisions.begin(); ic != collisions.end(); ++ic)
      {
        int collisionid = *ic;
        std::map<int, Side*>::const_iterator is = shadowsides.find(collisionid);
        if (is != shadowsides.end())
        {
          Side* s = is->second;
          if (s == cutside) continue;

          if (cutside->AllPointsCommon(*s)) skip_cutside = true;
        }
        if (skip_cutside) break;
      }
      if (skip_cutside) continue;
      // add the cutside to all the elements which have been found
      for (std::set<int>::iterator ic = collisions.begin(); ic != collisions.end(); ++ic)
      {
        int collisionid = *ic;
        std::map<int, Side*>::const_iterator is = shadowsides.find(collisionid);
        if (is != shadowsides.end())
        {
          Side* s = is->second;
          if (s == cutside) continue;
          if (not cutside->HaveCommonNode(*s))
          {
            // if "side" and "cutside" has two nodes at the same position, we merge them
            // if these two sides overlap within BB only because of coinciding nodes,
            // we do not take this as self-cut
            if (MergeCoincidingNodes(cutside, s))
            {
              // even after merging nodes, if two sides not have a common node --> self-cut possible
              if (not cutside->HaveCommonNode(*s))
              {
                cutside->GetCuttingSide(s);
              }
            }
            else
            {
              cutside->GetCuttingSide(s);
            }
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * If two nodes are at the same position, then this function will carry out
 * equivalent operations such that only one node is available at this place
 * Connectivity informations are incorporated accordingly                       sudhakar 10/13
 *
 * In the following, nodes 4-5 and 3-6 are at the same position, and the two
 * Quads are not connected. In order to avoid many complicated routines of handling such
 * coinciding nodes in self-cut libraries, we delete one corresponding node and ensure
 * connectivity betweent the Quads. The cut library can easily handle this modified confi.
 *
 *   1     4,5     8               1      4      8
 *    +-----+-----+                 +-----+-----+
 *    |     |     |                 |     |     |
 *    |     |     |     ===>        |     |     |
 *    |     |     |                 |     |     |
 *    +-----+-----+                 +-----+-----+
 *   2     3,6     7               2      6     7
 *
 *   Rememeber : Which one of the nodes is considered among the coinciding nodes is
 *   not reproducible (just like many operations in cut libraries)
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::SelfCut::MergeCoincidingNodes(Side* keep, Side* replace)
{
  bool merge = false;

  const std::vector<Node*>& cutnodes = keep->Nodes();

  const std::vector<Node*>& scnodes = replace->Nodes();

  // store the nodes that should be replaced with other ids
  std::vector<std::pair<const Node*, const Node*>> ids_replace;

  for (std::vector<Node*>::const_iterator noit = scnodes.begin(); noit != scnodes.end(); noit++)
  {
    const Node* scnod = *noit;

    for (std::vector<Node*>::const_iterator cutit = cutnodes.begin(); cutit != cutnodes.end();
         cutit++)
    {
      const Node* cutnod = *cutit;

      // if both nodes are the same - no need to do anything
      if (scnod->Id() == cutnod->Id()) continue;

      // Two nodes are at the same location in the initial mesh
      if (scnod->point()->Id() == cutnod->point()->Id())
      {
        std::pair<const Node*, const Node*> nodpair;
        nodpair.first = scnod;
        nodpair.second = cutnod;

        ids_replace.push_back(nodpair);
      }

      // In moving boundary simulation, two nodes from different discretizations
      // move such that two nodes reached to the same position
      // if we implement this, the previous loop (two nodes at the same position in initial mesh)
      // can be erased
      // I am not even sure whether this special check is necessary (Sudhakar)
      // because we always works with a copy of cut mesh, and only one point is
      // created for two coinciding nodes ---> need to be checked
      else if (cutnod->isAtSameLocation(scnod))
      {
        dserror("not implemented yet");
      }

      // now that we have all the informations about coinciding nodes
      // Do the necessary operations in the mesh to consider only one node
      if (ids_replace.size() > 0)
      {
        operationsForNodeMerging(ids_replace, true);
        merge = true;
      }
    }
  }

  return merge;
}

/*-----------------------------------------------------------------------------------------*
 * After all informations about coinciding nodes are available, perform the       sudhakar 10/13
 * operations in our data structure such that only one node is available and
 * ensure connectivity between the corresponding elements
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::operationsForNodeMerging(
    std::vector<std::pair<const Node*, const Node*>> repl, bool initial)
{
  // The following "edge" and "side" data structure has to be modified accodingly
  // the node numbers has to be set accordingly
  std::map<plain_int_set, Teuchos::RCP<Edge>>& edges =
      const_cast<std::map<plain_int_set, Teuchos::RCP<Edge>>&>(mesh_.Edges());
  std::map<plain_int_set, Teuchos::RCP<Side>>& sides =
      const_cast<std::map<plain_int_set, Teuchos::RCP<Side>>&>(mesh_.Sides());

  /*for (std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator it =
  edges.begin();it!=edges.end();it++)
  {
    plain_int_set pis = it->first;
    std::cout<<"edge ";
    for( unsigned i=0;i<pis.size();i++ )
      std::cout<<pis[i]<<" ";
    std::cout<<"\n";
  }*/

  // consider each node that needs to be replaced
  for (std::vector<std::pair<const Node*, const Node*>>::iterator it = repl.begin();
       it != repl.end(); it++)
  {
    std::pair<const Node*, const Node*>& pai = *it;
    Node* nod = const_cast<Node*>(pai.first);
    Node* replwith = const_cast<Node*>(pai.second);

    const plain_side_set& sideset = nod->Sides();

    // consider each side that are associated with this node
    // and replace it with correct node
    // Within this, internally edge information gets modified correspondingly
    for (plain_side_set::const_iterator i = sideset.begin(); i != sideset.end(); ++i)
    {
      Side* s = *i;
      s->replaceNodes(nod, replwith);
    }

    // modify the "edge" and "side data structure in the mesh
    // delete node "nod" by "replwith" in both edges and sides data-structure
    ModifyEdgeOrSideMap<plain_int_set, Teuchos::RCP<Edge>>(edges, nod->Id(), replwith->Id());
    ModifyEdgeOrSideMap<plain_int_set, Teuchos::RCP<Side>>(sides, nod->Id(), replwith->Id());
  }
}
/*----------------------------------------------------------------------------------------------------*
 * Templated function to delete "nod" by "replwith" in both edge and side data-structures sudhakar
 *10/13 This is a bit tricky because we need to modify the "key" of this std::map "Key" of a
 *std::map cannot be modified directly. So, A new element with correct information is added, and old
 *element is deleted
 *----------------------------------------------------------------------------------------------------*/
template <typename A, typename B>
void GEO::CUT::SelfCut::ModifyEdgeOrSideMap(std::map<A, B>& data, int nod, int replwith)
{
  typedef typename std::map<A, B>::iterator ittype;
  typedef typename A::iterator A_ittype;

  std::vector<ittype> eraseThese;

  // -----
  // loop over each element in data to check whether "nod" is available
  // if found, then we create a new element for data with correct "key" in the map
  // we add this and delete the element with wrong key
  // -----
  for (ittype itmap = data.begin(); itmap != data.end(); itmap++)
  {
    A idset = itmap->first;
    A newids;
    bool modi = false;

    for (A_ittype ii = idset.begin(); ii != idset.end(); ii++)
    {
      if (*ii == nod)
      {
        newids.insert(replwith);
        modi = true;
      }
      else
        newids.insert(*ii);
    }

    // insert correct elements
    if (modi)
    {
      data[newids] = data[idset];
      eraseThese.push_back(itmap);
    }
  }

  typedef typename std::vector<ittype>::iterator itc;

  // delete the elements with "wrong" key
  for (itc ite = eraseThese.begin(); ite != eraseThese.end(); ite++)
  {
    ittype era = *ite;
    data.erase(era);
  }
}

/*-------------------------------------------------------------------------------------*
 * finds cutpoints by searching the intersection between the edges
 * of one side with the other side and vice versa                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::FindSelfCutPoints()
{
  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const plain_side_set& possiblecuttingsides = cutside->CuttingSides();
    plain_side_set nocuttingsides;
    for (plain_side_set::const_iterator i = possiblecuttingsides.begin();
         i != possiblecuttingsides.end(); ++i)
    {
      Side* possiblecuttingside = *i;

      PointSet selfcutpoints;
      PerformSelfCut(*cutside, *possiblecuttingside, selfcutpoints);
      PerformSelfCut(*possiblecuttingside, *cutside, selfcutpoints);
      cutside->GetSelfCutPoints(selfcutpoints);

      // possible cut sides that are intersecting the bounding box of "this" side
      // but not actually cutting, are stored for deletion
      if (selfcutpoints.size() == 0)
      {
        nocuttingsides.insert(possiblecuttingside);
        possiblecuttingside->EraseCuttingSide(cutside);
      }
    }
    for (plain_side_set::iterator i = nocuttingsides.begin(); i != nocuttingsides.end(); ++i)
    {
      Side* nocuttingside = *i;
      cutside->EraseCuttingSide(nocuttingside);
    }
  }
#ifdef DEBUG_SELFCUT
  Debug("FindSelfCutPoints");
#endif
}

/*-------------------------------------------------------------------------------------*
 * gets all cutted sides and their nodes and edges to store them as
 * private variables                                                         wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::GetSelfCutObjects()
{
  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Teuchos::RCP<Side> cutside = i->second;
    if (cutside->CuttingSides().size() != 0)
    {
      plain_int_set cutsidenodeids = i->first;
      selfcut_sides_[cutsidenodeids] = cutside;
    }
  }

  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    GetSelfCutEdges(*cutside);
  }
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<Node*>& cutsidenodes = cutside->Nodes();
    for (std::vector<Node*>::const_iterator i = cutsidenodes.begin(); i != cutsidenodes.end(); ++i)
    {
      int cutsidenodeid = (*i)->Id();
      const std::map<int, Teuchos::RCP<Node>>& cutsidenodercp = mesh_.Nodes();
      std::map<int, Teuchos::RCP<Node>>::const_iterator nodeiterator =
          cutsidenodercp.find(cutsidenodeid);
      selfcut_nodes_[cutsidenodeid] = nodeiterator->second;
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * creates a node at every selfcutpoint with respect to the
 * corresponding cutsides                                                   wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CreateSelfCutNodes()
{
  std::vector<Node*> cuttedsidesnodes;
  for (std::map<int, Teuchos::RCP<Node>>::iterator i = selfcut_nodes_.begin();
       i != selfcut_nodes_.end(); ++i)
  {
    Node* cuttedsidesnode = &*i->second;
    cuttedsidesnodes.push_back(cuttedsidesnode);
  }

  // we take all selfcut points of each side
  // 1.if these points are at any of the nodes of cut side, then no need to create a new node
  //   just set their position as "oncutsurface"
  // 2.If not we create a new node at this location and add it to the mesh, and corresponding cut
  // side
  //   set their position as "oncutsurface"
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    const PointSet& cutsideselfcutpoints = cutside->SelfCutPoints();
    for (PointSet::const_iterator i = cutsideselfcutpoints.begin(); i != cutsideselfcutpoints.end();
         ++i)
    {
      Point* cutsideselfcutpoint = *i;
      if (cutsideselfcutpoint->NodalPoint(cuttedsidesnodes))
      {
        for (std::vector<Node*>::iterator i = cuttedsidesnodes.begin(); i != cuttedsidesnodes.end();
             ++i)
        {
          Node* cuttedsidesnode = *i;
          if (cuttedsidesnode->point() == cutsideselfcutpoint)
          {
            cuttedsidesnode->SelfCutPosition(Point::oncutsurface);
            cutside->GetSelfCutNode(cuttedsidesnode);
            break;
          }
        }
      }
      else
      {
        const std::map<int, Teuchos::RCP<Node>>& cutsidenodes = mesh_.Nodes();
        int selfcutnodeid = cutsidenodes.begin()->first - 1;
        Node* selfcutnode = new Node(selfcutnodeid, cutsideselfcutpoint, 0.0);
        selfcutnode->SelfCutPosition(Point::oncutsurface);
        mesh_.GetNode(selfcutnodeid, selfcutnode);
        const std::map<int, Teuchos::RCP<Node>>& selfcutnodercp = mesh_.Nodes();
        std::map<int, Teuchos::RCP<Node>>::const_iterator nodeiterator =
            selfcutnodercp.find(selfcutnodeid);
        selfcut_nodes_[selfcutnodeid] = nodeiterator->second;
        cutside->GetSelfCutNode(selfcutnode);
        cuttedsidesnodes.push_back(selfcutnode);
      }
    }
  }
#ifdef DEBUG_SELFCUT
  Debug("CreateSelfCutNodes");
#endif
}

/*-------------------------------------------------------------------------------------*
 * creates a edge between two selfcutnodes with respect to the
 * corresponding cutsides                                                   wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CreateSelfCutEdges()
{
  // Find common nodes between a cut side and its self-cutting sides
  // create edges between these common nodes
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    const plain_side_set& cuttingsides = cutside->CuttingSides();
    const plain_node_set& cutside_nodes = cutside->SelfCutNodes();
    for (plain_side_set::const_iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
    {
      Side* cuttingside = *i;
      const plain_node_set& cuttingside_nodes = cuttingside->SelfCutNodes();
      std::vector<Node*> commonselfcutnodes;
      plain_int_set commonselfcutnodeids;
      for (plain_node_set::const_iterator i = cutside_nodes.begin(); i != cutside_nodes.end(); ++i)
      {
        Node* cutsideselfcutnode = *i;
        plain_node_set::const_iterator nodeiterator = cuttingside_nodes.find(cutsideselfcutnode);
        if (nodeiterator != cuttingside_nodes.end())
        {
          Node* commonselfcutnode = *nodeiterator;
          commonselfcutnodes.push_back(commonselfcutnode);
          commonselfcutnodeids.insert(commonselfcutnode->Id());
        }
      }
      if (commonselfcutnodeids.size() == 2)
      {
        std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator edgeiterator =
            selfcut_edges_.find(commonselfcutnodeids);
        if (edgeiterator !=
            selfcut_edges_.end())  // this edge is already formed when processing other side
        {
          Edge* commonedge = &*edgeiterator->second;
          commonedge->SelfCutPosition(Point::oncutsurface);
          cutside->GetSelfCutEdge(commonedge);
        }
        else
        {
          Teuchos::RCP<Edge> selfcutedge =
              GEO::CUT::Edge::Create(DRT::Element::line2, commonselfcutnodes);
          selfcutedge->SelfCutPosition(Point::oncutsurface);
          mesh_.GetEdge(commonselfcutnodeids, selfcutedge);
          const std::map<plain_int_set, Teuchos::RCP<Edge>>& cutsideedgercp = mesh_.Edges();
          std::map<plain_int_set, Teuchos::RCP<Edge>>::const_iterator edgeiterator =
              cutsideedgercp.find(commonselfcutnodeids);
          selfcut_edges_[commonselfcutnodeids] = edgeiterator->second;
          cutside->GetSelfCutEdge(selfcutedge.get());
        }
      }
    }
  }
#ifdef DEBUG_SELFCUT
  Debug("CreateSelfCutEdges");
#endif
}

/*-------------------------------------------------------------------------------------*
 * finds triangles with respect to the self cut by using
 * a pointgraph and a triangulation method                                  wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::FindSelfCutTriangulation()
{
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = i->second.get();
    const std::vector<Node*>& cutsidenodes = cutside->Nodes();

    // find the equation of plane of this SIDE
    // this eqn is used to make sure that the axillary sides that resulting from
    // the triangulation preserve the normal direction of SIDE
    std::vector<Point*> cutsidepoints(cutsidenodes.size());
    for (unsigned ic = 0; ic < cutsidenodes.size(); ic++)
      cutsidepoints[ic] = cutsidenodes[ic]->point();

    cutsidepoints[0] = cutsidenodes[0]->point();
    cutsidepoints[1] = cutsidenodes[1]->point();
    cutsidepoints[2] = cutsidenodes[2]->point();
    std::vector<double> cutsideplane = KERNEL::EqnPlaneOfPolygon(cutsidepoints);

    // -----
    // STEP 1 : Create facets on cut side by taking into account the self-cut points
    // -----
    IMPL::PointGraph pointgraph(cutside);
    Cycle* maincycle = &*(pointgraph.fbegin());
    bool IsCorrectNormal = CheckNormal(cutsideplane, *maincycle);
    std::vector<Cycle*> maincycles;
    std::vector<Teuchos::RCP<Cycle>> newmaincyclercp;
    for (IMPL::PointGraph::facet_iterator i = pointgraph.fbegin(); i != pointgraph.fend(); ++i)
    {
      Cycle* maincycle = &*i;
      if (not IsCorrectNormal)
      {
        // normal is not in correct direction
        // we reverse the ordering of points to correct this
        maincycle->reverse();
      }
      maincycles.push_back(maincycle);
    }

    std::vector<Cycle*> holecycles;
    bool holenormalconservation = true;
    for (IMPL::PointGraph::hole_iterator i = pointgraph.hbegin(); i != pointgraph.hend(); ++i)
    {
      std::vector<Cycle> oldholecycles = *i;
      Cycle* holecycle = &oldholecycles[0];
      holenormalconservation = CheckNormal(cutsideplane, *holecycle);
      break;
    }
    bool firstinnerloop = true;
    for (IMPL::PointGraph::hole_iterator i = pointgraph.hbegin(); i != pointgraph.hend(); ++i)
    {
      std::vector<Cycle>* holecyclesold = &*i;
      firstinnerloop = true;
      for (std::vector<Cycle>::iterator i = holecyclesold->begin(); i != holecyclesold->end(); ++i)
      {
        Cycle* holecycle = &*i;
        if ((holenormalconservation and firstinnerloop) or
            (not holenormalconservation and not firstinnerloop))
        {
          maincycles.push_back(holecycle);
        }
        if ((holenormalconservation and not firstinnerloop) or
            (not holenormalconservation and firstinnerloop))
        {
          holecycles.push_back(holecycle);
        }
        firstinnerloop = false;
      }
    }
    std::vector<std::vector<Point*>> mainholecyclepoints;
    for (std::vector<Cycle*>::iterator i = holecycles.begin(); i != holecycles.end(); ++i)
    {
      Cycle* holecycle = *i;
      std::vector<Point*> holecyclepoints = (*holecycle)();
      mainholecyclepoints.push_back(holecyclepoints);
    }

    // -----
    // STEP 2 : Triangulate each facet using Earclipping algorithm
    // -----
    std::vector<std::vector<Point*>> allmaincycletriangles;
    for (std::vector<Cycle*>::iterator i = maincycles.begin(); i != maincycles.end(); ++i)
    {
      Cycle* maincycle = *i;
      std::vector<Point*> maincyclepoints = (*maincycle)();
      TriangulateFacet triangulatefacet(maincyclepoints, mainholecyclepoints);
      triangulatefacet.EarClippingWithHoles(cutside);
      std::vector<std::vector<Point*>> maincycletriangles = triangulatefacet.GetSplitCells();
      for (std::vector<std::vector<Point*>>::iterator i = maincycletriangles.begin();
           i != maincycletriangles.end(); ++i)
      {
        std::vector<Point*> maincycletriangle = *i;
        if (KERNEL::IsOnLine(maincycletriangle[0], maincycletriangle[1], maincycletriangle[2]))
        {
          ErrorStatusText(*cutside);
          ErrorGmsh("triangle_with_collinear_points.pos", *cutside);
          //          Output(*cutside, "FindSelfCutTriangulation");
          std::cout
              << "WARNING:::selfcut algorithm produced a triangle with all points on a line\n";
        }
        if (maincycletriangle.size() != 3)
        {
          //          Output(*cutside, "FindSelfCutTriangulation");
          throw std::runtime_error(
              "SelfCut: triangulation unsuccessful; triangle without 3 points");
        }
        if (maincycletriangle[0]->Id() == maincycletriangle[1]->Id() or
            maincycletriangle[0]->Id() == maincycletriangle[2]->Id() or
            maincycletriangle[1]->Id() == maincycletriangle[2]->Id())
        {
          throw std::runtime_error(
              "SelfCut: triangulation unsuccessful; triangle with two identical points");
        }
        cutside->GetSelfCutTriangle(maincycletriangle);
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * creates triangular sides out of the self cut triangles                   wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CreateSelfCutSides()
{
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    plain_node_set allcutsidenodes = cutside->SelfCutNodes();
    const std::vector<Node*>& cutsidenodes = cutside->Nodes();

    // store all the nodes (self-cut nodes + side nodes)
    // this is used to check we really have a node at each point of triangle
    for (std::vector<GEO::CUT::Node*>::const_iterator i = cutsidenodes.begin();
         i != cutsidenodes.end(); ++i)
    {
      Node* cutsidenode = *i;
      allcutsidenodes.insert(cutsidenode);
    }

    const std::vector<std::vector<Point*>>& selfcuttriangles = cutside->SelfCutTriangles();

    for (std::vector<std::vector<Point*>>::const_iterator i = selfcuttriangles.begin();
         i != selfcuttriangles.end(); ++i)
    {
      std::vector<Point*> selfcuttriangle = *i;
      std::vector<int> selfcutsidenodeids;
      plain_int_set selfcutsidenodeidsset;
      for (std::vector<GEO::CUT::Point*>::iterator i = selfcuttriangle.begin();
           i != selfcuttriangle.end(); ++i)
      {
        Point* selfcuttrianglepoint = *i;
        for (plain_node_set::iterator i = allcutsidenodes.begin(); i != allcutsidenodes.end(); ++i)
        {
          Node* allcutsidenode = *i;
          if (allcutsidenode->point() == selfcuttrianglepoint)
          {
            int allcutsidenodeid = allcutsidenode->Id();
            selfcutsidenodeids.push_back(allcutsidenodeid);
            selfcutsidenodeidsset.insert(allcutsidenodeid);
          }
        }
      }
      if (selfcutsidenodeidsset.size() != 3)
      {
        std::cout << "selfcutsidenodeidsset.size(): " << selfcutsidenodeidsset.size() << "\n";
        const std::vector<Node*>& cutsidenodes = cutside->Nodes();
        for (std::vector<Node*>::const_iterator i = cutsidenodes.begin(); i != cutsidenodes.end();
             ++i)
        {
          Node* cutsidenode = *i;
          std::cout << "Ids: " << cutsidenode->Id() << "\n";
        }
        //        Output(*cutside, "CreateSelfCutSides");
        throw std::runtime_error("SelfCut: creating sides unsuccessful; didn't find 3-id-set");
      }
      if (selfcutsidenodeids.size() != 3)
      {
        std::cout << "selfcutsidenodeids.size(): " << selfcutsidenodeids.size() << "\n";
        //        Output(*cutside, "CreateSelfCutSides");
        throw std::runtime_error("SelfCut: creating sides unsuccessful; didn't find 3 ids");
      }
      if (selfcutsidenodeids[0] == selfcutsidenodeids[1] or
          selfcutsidenodeids[0] == selfcutsidenodeids[2] or
          selfcutsidenodeids[1] == selfcutsidenodeids[2])
      {
        //        Output(*cutside, "CreateSelfCutSides");
        throw std::runtime_error(
            "SelfCut: triangulation unsuccessful; triangle with two identical points");
      }
      Side* selfcutside = mesh_.CreateSide(cutside->Id(), selfcutsidenodeids, DRT::Element::tri3);
      meshhandle_.AddSubSide(selfcutside);
      const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsideedgercp = mesh_.Sides();
      std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator sideiterator =
          cutsideedgercp.find(selfcutsidenodeidsset);
      selfcut_sides_[selfcutsidenodeidsset] = sideiterator->second;
      GetSelfCutEdges(*selfcutside);
    }
  }
#ifdef DEBUG_SELFCUT
  Debug("CreateSelfCutSides");
#endif
}

/*-------------------------------------------------------------------------------------*
 * erases all cutsides which are cut by another side                        wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseCuttedSides()
{
  std::vector<plain_int_set> cutsideids;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    if (cutside->CuttingSides().size() != 0 and cutside->SelfCutTriangles().size() != 1)
    {
      cutsideids.push_back(i->first);
      EraseSidePointer(*cutside, true);
    }
  }
  EraseSide(cutsideids);
}

/*-------------------------------------------------------------------------------------*
 * erases all edges which are cut by a cutside                              wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseCuttedEdges()
{
  std::vector<plain_int_set> cutsideedgeids;
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator i = selfcut_edges_.begin();
       i != selfcut_edges_.end(); ++i)
  {
    Edge* cutsideedge = &*i->second;
    if (cutsideedge->CutPoints().size() != 2 or cutsideedge->Sides().size() == 0)
    {
      if (!ConnectedtoBackground(cutsideedge))
      {
        cutsideedgeids.push_back(i->first);
        EraseEdgePointer(*cutsideedge);
      }
    }
  }
  EraseEdge(cutsideedgeids);
}

bool GEO::CUT::SelfCut::ConnectedtoBackground(Edge* edge)
{
  for (plain_side_set::const_iterator sit = edge->Sides().begin(); sit != edge->Sides().end();
       ++sit)
  {
    if ((*sit)->Id() < 0) return true;
  }
  return false;
}

bool GEO::CUT::SelfCut::ConnectedtoBackground(Node* node)
{
  for (plain_side_set::const_iterator sit = node->Sides().begin(); sit != node->Sides().end();
       ++sit)
  {
    if ((*sit)->Id() < 0) return true;
  }
  return false;
}

/*-------------------------------------------------------------------------------------*
 * locates the position of nodes, edges and sides of
 * a structure body with respect to the other bodys                         wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::DetermineSelfCutPosition()
{
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator i = selfcut_edges_.begin();
       i != selfcut_edges_.end(); ++i)
  {
    Edge* cutsideedge = &*i->second;
    std::vector<Side*> selfcutsides;
    plain_side_set cutsides;
    if (cutsideedge->SelfCutPosition() == Point::oncutsurface)
    {
      cutsides = cutsideedge->Sides();  // always 4 sides
      for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
      {
        Side* cutside = *i;
        const std::vector<Node*>& cutsidenodes = cutside->Nodes();
        for (std::vector<Node*>::const_iterator i = cutsidenodes.begin(); i != cutsidenodes.end();
             ++i)
        {
          Node* cutsidenode = *i;
          if (cutsidenode->SelfCutPosition() == Point::undecided)
          {
            selfcutsides.push_back(cutside);
            continue;
          }
        }
      }
    }

    for (std::vector<Side*>::iterator i = selfcutsides.begin(); i != selfcutsides.end(); ++i)
    {
      Side* selfcutside = *i;
      const std::vector<Node*>& selfcutsidenodes = selfcutside->Nodes();
      Node* undecidednode = NULL;
      Node* onselfcutedgenode = NULL;
      for (std::vector<Node*>::const_iterator i = selfcutsidenodes.begin();
           i != selfcutsidenodes.end(); ++i)
      {
        Node* selfcutsidenode = *i;
        if ((selfcutsidenode->SelfCutPosition() == Point::oncutsurface) &&
            (onselfcutedgenode == NULL))
        {
          onselfcutedgenode = selfcutsidenode;
        }
        else if (selfcutsidenode->SelfCutPosition() == Point::undecided)
        {
          undecidednode = selfcutsidenode;
        }
      }
      Side* otherselfcutside = NULL;
      for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
      {
        Side* anotherselfcutside = *i;
        if (anotherselfcutside->Id() != selfcutside->Id())
        {
          LINALG::Matrix<3, 1> normal;
          LINALG::Matrix<2, 1> center(true);
          anotherselfcutside->Normal(center, normal, false);
          double norm = normal.Norm2();
          if (norm > 1e-8)
          {
            otherselfcutside = anotherselfcutside;
            break;
          }
          else
          {
            std::cout << "==| WARNING: Avoided anotherside with normal norm of " << norm << " ( "
                      << myrank_ << ")! |==" << std::endl;
          }
        }
      }
      if (otherselfcutside == NULL or undecidednode == NULL)
      {
        continue;
      }

      // this means that we are deling with self-cut with coinciding nodes
      // we already merged the nodes appropriately
      // hence in this case, all nodes associated are always "outside" --> should be included
      /*if( selfcutside->HaveCommonNode(*otherselfcutside) )
      {
        undecidednode->SelfCutPosition(Point::oncutsurface);
        continue;
      }*/

      LINALG::Matrix<3, 1> oncut_cord;
      LINALG::Matrix<3, 1> oncut_cord_loc;
      LINALG::Matrix<2, 1> oncut_cord_loc2;
      LINALG::Matrix<3, 1> otherSideNormal;
      onselfcutedgenode->point()->Coordinates(oncut_cord.A());
      otherselfcutside->LocalCoordinates(oncut_cord, oncut_cord_loc, false);
      oncut_cord_loc2(0) = oncut_cord_loc(0);
      oncut_cord_loc2(1) = oncut_cord_loc(1);
      otherselfcutside->Normal(oncut_cord_loc2, otherSideNormal);
      LINALG::Matrix<3, 1> undecidedpointcoordinates;
      LINALG::Matrix<3, 1> differencebetweenpoints;
      undecidednode->point()->Coordinates(undecidedpointcoordinates.A());
      differencebetweenpoints.Update(1.0, oncut_cord, -1.0, undecidedpointcoordinates);
      double norm = differencebetweenpoints.Norm2();
      if (norm < SELF_CUT_POS_TOL)
      {
        std::cout << "==| WARNING: Avoided this point with distance norm of " << norm << " ( "
                  << myrank_ << ")! |==" << std::endl;
        continue;
      }
      else
        differencebetweenpoints.Scale(1. / norm);
      LINALG::Matrix<1, 1> innerproduct;
      innerproduct.MultiplyTN(differencebetweenpoints, otherSideNormal);
      if (innerproduct(0) > 0 and fabs(innerproduct(0)) > SELF_CUT_POS_TOL)
      {
        undecidednode->SelfCutPosition(Point::inside);
      }
      else if (innerproduct(0) < 0 and fabs(innerproduct(0)) > SELF_CUT_POS_TOL)
      {
        undecidednode->SelfCutPosition(Point::outside);
      }
      else if (innerproduct(0) < 0)
      {
        std::cout << "==| WARNING: I assign outside for an innerproduct of " << innerproduct(0)
                  << " ( " << myrank_ << ")! |==" << std::endl;
        undecidednode->SelfCutPosition(Point::outside);
      }
      else
      {
        std::cout << "==| WARNING: I assign inside for an innerproduct of " << innerproduct(0)
                  << " ( " << myrank_ << ")! |==" << std::endl;
        undecidednode->SelfCutPosition(Point::inside);
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * creates triangular sides with respect to the corresponding
 * cutsides by using a pointgraph and a triangulation method                wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::PropagateSelfCutPosition()
{
  plain_side_set undecidedsides;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    if (cutside->SelfCutPosition() == Point::undecided)
    {
      const std::vector<Edge*>& cutsideedges = cutside->Edges();
      Point::PointPosition cutsideedgeposition = Point::undecided;
      for (std::vector<Edge*>::const_iterator i = cutsideedges.begin(); i != cutsideedges.end();
           ++i)
      {
        Edge* cutsideedge = *i;
        Point::PointPosition cutsideedgepos = cutsideedge->SelfCutPosition();
        switch (cutsideedgepos)
        {
          case Point::undecided:
          case Point::oncutsurface:
            break;
          case Point::inside:
          case Point::outside:
            if ((cutsideedgeposition != Point::undecided) and
                (cutsideedgepos != cutsideedgeposition))
            {
              throw std::runtime_error("SelfCut: mixed PointPosition of cutsideedges");
            }
            cutsideedgeposition = cutsideedgepos;
            break;
        }
      }
      if (cutsideedgeposition != Point::undecided)
      {
        cutside->GetSelfCutPosition(cutsideedgeposition);
      }
      else
      {
        undecidedsides.insert(cutside);
      }
    }
  }
  if (selfcut_sides_.size() == undecidedsides.size() and selfcut_sides_.size() > 0)
  {
    throw std::runtime_error("SelfCut: all cutsides undecided");
  }
  while (undecidedsides.size() > 0)
  {
    unsigned undecidedsidesize = undecidedsides.size();
    for (plain_side_set::iterator ui = undecidedsides.begin(); ui != undecidedsides.end();)
    {
      Side* undecidedside = *ui;
      if (undecidedside->SelfCutPosition() == Point::undecided)
      {
        bool done = false;
        const std::vector<Edge*>& undecidedcutsideedges = undecidedside->Edges();
        for (std::vector<Edge*>::const_iterator i = undecidedcutsideedges.begin();
             i != undecidedcutsideedges.end(); ++i)
        {
          Edge* undecidedcutsideedge = *i;
          const plain_side_set& undecidedcutsideedgesides = undecidedcutsideedge->Sides();
          Side* siblingside = NULL;
          for (plain_side_set::const_iterator i = undecidedcutsideedgesides.begin();
               i != undecidedcutsideedgesides.end(); ++i)
          {
            Side* undecidedcutsideedgeside = *i;
            if (undecidedcutsideedgeside->Id() == undecidedside->Id() and
                undecidedcutsideedgeside->SelfCutPosition() != Point::undecided and
                undecidedcutsideedgeside != undecidedside)
            {
              siblingside = undecidedcutsideedgeside;
              break;
            }
          }
          if (siblingside != NULL)
          {
            Point::PointPosition siblingsideposition = siblingside->SelfCutPosition();
            switch (siblingsideposition)
            {
              case Point::undecided:
                if (undecidedsides.count(siblingside))
                {
                  throw std::runtime_error("SelfCut: uncomplete set of undecided cutsides");
                }
                break;
              case Point::oncutsurface:
                throw std::runtime_error("SelfCut: illegal side position");
                break;
              case Point::inside:
              case Point::outside:
                if (undecidedcutsideedge->SelfCutPosition() == Point::oncutsurface)
                {
                  undecidedside->GetSelfCutPosition(
                      siblingsideposition == Point::inside ? Point::outside : Point::inside);
                }
                else
                {
                  undecidedside->GetSelfCutPosition(siblingsideposition);
                }
                done = true;
                break;
            }
            if (done)
            {
              break;
            }
          }
        }
        if (done)
        {
          set_erase(undecidedsides, ui);
        }
        else
        {
          ui++;
        }
      }
      else
      {
        set_erase(undecidedsides, ui);
      }
    }
    if (undecidedsidesize == undecidedsides.size())
      throw std::runtime_error("SelfCut: no progress in cutside position");
  }
}

/*-------------------------------------------------------------------------------------*
 * erases sides which lies inside a structure body by locating
 * their position                                                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseInsideSides()
{
  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();
  std::vector<plain_int_set> cutsideids;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    if (cutside->SelfCutPosition() == Point::inside)
    {
      cutsideids.push_back(i->first);
      EraseSidePointer(*cutside, false);
    }
  }
  EraseInsideSide(cutsideids);
  if (mesh_.Sides().size() == 0) dserror("All self-cut positions are undecided\n");
}

/*-------------------------------------------------------------------------------------*
 * erases edges which lies inside a structure body by locating
 * their position                                                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseInsideEdges()
{
  const std::map<plain_int_set, Teuchos::RCP<Edge>>& cutsideedges = mesh_.Edges();
  std::vector<plain_int_set> cutsideedgeids;
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::const_iterator i = cutsideedges.begin();
       i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = &*i->second;
    if (cutsideedge->SelfCutPosition() == Point::inside)
    {
      if (!ConnectedtoBackground(cutsideedge))
      {
        cutsideedgeids.push_back(i->first);
        EraseEdgePointer(*cutsideedge);
      }
    }
  }
  EraseEdge(cutsideedgeids);
}

/*-------------------------------------------------------------------------------------*
 * erases nodes which lies inside a structure body by locating
 * there position                                                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseInsideNodes()
{
  const std::map<int, Teuchos::RCP<Node>>& cutsidenodes = mesh_.Nodes();
  std::vector<int> cutsidenodeids;
  for (std::map<int, Teuchos::RCP<Node>>::const_iterator i = cutsidenodes.begin();
       i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = &*i->second;
    if (cutsidenode->SelfCutPosition() == Point::inside)
    {
      if (!ConnectedtoBackground(cutsidenode))
      {
        cutsidenodeids.push_back(i->first);
      }
    }
  }
  for (std::vector<int>::iterator i = cutsidenodeids.begin(); i != cutsidenodeids.end(); ++i)
  {
    int cutsidenodeid = *i;
    selfcut_nodes_.erase(cutsidenodeid);
    mesh_.MoveNodetoStorage(cutsidenodeid);
  }
}

/*-------------------------------------------------------------------------------------*
 * construct the connectivity of the nodes to find potential islands in the cut mesh
 *                                                                          wirtz 07/16
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::ConstructConnectivity()
{
  const std::map<int, Teuchos::RCP<Node>>& cutsidenodes = mesh_.Nodes();
  plain_node_set remainingnodes;
  int count = 0;

  for (std::map<int, Teuchos::RCP<Node>>::const_iterator i = cutsidenodes.begin();
       i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = &*i->second;
    remainingnodes.insert(cutsidenode);
  }

  while (remainingnodes.size() > 0)
  {
    Node* node = *remainingnodes.begin();
    selfcut_connectivity_[count].insert(node);
    remainingnodes.erase(node);
    NextNode(node, remainingnodes, count);
    count++;
  }
}

/*-------------------------------------------------------------------------------------*
 * find the next node for the construction of the connectivity              wirtz 07/16
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::NextNode(Node* node, plain_node_set& remainingnodes, int count)
{
  const plain_edge_set& nodeedges = node->Edges();
  for (plain_edge_set::const_iterator i = nodeedges.begin(); i != nodeedges.end(); ++i)
  {
    Edge* nodeedge = *i;
    const std::vector<Node*>& nodeedgenodes = nodeedge->Nodes();
    for (std::vector<Node*>::const_iterator i = nodeedgenodes.begin(); i != nodeedgenodes.end();
         ++i)
    {
      Node* nodeedgenode = *i;
      plain_node_set::iterator remainingnodesiter = remainingnodes.find(nodeedgenode);
      if (remainingnodesiter != remainingnodes.end())
      {
        selfcut_connectivity_[count].insert(nodeedgenode);
        remainingnodes.erase(nodeedgenode);
        NextNode(nodeedgenode, remainingnodes, count);
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * identify the islands                                                     wirtz 07/16
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::FindIslands()
{
  plain_side_set selfcutsides;

  const std::map<plain_int_set, Teuchos::RCP<Side>>& fullcut_sides = mesh_.Sides();

  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = fullcut_sides.begin();
       i != fullcut_sides.end(); ++i)
  {
    Side* cutside = &*i->second;
    if (cutside->SelfCutPosition() == Point::outside) selfcutsides.insert(cutside);
  }

  std::vector<plain_edge_set> nonconnected_selfcutedges;
  while (selfcutsides.size() > 0)
  {
    plain_side_set islandsides;
    bool IsIsland = false;
    Side* cutside = *selfcutsides.begin();
    islandsides.insert(cutside);
    selfcutsides.erase(cutside);

    Teuchos::RCP<BoundingBox> tmp_bb = Teuchos::rcp(BoundingBox::Create(*cutside));
    NextSides(cutside, tmp_bb, selfcutsides, islandsides, IsIsland);

    if (tmp_bb->Diagonal() <=
            mesh_.GetOptions().SelfCut_IslandGeomMultiplicator() * meshsizeparam_ &&
        IsIsland)
    {
      for (plain_side_set::iterator i = islandsides.begin(); i != islandsides.end(); ++i)
      {
        Side* cutside = *i;
        cutside->GetSelfCutPosition(Point::inside);  // overwrite?
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * finds the next sides in a recursive way
 *                                                                          wirtz 02/15
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::NextSides(Side* cutside, Teuchos::RCP<GEO::CUT::BoundingBox>& tmpbb,
    plain_side_set& selfcutsides, plain_side_set& islandsides, bool& IsIsland)
{
  const std::vector<Edge*>& cutsideedges = cutside->Edges();
  for (std::vector<Edge*>::const_iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    {
      const plain_side_set& cutsidenodesides = cutsideedge->Sides();
      for (plain_side_set::const_iterator i = cutsidenodesides.begin(); i != cutsidenodesides.end();
           ++i)
      {
        Side* cutsidenodeside = *i;
        if (cutsidenodeside->SelfCutPosition() == Point::outside)
        {
          plain_side_set::iterator islandsideiter = islandsides.find(cutsidenodeside);
          if (islandsideiter == islandsides.end())  // side has not been treated jet
          {
            if (tmpbb->Diagonal() <=
                mesh_.GetOptions().SelfCut_IslandGeomMultiplicator() * meshsizeparam_)
            {
              IsIsland = true;
              tmpbb->AddPoints(cutsidenodeside->Nodes());
            }
            else
              IsIsland = false;
            islandsides.insert(cutsidenodeside);
            selfcutsides.erase(cutsidenodeside);
            NextSides(cutsidenodeside, tmpbb, selfcutsides, islandsides, IsIsland);
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * Status of the cutted sides for text viewer                               wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CuttedSideStatusText()
{
  int selfcutsidesize = selfcut_sides_.size();
  int j = 1;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    PlotAcross();
    std::cout << "\nCutside " << j << "/" << selfcutsidesize << ":\n";
    SidePlotHead();
    SidePlot(*cutside);

    std::cout << "\nCutting Sides:\n";
    SidePlotHead();
    const plain_side_set& cuttingsides = cutside->CuttingSides();
    for (plain_side_set::const_iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
    {
      Side* cuttingside = *i;
      SidePlot(*cuttingside);
    }

    const PointSet& selfcutpoints = cutside->SelfCutPoints();
    std::cout << "\nSelfcutpoints:\n";
    PointPlotHead();
    for (PointSet::const_iterator i = selfcutpoints.begin(); i != selfcutpoints.end(); ++i)
    {
      Point* selfcutpoint = *i;
      PointPlot(*selfcutpoint);
      std::cout << "\n";
    }

    const plain_node_set& selfcutnodes = cutside->SelfCutNodes();
    std::cout << "\nSelfcutnodes:\n";
    NodePlotHead();
    for (plain_node_set::const_iterator i = selfcutnodes.begin(); i != selfcutnodes.end(); ++i)
    {
      Node* selfcutnode = *i;
      NodePlot(*selfcutnode);
    }

    const plain_edge_set& selfcutedges = cutside->SelfCutEdges();
    std::cout << "\nSelfcutedges:\n";
    EdgePlotHead();
    for (plain_edge_set::const_iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
    {
      Edge* selfcutedge = *i;
      EdgePlot(*selfcutedge);
    }

    const std::vector<std::vector<GEO::CUT::Point*>>& selfcuttriangles =
        cutside->SelfCutTriangles();
    std::cout << "\nSelfcuttriangles:\n";
    for (std::vector<std::vector<GEO::CUT::Point*>>::const_iterator i = selfcuttriangles.begin();
         i != selfcuttriangles.end(); ++i)
    {
      std::vector<GEO::CUT::Point*> selfcuttriangle = *i;
      PointPlotHead();
      for (std::vector<GEO::CUT::Point*>::iterator i = selfcuttriangle.begin();
           i != selfcuttriangle.end(); ++i)
      {
        Point* selfcuttrianglepoint = *i;
        PointPlot(*selfcuttrianglepoint);
        std::cout << "\n";
      }
    }
    PlotAcross();
    j++;
  }
}

/*-------------------------------------------------------------------------------------*
 * Status of the cutmesh for text viewer                                    wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CutMeshStatusText()
{
  PlotAcross();
  std::cout << "\n" << selfcut_nodes_.size() << " Nodes: \n";
  NodePlotHead();
  for (std::map<int, Teuchos::RCP<Node>>::iterator i = selfcut_nodes_.begin();
       i != selfcut_nodes_.end(); ++i)
  {
    Node* node = &*i->second;
    NodePlot(*node);
  }
  PlotAcross();
  std::cout << "\n" << selfcut_edges_.size() << " Edges: \n";
  EdgePlotHead();
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator i = selfcut_edges_.begin();
       i != selfcut_edges_.end(); ++i)
  {
    Edge* edge = &*i->second;
    EdgePlot(*edge);
  }
  PlotAcross();
  std::cout << "\n" << selfcut_sides_.size() << " Sides: \n";
  SidePlotHead();
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* side = &*i->second;
    SidePlot(*side);
  }
  PlotAcross();
}

/*-------------------------------------------------------------------------------------*
 * Status of one problematic side for text viewer                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::ErrorStatusText(Side& cutside)
{
  int selfcutsidesize = 1;
  int j = 1;
  PlotAcross();
  std::cout << "\nCutside " << j << "/" << selfcutsidesize << ":\n";
  SidePlotHead();
  SidePlot(cutside);
  SideNodePlot(cutside);

  std::cout << "\nCutting Sides:\n";
  SidePlotHead();
  const plain_side_set& cuttingsides = cutside.CuttingSides();
  for (plain_side_set::const_iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
  {
    Side* cuttingside = *i;
    SidePlot(*cuttingside);
    SideNodePlot(*cuttingside);
  }

  const PointSet& selfcutpoints = cutside.SelfCutPoints();
  std::cout << "\nSelfcutpoints:\n";
  PointPlotHead();
  for (PointSet::const_iterator i = selfcutpoints.begin(); i != selfcutpoints.end(); ++i)
  {
    Point* selfcutpoint = *i;
    PointPlot(*selfcutpoint);
    std::cout << "\n";
  }

  const plain_node_set& selfcutnodes = cutside.SelfCutNodes();
  std::cout << "\nSelfcutnodes:\n";
  NodePlotHead();
  for (plain_node_set::const_iterator i = selfcutnodes.begin(); i != selfcutnodes.end(); ++i)
  {
    Node* selfcutnode = *i;
    NodePlot(*selfcutnode);
  }

  const plain_edge_set& selfcutedges = cutside.SelfCutEdges();
  std::cout << "\nSelfcutedges:\n";
  EdgePlotHead();
  for (plain_edge_set::const_iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
  {
    Edge* selfcutedge = *i;
    EdgePlot(*selfcutedge);
  }

  const std::vector<std::vector<GEO::CUT::Point*>>& selfcuttriangles = cutside.SelfCutTriangles();
  std::cout << "\nSelfcuttriangles:\n";
  for (std::vector<std::vector<GEO::CUT::Point*>>::const_iterator i = selfcuttriangles.begin();
       i != selfcuttriangles.end(); ++i)
  {
    std::vector<GEO::CUT::Point*> selfcuttriangle = *i;
    PointPlotHead();
    for (std::vector<GEO::CUT::Point*>::iterator i = selfcuttriangle.begin();
         i != selfcuttriangle.end(); ++i)
    {
      Point* selfcuttrianglepoint = *i;
      PointPlot(*selfcuttrianglepoint);
      std::cout << "\n";
    }
  }
  PlotAcross();
  j++;
}

/*-------------------------------------------------------------------------------------*
 * Status of the cutted sides for gmsh                                      wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::CuttedSideStatusGmsh(const std::string& name)
{
  std::ofstream file(name.c_str());
  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();
  int cutsidessize = cutsides.size();
  file << "View \"" << cutsidessize << " Sides\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    file.precision(16);
    int cutsidetype = i->first.size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      throw std::runtime_error("SelfCut: irregular side");
    }

    file << "){";
    Point::PointPosition cutsideselfcutposition = cutside->SelfCutPosition();
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  const std::map<plain_int_set, Teuchos::RCP<Edge>>& cutsideedges = mesh_.Edges();
  int cutedgessize = cutsideedges.size();
  file << "View \"" << cutedgessize << " Edges\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::const_iterator i = cutsideedges.begin();
       i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = &*i->second;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    cutsideedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsideedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  const std::map<int, Teuchos::RCP<Node>>& cutsidenodes = mesh_.Nodes();
  int cutnodessize = cutsidenodes.size();
  file << "View \"" << cutnodessize << " Nodes\" {\n";
  for (std::map<int, Teuchos::RCP<Node>>::const_iterator i = cutsidenodes.begin();
       i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = &*i->second;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    cutsidenode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsidenode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";
}

/*-------------------------------------------------------------------------------------*
 * Status of the cutmesh in gmsh                                            wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::WallGmsh(const std::string& name)
{
  std::stringstream str;
  str << name << "_" << myrank_ << ".pos";
  std::ofstream file(str.str().c_str());

  /*  BoundingBox elementbox;
   double x[3];
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   plain_side_set cutsides;
   for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=selfcut_sides_.begin();
   i!=selfcut_sides_.end();
   ++i )
   {
   Side * cutside = &*i->second;
   BoundingBox sidebox;
   sidebox.Assign( *cutside );
   if( sidebox.Within(1.0, elementbox) )
   {
   cutsides.insert( cutside );
   }
   }*/

  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();

  int cutsidessize = cutsides.size();
  file << "View \"" << cutsidessize << " Sides\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    file.precision(16);
    int cutsidetype = cutside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      throw std::runtime_error("SelfCut: irregular side");
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutside->SelfCutPosition();
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set cutsideedges;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<Edge*>& cutedges = cutside->Edges();
    for (std::vector<Edge*>::const_iterator i = cutedges.begin(); i != cutedges.end(); ++i)
    {
      Edge* cutsideedge = *i;
      cutsideedges.insert(cutsideedge);
    }
  }
  int cutsideedgessize = cutsideedges.size();
  file << "View \"" << cutsideedgessize << " Edges\" {\n";
  for (plain_edge_set::iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    cutsideedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsideedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set cutsidenodes;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<Node*>& cutnodes = cutside->Nodes();
    for (std::vector<Node*>::const_iterator i = cutnodes.begin(); i != cutnodes.end(); ++i)
    {
      Node* cutsidenode = *i;
      cutsidenodes.insert(cutsidenode);
    }
  }
  int cutnodessize = cutsidenodes.size();
  file << "View \"" << cutnodessize << " Nodes\" {\n";
  for (plain_node_set::iterator i = cutsidenodes.begin(); i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    cutsidenode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsidenode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_side_set cuttingsides;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const plain_side_set& selfcutsides = cutside->CuttingSides();
    for (plain_side_set::const_iterator i = selfcutsides.begin(); i != selfcutsides.end(); ++i)
    {
      Side* cuttingside = *i;
      cuttingsides.insert(cuttingside);
    }
  }
  int cuttingsidessize = cuttingsides.size();
  file << "View \"" << cuttingsidessize << " Cuttingsides\" {\n";
  for (plain_side_set::iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
  {
    Side* cuttingside = *i;
    file.precision(16);
    int cutsidetype = cuttingside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      throw std::runtime_error("SelfCut: irregular side");
    }
    file << "){";
    double cutsideselfcutposition = -0.5;
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set selfcutnodes;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const plain_node_set& selfnodes = cutside->SelfCutNodes();
    for (plain_node_set::const_iterator i = selfnodes.begin(); i != selfnodes.end(); ++i)
    {
      Node* selfcutnode = *i;
      selfcutnodes.insert(selfcutnode);
    }
  }
  int selfcutnodessize = selfcutnodes.size();
  file << "View \"" << selfcutnodessize << " Selfcutnodes\" {\n";
  for (plain_node_set::iterator i = selfcutnodes.begin(); i != selfcutnodes.end(); ++i)
  {
    Node* selfcutnode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    selfcutnode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutnode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set selfcutedges;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const plain_edge_set& selfedges = cutside->SelfCutEdges();
    for (plain_edge_set::const_iterator i = selfedges.begin(); i != selfedges.end(); ++i)
    {
      Edge* selfcutedge = *i;
      selfcutedges.insert(selfcutedge);
    }
  }
  int selfcutedgessize = selfcutedges.size();
  file << "View \"" << selfcutedgessize << " Selfcutedges\" {\n";
  for (plain_edge_set::iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
  {
    Edge* selfcutedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    selfcutedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  std::vector<std::vector<Point*>> selfcuttriangles;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<std::vector<Point*>>& selftriangles = cutside->SelfCutTriangles();
    for (std::vector<std::vector<Point*>>::const_iterator i = selftriangles.begin();
         i != selftriangles.end(); ++i)
    {
      std::vector<Point*> selfcuttriangle = *i;
      selfcuttriangles.push_back(selfcuttriangle);
    }
  }
  int selfcuttrianglessize = selfcuttriangles.size();
  file << "View \"" << selfcuttrianglessize << " Selfcuttriangles\" {\n";
  for (std::vector<std::vector<Point*>>::iterator i = selfcuttriangles.begin();
       i != selfcuttriangles.end(); ++i)
  {
    std::vector<Point*> selfcuttriangle = *i;
    int selfcuttriangleposition = 0;
    int elementfacetpointssize = selfcuttriangle.size();
    for (std::vector<Point*>::iterator i = selfcuttriangle.begin(); i != selfcuttriangle.end(); ++i)
    {
      Point* elementfacetpoint1 = *i;
      Point* elementfacetpoint2 = NULL;
      if (i + 1 != selfcuttriangle.end())
      {
        elementfacetpoint2 = *(i + 1);
      }
      else
      {
        elementfacetpoint2 = *(i + 1 - elementfacetpointssize);
      }
      file << "SL (";
      LINALG::Matrix<3, 1> elementfacetpoint1coordinates;
      elementfacetpoint1->Coordinates(elementfacetpoint1coordinates.A());
      LINALG::Matrix<3, 1> elementfacetpoint2coordinates;
      elementfacetpoint2->Coordinates(elementfacetpoint2coordinates.A());
      file << elementfacetpoint1coordinates(0, 0) << "," << elementfacetpoint1coordinates(1, 0)
           << "," << elementfacetpoint1coordinates(2, 0) << ","
           << elementfacetpoint2coordinates(0, 0) << "," << elementfacetpoint2coordinates(1, 0)
           << "," << elementfacetpoint2coordinates(2, 0);
      file << "){";
      for (int i = 0; i < 2; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << selfcuttriangleposition;
      }
      file << "};\n";
    }
  }
  file << "};\n";
}


/*-------------------------------------------------------------------------------------*
 * Status of the selfcut objects in gmsh                                    wirtz 07/16
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::SCObjectsGmsh(const std::string& name)
{
  std::ofstream file(name.c_str());

  /*  BoundingBox elementbox;
   double x[3];
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   plain_side_set cutsides;
   for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=selfcut_sides_.begin();
   i!=selfcut_sides_.end();
   ++i )
   {
   Side * cutside = &*i->second;
   BoundingBox sidebox;
   sidebox.Assign( *cutside );
   if( sidebox.Within(1.0, elementbox) )
   {
   cutsides.insert( cutside );
   }
   }*/

  const std::map<plain_int_set, Teuchos::RCP<Side>>& cutsides = mesh_.Sides();

  int cutsidessize = cutsides.size();
  file << "View \"" << cutsidessize << " Sides\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    file.precision(16);
    int cutsidetype = cutside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 4 |==" << std::endl;
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutside->SelfCutPosition();
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set cutsideedges;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<Edge*>& cutedges = cutside->Edges();
    for (std::vector<Edge*>::const_iterator i = cutedges.begin(); i != cutedges.end(); ++i)
    {
      Edge* cutsideedge = *i;
      cutsideedges.insert(cutsideedge);
    }
  }
  int cutsideedgessize = cutsideedges.size();
  file << "View \"" << cutsideedgessize << " Edges\" {\n";
  for (plain_edge_set::iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    cutsideedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsideedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set cutsidenodes;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = cutsides.begin();
       i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    const std::vector<Node*>& cutnodes = cutside->Nodes();
    for (std::vector<Node*>::const_iterator i = cutnodes.begin(); i != cutnodes.end(); ++i)
    {
      Node* cutsidenode = *i;
      cutsidenodes.insert(cutsidenode);
    }
  }
  int cutnodessize = cutsidenodes.size();
  file << "View \"" << cutnodessize << " Nodes\" {\n";
  for (plain_node_set::iterator i = cutsidenodes.begin(); i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    cutsidenode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsidenode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";



  int selfcutsidessize = selfcut_sides_.size();
  file << "View \"" << selfcutsidessize << " Sides\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::const_iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* selfcutside = &*i->second;
    file.precision(16);
    int selfcutsidetype = selfcutside->Nodes().size();
    if (selfcutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> selfcutsidecoordinates;
      selfcutside->Coordinates(selfcutsidecoordinates.A());
      for (int i = 0; i < selfcutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << selfcutsidecoordinates(0, i) << "," << selfcutsidecoordinates(1, i) << ","
             << selfcutsidecoordinates(2, i);
      }
    }
    else if (selfcutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> selfcutsidecoordinates;
      selfcutside->Coordinates(selfcutsidecoordinates.A());
      for (int i = 0; i < selfcutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << selfcutsidecoordinates(0, i) << "," << selfcutsidecoordinates(1, i) << ","
             << selfcutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 5 |==" << std::endl;
    }
    file << "){";
    Point::PointPosition selfcutsideselfcutposition = selfcutside->SelfCutPosition();
    for (int i = 0; i < selfcutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << selfcutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  int selfcutedgessize = selfcut_edges_.size();
  file << "View \"" << selfcutedgessize << " Edges\" {\n";
  for (std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator i = selfcut_edges_.begin();
       i != selfcut_edges_.end(); ++i)
  {
    Edge* selfcutedge = &*i->second;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> selfcutedgecoordinates;
    selfcutedge->Coordinates(selfcutedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << selfcutedgecoordinates(0, i) << "," << selfcutedgecoordinates(1, i) << ","
           << selfcutedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition selfcutedgeposition = selfcutedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << selfcutedgeposition;
    }
    file << "};\n";
  }
  file << "};\n";

  int selfcutnodessize = selfcut_nodes_.size();
  file << "View \"" << selfcutnodessize << " Nodes\" {\n";
  for (std::map<int, Teuchos::RCP<Node>>::iterator i = selfcut_nodes_.begin();
       i != selfcut_nodes_.end(); ++i)
  {
    Node* selfcutnode = &*i->second;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> selfcutnodecoordinates;
    selfcutnode->Coordinates(selfcutnodecoordinates.A());
    file << selfcutnodecoordinates(0, 0) << "," << selfcutnodecoordinates(1, 0) << ","
         << selfcutnodecoordinates(2, 0);
    file << "){";
    Point::PointPosition selfcutnodeposition = selfcutnode->SelfCutPosition();
    file << selfcutnodeposition;
    file << "};\n";
  }
  file << "};\n";
}


/*-------------------------------------------------------------------------------------*
 * Status of the cutmesh in gmsh for my CMGM (old)                          wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::SCmgmGmsh(const std::string& name)
{
  std::ofstream file(name.c_str());

  /*  BoundingBox elementbox;
   double x[3];
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   plain_side_set cutsides;
   for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=selfcut_sides_.begin();
   i!=selfcut_sides_.end();
   ++i )
   {
   Side * cutside = &*i->second;
   BoundingBox sidebox;
   sidebox.Assign( *cutside );
   if( sidebox.Within(1.0, elementbox) )
   {
   cutsides.insert( cutside );
   }
   }*/
  plain_side_set cutsides;
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    cutsides.insert(cutside);
  }
  int cutsidessize = cutsides.size();
  file << "View \"" << cutsidessize << " Sides\" {\n";
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    file.precision(16);
    int cutsidetype = cutside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 6 |==" << std::endl;
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutside->SelfCutPosition();
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set cutsideedges;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<Edge*>& cutedges = cutside->Edges();
    for (std::vector<Edge*>::const_iterator i = cutedges.begin(); i != cutedges.end(); ++i)
    {
      Edge* cutsideedge = *i;
      cutsideedges.insert(cutsideedge);
    }
  }
  int cutsideedgessize = cutsideedges.size();
  file << "View \"" << cutsideedgessize << " Edges\" {\n";
  for (plain_edge_set::iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    cutsideedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsideedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set cutsidenodes;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<Node*>& cutnodes = cutside->Nodes();
    for (std::vector<Node*>::const_iterator i = cutnodes.begin(); i != cutnodes.end(); ++i)
    {
      Node* cutsidenode = *i;
      cutsidenodes.insert(cutsidenode);
    }
  }
  int cutnodessize = cutsidenodes.size();
  file << "View \"" << cutnodessize << " Nodes\" {\n";
  for (plain_node_set::iterator i = cutsidenodes.begin(); i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    cutsidenode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsidenode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_side_set cuttingsides;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_side_set& selfcutsides = cutside->CuttingSides();
    for (plain_side_set::const_iterator i = selfcutsides.begin(); i != selfcutsides.end(); ++i)
    {
      Side* cuttingside = *i;
      cuttingsides.insert(cuttingside);
    }
  }
  int cuttingsidessize = cuttingsides.size();
  file << "View \"" << cuttingsidessize << " Cuttingsides\" {\n";
  for (plain_side_set::iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
  {
    Side* cuttingside = *i;
    file.precision(16);
    int cutsidetype = cuttingside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 7 |==" << std::endl;
    }
    file << "){";
    double cutsideselfcutposition = -0.5;
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set selfcutnodes;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_node_set& selfnodes = cutside->SelfCutNodes();
    for (plain_node_set::const_iterator i = selfnodes.begin(); i != selfnodes.end(); ++i)
    {
      Node* selfcutnode = *i;
      selfcutnodes.insert(selfcutnode);
    }
  }
  int selfcutnodessize = selfcutnodes.size();
  file << "View \"" << selfcutnodessize << " Selfcutnodes\" {\n";
  for (plain_node_set::iterator i = selfcutnodes.begin(); i != selfcutnodes.end(); ++i)
  {
    Node* selfcutnode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    selfcutnode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutnode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set selfcutedges;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_edge_set& selfedges = cutside->SelfCutEdges();
    for (plain_edge_set::const_iterator i = selfedges.begin(); i != selfedges.end(); ++i)
    {
      Edge* selfcutedge = *i;
      selfcutedges.insert(selfcutedge);
    }
  }
  int selfcutedgessize = selfcutedges.size();
  file << "View \"" << selfcutedgessize << " Selfcutedges\" {\n";
  for (plain_edge_set::iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
  {
    Edge* selfcutedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    selfcutedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  std::vector<std::vector<Point*>> selfcuttriangles;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<std::vector<Point*>>& selftriangles = cutside->SelfCutTriangles();
    for (std::vector<std::vector<Point*>>::const_iterator i = selftriangles.begin();
         i != selftriangles.end(); ++i)
    {
      std::vector<Point*> selfcuttriangle = *i;
      selfcuttriangles.push_back(selfcuttriangle);
    }
  }
  int selfcuttrianglessize = selfcuttriangles.size();
  file << "View \"" << selfcuttrianglessize << " Selfcuttriangles\" {\n";
  for (std::vector<std::vector<Point*>>::iterator i = selfcuttriangles.begin();
       i != selfcuttriangles.end(); ++i)
  {
    std::vector<Point*> selfcuttriangle = *i;
    int selfcuttriangleposition = 0;
    int elementfacetpointssize = selfcuttriangle.size();
    for (std::vector<Point*>::iterator i = selfcuttriangle.begin(); i != selfcuttriangle.end(); ++i)
    {
      Point* elementfacetpoint1 = *i;
      Point* elementfacetpoint2 = NULL;
      if (i + 1 != selfcuttriangle.end())
      {
        elementfacetpoint2 = *(i + 1);
      }
      else
      {
        elementfacetpoint2 = *(i + 1 - elementfacetpointssize);
      }
      file << "SL (";
      LINALG::Matrix<3, 1> elementfacetpoint1coordinates;
      elementfacetpoint1->Coordinates(elementfacetpoint1coordinates.A());
      LINALG::Matrix<3, 1> elementfacetpoint2coordinates;
      elementfacetpoint2->Coordinates(elementfacetpoint2coordinates.A());
      file << elementfacetpoint1coordinates(0, 0) << "," << elementfacetpoint1coordinates(1, 0)
           << "," << elementfacetpoint1coordinates(2, 0) << ","
           << elementfacetpoint2coordinates(0, 0) << "," << elementfacetpoint2coordinates(1, 0)
           << "," << elementfacetpoint2coordinates(2, 0);
      file << "){";
      for (int i = 0; i < 2; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << selfcuttriangleposition;
      }
      file << "};\n";
    }
  }
  file << "};\n";
}

/*-------------------------------------------------------------------------------------*
 * Status of all sides in separates files for gmsh                          wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::AllSingleGmsh(const std::string& location)
{
  for (std::map<plain_int_set, Teuchos::RCP<Side>>::iterator i = selfcut_sides_.begin();
       i != selfcut_sides_.end(); ++i)
  {
    Side* cutside = &*i->second;
    std::stringstream filenamegmsh;

    plain_int_set::const_iterator it = (i->first).begin();
    const int first = *it;
    it++;
    const int second = *it;
    it++;
    const int third = *it;

    filenamegmsh << location << "_SelfCutSide_" << first << "_" << second << "_" << third << ".pos";
    ErrorGmsh(filenamegmsh.str(), *cutside);
  }
}

/*-------------------------------------------------------------------------------------*
 * Status of one problematic side for gmsh                                  wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::ErrorGmsh(const std::string& name, Side& cutside)
{
  std::ofstream file(name.c_str());

  /*  BoundingBox elementbox;
   double x[3];
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.52;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.035;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = -0.005;
   elementbox.AddPoint( x );
   x[0] = 0.53;
   x[1] = 0.045;
   x[2] = 0.005;
   elementbox.AddPoint( x );
   plain_side_set cutsides;
   for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=selfcut_sides_.begin();
   i!=selfcut_sides_.end();
   ++i )
   {
   Side * cutside = &*i->second;
   BoundingBox sidebox;
   sidebox.Assign( *cutside );
   if( sidebox.Within(1.0, elementbox) )
   {
   cutsides.insert( cutside );
   }
   }*/
  plain_side_set cutsides;
  /*  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=selfcut_sides_.begin();
   i!=selfcut_sides_.end();
   ++i )
   {
   Side * cutside = &*i->second;
   cutsides.insert( cutside );
   }*/
  cutsides.insert(&cutside);

  int cutsidessize = cutsides.size();
  file << "View \"" << cutsidessize << " Sides\" {\n";
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    file.precision(16);
    int cutsidetype = cutside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cutside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 8 |==" << std::endl;
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutside->SelfCutPosition();
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set cutsideedges;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<Edge*>& cutedges = cutside->Edges();
    for (std::vector<Edge*>::const_iterator i = cutedges.begin(); i != cutedges.end(); ++i)
    {
      Edge* cutsideedge = *i;
      cutsideedges.insert(cutsideedge);
    }
  }
  int cutsideedgessize = cutsideedges.size();
  file << "View \"" << cutsideedgessize << " Edges\" {\n";
  for (plain_edge_set::iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    cutsideedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsideedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set cutsidenodes;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<Node*>& cutnodes = cutside->Nodes();
    for (std::vector<Node*>::const_iterator i = cutnodes.begin(); i != cutnodes.end(); ++i)
    {
      Node* cutsidenode = *i;
      cutsidenodes.insert(cutsidenode);
    }
  }
  int cutnodessize = cutsidenodes.size();
  file << "View \"" << cutnodessize << " Nodes\" {\n";
  for (plain_node_set::iterator i = cutsidenodes.begin(); i != cutsidenodes.end(); ++i)
  {
    Node* cutsidenode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    cutsidenode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = cutsidenode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_side_set cuttingsides;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_side_set& selfcutsides = cutside->CuttingSides();
    for (plain_side_set::const_iterator i = selfcutsides.begin(); i != selfcutsides.end(); ++i)
    {
      Side* cuttingside = *i;
      cuttingsides.insert(cuttingside);
    }
  }
  int cuttingsidessize = cuttingsides.size();
  file << "View \"" << cuttingsidessize << " Cuttingsides\" {\n";
  for (plain_side_set::iterator i = cuttingsides.begin(); i != cuttingsides.end(); ++i)
  {
    Side* cuttingside = *i;
    file.precision(16);
    int cutsidetype = cuttingside->Nodes().size();
    if (cutsidetype == 3)
    {
      file << "ST (";
      LINALG::Matrix<3, 3> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else if (cutsidetype == 4)
    {
      file << "SQ (";
      LINALG::Matrix<3, 4> cutsidecoordinates;
      cuttingside->Coordinates(cutsidecoordinates.A());
      for (int i = 0; i < cutsidetype; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << cutsidecoordinates(0, i) << "," << cutsidecoordinates(1, i) << ","
             << cutsidecoordinates(2, i);
      }
    }
    else
    {
      // throw std::runtime_error("SelfCut: irregular side");
      std::cout << "==| WARNING: SelfCut: irregular side 9 |==" << std::endl;
    }
    file << "){";
    double cutsideselfcutposition = -0.5;
    for (int i = 0; i < cutsidetype; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  plain_node_set selfcutnodes;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_node_set& selfnodes = cutside->SelfCutNodes();
    for (plain_node_set::const_iterator i = selfnodes.begin(); i != selfnodes.end(); ++i)
    {
      Node* selfcutnode = *i;
      selfcutnodes.insert(selfcutnode);
    }
  }
  int selfcutnodessize = selfcutnodes.size();
  file << "View \"" << selfcutnodessize << " Selfcutnodes\" {\n";
  for (plain_node_set::iterator i = selfcutnodes.begin(); i != selfcutnodes.end(); ++i)
  {
    Node* selfcutnode = *i;
    file.precision(16);
    file << "SP (";
    LINALG::Matrix<3, 1> cutsidenodecoordinates;
    selfcutnode->Coordinates(cutsidenodecoordinates.A());
    file << cutsidenodecoordinates(0, 0) << "," << cutsidenodecoordinates(1, 0) << ","
         << cutsidenodecoordinates(2, 0);
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutnode->SelfCutPosition();
    file << cutsideselfcutposition;
    file << "};\n";
  }
  file << "};\n";

  plain_edge_set selfcutedges;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const plain_edge_set& selfedges = cutside->SelfCutEdges();
    for (plain_edge_set::const_iterator i = selfedges.begin(); i != selfedges.end(); ++i)
    {
      Edge* selfcutedge = *i;
      selfcutedges.insert(selfcutedge);
    }
  }
  int selfcutedgessize = selfcutedges.size();
  file << "View \"" << selfcutedgessize << " Selfcutedges\" {\n";
  for (plain_edge_set::iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
  {
    Edge* selfcutedge = *i;
    file.precision(16);
    file << "SL (";
    LINALG::Matrix<3, 2> cutsideedgecoordinates;
    selfcutedge->Coordinates(cutsideedgecoordinates.A());
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideedgecoordinates(0, i) << "," << cutsideedgecoordinates(1, i) << ","
           << cutsideedgecoordinates(2, i);
    }
    file << "){";
    Point::PointPosition cutsideselfcutposition = selfcutedge->SelfCutPosition();
    for (int i = 0; i < 2; ++i)
    {
      if (i > 0)
      {
        file << ", ";
      }
      file << cutsideselfcutposition;
    }
    file << "};\n";
  }
  file << "};\n";

  std::vector<std::vector<Point*>> selfcuttriangles;
  for (plain_side_set::iterator i = cutsides.begin(); i != cutsides.end(); ++i)
  {
    Side* cutside = *i;
    const std::vector<std::vector<Point*>>& selftriangles = cutside->SelfCutTriangles();
    for (std::vector<std::vector<Point*>>::const_iterator i = selftriangles.begin();
         i != selftriangles.end(); ++i)
    {
      std::vector<Point*> selfcuttriangle = *i;
      selfcuttriangles.push_back(selfcuttriangle);
    }
  }
  int selfcuttrianglessize = selfcuttriangles.size();
  file << "View \"" << selfcuttrianglessize << " Selfcuttriangles\" {\n";
  for (std::vector<std::vector<Point*>>::iterator i = selfcuttriangles.begin();
       i != selfcuttriangles.end(); ++i)
  {
    std::vector<Point*> selfcuttriangle = *i;
    int selfcuttriangleposition = 0;
    int elementfacetpointssize = selfcuttriangle.size();
    for (std::vector<Point*>::iterator i = selfcuttriangle.begin(); i != selfcuttriangle.end(); ++i)
    {
      Point* elementfacetpoint1 = *i;
      Point* elementfacetpoint2 = NULL;
      if (i + 1 != selfcuttriangle.end())
      {
        elementfacetpoint2 = *(i + 1);
      }
      else
      {
        elementfacetpoint2 = *(i + 1 - elementfacetpointssize);
      }
      file << "SL (";
      LINALG::Matrix<3, 1> elementfacetpoint1coordinates;
      elementfacetpoint1->Coordinates(elementfacetpoint1coordinates.A());
      LINALG::Matrix<3, 1> elementfacetpoint2coordinates;
      elementfacetpoint2->Coordinates(elementfacetpoint2coordinates.A());
      file << elementfacetpoint1coordinates(0, 0) << "," << elementfacetpoint1coordinates(1, 0)
           << "," << elementfacetpoint1coordinates(2, 0) << ","
           << elementfacetpoint2coordinates(0, 0) << "," << elementfacetpoint2coordinates(1, 0)
           << "," << elementfacetpoint2coordinates(2, 0);
      file << "){";
      for (int i = 0; i < 2; ++i)
      {
        if (i > 0)
        {
          file << ", ";
        }
        file << selfcuttriangleposition;
      }
      file << "};\n";
    }
  }
  file << "};\n";
}

/*-------------------------------------------------------------------------------------*
 * cuts the edges of the first cutside with the other cutside               wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::PerformSelfCut(Side& cutside, Side& otherside, PointSet& selfcutpoints)
{
  const std::vector<Edge*>& cutsideedges = cutside.Edges();
  for (std::vector<Edge*>::const_iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    PointSet sidepairselfcutpoints;
    cutsideedge->FindCutPointsMeshCut(mesh_, NULL, cutside, otherside, &sidepairselfcutpoints);
    for (PointSet::iterator i = sidepairselfcutpoints.begin(); i != sidepairselfcutpoints.end();
         ++i)
    {
      Point* edgeselfcutpoint = *i;
      cutsideedge->AddPoint(edgeselfcutpoint);
      selfcutpoints.insert(edgeselfcutpoint);
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * gets the edges of the new created sides                                  wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::GetSelfCutEdges(Side& cutside)
{
  const std::vector<Edge*>& cutsideedges = cutside.Edges();
  for (std::vector<Edge*>::const_iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    plain_int_set cutsideedgenodeids;
    const std::vector<Node*>& cutsideedgenodes = cutsideedge->Nodes();
    for (std::vector<Node*>::const_iterator i = cutsideedgenodes.begin();
         i != cutsideedgenodes.end(); ++i)
    {
      Node* cutsideedgenode = *i;
      int cutsideedgenodeid = cutsideedgenode->Id();
      cutsideedgenodeids.insert(cutsideedgenodeid);
    }

    std::map<plain_int_set, Teuchos::RCP<Edge>>::iterator edgeiterator =
        selfcut_edges_.find(cutsideedgenodeids);
    if (edgeiterator == selfcut_edges_.end())
    {
      const std::map<plain_int_set, Teuchos::RCP<Edge>>& cutsideedgercp = mesh_.Edges();
      std::map<plain_int_set, Teuchos::RCP<Edge>>::const_iterator edgeiterator =
          cutsideedgercp.find(cutsideedgenodeids);
      selfcut_edges_[cutsideedgenodeids] = edgeiterator->second;
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * checks if the direction of rotation of the cycle is correct              wirtz 05/13
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::SelfCut::CheckNormal(std::vector<double> cutsideplane, Cycle& maincycle)
{
  std::vector<Point*> maincyclepoints = maincycle();
  std::vector<double> maincycleplane = KERNEL::EqnPlaneOfPolygon(maincyclepoints);
  if (fabs(cutsideplane[0]) > TOL_EQN_PLANE and cutsideplane[0] * maincycleplane[0] < 0.0)
  {
    return false;
  }
  else if (fabs(cutsideplane[1]) > TOL_EQN_PLANE and cutsideplane[1] * maincycleplane[1] < 0.0)
  {
    return false;
  }
  else if (fabs(cutsideplane[2]) > TOL_EQN_PLANE and cutsideplane[2] * maincycleplane[2] < 0.0)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/*-------------------------------------------------------------------------------------*
 * deletes all pointers which are pointing to a soon to be erased side      wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseSidePointer(Side& cutside, bool erase_nodepointers)
{
  const std::vector<Edge*>& cutsideedges = cutside.Edges();
  for (std::vector<Edge*>::const_iterator i = cutsideedges.begin(); i != cutsideedges.end(); ++i)
  {
    Edge* cutsideedge = *i;
    cutsideedge->EraseCutSide(&cutside);
  }
  const PointSet& cutsidepoints = cutside.CutPoints();
  for (PointSet::const_iterator i = cutsidepoints.begin(); i != cutsidepoints.end(); ++i)
  {
    Point* cutsidepoint = *i;
    if (erase_nodepointers) cutsidepoint->EraseCutSide(&cutside);
    cutsidepoint->ErasedContainingCutPairs(&cutside);
  }
}

/*-------------------------------------------------------------------------------------*
 * erases a side                                                            wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseSide(std::vector<plain_int_set>& cutsideids)
{
  for (std::vector<plain_int_set>::iterator i = cutsideids.begin(); i != cutsideids.end(); ++i)
  {
    const plain_int_set& cutsideid = *i;
    selfcut_sides_.erase(cutsideid);
    meshhandle_.RemoveSubSide(&(*mesh_.Sides().at(cutsideid)));
    mesh_.EraseSide(cutsideid);
  }
}

void GEO::CUT::SelfCut::EraseInsideSide(std::vector<plain_int_set>& cutsideids)
{
  for (std::vector<plain_int_set>::iterator i = cutsideids.begin(); i != cutsideids.end(); ++i)
  {
    const plain_int_set& cutsideid = *i;
    selfcut_sides_.erase(cutsideid);
    meshhandle_.MarkSubSideasUnphysical(&(*mesh_.Sides().at(cutsideid)));
    mesh_.MoveSidetoStorage(cutsideid);
  }
}

/*-------------------------------------------------------------------------------------*
 * deletes all pointers which are pointing to a soon to be erased edge      wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseEdgePointer(Edge& cutsideedge)
{
  const std::vector<Node*>& cutsideedgenodes = cutsideedge.Nodes();
  for (std::vector<Node*>::const_iterator i = cutsideedgenodes.begin(); i != cutsideedgenodes.end();
       ++i)
  {
    Node* cutsideedgenode = *i;
    cutsideedgenode->EraseCutSideEdge(&cutsideedge);
  }
  const PointPositionSet& cutsideedgepoints = cutsideedge.CutPoints();
  for (PointPositionSet::const_iterator i = cutsideedgepoints.begin(); i != cutsideedgepoints.end();
       ++i)
  {
    Point* cutsideedgepoint = *i;
    cutsideedgepoint->EraseCutSideEdge(&cutsideedge);
    cutsideedgepoint->ErasedContainingCutPairs(&cutsideedge);
  }
}

/*-------------------------------------------------------------------------------------*
 * erases an edge                                                           wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::SelfCut::EraseEdge(std::vector<plain_int_set>& cutsideedgeids)
{
  for (std::vector<plain_int_set>::iterator i = cutsideedgeids.begin(); i != cutsideedgeids.end();
       ++i)
  {
    const plain_int_set& cutsideedgeid = *i;
    selfcut_edges_.erase(cutsideedgeid);
    mesh_.EraseEdge(cutsideedgeid);
  }
}

void GEO::CUT::SelfCut::PlotAcross()
{
  std::cout << "\n================================================================================="
               "===================================\n";
}

void GEO::CUT::SelfCut::PointPlotHead()
{
  std::cout << "x\ty\tz\t# pid\n"
            << "-----------------------------\n";
}

void GEO::CUT::SelfCut::NodePlotHead()
{
  std::cout << "x\ty\tz\t# pid\t# nid\tSelfCutPosition\n"
            << "-------------------------------------------------------\n";
}

void GEO::CUT::SelfCut::EdgePlotHead()
{
  for (int i = 1; i < 3; ++i)
  {
    std::cout << "# nid_" << i << "\t";
  }
  std::cout << "SelfCutPosition\n"
            << "-------------------------------\n";
}

void GEO::CUT::SelfCut::SidePlotHead()
{
  for (int i = 1; i < 4; ++i)
  {
    std::cout << "# nid_" << i << "\t";
  }
  std::cout << "# sid\tSelfCutPosition\n"
            << "-----------------------------------------------\n";
}

void GEO::CUT::SelfCut::PointPlot(Point& point)
{
  LINALG::Matrix<3, 1> pointcoordinates;
  point.Coordinates(pointcoordinates.A());
  std::cout << std::setprecision(16) << pointcoordinates(0) << "\t" << std::setprecision(16)
            << pointcoordinates(1) << "\t" << std::setprecision(16) << pointcoordinates(2) << "\t# "
            << point.Id();
}

void GEO::CUT::SelfCut::NodePlot(Node& node)
{
  PointPlot(*node.point());
  std::cout << "\t# " << node.Id() << "\t" << node.SelfCutPosition() << "\n";
}

void GEO::CUT::SelfCut::EdgePlot(Edge& edge)
{
  const std::vector<Node*>& edgenodes = edge.Nodes();
  for (unsigned int i = 0; i < edgenodes.size(); ++i)
  {
    std::cout << "# " << edgenodes[i]->Id() << "\t";
  }
  std::cout << edge.SelfCutPosition() << "\n";
}

void GEO::CUT::SelfCut::SidePlot(Side& side)
{
  const std::vector<Node*>& sidenodes = side.Nodes();
  for (unsigned int i = 0; i < sidenodes.size(); ++i)
  {
    std::cout << "# " << sidenodes[i]->Id() << "\t";
  }
  std::cout << "# " << side.Id() << "\t" << side.SelfCutPosition() << "\n";
}

void GEO::CUT::SelfCut::SideNodePlot(Side& side)
{
  const std::vector<Node*>& sidenodes = side.Nodes();
  NodePlotHead();
  for (std::vector<Node*>::const_iterator i = sidenodes.begin(); i != sidenodes.end(); ++i)
  {
    Node* sidenode = *i;
    NodePlot(*sidenode);
  }
}
