/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of all sorts of conditions on discretization

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_legacy_enum_definitions_problem_type.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Build boundary condition geometries (public)             mwgee 01/07|
 *----------------------------------------------------------------------*/
void Core::FE::Discretization::boundary_conditions_geometry()
{
  // As a first step we delete ALL references to any conditions
  // in the discretization
  for (int i = 0; i < num_my_col_nodes(); ++i) l_col_node(i)->clear_conditions();
  for (int i = 0; i < num_my_col_elements(); ++i) l_col_element(i)->clear_conditions();

  // now we delete all old geometries that are attached to any conditions
  // and set a communicator to the condition
  for (auto& [name, condition] : condition_)
  {
    condition->clear_geometry();
  }

  // for all conditions, we set a ptr in the nodes to the condition
  for (auto& [name, condition] : condition_)
  {
    const std::vector<int>* nodes = condition->get_nodes();
    // There might be conditions that do not have a nodal cloud
    if (!nodes) continue;
    for (int node : *nodes)
    {
      if (!node_col_map()->MyGID(node)) continue;
      Core::Nodes::Node* actnode = g_node(node);
      if (!actnode) FOUR_C_THROW("Cannot find global node");
      actnode->set_condition(name, condition);
    }
  }

  // create a map that holds the overall number of created elements
  // associated with a specific condition type
  std::map<std::string, int> numele;

  // Loop all conditions and build geometry description if desired
  for (auto& [name, condition] : condition_)
  {
    // flag whether new elements have been created for this condition
    bool havenewelements = false;

    // do not build geometry description for this condition
    if (!condition->geometry_description()) continue;
    // do not build geometry description for this condition
    else if (condition->g_type() == Core::Conditions::geometry_type_no_geom)
      continue;
    // do not build anything for point wise conditions
    else if (condition->g_type() == Core::Conditions::geometry_type_point)
      continue;
    // build element geometry description without creating new elements if the
    // condition is not a boundary condition; this would be:
    //  - line conditions in 1D
    //  - surface conditions in 2D
    //  - volume conditions in 3D
    else if (Core::Conditions::geometry_type_to_dim.at(condition->g_type()) == n_dim_)
      havenewelements = build_volumesin_condition(name, condition);
    // dimension of condition must not larger than the one of the problem itself
    else if (Core::Conditions::geometry_type_to_dim.at(condition->g_type()) > n_dim_)
      FOUR_C_THROW("Dimension of condition is larger than the problem dimension.");
    // build a line element geometry description
    else if (condition->g_type() == Core::Conditions::geometry_type_line)
      havenewelements = build_linesin_condition(name, condition);
    // build a surface element geometry description
    else if (condition->g_type() == Core::Conditions::geometry_type_surface)
      havenewelements = build_surfacesin_condition(name, condition);
    // this should be it. if not: FOUR_C_THROW.
    else
      FOUR_C_THROW("Somehow the condition geometry does not fit to the problem dimension.");

    if (havenewelements)
    {
      // determine the local number of created elements associated with
      // the active condition
      int localcount = 0;
      for (const auto& [ele_id, ele] : condition->geometry())
      {
        // do not count ghosted elements
        if (ele->owner() == get_comm().MyPID())
        {
          localcount += 1;
        }
      }

      // determine the global number of created elements associated with
      // the active condition
      int count;
      get_comm().SumAll(&localcount, &count, 1);

      if (numele.find(name) == numele.end())
      {
        numele[name] = 0;
      }

      // adjust the IDs of the elements associated with the active
      // condition in order to obtain unique IDs within one condition type
      condition->adjust_id(numele[name]);

      // adjust the number of elements associated with the current
      // condition type
      numele[name] += count;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::Discretization::assign_global_i_ds(const Epetra_Comm& comm,
    const std::map<std::vector<int>, Teuchos::RCP<Core::Elements::Element>>& elementmap,
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& finalgeometry)
{
  // pack elements on all processors

  int size = 0;
  for (const auto& [nodes, ele] : elementmap)
  {
    size += nodes.size() + 1;
  }
  std::vector<int> sendblock;
  sendblock.reserve(size);
  for (const auto& [nodes, ele] : elementmap)
  {
    sendblock.push_back(nodes.size());
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(sendblock));
  }

  // communicate elements to processor 0

  int mysize = sendblock.size();
  comm.SumAll(&mysize, &size, 1);
  int mypos = Core::LinAlg::FindMyPos(sendblock.size(), comm);

  std::vector<int> send(size);
  std::fill(send.begin(), send.end(), 0);
  std::copy(sendblock.begin(), sendblock.end(), &send[mypos]);
  sendblock.clear();
  std::vector<int> recv(size);
  comm.SumAll(send.data(), recv.data(), size);

  send.clear();

  // unpack, unify and sort elements on processor 0

  if (comm.MyPID() == 0)
  {
    std::set<std::vector<int>> elements;
    int index = 0;
    while (index < static_cast<int>(recv.size()))
    {
      int esize = recv[index];
      index += 1;
      std::vector<int> element;
      element.reserve(esize);
      std::copy(&recv[index], &recv[index + esize], std::back_inserter(element));
      index += esize;
      elements.insert(element);
    }
    recv.clear();

    // pack again to distribute pack to all processors

    send.reserve(index);
    for (const auto& ele : elements)
    {
      send.push_back(ele.size());
      std::copy(ele.begin(), ele.end(), std::back_inserter(send));
    }
    size = send.size();
  }
  else
  {
    recv.clear();
  }

  // broadcast sorted elements to all processors

  comm.Broadcast(&size, 1, 0);
  send.resize(size);
  comm.Broadcast(send.data(), send.size(), 0);

  // Unpack sorted elements. Take element position for gid.

  int index = 0;
  int gid = 0;
  while (index < static_cast<int>(send.size()))
  {
    int esize = send[index];
    index += 1;
    std::vector<int> element;
    element.reserve(esize);
    std::copy(&send[index], &send[index + esize], std::back_inserter(element));
    index += esize;

    // set gid to my elements
    auto iter = elementmap.find(element);
    if (iter != elementmap.end())
    {
      iter->second->set_id(gid);
      finalgeometry[gid] = iter->second;
    }

    gid += 1;
  }
}  // assign_global_i_ds


/*----------------------------------------------------------------------*
 |  Build line geometry in a condition (public)              mwgee 01/07|
 *----------------------------------------------------------------------*/
/* Hopefully improved by Heiner (h.kue 09/07) */
bool Core::FE::Discretization::build_linesin_condition(
    const std::string& name, Teuchos::RCP<Core::Conditions::Condition> cond)
{
  /* First: Create the line objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->get_nodes();
  if (!nodeids) FOUR_C_THROW("Cannot find array 'Node Ids' in condition");

  // ptrs to my row/column nodes of those
  std::map<int, Core::Nodes::Node*> colnodes;

  for (const auto& nodeid : *nodeids)
  {
    if (node_col_map()->MyGID(nodeid))
    {
      Core::Nodes::Node* actnode = g_node(nodeid);
      if (!actnode) FOUR_C_THROW("Cannot find global node");
      colnodes[actnode->id()] = actnode;
    }
  }

  // map of lines in our cloud: (node_ids) -> line
  std::map<std::vector<int>, Teuchos::RCP<Core::Elements::Element>> linemap;
  // loop these nodes and build all lines attached to them
  for (const auto& [node_id, actnode] : colnodes)
  {
    // loop all elements attached to actnode
    Core::Elements::Element** elements = actnode->elements();
    for (int i = 0; i < actnode->num_element(); ++i)
    {
      // loop all lines of all elements attached to actnode
      const int numlines = elements[i]->num_line();
      if (!numlines) continue;
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines = elements[i]->lines();
      if (lines.size() == 0) FOUR_C_THROW("Element returned no lines");
      for (int j = 0; j < numlines; ++j)
      {
        Teuchos::RCP<Core::Elements::Element> actline = lines[j];
        // find lines that are attached to actnode
        const int nnodeperline = actline->num_node();
        Core::Nodes::Node** nodesperline = actline->nodes();
        if (!nodesperline) FOUR_C_THROW("Line returned no nodes");
        for (int k = 0; k < nnodeperline; ++k)
        {
          if (nodesperline[k]->id() == actnode->id())
          {
            // line is attached to actnode
            // see whether all nodes on the line are in our nodal cloud
            bool allin = true;
            for (int l = 0; l < nnodeperline; ++l)
            {
              if (colnodes.find(nodesperline[l]->id()) == colnodes.end())
              {
                allin = false;
                break;
              }
            }  // for (int l=0; l<nnodeperline; ++l)
            // if all nodes on line are in our cloud, add line
            if (allin)
            {
              std::vector<int> nodes(actline->num_node());
              transform(actline->nodes(), actline->nodes() + actline->num_node(), nodes.begin(),
                  std::mem_fn(&Core::Nodes::Node::id));
              sort(nodes.begin(), nodes.end());

              if (linemap.find(nodes) == linemap.end())
              {
                Teuchos::RCP<Core::Elements::Element> line = Teuchos::rcp(actline->clone());
                // Set owning process of line to node with smallest gid.
                line->set_owner(elements[i]->owner());
                linemap[nodes] = line;
              }
            }
            break;
          }  // if (nodesperline[k] == actnode)
        }
      }  // for (int j=0; j<numlines; ++j)
    }    // for (int i=0; i<actnode->NumElement(); ++i)
  }      // for (fool=nodes.begin(); fool != nodes.end(); ++fool)


  // Lines be added to the condition: (line_id) -> (line).
  auto finallines = Teuchos::rcp(new std::map<int, Teuchos::RCP<Core::Elements::Element>>());

  assign_global_i_ds(get_comm(), linemap, *finallines);

  cond->add_geometry(finallines);

  // elements where created that need new unique ids
  bool havenewelements = true;
  // note: this seems useless since this build_linesin_condition always creates
  //       elements. However, this function is overloaded in
  //       MeshfreeDiscretization where it does not necessarily build elements.
  return havenewelements;

}  // Core::FE::Discretization::build_linesin_condition


/*----------------------------------------------------------------------*
 |  Build surface geometry in a condition (public)          rauch 10/16 |
 *----------------------------------------------------------------------*/
bool Core::FE::Discretization::build_surfacesin_condition(
    const std::string& name, Teuchos::RCP<Core::Conditions::Condition> cond)
{
  // these conditions are special since associated volume conditions also need
  // to be considered.
  // we want to allow building surfaces where two volume elements which belong to
  // the same discretization share a common surface. the condition surface element however,
  // is associated to only one of the two volume elements.
  // these volume elements are inserted into VolEleIDs via the method find_associated_ele_i_ds.
  std::set<int> VolEleIDs;
  if (cond->type() == Core::Conditions::StructFluidSurfCoupling)
  {
    if ((cond->parameters().get<std::string>("field")) == "structure")
    {
      find_associated_ele_i_ds(cond, VolEleIDs, "StructFluidVolCoupling");
    }
  }
  else if (cond->type() == Core::Conditions::RedAirwayTissue)
    find_associated_ele_i_ds(cond, VolEleIDs, "StructFluidVolCoupling");

  /* First: Create the surface objects that belong to the condition. */

  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->get_nodes();
  if (!nodeids) FOUR_C_THROW("Cannot find array 'Node Ids' in condition");

  // ptrs to my row/column nodes of those
  std::map<int, Core::Nodes::Node*> myrownodes;
  std::map<int, Core::Nodes::Node*> mycolnodes;
  for (const auto& nodeid : *nodeids)
  {
    if (node_col_map()->MyGID(nodeid))
    {
      Core::Nodes::Node* actnode = g_node(nodeid);
      if (!actnode) FOUR_C_THROW("Cannot find global node");
      mycolnodes[actnode->id()] = actnode;
    }
    if (node_row_map()->MyGID(nodeid))
    {
      Core::Nodes::Node* actnode = g_node(nodeid);
      if (!actnode) FOUR_C_THROW("Cannot find global node");
      myrownodes[actnode->id()] = actnode;
    }
  }

  // map of surfaces in this cloud: (node_ids) -> (surface)
  std::map<std::vector<int>, Teuchos::RCP<Core::Elements::Element>> surfmap;

  // loop these column nodes and build all surfs attached to them.
  // we have to loop over column nodes because it can happen that
  // we want to create a surface on an element face which has only
  // ghosted nodes. if the resulting condition geometry is then used for
  // cloning a discretization from it, we copy the condition to the
  // cloned surface discretization, and if we build the geometry of
  // this surface discretization, this way we make sure that we do not
  // miss a surface element. otherwise, we would miss the surface element
  // of the face which has now only ghosted nodes because we would only
  // look at row nodes but the considered face has no row node.
  for (const auto& [node_id, actnode] : mycolnodes)
  {
    Core::Elements::Element** elements = actnode->elements();

    // loop all elements attached to actnode
    for (int i = 0; i < actnode->num_element(); ++i)
    {
      // special treatment of RedAirwayTissue and StructFluidVolCoupling
      if (VolEleIDs.size())
        if (VolEleIDs.find(elements[i]->id()) == VolEleIDs.end()) continue;

      // loop all surfaces of all elements attached to actnode
      const int numsurfs = elements[i]->num_surface();
      if (!numsurfs) continue;
      std::vector<Teuchos::RCP<Core::Elements::Element>> surfs = elements[i]->surfaces();
      if (surfs.size() == 0) FOUR_C_THROW("Element does not return any surfaces");

      // loop all surfaces of all elements attached to actnode
      for (int j = 0; j < numsurfs; ++j)
      {
        Teuchos::RCP<Core::Elements::Element> actsurf = surfs[j];
        // find surfs attached to actnode
        const int nnodepersurf = actsurf->num_node();
        Core::Nodes::Node** nodespersurf = actsurf->nodes();
        if (!nodespersurf) FOUR_C_THROW("Surface returned no nodes");
        for (int k = 0; k < nnodepersurf; ++k)
        {
          if (nodespersurf[k]->id() == actnode->id())
          {
            // surface is attached to actnode
            // see whether all  nodes on the surface are in our cloud
            bool is_conditioned_surface = true;
            for (int l = 0; l < nnodepersurf; ++l)
            {
              if (mycolnodes.find(nodespersurf[l]->id()) == mycolnodes.end())
              {
                is_conditioned_surface = false;
                // continue with next element surface
                break;
              }
            }
            // if all nodes are in our cloud, add surface
            if (is_conditioned_surface)
            {
              // remove internal surfaces that are connected to two volume elements
              Core::Nodes::Node** actsurfnodes = actsurf->nodes();

              // get sorted vector of node ids
              std::vector<int> nodes(actsurf->num_node());
              transform(actsurf->nodes(), actsurf->nodes() + actsurf->num_node(), nodes.begin(),
                  std::mem_fn(&Core::Nodes::Node::id));
              sort(nodes.begin(), nodes.end());

              // special treatment of RedAirwayTissue and StructFluidSurfCoupling
              if ((cond->type() == Core::Conditions::StructFluidSurfCoupling) or
                  (cond->type() == Core::Conditions::RedAirwayTissue))
              {
                // fill map with 'surfmap' std::vector<int> -> Teuchos::RCP<Core::Elements::Element>
                {
                  // indicator for volume element underlying the surface element
                  int identical = 0;
                  // owner of underlying volume element
                  int voleleowner = -1;
                  // number of considered nodes
                  const int numnode = (int)nodes.size();

                  // set of adjacent elements
                  std::set<Core::Elements::Element*> adjacentvoleles;

                  // get all volume elements connected to the nodes of this surface element
                  for (int n = 0; n < numnode; ++n)
                  {
                    Core::Nodes::Node* actsurfnode = actsurfnodes[n];
                    Core::Elements::Element** eles = actsurfnode->elements();
                    int numeles = actsurfnode->num_element();
                    for (int e = 0; e < numeles; ++e)
                    {
                      // do not consider volume elements that do not belong to VolEleIDs
                      if (VolEleIDs.size())
                        if (VolEleIDs.find(eles[e]->id()) == VolEleIDs.end()) continue;
                      adjacentvoleles.insert(eles[e]);
                    }
                  }

                  // get surfaces of all adjacent vol eles and check how often actsurf is included
                  // via comparison of node ids
                  for (const auto& adjacent_ele : adjacentvoleles)
                  {
                    std::vector<Teuchos::RCP<Core::Elements::Element>> adjacentvolelesurfs =
                        adjacent_ele->surfaces();
                    const int adjacentvolelenumsurfs = adjacent_ele->num_surface();
                    for (int n = 0; n < adjacentvolelenumsurfs; ++n)
                    {
                      // current surf of adjacent vol ele
                      std::vector<int> nodesadj(adjacentvolelesurfs[n]->num_node());
                      transform(adjacentvolelesurfs[n]->nodes(),
                          adjacentvolelesurfs[n]->nodes() + adjacentvolelesurfs[n]->num_node(),
                          nodesadj.begin(), std::mem_fn(&Core::Nodes::Node::id));
                      sort(nodesadj.begin(), nodesadj.end());

                      if (nodes.size() == nodesadj.size())
                      {
                        if (std::equal(nodes.begin(), nodes.end(), nodesadj.begin()))
                        {
                          identical++;
                          voleleowner = adjacent_ele->owner();
                        }  // if node ids are identical
                      }    // if number of nodes matches
                    }      // loop over all element surfaces
                  }        // loop over all adjacent volume elements

                  if (identical == 0)
                    FOUR_C_THROW("surface found with missing underlying volume element");
                  else if (identical > 1)
                  {
                    // do nothing
                  }
                  else
                  {
                    // now we can add the surface
                    if (surfmap.find(nodes) == surfmap.end())
                    {
                      Teuchos::RCP<Core::Elements::Element> surf = Teuchos::rcp(actsurf->clone());
                      // Set owning processor of surface owner of underlying volume element.
                      surf->set_owner(voleleowner);
                      surfmap[nodes] = surf;
                    }  // if surface not yet in map
                  }    // end if unique underlying vol ele was found
                }      // map 'surfmap' is now filled

                // continue with next element surface
                break;
              }  // if special treatment of RedAirwayTissue and StructFluidSurfCoupling
              else
              {
                // now we can add the surface
                if (surfmap.find(nodes) == surfmap.end())
                {
                  Teuchos::RCP<Core::Elements::Element> surf = Teuchos::rcp(actsurf->clone());
                  // Set owning processor of surface owner of underlying volume element.
                  surf->set_owner(elements[i]->owner());
                  surfmap[nodes] = surf;
                }  // if surface not yet in map
              }    // else if standard case
            }      // if all nodes of surface belong to condition (is_conditioned_surface == true)
          }        // if surface contains conditioned row node
        }          // loop over all nodes of element surface
      }            // loop over all element surfaces
    }              // loop over all adjacent elements of conditioned row node
  }                // loop over all conditioned row nodes

  // surfaces be added to the condition: (surf_id) -> (surface).
  Teuchos::RCP<std::map<int, Teuchos::RCP<Core::Elements::Element>>> final_geometry =
      Teuchos::rcp(new std::map<int, Teuchos::RCP<Core::Elements::Element>>());

  assign_global_i_ds(get_comm(), surfmap, *final_geometry);
  cond->add_geometry(final_geometry);

  // elements were created that need new unique ids
  bool havenewelements = true;
  // note: this seems useless since this build_surfacesin_condition always
  //       creates elements. However, this function is overloaded in
  //       MeshfreeDiscretization where it does not necessarily build elements.
  return havenewelements;

}  // Core::FE::Discretization::build_surfacesin_condition


/*----------------------------------------------------------------------*
 |  Build volume geometry in a condition (public)            mwgee 01/07|
 *----------------------------------------------------------------------*/
bool Core::FE::Discretization::build_volumesin_condition(
    const std::string& name, Teuchos::RCP<Core::Conditions::Condition> cond)
{
  // get ptrs to all node ids that have this condition
  const std::vector<int>* nodeids = cond->get_nodes();
  if (!nodeids) FOUR_C_THROW("Cannot find array 'Node Ids' in condition");

  // extract colnodes on this proc from condition
  const Epetra_Map* colmap = node_col_map();
  std::set<int> mynodes;

  std::remove_copy_if(nodeids->begin(), nodeids->end(), std::inserter(mynodes, mynodes.begin()),
      std::not_fn(Core::Conditions::MyGID(colmap)));

  // this is the map we want to construct
  auto geom = Teuchos::rcp(new std::map<int, Teuchos::RCP<Core::Elements::Element>>());

  for (const auto& [ele_id, actele] : element_)
  {
    std::vector<int> myelenodes(actele->node_ids(), actele->node_ids() + actele->num_node());

    // check whether all node ids of the element are nodes belonging
    // to the condition and stored on this proc
    bool allin = true;
    for (const auto& myid : myelenodes)
    {
      if (mynodes.find(myid) == mynodes.end())
      {
        // myid is not in the condition
        allin = false;
        break;
      }
    }

    if (allin)
    {
      (*geom)[ele_id] = actele;
    }
  }

  cond->add_geometry(geom);

  // no elements where created to assign new unique ids to
  return false;
}  // Core::FE::Discretization::build_volumesin_condition



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::FE::Discretization::find_associated_ele_i_ds(
    Teuchos::RCP<Core::Conditions::Condition> cond, std::set<int>& VolEleIDs,
    const std::string& name)
{
  // determine constraint number
  int condID = cond->parameters().get<int>("coupling id");

  std::vector<Core::Conditions::Condition*> volconds;
  get_condition(name, volconds);

  for (auto& actvolcond : volconds)
  {
    if (actvolcond->parameters().get<int>("coupling id") == condID)
    {
      // get ptrs to all node ids that have this condition
      const std::vector<int>* nodeids = actvolcond->get_nodes();
      if (!nodeids) FOUR_C_THROW("Cannot find array 'Node Ids' in condition");

      // extract colnodes on this proc from condition
      const Epetra_Map* colmap = node_col_map();
      std::set<int> mynodes;

      std::remove_copy_if(nodeids->begin(), nodeids->end(), std::inserter(mynodes, mynodes.begin()),
          std::not_fn(Core::Conditions::MyGID(colmap)));

      for (const auto& [ele_id, actele] : element_)
      {
        std::vector<int> myelenodes(actele->node_ids(), actele->node_ids() + actele->num_node());

        // check whether all node ids of the element are nodes belonging
        // to the condition and stored on this proc
        bool allin = true;
        for (const auto& myid : myelenodes)
        {
          if (mynodes.find(myid) == mynodes.end())
          {
            // myid is not in the condition
            allin = false;
            break;
          }
        }

        if (allin)
        {
          VolEleIDs.insert(actele->id());
        }
      }
    }
  }
}  // Core::FE::Discretization::find_associated_ele_i_ds

FOUR_C_NAMESPACE_CLOSE
