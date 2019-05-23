/*----------------------------------------------------------------------*/
/*!

\brief solid-shell body creation by extruding surface

\maintainer Martin Kronbichler

\level 1

Here everything related with solid-shell body extrusion
*/
/*----------------------------------------------------------------------*/
#include "pre_exodus_soshextrusion.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_validate.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "pre_exodus_centerline.H"



/* Method to extrude a surface to become a volumetric body */
EXODUS::Mesh EXODUS::SolidShellExtrusion(EXODUS::Mesh& basemesh, double thickness, int layers,
    int seedid, int gmsh, int concat2loose, int diveblocks, const std::string cline,
    const std::vector<double> coordcorr)
{
  int highestnid =
      basemesh
          .GetNumNodes();  // this is the currently highest id, a new node must become highestnid+1
  // map<int,std::vector<double> > newnodes;
  Teuchos::RCP<std::map<int, std::vector<double>>> newnodes =
      Teuchos::rcp(new std::map<int, std::vector<double>>);      // here the new nodes are stored
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> neweblocks;  // here the new EBlocks are stored
  std::map<int, EXODUS::NodeSet> newnodesets;                    // here the new NS are stored
  std::map<int, EXODUS::SideSet> newsidesets;                    // here the new SS are stored
  int highestblock = basemesh.GetNumElementBlocks();  // check whether there are ebs at all
  int highestns = basemesh.GetNumNodeSets();          // check whether there are nss at all
  int highestss = basemesh.GetNumSideSets();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = basemesh.GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;
  if (highestblock != 0)
    highestblock = ebs.rbegin()->first +
                   1;  // if there are ebs get the highest number, not necessarily consecutive
  else
    highestblock = 1;  // case of no eblocks at all -> starts with 1 due to exodus format

  std::map<int, EXODUS::NodeSet> nss = basemesh.GetNodeSets();
  std::map<int, EXODUS::NodeSet>::const_iterator i_nss;
  if (highestns != 0)
    highestns = nss.rbegin()->first +
                1;  // if there are nss get the highest number, not necessarily consecutive
  else
    highestns = 1;  // case of no nodesets at all -> starts with 1 due to exodus format

  std::map<int, EXODUS::SideSet> sss = basemesh.GetSideSets();
  std::map<int, EXODUS::SideSet>::const_iterator i_sss;
  if (highestss != 0)
    highestss = sss.rbegin()->first +
                1;  // if there are ss get the highest number, not necessarily consecutive
  else
    highestss = 1;  // case of no sidesets at all -> starts with 1 due to exodus format

  /* Extrusion is always based on a connectivity map*/
  // map of connectivity maps to be extruded
  std::map<int, std::map<int, std::vector<int>>> extrusion_conns;
  std::map<int, std::map<int, std::vector<int>>>::const_iterator i_extr;
  int extrusioncounter = 0;
  std::map<int, ExtrusionType> extrusion_types;
  std::map<int, ExtrusionType>::const_iterator i_exty;

  // counter for concatenated nodes during extrusion
  int concat_counter = 0;

  std::map<int, std::vector<int>> node_pair;  // stores new node id with base node id
  std::map<int, int> inv_node_pair;  // inverse map of above new node (first of layers) -> base node

  // store average node normals for each basemesh node
  std::map<int, std::vector<double>> node_normals;
  // store average node normals for base node of new extrusion
  std::map<int, std::vector<double>> ext_node_normals;

  std::set<int> nodes_from_sideset;
  int extrusion_sideset_id = -1;
  // loop through all SideSets to check for extrusion
  for (i_sss = sss.begin(); i_sss != sss.end(); ++i_sss)
  {
    bool toextrude = CheckExtrusion(
        i_sss->second);  // currently no rule is applied, the one (and only) sideset is extruded
    extrusion_sideset_id = i_sss->first;
    if (toextrude)
    {
      std::map<int, std::vector<int>> sidesetconn = basemesh.GetSideSetConn(i_sss->second, true);
      extrusion_conns.insert(
          std::pair<int, std::map<int, std::vector<int>>>(extrusioncounter, sidesetconn));
      // extrusion_conns.insert(std::pair<int,std::map<int,std::vector<int> >
      // >(extrusioncounter,basemesh.GetSideSetConn(i_sss->second)));
      extrusion_types.insert(std::pair<int, ExtrusionType>(extrusioncounter, sideset));
      extrusioncounter++;
      // create NodeSet from SideSet to apply bc e.g. pressure on extruded surface
      nodes_from_sideset = basemesh.GetSideSetNodes(i_sss->second, sidesetconn);
      // testing: also write out EBlocks of extruding sideset surface
      // vector<EXODUS::ElementBlock> sideblocks =
      // basemesh.SideSetToEBlocks(i_sss->second,sidesetconn);
      // neweblocks.insert(std::pair<int,EXODUS::ElementBlock>(highestblock,sideblocks[0]));
      // highestblock ++;
      // neweblocks.insert(std::pair<int,EXODUS::ElementBlock>(highestblock,sideblocks[1]));
      // highestblock ++;
    }
  }

  if (extrusioncounter == 0)
  {
    // loop through all EBlocks to check for extrusion blocks
    for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
    {
      bool toextrude = CheckExtrusion(*i_ebs->second);
      if (toextrude)
      {
        extrusion_conns.insert(std::pair<int, std::map<int, std::vector<int>>>(
            extrusioncounter, *(i_ebs->second->GetEleConn())));
        extrusion_types.insert(std::pair<int, ExtrusionType>(extrusioncounter, eblock));
        extrusioncounter++;
      }
    }
  }
  // obtain variable extrusion thicknesses (at nodes) from centerline
  std::map<int, double> nd_ts;
  if (cline != "")
  {
    // pass extrusion connectivity for gmsh-debugplot, use first as we up to now use only ONE
    // sideset-extrusionconn
    nd_ts = EXODUS::NdCenterlineThickness(
        cline, nodes_from_sideset, extrusion_conns.begin()->second, basemesh, thickness, coordcorr);
  }


  // loop all existing extrude connectivities
  std::cout << "Extruding surfaces..." << std::endl;
  for (i_extr = extrusion_conns.begin(); i_extr != extrusion_conns.end(); ++i_extr)
  {
    // get connectivity
    const std::map<int, std::vector<int>> ele_conn = i_extr->second;

    // Create Node to Element Connectivity
    const std::map<int, std::set<int>> node_conn = NodeToEleConn(ele_conn);

    // Create Element to Element Connectivity (sharing an edge)
    const std::map<int, std::vector<int>> ele_neighbor = EleNeighbors(ele_conn, node_conn);

    // fill set of "free" edge nodes, i.e. nodes with no ele edge neighbor
    const std::set<int> free_edge_nodes = FreeEdgeNodes(ele_conn, ele_neighbor);

    // loop through all its elements to create new connectivity  ***************
    Teuchos::RCP<std::map<int, std::vector<int>>> newconn =
        Teuchos::rcp(new std::map<int, std::vector<int>>);
    int newele = 0;

    std::map<int, std::vector<int>>
        newsideset;  // this sideset will become the new extruded out-'side'
    std::map<int, std::vector<int>>::iterator i_ss;

    // another map is needed for twistfixing: in case of layered extrusion it
    // consists of elements enclosing the layers, i.e. from most inner to most outer node
    std::map<int, std::vector<int>> encloseconn;
    std::map<int, std::vector<int>> layerstack;

    // set of elements already done
    std::set<int> doneles;
    std::set<int>::iterator i_doneles;

    // here we store a potential set AND vector of eles still to extrude
    std::set<int> todo_eleset;
    std::map<int, std::vector<int>>::const_iterator i_ele;
    for (i_ele = ele_conn.begin(); i_ele != ele_conn.end(); ++i_ele)
      todo_eleset.insert(i_ele->first);
    std::set<int>::iterator i_todo_eleset;
    std::vector<int> todo_eles;
    int todo_counter = 0;

    // define starte ele
    int startele = 0;

    // do check for normal direction at nodes, so far by scalar product
    // if false just the node numbering decides which is the extruding direction
    bool normcheck = true;
    if (seedid == -1)
      normcheck = false;
    else if (unsigned(seedid) > ele_conn.size())
      dserror("SeedID out of range!");
    else
      startele = seedid;

    // todo_eleset.insert(startele); // lets insert the very first one for a start
    todo_eles.push_back(startele);
    todo_eleset.erase(startele);

    // for the first element we set up everything *****************************
    std::vector<int> actelenodes = ele_conn.find(startele)->second;
    std::vector<int>::const_iterator i_node;
    int newid;
    bool havehexes = false;
    if (actelenodes.size() == 4) havehexes = true;

    // calculate the normal at the first node which defines the "outside" dir
    std::vector<int> myNodeNbrs = FindNodeNeighbors(actelenodes, actelenodes.front());
    std::vector<double> first_normal =
        Normal(myNodeNbrs[1], actelenodes.front(), myNodeNbrs[0], basemesh);

    // EXODUS::PlotStartEleGmsh(startele,actelenodes,basemesh,actelenodes.front(),first_normal);

    // create a new element
    std::vector<int> newelenodes;
    // create a vector of nodes for each layer which will form layered eles
    std::vector<std::vector<int>> layer_nodes(layers + 1);

    for (i_node = actelenodes.begin(); i_node < actelenodes.end(); ++i_node)
    {
      // place new node at new position
      std::vector<double> actcoords = basemesh.GetNode(*i_node);  // curr position
      // new base position equals actele (matching mesh!)
      const std::vector<double> newcoords = actcoords;

      // concatenating or node-merging feature
      newid = *i_node;
      concat_counter++;

      // create new node = loose concatenating
      if (concat_counter >= concat2loose)
      {
        newid = highestnid + 1;  // here just raise for each basenode
        highestnid++;
        concat_counter = 0;
      }

      // put new coords into newnode map
      newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));

      // put new node into map of OldNodeToNewNode
      std::vector<int> newids(1, newid);
      node_pair.insert(std::pair<int, std::vector<int>>(*i_node, newids));
      inv_node_pair.insert(std::pair<int, int>(newid, *i_node));
      // insert node into base layer
      layer_nodes[0].push_back(newid);

      // calculate extruding direction
      std::vector<double> normal = NodeToAvgNormal(
          *i_node, actelenodes, first_normal, node_conn, ele_conn, basemesh, normcheck);
      node_normals.insert(std::pair<int, std::vector<double>>(*i_node, normal));
      ext_node_normals.insert(std::pair<int, std::vector<double>>(newid, normal));
      // ele_nodenormals.push_back(normal);

      // create layers at this node location
      for (int i_layer = 1; i_layer <= layers; ++i_layer)
      {
        if (nd_ts.size() != 0)
        {
          thickness = nd_ts.find(*i_node)->second;  // get current node thickness
          if (thickness < 1E-12)
            std::cout << "Node: " << (*i_node) << " has zero thickness!" << std::endl;
        }
        const std::vector<double> newcoords =
            ExtrudeNodeCoords(actcoords, thickness, i_layer, layers, normal);
        // numbering of new ids nodewise not layerwise as may be expected
        newid = highestnid + 1;
        ++highestnid;
        // put new coords into newnode map
        newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));
        // put new node into map of OldNodeToNewNode
        node_pair[*i_node].push_back(newid);
        // finally store this node where it will be connected to an ele
        layer_nodes[i_layer].push_back(newid);
      }
    }

    /* find out initial sign of jacobian for this first extruded element as reference
     * The issue is that we can extrude into inside and outside therefore
     * maybe produce all negative elements however not twisted
     * Just make sure your initial element is NOT twisted (choose appropriate seed) */
    std::vector<int> firstele(layer_nodes[0]);
    firstele.insert(firstele.end(), layer_nodes[layers].begin(), layer_nodes[layers].end());
    int nnodes = firstele.size();  // could be either 6 or 8
    std::map<int, std::vector<double>> coords;
    for (int i = 0; i < nnodes; ++i)
    {
      coords.insert(
          std::pair<int, std::vector<double>>(firstele[i], newnodes->find(firstele[i])->second));
    }
    int initelesign = EleSaneSign(firstele, coords);
    encloseconn.insert(std::pair<int, std::vector<int>>(newele, firstele));
    int stackbase = newele;
    std::vector<int> empty;
    layerstack.insert(std::pair<int, std::vector<int>>(stackbase, empty));  // **************

    doneles.insert(startele);  // the first element is done ************************
    // form every new layer element
    for (int i_layer = 0; i_layer < layers; ++i_layer)
    {
      std::vector<int> basenodes = layer_nodes[i_layer];
      std::vector<int> ceilnodes = layer_nodes[i_layer + 1];
      newelenodes.insert(newelenodes.begin(), basenodes.begin(), basenodes.end());
      newelenodes.insert(newelenodes.end(), ceilnodes.begin(), ceilnodes.end());
      // insert new element into new connectivity
      newconn->insert(std::pair<int, std::vector<int>>(newele, newelenodes));
      layerstack.find(stackbase)->second.push_back(newele);

      // if we are at top layer put the corresponding side into newsideset
      if (i_layer == layers - 1)
      {
        std::vector<int> ss(3);  // third pos for later eblockid
        ss.at(0) = newele;       // first entry is element id
        if (newelenodes.size() == 8)
          ss.at(1) = 6;  // hexcase: top face id = 6 // bottom face ID = 5 //TODO
        else if (newelenodes.size() == 6)
          ss.at(1) = 5;  // wedgecase: top face id
        else
          dserror("wrong number of elenodes!");
        newsideset.insert(std::pair<int, std::vector<int>>(newele, ss));
      }

      ++newele;
      newelenodes.clear();
    }
    layer_nodes.clear();


    // while (todo_eles.size() > 0){ // start of going through ele neighbors///////
    while (doneles.size() < ele_conn.size())
    {
      // find an actele still to do
      int actele = todo_eles[todo_counter];
      // fall back if vector is dubious
      if (unsigned(actele) > ele_conn.size())
      {
        i_todo_eleset = todo_eleset.begin();
        actele = *i_todo_eleset;
        todo_eleset.erase(actele);
      }
      // delete actele from todo_eles list
      else
        todo_eleset.erase(actele);

      // get nodes of actual element
      std::vector<int> actelenodes = ele_conn.find(actele)->second;

      // get edge neighbors
      std::vector<int> actneighbors = ele_neighbor.at(actele);

      std::vector<int>::const_iterator i_nbr;
      unsigned int edge = 0;  // edge counter
      for (i_nbr = actneighbors.begin(); i_nbr < actneighbors.end(); ++i_nbr)
      {
        int actneighbor = *i_nbr;
        // check for undone neighbor
        i_doneles = doneles.find(actneighbor);
        if ((actneighbor != -1) && (i_doneles == doneles.end()))
        {
          std::vector<int> actneighbornodes = ele_conn.find(actneighbor)->second;
          // lets extrude *****************************************************

          // get nodepair of actual edge
          int firstedgenode = actelenodes.at(edge);
          int secedgenode;
          // switch in case of last to first node edge
          if (edge == actelenodes.size() - 1)
            secedgenode = actelenodes.at(0);
          else
            secedgenode = actelenodes.at(edge + 1);

          // create a vector of nodes for each layer which will form layered eles
          std::vector<std::vector<int>> layer_nodes(layers + 1);

          /* the new elements orientation is opposite the current one
           * therefore the FIRST node is secedgenode and it does exist */
          for (int i_layer = 0; i_layer <= layers; ++i_layer)
          {
            if (node_pair.find(secedgenode) == node_pair.end())
            {
              std::map<int, std::vector<int>> leftovers = EXODUS::ExtrusionErrorOutput(
                  secedgenode, todo_counter, doneles, ele_conn, todo_eleset);
              // EXODUS::PlotEleConnGmsh(leftovers,*basemesh.GetNodes());
              std::cout << "Check the gmsh-file 'extrusionproblems.gmsh' " << std::endl;
              EXODUS::PlotEleConnGmsh(ele_conn, *basemesh.GetNodes(), leftovers);
              std::cout << "Check the gmsh-file 'neighbors.gmsh' " << std::endl;
              std::vector<double> normal(3, 1.0);  // = node_normals.find(firstedgenode)->second; //
                                                   // alternatively use normal at firstedgenode
              EXODUS::PlotEleNbrs(actelenodes, actneighbors, ele_conn, basemesh, secedgenode,
                  normal, node_conn, node_normals);
              dserror("Mesh problem!");
            }
            else
            {
              newid = node_pair.find(secedgenode)->second[i_layer];
              // finally store this node where it will be connected to an ele
              layer_nodes[i_layer].push_back(newid);
            }
          }


          /* the new elements orientation is opposite the current one
           * therefore the SECOND node is firstedgenode and it also does exist */
          for (int i_layer = 0; i_layer <= layers; ++i_layer)
          {
            newid = node_pair.find(firstedgenode)->second[i_layer];
            // finally store this node where it will be connected to an ele
            layer_nodes[i_layer].push_back(newid);
          }


          /* the new elements orientation is opposite the current one
           * therefore the THIRD node is gained from the neighbor ele */

          std::vector<int> actnbrnodes = ele_conn.find(actneighbor)->second;
          int thirdnode = FindEdgeNeighbor(actnbrnodes, firstedgenode, secedgenode);

          // check if new node already exists
          if (node_pair.find(thirdnode) == node_pair.end())
          {
            // place node at new position
            std::vector<double> actcoords = basemesh.GetNode(thirdnode);  // curr position
            // new base position equals actele (matching mesh!)
            const std::vector<double> newcoords = actcoords;

            // concatenating or node-merging feature
            newid = thirdnode;
            concat_counter++;

            // create new node = loose concatenating
            if (concat_counter >= concat2loose)
            {
              newid = highestnid + 1;  // here just raise for each basenode
              highestnid++;
              concat_counter = 0;
            }

            // put new coords into newnode map
            newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));

            // put new node into map of OldNodeToNewNode
            std::vector<int> newids(1, newid);
            node_pair.insert(std::pair<int, std::vector<int>>(thirdnode, newids));
            inv_node_pair.insert(std::pair<int, int>(newid, thirdnode));

            // insert node into base layer
            layer_nodes[0].push_back(newid);

            // get reference direction for new avg normal here firstedgenode
            std::vector<double> refnormal = node_normals.find(firstedgenode)->second;
            // calculate extruding direction
            std::vector<double> normal = NodeToAvgNormal(
                thirdnode, actneighbornodes, refnormal, node_conn, ele_conn, basemesh, normcheck);
            node_normals.insert(std::pair<int, std::vector<double>>(thirdnode, normal));
            ext_node_normals.insert(std::pair<int, std::vector<double>>(newid, normal));

            // create layers at this node location
            for (int i_layer = 1; i_layer <= layers; ++i_layer)
            {
              if (nd_ts.size() != 0)
              {
                thickness = nd_ts.find(thirdnode)->second;  // get current node thickness
                if (thickness < 1E-12)
                  std::cout << "Node: " << (thirdnode) << " has zero thickness!" << std::endl;
              }
              const std::vector<double> newcoords =
                  ExtrudeNodeCoords(actcoords, thickness, i_layer, layers, normal);
              newid = highestnid + 1;
              ++highestnid;
              // put new coords into newnode map
              newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));

              // put new node into map of OldNodeToNewNode
              node_pair[thirdnode].push_back(newid);
              // finally store this node where it will be connected to an ele
              layer_nodes[i_layer].push_back(newid);
            }
            // EXODUS::PlotEleNbrs(actelenodes,actneighbors,ele_conn,basemesh,thirdnode,normal,node_conn,node_normals);
          }
          else
          {
            for (int i_layer = 0; i_layer <= layers; ++i_layer)
            {
              newid = node_pair.find(thirdnode)->second[i_layer];
              // finally store this node where it will be connected to an ele
              layer_nodes[i_layer].push_back(newid);
            }
          }

          if (actnbrnodes.size() > 3)
          {  // in case of neighbor not being a tri3
            havehexes = true;

            /* the new elements orientation is opposite the current one
             * therefore the FOURTH node is gained from the neighbor ele */

            int fourthnode = FindEdgeNeighbor(actnbrnodes, thirdnode, firstedgenode);

            // check if new node already exists
            if (node_pair.find(fourthnode) == node_pair.end())
            {
              // place new node at new position
              std::vector<double> actcoords = basemesh.GetNode(fourthnode);  // curr position
              // new base position equals actele (matching mesh!)
              const std::vector<double> newcoords = actcoords;

              // concatenating or node-merging feature
              newid = fourthnode;
              concat_counter++;

              // create new node = loose concatenating
              if (concat_counter >= concat2loose)
              {
                newid = highestnid + 1;  // here just raise for each basenode
                highestnid++;
                concat_counter = 0;
              }

              // put new coords into newnode map
              newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));

              // put new node into map of OldNodeToNewNode
              std::vector<int> newids(1, newid);
              node_pair.insert(std::pair<int, std::vector<int>>(fourthnode, newids));
              inv_node_pair.insert(std::pair<int, int>(newid, fourthnode));

              // insert node into base layer
              layer_nodes[0].push_back(newid);

              // get reference direction for new avg normal here firstedgenode
              std::vector<double> refnormal = node_normals.find(secedgenode)->second;
              // calculate extruding direction
              std::vector<double> normal = NodeToAvgNormal(fourthnode, actneighbornodes, refnormal,
                  node_conn, ele_conn, basemesh, normcheck);
              node_normals.insert(std::pair<int, std::vector<double>>(fourthnode, normal));
              ext_node_normals.insert(std::pair<int, std::vector<double>>(newid, normal));

              // create layers at this node location
              for (int i_layer = 1; i_layer <= layers; ++i_layer)
              {
                if (nd_ts.size() != 0)
                {
                  thickness = nd_ts.find(fourthnode)->second;  // get current node thickness
                  if (thickness < 1E-12)
                    std::cout << "Node: " << (fourthnode) << " has zero thickness!" << std::endl;
                }
                const std::vector<double> newcoords =
                    ExtrudeNodeCoords(actcoords, thickness, i_layer, layers, normal);
                newid = highestnid + 1;
                ++highestnid;
                // put new coords into newnode map
                newnodes->insert(std::pair<int, std::vector<double>>(newid, newcoords));

                // put new node into map of OldNodeToNewNode
                node_pair[fourthnode].push_back(newid);
                // finally store this node where it will be connected to an ele
                layer_nodes[i_layer].push_back(newid);
              }
            }
            else
            {
              for (int i_layer = 0; i_layer <= layers; ++i_layer)
              {
                newid = node_pair.find(fourthnode)->second[i_layer];
                // finally store this node where it will be connected to an ele
                layer_nodes[i_layer].push_back(newid);
              }
            }
          }  // end of 4-node case

          // create enclose-element and insert into map
          std::vector<int> ele(layer_nodes[0]);
          ele.insert(ele.end(), layer_nodes[layers].begin(), layer_nodes[layers].end());
          encloseconn.insert(std::pair<int, std::vector<int>>(newele, ele));
          int stackbase = newele;
          std::vector<int> empty;
          layerstack.insert(std::pair<int, std::vector<int>>(stackbase, empty));  // **************

          // form every new layer element
          for (int i_layer = 0; i_layer < layers; ++i_layer)
          {
            std::vector<int> basenodes = layer_nodes[i_layer];
            std::vector<int> ceilnodes = layer_nodes[i_layer + 1];
            newelenodes.insert(newelenodes.begin(), basenodes.begin(), basenodes.end());
            newelenodes.insert(newelenodes.end(), ceilnodes.begin(), ceilnodes.end());
            // insert new element into new connectivity
            newconn->insert(std::pair<int, std::vector<int>>(newele, newelenodes));
            layerstack.find(stackbase)->second.push_back(newele);

            // if we are at top layer put the corresponding side into newsideset
            if (i_layer == layers - 1)
            {
              std::vector<int> ss(3);  // third pos for later eblockid
              ss.at(0) = newele;       // first entry is element id
              if (newelenodes.size() == 8)
                ss.at(1) = 6;  // hexcase: top face id = 6, bottom = 5 //TODO
              else if (newelenodes.size() == 6)
                ss.at(1) = 5;  // wedgecase: top face id
              else
                dserror("wrong number of elenodes!");
              newsideset.insert(std::pair<int, std::vector<int>>(newele, ss));
            }

            ++newele;
            newelenodes.clear();
          }
          layer_nodes.clear();

          // insert actneighbor into done elements
          doneles.insert(actneighbor);

          // neighbor eles are possible next "center" eles
          // todo_eleset.insert(actneighbor);
          todo_eles.push_back(actneighbor);

        }        // end of if undone->extrude this neighbor, next neighbor *************
        ++edge;  // next element edge

      }  // end of this "center" - element ///////////////////////////////////////
      todo_counter++;

      if (gmsh == todo_counter) PlotEleConnGmsh(*newconn, *newnodes, todo_counter);

    }  // end of extruding all elements in connectivity

    // repair twisted extrusion elements
    std::cout << "Repairing twisted elements..." << std::endl;
    std::map<int, std::set<int>> ext_node_conn = NodeToEleConn(*newconn);
    int twistcounter = RepairTwistedExtrusion(thickness, nd_ts, layers, initelesign, highestnid,
        encloseconn, layerstack, *newconn, *newnodes, ext_node_conn, ext_node_normals, node_pair,
        inv_node_pair);
    std::cout << twistcounter << " makro-elements and all their layer-elements repaired."
              << std::endl;

    // Gmsh Debug-Output of extrusion Mesh
    if (gmsh == 0) PlotEleConnGmsh(*newconn, *newnodes);

    // create new Element Blocks
    std::ostringstream blockname;
    blockname << "ext";
    // std::string blockname = "extr";
    switch (extrusion_types.find(i_extr->first)->second)
    {
      case eblock:
      {  // Eblocks have only one type of eles
         //      const int numnodes = newconn->find(0)->second.size();
         //      ElementBlock::Shape newshape = ElementBlock::dis_none;
         //      if (numnodes == 6) newshape = ElementBlock::wedge6;
         //      else if (numnodes == 8) newshape = ElementBlock::hex8;
         //      else dserror("Number of basenodes for extrusion not supported");
         //      blockname << highestblock;
         //      Teuchos::RCP<EXODUS::ElementBlock> neweblock = Teuchos::rcp(new
         //      ElementBlock(newshape,newconn,blockname.str()));
         //      neweblocks.insert(std::pair<int,Teuchos::RCP<EXODUS::ElementBlock>
         //      >(highestblock,neweblock)); highestblock ++; break;
      }
      case sideset:
      {  // SideSets can have different types of eles and have to be checked individually
        std::map<int, std::vector<int>>::const_iterator i_ele;
        Teuchos::RCP<std::map<int, std::vector<int>>> hexconn =
            Teuchos::rcp(new std::map<int, std::vector<int>>);
        int hexcounter = 0;
        Teuchos::RCP<std::map<int, std::vector<int>>> wegconn =
            Teuchos::rcp(new std::map<int, std::vector<int>>);
        int wegcounter = 0;

        // case of 2 eblocks over extrusion (diveblocks > 0) ===================
        Teuchos::RCP<std::map<int, std::vector<int>>> hexconn2 =
            Teuchos::rcp(new std::map<int, std::vector<int>>);
        int hexcounter2 = 0;
        Teuchos::RCP<std::map<int, std::vector<int>>> wegconn2 =
            Teuchos::rcp(new std::map<int, std::vector<int>>);
        int wegcounter2 = 0;
        int innerelecounter = -1;
        // =====================================================================

        for (i_ele = newconn->begin(); i_ele != newconn->end(); ++i_ele)
        {
          int numnodes = i_ele->second.size();
          if (numnodes == 8)
          {
            if (diveblocks == 0)
            {  // standard case: just one eblock over extrusion
              hexconn->insert(std::pair<int, std::vector<int>>(hexcounter, i_ele->second));
              if (newsideset.find(i_ele->first) != newsideset.end())
              {
                newsideset.find(i_ele->first)->second.at(0) =
                    hexcounter;  // update sideset with new element ids
                // store eblockid the sideset is related to at 3rd position, in hexcase it will be
                // highestblock
                newsideset.find(i_ele->first)->second.at(2) = highestblock;
              }
              hexcounter++;
            }
            else
            {  // special case of 2 eblocks over extrusion
              // initialize counter if we have a base-ele
              if (layerstack.find(i_ele->first) != layerstack.end()) innerelecounter = 0;

              if (innerelecounter < diveblocks)
              {  // insert ele into inner layer eblock
                hexconn->insert(std::pair<int, std::vector<int>>(hexcounter, i_ele->second));
                hexcounter++;
                innerelecounter++;
              }
              else
              {  // insert ele into outer layer eblock
                hexconn2->insert(std::pair<int, std::vector<int>>(hexcounter2, i_ele->second));
                if (newsideset.find(i_ele->first) != newsideset.end())
                {
                  newsideset.find(i_ele->first)->second.at(0) =
                      hexcounter2;  // update sideset with new element ids
                  // store eblockid the sideset is related to at 3rd position, it will be
                  // highestblock + 1
                  newsideset.find(i_ele->first)->second.at(2) = highestblock + int(diveblocks > 0);
                }
                hexcounter2++;
              }
            }
          }
          else if (numnodes == 6)
          {
            if (diveblocks == 0)
            {  // standard case: just one eblock over extrusion
              wegconn->insert(std::pair<int, std::vector<int>>(wegcounter, i_ele->second));
              if (newsideset.find(i_ele->first) != newsideset.end())
              {
                newsideset.find(i_ele->first)->second.at(0) =
                    wegcounter;  // update sideset with new element ids
                // store eblockid the sideset is related to at 3rd position, in wegcase it will be
                // highestblock or if we have hexes highestblock+1
                newsideset.find(i_ele->first)->second.at(2) = highestblock + int(havehexes);
              }
              wegcounter++;
            }
            else
            {  // special case of 2 eblocks over extrusion
              // initialize counter if we have a base-ele
              if (layerstack.find(i_ele->first) != layerstack.end()) innerelecounter = 0;

              if (innerelecounter < diveblocks)
              {  // insert ele into inner layer eblock
                wegconn->insert(std::pair<int, std::vector<int>>(wegcounter, i_ele->second));
                wegcounter++;
                innerelecounter++;
              }
              else
              {  // insert ele into outer layer eblock
                wegconn2->insert(std::pair<int, std::vector<int>>(wegcounter2, i_ele->second));
                if (newsideset.find(i_ele->first) != newsideset.end())
                {
                  newsideset.find(i_ele->first)->second.at(0) =
                      wegcounter2;  // update sideset with new element ids
                  // store eblockid the sideset is related to at 3rd position, in wegcase it will be
                  // highestblock or if we have hexes highestblock+2+1
                  newsideset.find(i_ele->first)->second.at(2) =
                      highestblock + 2 * int(havehexes) + int(diveblocks > 0);
                }
                wegcounter2++;
              }
            }
          }
          else
            dserror("Number of basenodes for extrusion not supported");
        }


        if (hexcounter > 0)
        {
          std::ostringstream hexblockname;
          hexblockname << blockname.str() << "h" << highestblock;
          const std::string hexname = hexblockname.str();
          Teuchos::RCP<EXODUS::ElementBlock> neweblock =
              Teuchos::rcp(new ElementBlock(ElementBlock::hex8, hexconn, hexname));
          neweblocks.insert(
              std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(highestblock, neweblock));
          highestblock++;
        }
        if (hexcounter2 > 0)
        {
          std::ostringstream hexblockname;
          hexblockname << blockname.str() << "h" << highestblock;
          const std::string hexname = hexblockname.str();
          Teuchos::RCP<EXODUS::ElementBlock> neweblock =
              Teuchos::rcp(new ElementBlock(ElementBlock::hex8, hexconn2, hexname));
          neweblocks.insert(
              std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(highestblock, neweblock));
          highestblock++;
        }

        if (wegcounter > 0)
        {
          std::ostringstream wegblockname;
          wegblockname << blockname.str() << "w" << highestblock;
          Teuchos::RCP<EXODUS::ElementBlock> neweblock =
              Teuchos::rcp(new ElementBlock(ElementBlock::wedge6, wegconn, wegblockname.str()));
          neweblocks.insert(
              std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(highestblock, neweblock));
          highestblock++;
        }
        if (wegcounter2 > 0)
        {
          std::ostringstream wegblockname;
          wegblockname << blockname.str() << "w" << highestblock;
          Teuchos::RCP<EXODUS::ElementBlock> neweblock =
              Teuchos::rcp(new ElementBlock(ElementBlock::wedge6, wegconn2, wegblockname.str()));
          neweblocks.insert(
              std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(highestblock, neweblock));
          highestblock++;
        }

        // put new sideset into map
        const std::string sidesetname = "extsideset";
        EXODUS::SideSet newSideSet = EXODUS::SideSet(newsideset, sidesetname);
        newsidesets.insert(std::pair<int, EXODUS::SideSet>(highestss, newSideSet));
        break;
      }
      default:
        dserror("unrecognized extrude type");
    }

    /* Create new "extruded" NodeSet of existing NodeSet,
     * e.g. a new surface NodeSet out of an existing line NodeSet */
    for (i_nss = nss.begin(); i_nss != nss.end(); ++i_nss)
    {
      EXODUS::NodeSet existing_ns = i_nss->second;
      const std::set<int> extruded_nodes =
          FindExtrudedNodes(free_edge_nodes, node_pair, existing_ns.GetNodeSet());
      if (extruded_nodes.size() != 0)
      {
        std::ostringstream nodesetname;
        nodesetname << "ext" << existing_ns.GetName();
        EXODUS::NodeSet newnodeset(extruded_nodes, nodesetname.str(), nodesetname.str());
        newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, newnodeset));
        highestns++;
      }
    }



    // additionally create new NodeSet with ALL nodes at newly created "free" faces
    // this is the sum of above, but could be handy for applying just one BC
    std::set<int> free_nodes = FreeFaceNodes(free_edge_nodes, node_pair);
    std::ostringstream nodesetname;
    nodesetname << "ext_free_bnd" << highestns;
    ;
    EXODUS::NodeSet newnodeset(free_nodes, nodesetname.str(), nodesetname.str());
    newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, newnodeset));
    highestns++;

  }  // end of extruding

  // Won't work a second time if extruded in "wrong direction"! //TODO
  // create 2 NodeSets consisting of the inner extrusion face and the outer extrusion face
  // respectively
  std::set<int> nodes_extrusion_base;
  std::set<int> nodes_extrusion_roof;
  std::map<int, std::vector<double>>::const_iterator it;
  for (it = node_normals.begin(); it != node_normals.end(); ++it)
  {
    nodes_extrusion_base.insert(node_pair.find(it->first)->second.front());
    nodes_extrusion_roof.insert(node_pair.find(it->first)->second.back());
  }
  std::ostringstream nodesetname;
  nodesetname << "base_nodes" << highestns;  // sideset.GetName() << "nodes";
  std::string propname = "";
  EXODUS::NodeSet nodeset_extrusion_base(nodes_extrusion_base, nodesetname.str(), propname);
  newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, nodeset_extrusion_base));
  highestns++;

  std::ostringstream nodesetnamer;
  nodesetnamer << "roof_nodes" << highestns;  // sideset.GetName() << "nodes";
  EXODUS::NodeSet nodeset_extrusion_roof(nodes_extrusion_roof, nodesetnamer.str(), propname);
  newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, nodeset_extrusion_roof));
  highestns++;

  // extrude NodeSets which transfers a marked NodeSet to its extruded base- and roof-NodeSet
  // loop through all NodeSets to check for extrusion ones
  for (i_nss = nss.begin(); i_nss != nss.end(); ++i_nss)
  {
    bool toextrude = CheckExtrusion(i_nss->second);
    if (toextrude)
    {
      std::set<int> nodes_from_nodeset = (i_nss->second).GetNodeSet();
      if (node_pair.find(*(nodes_from_nodeset.begin())) != node_pair.end())
      {
        std::set<int> nodes_extrusion_base;
        std::set<int> nodes_extrusion_roof;
        std::set<int>::iterator it;
        for (it = nodes_from_nodeset.begin(); it != nodes_from_nodeset.end(); ++it)
        {
          nodes_extrusion_base.insert(node_pair.find(*it)->second.front());
          nodes_extrusion_roof.insert(node_pair.find(*it)->second.back());
        }
        std::ostringstream nodesetname;
        std::string propname = "";
        nodesetname << "b_" << (i_nss->second).GetName();
        EXODUS::NodeSet nodeset_extrusion_base(nodes_extrusion_base, nodesetname.str(), propname);
        newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, nodeset_extrusion_base));
        highestns++;

        std::ostringstream nodesetnamer;
        nodesetnamer << "r_" << (i_nss->second).GetName();
        EXODUS::NodeSet nodeset_extrusion_roof(nodes_extrusion_roof, nodesetnamer.str(), propname);
        newnodesets.insert(std::pair<int, EXODUS::NodeSet>(highestns, nodeset_extrusion_roof));
        highestns++;
      }
    }
  }
  // rest works fine! */ //TODO
  // check and make flat extrusion surfaces for symmetry or other boundary conditions
  // loop through all NodeSets to check for flat extrusion
  for (i_nss = nss.begin(); i_nss != nss.end(); ++i_nss)
  {
    bool FlatEx = CheckFlatEx(i_nss->second);
    if (FlatEx)
    {
      std::cout << "Flattening Nodeset " << (i_nss->second).GetName()
                << " to its Normal: ";  //"..." <<std::endl;
      std::set<int> nodes_from_nodeset = (i_nss->second).GetNodeSet();
      std::set<int>::iterator it;

      // ***** compute flat surface from nodeset *****
      // compute normal
      std::set<int>::iterator surfit = nodes_from_nodeset.begin();
      int origin = *surfit;  // get first set node
      ++surfit;
      int head1 = *surfit;  // get second set node
      ++surfit;

      // compute normal
      std::vector<double> facenormal;
      for (it = surfit; it != nodes_from_nodeset.end(); ++it)
      {
        // third node such that a proper normal can be computed
        std::set<int>::iterator thirdnode = it;
        facenormal = Normal(head1, origin, *thirdnode, basemesh);
        // leave for loop if valid normal is found
        if (facenormal.size() != 1) break;
      }

      PrintVec(std::cout, facenormal);

      // could a normal direction be computed?
      if (facenormal.size() == 1)
      {
        std::cout << "  Warning! No normal defined within flat nodeset '"
                  << (i_nss->second).GetName() << "', stop flattening" << std::endl;
      }
      else
      {
        for (it = nodes_from_nodeset.begin(); it != nodes_from_nodeset.end(); ++it)
        {
          if (node_pair.find(*it) != node_pair.end())
          {
            const std::vector<int> nodes2flatten = node_pair.find(*it)->second;
            std::vector<int>::const_iterator i_node;
            std::vector<double> normal = node_normals.find(*it)->second;

            // project orignormal into flat plane
            double projection =
                normal[0] * facenormal[0] + normal[1] * facenormal[1] + normal[2] * facenormal[2];
            double length = 0.0;
            for (int i = 0; i < 3; ++i)
            {
              normal[i] -= projection * facenormal[i];
              length += normal[i] * normal[i];
            }
            length = sqrt(length);
            for (int i = 0; i < 3; ++i) normal[i] /= length;

            // keep this wonderful flat normal
            node_normals.find(*it)->second = normal;

            // correct node coords for all layers
            int i_layer = 1;
            std::vector<double> basecoords = newnodes->find(nodes2flatten[0])->second;
            for (i_node = (nodes2flatten.begin() + 1); i_node != nodes2flatten.end(); ++i_node)
            {
              int node = *i_node;
              if (nd_ts.size() != 0)
                thickness = nd_ts.find(*it)->second;  // get current node thickness
              std::vector<double> correctcoords =
                  ExtrudeNodeCoords(basecoords, thickness, i_layer, layers, normal);
              newnodes->find(node)->second = correctcoords;
              ++i_layer;
            }
          }
        }
      }
    }
  }



  // get rid of original extrusion SideSet
  basemesh.EraseSideSet(extrusion_sideset_id);

  std::string newtitle = "extrusion";
  // map<int,EXODUS::SideSet> emptysideset;

  std::cout << "...done" << std::endl;

  EXODUS::Mesh extruded_mesh(basemesh, newnodes, neweblocks, newnodesets, newsidesets, newtitle);

  return extruded_mesh;
}


bool EXODUS::CheckExtrusion(const EXODUS::ElementBlock eblock)
{
  const EXODUS::ElementBlock::Shape myshape = eblock.GetShape();
  const std::string myname = eblock.GetName();
  if ((myname.find("toex") != std::string::npos) or (myname.find("extrude") != std::string::npos))
    if ((myshape == EXODUS::ElementBlock::shell4) || (myshape == EXODUS::ElementBlock::tri3))
      return true;
  return false;
}
bool EXODUS::CheckExtrusion(const EXODUS::SideSet sideset)
{
  const std::string myname = sideset.GetName();
  // if (myname.compare(0,7,"extrude") == 0)
  return true;
  // return false;
}
bool EXODUS::CheckExtrusion(const EXODUS::NodeSet nodeset)
{
  const std::string myname = nodeset.GetName();
  if ((myname.find("toex") != std::string::npos) or (myname.find("extrude") != std::string::npos))
    return true;
  return false;
}

bool EXODUS::CheckFlatEx(const EXODUS::NodeSet nodeset)
{
  const std::string myname = nodeset.GetName();
  const std::string mypropname = nodeset.GetPropName();
  if ((myname.find("flat") != std::string::npos) || (myname.find("FLAT") != std::string::npos) ||
      (mypropname.find("flat") != std::string::npos) ||
      (mypropname.find("FLAT") != std::string::npos))
  {
    return true;
  }
  return false;
}


int EXODUS::RepairTwistedExtrusion(const double thickness,  // extrusion thickness
    const std::map<int, double>& nd_ts,                     // map of variable thicknesses at node
    const int numlayers,                                    // number of layers in extrusion
    const int initelesign,  // the sign of det(J) of first extruded ele as reference for pos or neg
                            // extrusion (in)
    int& highestnid,        // current global highest node id (in/out)
    const std::map<int, std::vector<int>>&
        encloseconn,  // connectivity of enclosing elements to be tested for twist (in)
    const std::map<int, std::vector<int>>&
        layerstack,  // stack of layer elements above the base element
    std::map<int, std::vector<int>>&
        newconn,  // connectivity of created extrusion elements (in/out)
    std::map<int, std::vector<double>>& newnodes,  // created extrusion nodes (in/out)
    std::map<int, std::set<int>>&
        ext_node_conn,  // node to ele map of created extrusion connectivity (in/out)
    const std::map<int, std::vector<double>>&
        avgnode_normals,                         // averaged node normals at extrusion basenodes
    std::map<int, std::vector<int>>& node_pair,  // map basenode -> newnodes of extrusion
    const std::map<int, int>& inv_node_pair)     // inverse map of above
{
  std::map<int, std::vector<int>>::const_iterator i_encl;
  std::vector<int>::const_iterator it;
  int twistcounter = 0;
  int firstcheck = 0;
  int secondcheck = 0;
  int newnodesbyrepair = 0;
  std::map<int, std::vector<int>> repaired_conn;  // for gmsh plot
  double actthickness = thickness;

  for (i_encl = encloseconn.begin(); i_encl != encloseconn.end(); ++i_encl)
  {
    std::vector<int> actele = i_encl->second;
    int nnodes = actele.size();  // could be either 6 or 8
    std::map<int, std::vector<double>> coords;
    for (int i = 0; i < nnodes; ++i)
    {
      coords.insert(
          std::pair<int, std::vector<double>>(actele[i], newnodes.find(actele[i])->second));
    }
    int actelesign = EleSaneSign(actele, coords);

    if (actelesign != initelesign)
    {  // thus we have a twisted element
      ++twistcounter;
      // PlotEleGmsh(actele, newnodes);
      // std::cout << "twisted: ";PrintVec(std::cout,i_encl->second);

      // get extrusion base which is always the first half of actele
      std::vector<int> baseface;
      for (int i = 0; i < nnodes / 2; ++i) baseface.push_back(actele[i]);

      // store angular deviation between avgNormal to nodeNormal
      std::map<double, int>
          normaldeviations;  // carefull, we compare doubles here, should work though
      std::pair<std::map<double, int>::iterator, bool> doubleentry;
      const double epsilon = 1.0E-6;
      std::map<int, std::vector<double>> node_normals;
      // calculate node normals of baseface
      for (it = baseface.begin(); it != baseface.end(); ++it)
      {
        std::vector<int> nnbr = FindNodeNeighbors(baseface, *it);
        std::vector<double> normal = Normal(nnbr.front(), *it, nnbr.back(), coords);
        node_normals.insert(std::pair<int, std::vector<double>>(*it, normal));
        std::vector<double> avgnormal = avgnode_normals.find(*it)->second;
        // scalarproduct normal.avgNormal delivers deviation
        double dev =
            abs(normal[0] * avgnormal[0] + normal[1] * avgnormal[1] + normal[2] * avgnormal[2]);
        doubleentry = normaldeviations.insert(std::pair<double, int>(dev, *it));
        // in case of two exactly equal deviations we raise it a bit
        if (doubleentry.second == false)
          normaldeviations.insert(std::pair<double, int>(dev + epsilon, *it));
      }
      std::map<double, int>::iterator i_dev = normaldeviations.begin();

      while (actelesign != initelesign && i_dev != normaldeviations.end())
      {
        // get current repairnode and its normal and extrude to new coords
        int repairbasenode = i_dev->second;
        std::vector<double> repairnormal = node_normals.find(repairbasenode)->second;
        CheckNormDir(repairnormal, avgnode_normals.find(repairbasenode)->second);
        if (nd_ts.size() != 0)
        {
          if (nd_ts.size() != avgnode_normals.size()) std::cout << "Size mismatch!" << std::endl;
          if (nd_ts.find(inv_node_pair.find(repairbasenode)->second) == nd_ts.end())
            std::cout << "Node thickness not found!" << std::endl;
          actthickness = nd_ts.find(inv_node_pair.find(repairbasenode)->second)
                             ->second;  // get current node thickness
          if (actthickness < 1E-12)
            std::cout << "Node: " << (inv_node_pair.find(repairbasenode)->second)
                      << " has zero thickness: " << actthickness << std::endl;
        }
        std::vector<double> newcoords = ExtrudeNodeCoords(
            coords.find(repairbasenode)->second, actthickness, 1, 1, repairnormal);

        int repairnodepos = FindPosinVec(repairbasenode, baseface) + nnodes / 2;
        int repairnode = actele.at(repairnodepos);

        // create new node or move node if only one element is left at this node
        if (ext_node_conn.find(repairnode)->second.size() == 1)
        {
          // move existing node
          newnodes.find(repairnode)->second = newcoords;
          coords.find(repairnode)->second = newcoords;  // coords update

          // layer case rework layer elements within enclosing one
          if (numlayers > 1)
          {
            const std::vector<int> layereles = layerstack.find(i_encl->first)->second;
            std::vector<int>::const_iterator i_layerele;
            int i_layer = 1;
            for (i_layerele = layereles.begin(); i_layerele != (layereles.end() - 1); ++i_layerele)
            {
              if (nd_ts.size() != 0)
              {
                if (nd_ts.size() != avgnode_normals.size())
                  std::cout << "Size mismatch!" << std::endl;
                actthickness = nd_ts.find(inv_node_pair.find(repairbasenode)->second)
                                   ->second;  // get current node thickness
                if (actthickness < 1E-12)
                  std::cout << "Node: " << (inv_node_pair.find(repairbasenode)->second)
                            << " has zero thickness: " << actthickness << std::endl;
              }
              std::vector<double> newlayercoords =
                  ExtrudeNodeCoords(coords.find(repairbasenode)->second, actthickness, i_layer,
                      numlayers, repairnormal);
              ++i_layer;
              int layerrepairnode = newconn.find(*i_layerele)->second.at(repairnodepos);
              newnodes.find(layerrepairnode)->second = newlayercoords;
            }
          }
        }
        else
        {
          // we need a new node

          // layer case rework layer elements within enclosing one
          const std::vector<int> layereles = layerstack.find(i_encl->first)->second;
          std::vector<int>::const_iterator i_layerele;
          int i_layer = 1;
          for (i_layerele = layereles.begin(); i_layerele != layereles.end(); ++i_layerele)
          {
            if (nd_ts.size() != 0)
            {
              actthickness = nd_ts.find(inv_node_pair.find(repairbasenode)->second)
                                 ->second;  // get current node thickness
              if (actthickness < 1E-12)
                std::cout << "Node: " << (inv_node_pair.find(repairbasenode)->second)
                          << " has zero thickness!" << std::endl;
            }
            std::vector<double> newlayercoords =
                ExtrudeNodeCoords(coords.find(repairbasenode)->second, actthickness, i_layer,
                    numlayers, repairnormal);

            // replace also the base for inner inner
            if (i_layer > 1)
              newconn.find(*i_layerele)->second.at(repairnodepos - nnodes / 2) = highestnid;

            // create new node
            ++highestnid;
            ++newnodesbyrepair;
            // put new coords into newnode map
            newnodes.insert(std::pair<int, std::vector<double>>(highestnid, newlayercoords));
            coords.insert(
                std::pair<int, std::vector<double>>(highestnid, newlayercoords));  // coords update

            // replace previously connected node corresponding to repairnode with new node in
            // actlayerele
            newconn.find(*i_layerele)->second.at(repairnodepos) = highestnid;

            ++i_layer;
          }
          // update actual enclosing element for following double check
          actele.at(repairnodepos) = highestnid;
          // update ext_node_conn
          ext_node_conn.find(repairnode)->second.erase(i_encl->first);
        }

        // while loop update
        actelesign = EleSaneSign(actele, coords);
        ++i_dev;
        // PlotEleGmsh(actele, newnodes);
        ++firstcheck;
      }


      // check repaired element
      int repairedelesign = EleSaneSign(actele, coords);

      // double check
      if (repairedelesign != initelesign)
      {
        ++secondcheck;
        // PlotEleGmsh(actele, newnodes, i_encl->first);
        // still not sane!

        // extreme repair: align all normals in one direction
        it = baseface.begin();  // first node is the reference
        std::vector<double> refnormal = node_normals.find(*it)->second;
        CheckNormDir(refnormal, avgnode_normals.find(*it)->second);
        std::vector<double> repairnormal = node_normals.find(*it)->second;
        CheckNormDir(repairnormal, avgnode_normals.find(*it)->second);
        for (it = (baseface.begin()); it != baseface.end(); ++it)
        {
          int repairnodepos = FindPosinVec(*it, baseface) + nnodes / 2;
          int repairnode = actele.at(repairnodepos);
          // vector<double>repairnormal = node_normals.find(*it)->second;
          // CheckNormDir(repairnormal,refnormal);
          if (nd_ts.size() != 0)
          {
            actthickness =
                nd_ts.find(inv_node_pair.find(*it)->second)->second;  // get current node thickness
            if (actthickness < 1E-12)
              std::cout << "Node: " << (inv_node_pair.find(*it)->second) << " has zero thickness!"
                        << std::endl;
          }
          std::vector<double> newcoords =
              ExtrudeNodeCoords(coords.find(*it)->second, actthickness, 1, 1, repairnormal);


          // move node to this aligned position
          // actually due to repairing above all extrusion nodes must be already new and can just be
          // moved
          newnodes.find(repairnode)->second = newcoords;
          coords.find(repairnode)->second = newcoords;  // coords update

          // layer case rework layer elements within enclosing one
          if (numlayers > 1)
          {
            const std::vector<int> layereles = layerstack.find(i_encl->first)->second;
            std::vector<int>::const_iterator i_layerele;
            int i_layer = 1;
            for (i_layerele = layereles.begin(); i_layerele != (layereles.end() - 1); ++i_layerele)
            {
              if (nd_ts.size() != 0)
              {
                actthickness = nd_ts.find(inv_node_pair.find(*it)->second)
                                   ->second;  // get current node thickness
                if (actthickness < 1E-12)
                  std::cout << "Node: " << (inv_node_pair.find(*it)->second)
                            << " has zero thickness!" << std::endl;
              }
              std::vector<double> newlayercoords = ExtrudeNodeCoords(
                  coords.find(*it)->second, actthickness, i_layer, numlayers, repairnormal);
              ++i_layer;
              int layerrepairnode = newconn.find(*i_layerele)->second.at(repairnodepos);
              newnodes.find(layerrepairnode)->second = newlayercoords;
            }
          }
        }
        // PlotEleGmsh(actele, newnodes, i_encl->first);

        int doublerepairedelesign = EleSaneSign(actele, coords);
        if (doublerepairedelesign != initelesign)
          std::cout << "What?! Element still twisted, I give up! elesign=" << doublerepairedelesign
                    << ", reference elesign=" << initelesign << std::endl;

        // PlotEleConnGmsh(obstinates,newnodes);
      }

      repaired_conn.insert(std::pair<int, std::vector<int>>(twistcounter, actele));
    }
  }

  std::cout << "firstcheck: " << firstcheck << ", secondcheck: " << secondcheck << std::endl;
  std::cout << "During repair " << newnodesbyrepair << " new nodes have been created." << std::endl;
  if (repaired_conn.size() != 0) PlotEleConnGmsh(repaired_conn, newnodes, repaired_conn);

  return twistcounter;
}


std::vector<double> EXODUS::ExtrudeNodeCoords(const std::vector<double> basecoords,
    const double distance, const int layer, const int numlayers, const std::vector<double> normal)
{
  std::vector<double> newcoords(3);

  double actdistance = layer * distance / numlayers;

  newcoords[0] = basecoords[0] + actdistance * normal.at(0);
  newcoords[1] = basecoords[1] + actdistance * normal.at(1);
  newcoords[2] = basecoords[2] + actdistance * normal.at(2);
  return newcoords;
}

const std::map<int, std::set<int>> EXODUS::NodeToEleConn(
    const std::map<int, std::vector<int>> ele_conn)
{
  std::map<int, std::set<int>> node_conn;
  std::map<int, std::vector<int>>::const_iterator i_ele;

  // loop all elements for their nodes
  for (i_ele = ele_conn.begin(); i_ele != ele_conn.end(); ++i_ele)
  {
    std::vector<int> elenodes = i_ele->second;
    std::vector<int>::iterator i_node;
    // loop all nodes within element
    for (i_node = elenodes.begin(); i_node < elenodes.end(); ++i_node)
    {
      int nodeid = *i_node;
      // add this ele_id into set of nodeid
      node_conn[nodeid].insert(i_ele->first);
    }
  }
  return node_conn;
}

const std::map<int, std::vector<int>> EXODUS::EleNeighbors(
    const std::map<int, std::vector<int>> ele_conn, const std::map<int, std::set<int>>& node_conn)
{
  std::map<int, std::vector<int>> eleneighbors;
  std::map<int, std::vector<int>>::const_iterator i_ele;

  // loop all elements
  for (i_ele = ele_conn.begin(); i_ele != ele_conn.end(); ++i_ele)
  {
    int acteleid = i_ele->first;
    std::vector<int> actelenodes = i_ele->second;
    std::vector<int>::iterator i_node;
    std::map<int, std::set<int>> elepatch;  // consists of all elements connected by shared nodes
    // loop all nodes within element
    for (i_node = actelenodes.begin(); i_node < actelenodes.end(); ++i_node)
    {
      int nodeid = *i_node;
      // find all elements connected to this node
      const std::set<int> eles = node_conn.find(nodeid)->second;
      // add these eles into patch
      elepatch[nodeid].insert(eles.begin(), eles.end());
    }

    // default case: no neighbors at any edge
    std::vector<int> defaultnbrs(actelenodes.size(), -1);
    eleneighbors.insert(std::pair<int, std::vector<int>>(acteleid, defaultnbrs));
    int edgeid = 0;
    // now select those elements out of the patch which share an edge
    for (i_node = actelenodes.begin(); i_node < actelenodes.end(); ++i_node)
    {
      int firstedgenode = *i_node;
      int secedgenode;
      // edge direction according to order in elenodes, plus last to first
      if (i_node == (actelenodes.end() - 1))
        secedgenode = *actelenodes.begin();
      else
        secedgenode = *(i_node + 1);
      // find all elements connected to the first node
      const std::set<int> firsteles = elepatch.find(firstedgenode)->second;
      // loop over these elements to find the one sharing the secondedgenode
      std::set<int>::const_iterator it;
      for (it = firsteles.begin(); it != firsteles.end(); ++it)
      {
        const int trialele = *it;
        if (trialele != acteleid)
        {
          std::vector<int> neighbornodes = ele_conn.find(trialele)->second;
          bool found = FindinVec(secedgenode, neighbornodes);
          if (found)
          {
            eleneighbors[acteleid].at(edgeid) = trialele;
            break;
          }
        }
      }
      edgeid++;
    }  // end of loop i_nodes
  }    // end of loop all elements
  return eleneighbors;
}

const std::set<int> EXODUS::FreeEdgeNodes(const std::map<int, std::vector<int>>& ele_conn,
    const std::map<int, std::vector<int>>& ele_nbrs)
{
  std::set<int> freenodes;
  std::map<int, std::vector<int>>::const_iterator i_ele;
  for (i_ele = ele_nbrs.begin(); i_ele != ele_nbrs.end(); ++i_ele)
  {
    std::vector<int> actnbrs = i_ele->second;

    // a free edge is a -1 in ele_nbrs
    bool free_edge = FindinVec(-1, actnbrs);
    if (free_edge)
    {
      for (unsigned int i = 0; i < actnbrs.size(); ++i)
      {
        if (actnbrs[i] == -1)
        {
          int actele = i_ele->first;
          const std::vector<int> actelenodes = ele_conn.find(actele)->second;
          // insert pos-node (edgeID is related to first edge node)
          freenodes.insert(actelenodes.at(i));
          // insert second edge node (check if pos is last edge)
          if (i == actelenodes.size() - 1)
            freenodes.insert(actelenodes.at(0));
          else
            freenodes.insert(actelenodes.at(i + 1));
        }
      }
    }
  }

  return freenodes;
}

const std::set<int> EXODUS::FindExtrudedNodes(const std::set<int>& freedgenodes,
    const std::map<int, std::vector<int>>& nodepair, const std::set<int>& ns)
{
  std::set<int> extr_nodes;
  std::set<int>::const_iterator it;
  for (it = ns.begin(); it != ns.end(); ++it)
  {
    if (freedgenodes.find(*it) != freedgenodes.end())
    {
      // extr_nodes.insert(*it);   // insert also basenode
      std::vector<int> extnodes = nodepair.find(*it)->second;
      for (unsigned int i = 0; i < extnodes.size(); ++i) extr_nodes.insert(extnodes[i]);
    }
  }
  return extr_nodes;
}


std::set<int> EXODUS::FreeFaceNodes(
    const std::set<int>& freedgenodes, const std::map<int, std::vector<int>>& nodepair)
{
  std::set<int> freefacenodes;
  std::set<int>::const_iterator it;
  // std::cout << freedgenodes.size() << " zu " << nodepair.size() << std::endl;
  for (it = freedgenodes.begin(); it != freedgenodes.end(); ++it)
  {
    std::vector<int> facenodes = nodepair.find(*it)->second;
    for (unsigned int i = 0; i < facenodes.size(); ++i) freefacenodes.insert(facenodes[i]);
  }
  return freefacenodes;
}

std::vector<double> EXODUS::NodeToAvgNormal(const int node, const std::vector<int> elenodes,
    const std::vector<double> refnormdir, const std::map<int, std::set<int>>& nodetoele,
    const std::map<int, std::vector<int>>& ele_conn, const EXODUS::Mesh& basemesh,
    bool check_norm_scalarproduct)
{
  std::vector<int> myNodeNbrs = FindNodeNeighbors(elenodes, node);

  // calculate normal at node
  std::vector<double> normal = Normal(myNodeNbrs[1], node, myNodeNbrs[0], basemesh);
  if (check_norm_scalarproduct) CheckNormDir(normal, refnormdir);

  // look at neighbors
  std::vector<std::vector<double>> nbr_normals;
  const std::set<int> nbreles = nodetoele.find(node)->second;
  std::set<int>::const_iterator i_nbr;
  for (i_nbr = nbreles.begin(); i_nbr != nbreles.end(); ++i_nbr)
  {
    int nbr = *i_nbr;
    std::vector<int> nbrele = ele_conn.find(nbr)->second;
    if (nbrele != elenodes)
    {  // otherwise it's me, not a neighbor!
      std::vector<int> n_nbrs = FindNodeNeighbors(nbrele, node);
      std::vector<double> nbr_normal = Normal(n_nbrs[1], node, n_nbrs[0], basemesh);
      // check whether nbr_normal points into approx same dir
      // CheckNormDir(nbr_normal,refnormdir);
      nbr_normals.push_back(nbr_normal);
    }
  }

  // average node normal with all neighbors
  std::vector<double> myavgnormal = AverageNormal(normal, nbr_normals);
  if (check_norm_scalarproduct) CheckNormDir(myavgnormal, refnormdir);
  return myavgnormal;
}

std::vector<double> EXODUS::AverageNormal(
    const std::vector<double> n, const std::vector<std::vector<double>> nbr_ns)
{
  // if node has no neighbor avgnormal is normal
  if (nbr_ns.size() == 0) return n;

  // else do averaging
  std::vector<double> avgn = n;
  std::vector<std::vector<double>>::const_iterator i_nbr;

  //  // define lower bound for (nearly) parallel normals
  //  const double para = 1.0e-12;
  //
  //  for(i_nbr=nbr_ns.begin(); i_nbr < nbr_ns.end(); ++i_nbr){
  //    // cross-product with next neighbor normal
  //    std::vector<double> cross(3);
  //    std::vector<double> nbr_n = *i_nbr;
  //    cross[0] =    avgn[1]*nbr_n[2] - avgn[2]*nbr_n[1];
  //    cross[1] = - (avgn[0]*nbr_n[2] - avgn[2]*nbr_n[0]);
  //    cross[2] =    avgn[0]*nbr_n[1] - avgn[1]*nbr_n[0];
  //    double crosslength = cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];
  //
  //    if (crosslength<para){
  //    // if almost parallel do the easy way: average = mean
  //      avgn[0] = 0.5 * (avgn[0] + nbr_n[0]);
  //      avgn[1] = 0.5 * (avgn[1] + nbr_n[1]);
  //      avgn[2] = 0.5 * (avgn[2] + nbr_n[2]);
  //      avgn[0] += nbr_n[0];
  //      avgn[1] += nbr_n[1];
  //      avgn[2] += nbr_n[2];
  //
  //    } else {
  //    // do the Bischoff-Way:
  //      // left length
  //      double leftl = avgn[0]*avgn[0] + avgn[1]*avgn[1] + avgn[2]*avgn[2];
  //      // right length
  //      double rightl = nbr_n[0]*nbr_n[0] + nbr_n[1]*nbr_n[1] + nbr_n[2]*nbr_n[2];
  //      // mean
  //      avgn[0] = 0.5 * (avgn[0] + nbr_n[0]);
  //      avgn[1] = 0.5 * (avgn[1] + nbr_n[1]);
  //      avgn[2] = 0.5 * (avgn[2] + nbr_n[2]);
  //      // mean length
  //      double avgl = avgn[0]*avgn[0] + avgn[1]*avgn[1] + avgn[2]*avgn[2];
  //      // scale by mean of left and right normal
  //      avgn[0] = avgn[0] * 0.5*(leftl+rightl)/avgl;
  //      avgn[1] = avgn[1] * 0.5*(leftl+rightl)/avgl;
  //      avgn[2] = avgn[2] * 0.5*(leftl+rightl)/avgl;
  //    }
  //  } // average with next neighbor


  // new version: order-independent of normals, but without "parallel-check"
  double meanlength = avgn[0] * avgn[0] + avgn[1] * avgn[1] + avgn[2] * avgn[2];
  for (i_nbr = nbr_ns.begin(); i_nbr < nbr_ns.end(); ++i_nbr)
  {
    std::vector<double> nbr_n = *i_nbr;

    // sum a^i:
    avgn[0] += nbr_n[0];
    avgn[1] += nbr_n[1];
    avgn[2] += nbr_n[2];
    // sum |a^i|
    meanlength += nbr_n[0] * nbr_n[0] + nbr_n[1] * nbr_n[1] + nbr_n[2] * nbr_n[2];
  }

  // a^m = 1/n * sum a^i
  avgn[0] = avgn[0] / nbr_ns.size();
  avgn[1] = avgn[1] / nbr_ns.size();
  avgn[2] = avgn[2] / nbr_ns.size();

  // meanlength = 1/n * sum a^i
  meanlength = meanlength / nbr_ns.size();

  // |a^m|
  double am_length = avgn[0] * avgn[0] + avgn[1] * avgn[1] + avgn[2] * avgn[2];

  // a^m in Bischoff-Style (Diss. S. 129 Fig. 8.2(b)
  avgn[0] = avgn[0] * meanlength / am_length;
  avgn[1] = avgn[1] * meanlength / am_length;
  avgn[2] = avgn[2] * meanlength / am_length;

  // unit length:
  double length = sqrt(avgn[0] * avgn[0] + avgn[1] * avgn[1] + avgn[2] * avgn[2]);
  avgn[0] = avgn[0] / length;
  avgn[1] = avgn[1] / length;
  avgn[2] = avgn[2] / length;

  return avgn;
}

std::vector<double> EXODUS::Normal(int head1, int origin, int head2, const EXODUS::Mesh& basemesh)
{
  std::vector<double> normal(3);
  std::vector<double> h1 = basemesh.GetNode(head1);
  std::vector<double> h2 = basemesh.GetNode(head2);
  std::vector<double> o = basemesh.GetNode(origin);

  normal[0] = ((h1[1] - o[1]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[1] - o[1]));
  normal[1] = -((h1[0] - o[0]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[0] - o[0]));
  normal[2] = ((h1[0] - o[0]) * (h2[1] - o[1]) - (h1[1] - o[1]) * (h2[0] - o[0]));

  double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  const double epsilon = 1.E-4;
  if (length > epsilon)
  {
    normal[0] = normal[0] / length;
    normal[1] = normal[1] / length;
    normal[2] = normal[2] / length;
  }
  else
  {  // normal is undefined, vectors seem collinear
    normal.resize(1);
    normal[0] = 0.0;
  }

  return normal;
}

std::vector<double> EXODUS::Normal(
    int head1, const int origin, int head2, const std::map<int, std::vector<double>>& coords)
{
  std::vector<double> normal(3);
  std::vector<double> h1 = coords.find(head1)->second;
  std::vector<double> h2 = coords.find(head2)->second;
  std::vector<double> o = coords.find(origin)->second;

  normal[0] = ((h1[1] - o[1]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[1] - o[1]));
  normal[1] = -((h1[0] - o[0]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[0] - o[0]));
  normal[2] = ((h1[0] - o[0]) * (h2[1] - o[1]) - (h1[1] - o[1]) * (h2[0] - o[0]));

  double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  normal[0] = normal[0] / length;
  normal[1] = normal[1] / length;
  normal[2] = normal[2] / length;

  return normal;
}

void EXODUS::CheckNormDir(std::vector<double>& checkn, const std::vector<double> refn)
{
  // compute scalar product of the two normals
  double scp = checkn[0] * refn.at(0) + checkn[1] * refn.at(1) + checkn[2] * refn.at(2);
  if (scp < 0)
  {
    checkn[0] = -checkn[0];
    checkn[1] = -checkn[1];
    checkn[2] = -checkn[2];
  }
}

std::vector<double> EXODUS::MeanVec(const std::vector<std::vector<double>> baseVecs)
{
  if (baseVecs.size() == 0) dserror("baseVecs empty -> div by 0");
  std::vector<double> mean;
  std::vector<std::vector<double>>::const_iterator i_vec;
  // vector<double>::const_iterator i;
  for (i_vec = baseVecs.begin(); i_vec < baseVecs.end(); ++i_vec)
  {
    const std::vector<double> vec = *i_vec;
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      mean[i] += vec.at(i);
    }
  }
  for (unsigned int i = 0; i < mean.size(); ++i) mean[i] = mean[i] / baseVecs.size();
  return mean;
}

bool EXODUS::FindinVec(const int id, const std::vector<int> vec)
{
  std::vector<int>::const_iterator i;
  for (i = vec.begin(); i < vec.end(); ++i)
    if (*i == id) return true;
  return false;
}

int EXODUS::FindPosinVec(const int id, const std::vector<int> vec)
{
  for (unsigned int i = 0; i < vec.size(); ++i)
    if (vec.at(i) == id) return i;
  return -1;
}

int EXODUS::FindEdgeNeighbor(
    const std::vector<int> nodes, const int actnode, const int wrong_dir_node)
{
  // special case of very first node
  if (nodes.at(0) == actnode)
  {
    if (nodes.at(1) == wrong_dir_node)
      return nodes.back();
    else
      return nodes.at(1);
  }
  else if (nodes.back() == actnode)
  {
    // special case of very last node
    if (nodes.at(0) == wrong_dir_node)
      return nodes.at(nodes.size() - 2);
    else
      return nodes.at(0);
  }
  else
  {
    // case of somewhere in between
    std::vector<int>::const_iterator i;
    for (i = nodes.begin(); i < nodes.end(); ++i)
    {
      if (*i == actnode)
      {
        if (*(i + 1) == wrong_dir_node)
          return *(i - 1);
        else
          return *(i + 1);
      }
    }
  }
  return 0;  // never reached!
}

std::vector<int> EXODUS::FindNodeNeighbors(const std::vector<int> nodes, const int actnode)
{
  std::vector<int> neighbors(2);
  // special case of very first node
  if (nodes.at(0) == actnode)
  {
    neighbors[0] = nodes.back();
    neighbors[1] = nodes.at(1);
  }
  else if (nodes.back() == actnode)
  {
    neighbors[0] = nodes.at(nodes.size() - 2);
    neighbors[1] = nodes.front();
  }
  else
  {
    // case of somewhere in between
    std::vector<int>::const_iterator i;
    for (i = nodes.begin(); i < nodes.end(); ++i)
    {
      if (*i == actnode)
      {
        neighbors[0] = *(i - 1);
        neighbors[1] = *(i + 1);
      }
    }
  }
  return neighbors;
}

std::map<int, std::vector<int>> EXODUS::ExtrusionErrorOutput(const int secedgenode,
    const int todo_counter, const std::set<int>& doneles,
    const std::map<int, std::vector<int>>& ele_conn, const std::set<int>& todo_eleset)
{
  std::map<int, std::vector<int>> leftovers;
  std::cout << "There is a problem with neighbor nodes at node: " << secedgenode << std::endl;
  std::cout << "Current element to extrude is: " << todo_counter << std::endl;
  std::cout << "doneles.size: " << doneles.size() << ", ele_conn.size: " << ele_conn.size();
  std::cout << ", this means that " << ele_conn.size() - doneles.size() << " left to extrude."
            << std::endl;
  std::cout << "todo-eleset: ";
  PrintSet(std::cout, todo_eleset);

  std::set<int>::const_iterator it;
  for (it = todo_eleset.begin(); it != todo_eleset.end(); ++it)
  {
    leftovers.insert(
        std::pair<int, std::vector<int>>(ele_conn.find(*it)->first, ele_conn.find(*it)->second));
  }

  return leftovers;
}

void EXODUS::PlotEleGmsh(
    const std::vector<int> elenodes, const std::map<int, std::vector<double>>& nodes, const int id)
{
  std::stringstream filename;
  filename << "ele_" << id << ".gmsh";
  std::ofstream f_system(filename.str().c_str());
  // ofstream f_system("ele.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Element \" {" << std::endl;
  int numnodes = elenodes.size();
  if (numnodes == 6)
    gmshfilecontent << "SI(";
  else if (numnodes == 8)
    gmshfilecontent << "SH(";
  else if (numnodes == 3)
    gmshfilecontent << "ST(";
  else if (numnodes == 4)
    gmshfilecontent << "SQ(";
  for (unsigned int i = 0; i < elenodes.size(); ++i)
  {
    gmshfilecontent << nodes.find(elenodes.at(i))->second[0] << ",";
    gmshfilecontent << nodes.find(elenodes.at(i))->second[1] << ",";
    gmshfilecontent << nodes.find(elenodes.at(i))->second[2];
    if (i == (elenodes.size() - 1))
      gmshfilecontent << ")";
    else
      gmshfilecontent << ",";
  }
  gmshfilecontent << "{";
  for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << elenodes[i] << ",";
  gmshfilecontent << elenodes.back() << "};" << std::endl;
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}

void EXODUS::PlotStartEleGmsh(const int eleid, const std::vector<int> elenodes,
    const EXODUS::Mesh& basemesh, const int nodeid, const std::vector<double> normal)
{
  std::ofstream f_system("startele.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Start Element \" {" << std::endl;
  int numnodes = elenodes.size();
  if (numnodes == 3)
  {
    gmshfilecontent << "ST(" << basemesh.GetNode(elenodes.at(0))[0] << ","
                    << basemesh.GetNode(elenodes.at(0))[1] << ","
                    << basemesh.GetNode(elenodes.at(0))[2] << ","
                    << basemesh.GetNode(elenodes.at(1))[0] << ","
                    << basemesh.GetNode(elenodes.at(1))[1] << ","
                    << basemesh.GetNode(elenodes.at(1))[2] << ","
                    << basemesh.GetNode(elenodes.at(2))[0] << ","
                    << basemesh.GetNode(elenodes.at(2))[1] << ","
                    << basemesh.GetNode(elenodes.at(2))[2] << ")"
                    << "{" << eleid << "," << eleid << "," << eleid << "};" << std::endl;
  }
  else if (numnodes == 4)
  {
    gmshfilecontent << "SQ(" << basemesh.GetNode(elenodes.at(0))[0] << ","
                    << basemesh.GetNode(elenodes.at(0))[1] << ","
                    << basemesh.GetNode(elenodes.at(0))[2] << ","
                    << basemesh.GetNode(elenodes.at(1))[0] << ","
                    << basemesh.GetNode(elenodes.at(1))[1] << ","
                    << basemesh.GetNode(elenodes.at(1))[2] << ","
                    << basemesh.GetNode(elenodes.at(2))[0] << ","
                    << basemesh.GetNode(elenodes.at(2))[1] << ","
                    << basemesh.GetNode(elenodes.at(2))[2] << ","
                    << basemesh.GetNode(elenodes.at(3))[0] << ","
                    << basemesh.GetNode(elenodes.at(3))[1] << ","
                    << basemesh.GetNode(elenodes.at(3))[2] << ")"
                    << "{" << eleid << "," << eleid << "," << eleid << "," << eleid << "};"
                    << std::endl;
  }
  else
    dserror("numnodes not supported");
  gmshfilecontent << "};" << std::endl;
  gmshfilecontent << "View \" Normal \" {" << std::endl;
  gmshfilecontent << "VP(" << basemesh.GetNode(nodeid)[0] << "," << basemesh.GetNode(nodeid)[1]
                  << "," << basemesh.GetNode(nodeid)[2] << ")"
                  << "{" << normal.at(0) << "," << normal.at(1) << "," << normal.at(2) << "};"
                  << std::endl;
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}

void EXODUS::PlotEleNbrs(const std::vector<int> centerele, const std::vector<int> nbrs,
    const std::map<int, std::vector<int>>& baseconn, const EXODUS::Mesh& basemesh, const int nodeid,
    const std::vector<double> normal, const std::map<int, std::set<int>>& node_conn,
    const std::map<int, std::vector<double>> avg_nn)
{
  std::cout << "centerele: ";
  PrintVec(std::cout, centerele);
  std::cout << "neighboreles: ";
  PrintVec(std::cout, nbrs);
  std::set<int> patchnodes;
  std::ofstream f_system("neighbors.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Neighbors \" {" << std::endl;
  int numnodes = centerele.size();
  if (numnodes == 3)
    gmshfilecontent << "ST(";
  else if (numnodes == 4)
    gmshfilecontent << "SQ(";
  for (unsigned int i = 0; i < centerele.size(); ++i)
  {
    gmshfilecontent << basemesh.GetNode(centerele.at(i))[0] << ",";
    gmshfilecontent << basemesh.GetNode(centerele.at(i))[1] << ",";
    gmshfilecontent << basemesh.GetNode(centerele.at(i))[2];
    if (i == (centerele.size() - 1))
      gmshfilecontent << ")";
    else
      gmshfilecontent << ",";
  }
  gmshfilecontent << "{";
  for (unsigned int i = 0; i < (centerele.size() - 1); ++i) gmshfilecontent << 100 << ",";
  gmshfilecontent << 100 << "};" << std::endl;
  std::vector<int>::const_iterator i_nbr;
  for (unsigned int i_nbr = 0; i_nbr < nbrs.size(); ++i_nbr)
  {
    int eleid = nbrs.at(i_nbr);
    const std::vector<int> elenodes = baseconn.find(eleid)->second;
    std::cout << i_nbr << "th neighbor elenodes: ";
    PrintVec(std::cout, elenodes);
    int numnodes = elenodes.size();
    if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      patchnodes.insert(elenodes.at(i));
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[0] << ",";
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[1] << ",";
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << i_nbr << ",";
    gmshfilecontent << i_nbr << "};" << std::endl;
  }
  std::set<int>::iterator it;
  std::set<int>::iterator it2;
  std::set<int> patcheles;
  for (it = patchnodes.begin(); it != patchnodes.end(); ++it)
  {
    std::set<int> eles = node_conn.find(*it)->second;
    for (it2 = eles.begin(); it2 != eles.end(); ++it2) patcheles.insert(*it2);
  }
  for (it = patcheles.begin(); it != patcheles.end(); ++it)
  {
    const std::vector<int> elenodes = baseconn.find(*it)->second;
    int numnodes = elenodes.size();
    if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[0] << ",";
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[1] << ",";
      gmshfilecontent << basemesh.GetNode(elenodes.at(i))[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << -5 << ",";
    gmshfilecontent << -5 << "};" << std::endl;
  }


  gmshfilecontent << "};" << std::endl;
  gmshfilecontent << "View \" Avg Normals \" {" << std::endl;
  // plot avg node normals
  for (it = patchnodes.begin(); it != patchnodes.end(); ++it)
  {
    if (avg_nn.find(*it) != avg_nn.end())
    {
      gmshfilecontent << "VP(" << basemesh.GetNode(*it)[0] << "," << basemesh.GetNode(*it)[1] << ","
                      << basemesh.GetNode(*it)[2] << ")";
      std::vector<double> actn = avg_nn.find(*it)->second;
      gmshfilecontent << "{" << actn[0] << "," << actn[1] << "," << actn[2] << "};" << std::endl;
    }
  }
  gmshfilecontent << "};" << std::endl;

  gmshfilecontent << "View \" Normal \" {" << std::endl;
  gmshfilecontent << "VP(" << basemesh.GetNode(nodeid)[0] << "," << basemesh.GetNode(nodeid)[1]
                  << "," << basemesh.GetNode(nodeid)[2] << ")"
                  << "{" << normal.at(0) << "," << normal.at(1) << "," << normal.at(2) << "};"
                  << std::endl;
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}


void EXODUS::PlotEleConnGmsh(const std::map<int, std::vector<int>>& conn,
    const std::map<int, std::vector<double>>& nodes, const int plot_eleid)
{
  std::ofstream f_system("extrusionmesh.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Extrusion \" {" << std::endl;
  std::map<int, std::vector<int>>::const_iterator it;
  for (it = conn.begin(); it != conn.end(); ++it)
  {
    int eleid = it->first;
    // if(eleid >= plot_eleid-5){
    //  eleid = -1; // will turn plot_eleid into dark blue in gmsh
    //}
    const std::vector<int> elenodes = it->second;
    int numnodes = elenodes.size();
    if (numnodes == 6)
      gmshfilecontent << "SI(";
    else if (numnodes == 8)
      gmshfilecontent << "SH(";
    else if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      gmshfilecontent << nodes.find(elenodes.at(i))->second[0] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[1] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
    gmshfilecontent << eleid << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}

void EXODUS::PlotEleConnGmsh(const std::map<int, std::vector<int>>& conn,
    const std::map<int, std::vector<double>>& nodes,
    const std::map<int, std::vector<int>>& leftovers)
{
  std::ofstream f_system("extrusionproblems.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" base connectivity \" {" << std::endl;
  std::map<int, std::vector<int>>::const_iterator it;
  for (it = conn.begin(); it != conn.end(); ++it)
  {
    const std::vector<int> elenodes = it->second;
    int numnodes = elenodes.size();
    if (numnodes == 6)
      gmshfilecontent << "SI(";
    else if (numnodes == 8)
      gmshfilecontent << "SH(";
    else if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      gmshfilecontent << nodes.find(elenodes.at(i))->second[0] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[1] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << 0 << ",";
    gmshfilecontent << 0 << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;

  gmshfilecontent << "View \" leftovers \" {" << std::endl;
  for (it = leftovers.begin(); it != leftovers.end(); ++it)
  {
    const std::vector<int> elenodes = it->second;
    int eleid = it->first;
    int numnodes = elenodes.size();
    if (numnodes == 6)
      gmshfilecontent << "SI(";
    else if (numnodes == 8)
      gmshfilecontent << "SH(";
    else if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      // node map starts with 0 but exodus with 1!
      gmshfilecontent << nodes.find(elenodes.at(i))->second[0] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[1] << ",";
      gmshfilecontent << nodes.find(elenodes.at(i))->second[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
    gmshfilecontent << eleid << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;

  f_system << gmshfilecontent.str();
  f_system.close();
}
