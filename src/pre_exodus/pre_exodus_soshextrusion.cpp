/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_soshextrusion.cpp

\brief solid-shell body creation by extruding surface 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with solid-shell body extrusion
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#include "pre_exodus_soshextrusion.H"
#include "pre_exodus_reader.H"
#include "pre_node.H"
#include "pre_element.H"

using namespace std;
using namespace Teuchos;


/* Method to extrude a surface to become a volumetric body */
EXODUS::Mesh EXODUS::SolidShellExtrusion(EXODUS::Mesh basemesh, double thickness, int layers)
{
  int highestnid = basemesh.GetNumNodes() +1;
  map<int,vector<double> > newnodes;          // here the new nodes ar stored
  map<int,EXODUS::ElementBlock> neweblocks;   // here the new EBlocks are stores
  int highestblock = basemesh.GetNumElementBlocks();
  map<int,EXODUS::ElementBlock> ebs = basemesh.GetElementBlocks();
  map<int,EXODUS::ElementBlock>::const_iterator i_ebs;
  map<int,EXODUS::ElementBlock> extrudeblocks;
  
  // loop through all EBlocks to check for extrusion blocks
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs ){
    bool toextrude = CheckExtrusion(i_ebs->second);
    if (toextrude) extrudeblocks.insert(pair<int,ElementBlock>(i_ebs->first,i_ebs->second));
  }
  
  // loop all existing extrudeblocks
  for (i_ebs = extrudeblocks.begin(); i_ebs != extrudeblocks.end(); ++i_ebs){
    // get existing eblock
    EXODUS::ElementBlock extrudeblock = i_ebs->second;
    extrudeblock.Print(cout);
    
    // Create Node to Element Connectivity
    const map<int,set<int> > node_conn = NodeToEleConn(extrudeblock);
    
    // Create Element to Element Connectivity (sharing an edge)
    const map<int,vector<int> > ele_neighbor = EleNeighbors(extrudeblock,node_conn);
    PrintMap(cout,ele_neighbor);
    
    // loop through all its elements to create new connectivity  ***************
    map<int,vector<int> > newconn;
    int newele = 0;
    map<int,int> node_pair; // stores new node id with base node id
    
    // set of elements already done
    set<int> doneles;
    set<int>::iterator i_doneles;
    
   // here we store a potential set of eles still to extrude
    set<int> todo_eles;
    set<int>::iterator i_todo_eles;
    todo_eles.insert(0); // lets insert the very first one for a start
    
    // for the first element we set up everything *****************************
    vector<int> actelenodes = extrudeblock.GetEleNodes(0);
    vector<int>::const_iterator i_node;
    int newid;
    // create a new element
    vector<int> newelenodes;
    for (i_node=actelenodes.begin(); i_node < actelenodes.end(); ++i_node){
      newid = highestnid; ++ highestnid;
      
      // place new node at new position
      vector<double> actcoords = basemesh.GetNodeExo(*i_node); //curr position
      // calculate new position
      const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
      int newExoNid = ExoToStore(newid);
      // put new coords into newnode map
      newnodes.insert(pair<int,vector<double> >(newExoNid,newcoords));
      
      // put new node into map of OldNodeToNewNode
      node_pair.insert(std::pair<int,int>(*i_node,newid));
      // insert node into element
      newelenodes.push_back(newid);
    }    
    doneles.insert(0);    // the first element is done ************************
    // insert new element into new connectivity
    newconn.insert(pair<int,vector<int> >(newele,newelenodes));
    PrintVec(cout,newelenodes);
    newelenodes.clear();
    ++ newele;
    
    
    while (todo_eles.size() > 0){ // start of going through ele neighbors///////
      
      // find an actele still to do
      PrintSet(cout,todo_eles);
      i_todo_eles=todo_eles.begin();
      int actele = *i_todo_eles;
      // delete actele from todo_eles list
      todo_eles.erase(i_todo_eles);
      
      // get nodes of actual element
      vector<int> actelenodes = extrudeblock.GetEleNodes(actele);
      
      // get edge neighbors
      vector<int> actneighbors = ele_neighbor.at(actele);
      
      vector<int>::const_iterator i_nbr;
      int edge = 0; // edge counter
      for (i_nbr = actneighbors.begin(); i_nbr < actneighbors.end(); ++ i_nbr){
        int actneighbor = *i_nbr;
        // check for undone neighbor
        i_doneles = doneles.find(actneighbor);
        if ((actneighbor != -1) && (i_doneles == doneles.end())){
          // lets extrude *****************************************************
          
          // get nodepair of actual edge
          int firstedgenode = actelenodes.at(edge);
          int secedgenode;
          // switch in case of last to first node edge
          if (edge == signed(actelenodes.size())) secedgenode = actelenodes.at(0);
          else secedgenode = actelenodes.at(edge+1);
          
          // create a new element
          vector<int> newelenodes;
          
          /* the new elements orientation is opposite the current one
           * therefore the first node is secedgenode */
          
          // check if new node already exists
          if (node_pair.find(secedgenode)==node_pair.end()){
            // lets create a new node
            newid = highestnid; ++ highestnid;
            // place new node at new position
            vector<double> actcoords = basemesh.GetNodeExo(secedgenode); //curr position
            // calculate new position
            const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
            int newExoNid = ExoToStore(newid);
            // put new coords into newnode map
            newnodes.insert(pair<int,vector<double> >(newExoNid,newcoords));
            
            // put new node into map of OldNodeToNewNode
            node_pair.insert(std::pair<int,int>(secedgenode,newid));
            
          } else newid = node_pair.find(secedgenode)->second;
          
          // insert node into element
          newelenodes.push_back(newid);
          
          /* the new elements orientation is opposite the current one
           * therefore the second node is firstedgenode */
          
          // check if new node already exists
          if (node_pair.find(firstedgenode)==node_pair.end()){
            // lets create a new node
            newid = highestnid; ++ highestnid;
            // place new node at new position
            vector<double> actcoords = basemesh.GetNodeExo(firstedgenode); //curr position
            // calculate new position
            const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
            int newExoNid = ExoToStore(newid);
            // put new coords into newnode map
            newnodes.insert(pair<int,vector<double> >(newExoNid,newcoords));
            
            // put new node into map of OldNodeToNewNode
            node_pair.insert(std::pair<int,int>(firstedgenode,newid));
            
          } else newid = node_pair.find(firstedgenode)->second;
          
          // insert node into element
          newelenodes.push_back(newid);
          
          /* the new elements orientation is opposite the current one
           * therefore the remaining nodes are gained from the neighbor ele */
          
          vector<int> actnbrnodes = extrudeblock.GetEleNodes(actneighbor);
          int thirdnode = FindEdgeNeighbor(actnbrnodes,firstedgenode,secedgenode);
          
          // check if new node already exists
          if (node_pair.find(thirdnode)==node_pair.end()){
            // lets create a new node
            newid = highestnid; ++ highestnid;
            // place new node at new position
            vector<double> actcoords = basemesh.GetNodeExo(thirdnode); //curr position
            // calculate new position
            const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
            int newExoNid = ExoToStore(newid);
            // put new coords into newnode map
            newnodes.insert(pair<int,vector<double> >(newExoNid,newcoords));
            
            // put new node into map of OldNodeToNewNode
            node_pair.insert(std::pair<int,int>(thirdnode,newid));
            
          } else newid = node_pair.find(thirdnode)->second;
          
          // insert node into element
          newelenodes.push_back(newid);
         
          if (actelenodes.size() > 3){ // in case of not being a tri3
            
            /* the new elements orientation is opposite the current one
             * therefore the remaining nodes are gained from the neighbor ele */

            int fourthnode = FindEdgeNeighbor(actnbrnodes,thirdnode,firstedgenode);
            
            // check if new node already exists
            if (node_pair.find(fourthnode)==node_pair.end()){
              // lets create a new node
              newid = highestnid; ++ highestnid;
              // place new node at new position
              vector<double> actcoords = basemesh.GetNodeExo(fourthnode); //curr position
              // calculate new position
              const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
              int newExoNid = ExoToStore(newid);
              // put new coords into newnode map
              newnodes.insert(pair<int,vector<double> >(newExoNid,newcoords));
              
              // put new node into map of OldNodeToNewNode
              node_pair.insert(std::pair<int,int>(fourthnode,newid));
              
            } else newid = node_pair.find(fourthnode)->second;
            
            // insert node into element
            newelenodes.push_back(newid);
          }
          
          // insert actneighbor element into new connectivity map
          newconn.insert(pair<int,vector<int> >(newele,newelenodes));
          newelenodes.clear();
          ++ newele;
          
          // insert actneighbor into done elements
          doneles.insert(actneighbor);
          
          // neighbor eles are possible next "center" eles
          todo_eles.insert(actneighbor);  
          
        }// end of if undone->extrude this neighbor, next neighbor *************
        ++ edge; // next element edge
        
      }// end of this "center" - element ///////////////////////////////////////
      
    }// end of extruding all elements in block
    // debug
    PrintMap(cout,newconn);
   
    // create new Element Block
    std::ostringstream blockname;
    blockname << "extrude" << i_ebs->first;
    EXODUS::ElementBlock neweblock(ElementBlock::shell4,newconn,blockname.str());
    int newblockID = highestblock + i_ebs->first;
    neweblocks.insert(pair<int,EXODUS::ElementBlock>(newblockID,neweblock));
    
  } // end of extruding elementblocks
  
  string newtitle = "extrusion";
  map<int,EXODUS::NodeSet> emptynodeset;
  map<int,EXODUS::SideSet> emptysideset;

  EXODUS::Mesh extruded_mesh(basemesh,newnodes,neweblocks,emptynodeset,emptysideset,newtitle);
  
  return extruded_mesh;
}


bool EXODUS::CheckExtrusion(const EXODUS::ElementBlock eblock)
{
  const EXODUS::ElementBlock::Shape myshape = eblock.GetShape();
  const string myname = eblock.GetName();
  if (myshape == EXODUS::ElementBlock::shell4)
    if (myname.compare(0,7,"extrude") == 0) return true;
  return false;
}

vector<int> EXODUS::RiseNodeIds(int highestid, map<int,int> pair, const vector<int> elenodes)
{
  vector<int> newnodes(elenodes.size());
//  int newid;
//  map<int,int>::iterator it;
//  vector<int>::const_iterator i_node;
//  for (i_node=elenodes.begin(); i_node < elenodes.end(); ++i_node){
//    if (pair.find(*i_node)==pair.end()){ newid = highestid; ++ highestid;}
//    else newid = pair.find(*i_node)->second;
//    newnodes.push_back(newid);
//    pair.insert(std::pair<int,int>(*i_node,newid));
//  }
  return newnodes;
}

vector<double> EXODUS::ExtrudeNodeCoords(const vector<double> basecoords, double distance)
{
  vector<double> newcoords(3);
  newcoords[0] = basecoords[0];
  newcoords[1] = basecoords[1];
  newcoords[2] = basecoords[2] + distance;
  return newcoords;
}

const map<int,set<int> > EXODUS::NodeToEleConn(const EXODUS::ElementBlock eblock)
{
  map<int,set<int> > node_conn;
  const map<int,vector<int> > ele_conn = eblock.GetEleConn();
  map<int,vector<int> >::const_iterator i_ele;
  
  // loop all elements for their nodes
  for (i_ele = ele_conn.begin(); i_ele != ele_conn.end(); ++i_ele){
    vector<int> elenodes = i_ele->second;
    vector<int>::iterator i_node;
    // loop all nodes within element
    for (i_node = elenodes.begin(); i_node < elenodes.end(); ++i_node){
      int nodeid = *i_node;
      // add this ele_id into set of nodeid
      node_conn[nodeid].insert(i_ele->first);
    }
  }
  return node_conn;
}

const map<int,vector<int> > EXODUS::EleNeighbors(EXODUS::ElementBlock eblock, const map<int,set<int> > node_conn)
{
  map<int,vector<int> > eleneighbors;
  const map<int,vector<int> > ele_conn = eblock.GetEleConn();
  map<int,vector<int> >::const_iterator i_ele;

  // loop all elements
  for (i_ele = ele_conn.begin(); i_ele != ele_conn.end(); ++i_ele){
    int acteleid = i_ele->first;
    vector<int> actelenodes = i_ele->second;
    vector<int>::iterator i_node;
    map<int,set<int> > elepatch; // consists of all elements connected by shared nodes
    // loop all nodes within element
    for (i_node = actelenodes.begin(); i_node < actelenodes.end(); ++i_node){
      int nodeid = *i_node;
      // find all elements connected to this node
      const set<int> eles = node_conn.find(nodeid)->second;
      // add these eles into patch
      elepatch[nodeid].insert(eles.begin(),eles.end());
    }
    
    // default case: no neighbors at any edge
    vector<int> defaultnbrs(actelenodes.size(),-1);
    eleneighbors.insert(pair<int,vector<int> >(acteleid,defaultnbrs));
    int edgeid = 0;
    // now select those elements out of the patch which share an edge
    for (i_node = actelenodes.begin(); i_node < actelenodes.end(); ++i_node){
      int firstedgenode = *i_node;
      int secedgenode;
      // edge direction according to order in elenodes, plus last to first
      if (i_node == (actelenodes.end()-1)) secedgenode = *actelenodes.begin();
      else secedgenode = *(i_node + 1);
      // find all elements connected to the first node
      const set<int> firsteles = elepatch.find(firstedgenode)->second;
      // loop over these elements to find the one sharing the secondedgenode
      set<int>::const_iterator it;
      for(it = firsteles.begin(); it != firsteles.end(); ++it){
        const int trialele = *it;
        if (trialele != acteleid){
          vector<int> neighbornodes = ele_conn.find(trialele)->second;
          bool found = FindinVec(secedgenode,neighbornodes);
          if (found) {
            eleneighbors[acteleid].at(edgeid) = trialele;
            break;
          }
        }
      }
      edgeid ++;
    } // end of loop i_nodes
  } // end of loop all elements
  return eleneighbors;
}

bool EXODUS::FindinVec(const int id, const vector<int> vec)
{
  vector<int>::const_iterator i;
  for(i=vec.begin(); i<vec.end(); ++i) if (*i == id) return true;
  return false;
}

int EXODUS::FindEdgeNeighbor(const vector<int> nodes, const int actnode, const int wrong_dir_node)
{
  // special case of very first node
  if (nodes.at(0) == actnode){
    if (nodes.at(1) == wrong_dir_node) return nodes.at(nodes.size());
    else return nodes.at(1); 
  }
  // special case of very last node
  if (nodes.back() == actnode){
    if (nodes.at(0) == wrong_dir_node) return nodes.at(nodes.size()-2);
    else return nodes.at(0);
  }
  // case of somewhere in between
  vector<int>::const_iterator i;
  for(i=nodes.begin(); i<nodes.end(); ++i){
    if (*i == actnode){
      if (*(i+1) == wrong_dir_node) return *(i-1);
      else return *(i+1);
    }
  }
  return 0; // never reached!
}

void EXODUS::PrintMap(ostream& os,const map<int,vector<int> > mymap)
{
  map<int,vector<int> >::const_iterator iter;
  for(iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
      os << iter->first << ": ";
      vector<int> actvec = iter->second;
      vector<int>::iterator i;
      for (i=actvec.begin(); i<actvec.end(); ++i) {
        os << *i << ",";
      }
      os << endl;
  }
}

void EXODUS::PrintMap(ostream& os,const map<int,set<int> > mymap)
{
  map<int,set<int> >::const_iterator iter;
  for(iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
      os << iter->first << ": ";
      set<int> actset = iter->second;
      set<int>::iterator i;
      for (i=actset.begin(); i != actset.end(); ++i) {
        os << *i << ",";
      }
      os << endl;
  }
}

void EXODUS::PrintVec(ostream& os, const vector<int> actvec)
{
  vector<int>::const_iterator i;
  for (i=actvec.begin(); i<actvec.end(); ++i) {
    os << *i << ",";
  }
  os << endl;
}

void EXODUS::PrintSet(ostream& os, const set<int> actset)
{
  set<int>::iterator i;
  for (i=actset.begin(); i != actset.end(); ++i) {
    os << *i << ",";
  }
  os << endl;
  
}


#endif
