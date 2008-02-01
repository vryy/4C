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
  map<int,vector<double> > newnodes;
  int highestele = basemesh.GetNumEle();
  int highestblock = basemesh.GetNumElementBlocks();
  map<int,EXODUS::ElementBlock> ebs = basemesh.GetElementBlocks();
  map<int,EXODUS::ElementBlock>::const_iterator i_ebs;
  map<int,EXODUS::ElementBlock> extrudeblocks;
  
  // loop through all EBlocks to check for extrusion blocks
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs ){
    bool toextrude = CheckExtrusion(i_ebs->second);
    if (toextrude) extrudeblocks.insert(pair<int,ElementBlock>(i_ebs->first,i_ebs->second));
  }
  // create new Element Block map
  map<int,EXODUS::ElementBlock> neweblocks;
  
  // loop all existing extrudeblocks
  for (i_ebs = extrudeblocks.begin(); i_ebs != extrudeblocks.end(); ++i_ebs){
    // get existing eblock
    EXODUS::ElementBlock extrudeblock = i_ebs->second;
    int numele = extrudeblock.GetNumEle();
    extrudeblock.Print(cout);
    
    // Create Node to Element Connectivity
    const map<int,set<int> > node_conn = NodeToEleConn(extrudeblock);
    
    // Create Element to Element Connectivity (sharing an edge)
    const map<int,vector<int> > ele_neighbor = EleNeighbors(extrudeblock,node_conn);
    PrintMap(cout,ele_neighbor);
    
//    // Create new elements based on the basemesh connectivity
//    map<int,vector<int> > newelecon = CopyConnectedMesh(highestnode,
//                                highestele,basenodes,extrudeblock,ele_neighbor);
    
    // loop through all its elements to create new connectivity  ***************
    map<int,vector<int> > newconn;
    map<int,int> node_pair; // stores new node id with base node id
    
    // here we store a potential set of eles still to extrude
    set<int> todo_eles;
    set<int>::iterator i_todo_eles;
    todo_eles.insert(0); // lets insert the very first one for a start
    
    // for the first element we set up the whole
    vector<int> actelenodes = extrudeblock.GetEleNodes(0);
    vector<int>::const_iterator i_node;
    int newid;
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
    }    
    
    map<int,bool> doneles;
    doneles.insert(pair<int,bool>(0,true));
    
    while (todo_eles.size() > 0){
    //for (int iele = 0; iele < numele; ++iele) {
      
      // find an actele still to do
      i_todo_eles=todo_eles.begin();
      
      // start algo of finding edge neighbors
      // ***********************************
      // **********************************
      
      int actele = *i_todo_eles;
      // get nodes of actual element
      vector<int> actelenodes = extrudeblock.GetEleNodes(actele);
      
      /* get new nodes for "copy" element.
       * Here, checking of already copied nodes happens and the corresp map is filled */
      //vector<int> newelenodes_base = RiseNodeIds(highestnid,node_pair,actelenodes);
      vector<int> newelenodes(actelenodes.size());
      map<int,int>::iterator it;
      vector<int>::const_iterator i_node;
      for (i_node=actelenodes.begin(); i_node < actelenodes.end(); ++i_node){
        // check whether this node has not been done
        if (node_pair.find(*i_node)==node_pair.end()){
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
          
        } else newid = node_pair.find(*i_node)->second;
        // insert node into element
        newelenodes.push_back(newid);
      }
      
      
//      for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
//        vector<double> actcoords = basemesh.GetNodeExo(actelenodes[inode]);
//        const int newnodeid = ExoToStore(newelenodes_base.at(inode));
//        newnodes.insert(pair<int,vector<double> >(newnodeid,actcoords));
//      }
//      //highestnode += actelenodes.size();
//      //vector<int> newelenodes_ext = RiseNodeIds(highestnode,actelenodes);
//      for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
//        vector<double> actcoords = basemesh.GetNodeExo(actelenodes[inode]);
//        const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
//        const int newnodeid = ExoToStore(newelenodes_ext.at(inode));
//        newnodes.insert(pair<int,vector<double> >(newnodeid,newcoords));
//      }
//      highestnode += actelenodes.size();
//      vector<int> newhex8 = newelenodes_base;
//      newhex8.insert(newhex8.end(),newelenodes_ext.begin(),newelenodes_ext.end());
//      newconn.insert(pair<int,vector<int> >(iele+highestele, newhex8));
    }
    // *************************************************************************
    
    string newname = "anotherblock";
    EXODUS::ElementBlock neweblock(ElementBlock::hex8,newconn,newname);
    neweblock.Print(cout, true);
    int newblockID = highestblock + i_ebs->first;
    neweblocks.insert(pair<int,EXODUS::ElementBlock>(newblockID,neweblock));
  }
  
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
          if (found) {eleneighbors[acteleid].push_back(trialele); break;}
        }
      }
    }
  }
  return eleneighbors;
}

bool EXODUS::FindinVec(const int id, const vector<int> vec)
{
  vector<int>::const_iterator i;
  for(i=vec.begin(); i<vec.end(); ++i) if (*i == id) return true;
  return false;
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


#endif
