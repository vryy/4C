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
  int highestnode = basemesh.GetNumNodes();
  map<int,vector<double> > newnodes;
  
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
    extrudeblock.Print(cout,true);
    
    // loop through all its elements to create new connectivity
    map<int,vector<int> > newconn;
    for (int iele = 0; iele < numele; ++iele) {
      vector<int> actelenodes = extrudeblock.GetEleNodes(iele);
      vector<int> newelenodes_base = RiseNodeIds(highestnode,actelenodes);
      for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
        vector<double> actcoords = basemesh.GetNodeExo(actelenodes[inode]);
        const int newnodeid = ExoToStore(newelenodes_base.at(inode));
        newnodes.insert(pair<int,vector<double> >(newnodeid,actcoords));
      }
      highestnode += actelenodes.size();
      vector<int> newelenodes_ext = RiseNodeIds(highestnode,actelenodes);
      for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
        vector<double> actcoords = basemesh.GetNodeExo(actelenodes[inode]);
        const vector<double> newcoords = ExtrudeNodeCoords(actcoords, thickness);
        const int newnodeid = ExoToStore(newelenodes_ext.at(inode));
        newnodes.insert(pair<int,vector<double> >(newnodeid,newcoords));
      }
      highestnode += actelenodes.size();
      vector<int> newhex8 = newelenodes_base;
      newhex8.insert(newhex8.end(),newelenodes_ext.begin(),newelenodes_ext.end());
      newconn.insert(pair<int,vector<int> >(iele,newhex8));
    }
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

vector<int> EXODUS::RiseNodeIds(const int highestnode, const vector<int> elenodes)
{
  vector<int> newnodes(elenodes.size());
  int newid = highestnode +1;
  for (int i = 0; i < signed(elenodes.size()); ++i){
    newnodes[i] = newid;
    newid ++;
  }
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


#endif
