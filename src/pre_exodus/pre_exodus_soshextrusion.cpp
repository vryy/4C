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
#include "pre_node.H"
#include "pre_element.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::Soshextrusion::Soshextrusion(string exofilename,double thickness,int layers) :
Mesh(exofilename)
{
  for (int i = 1; i <= GetNumNodes(); ++i) {
    RCP<PreNode> actnode = GetNode(i);
    actnode->Print(cout, true);
  }
  
  int highestnode = GetNumNodes();
  
  map<int,ElementBlock> ebs = GetElementBlocks();
  map<int,ElementBlock>::const_iterator i_ebs;
  map<int,ElementBlock> extrudeblocks;
  
  // loop through all EBlocks to check for extrusion blocks
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs ){
    bool toextrude = CheckExtrusion(i_ebs->second);
    if (toextrude) extrudeblocks.insert(pair<int,ElementBlock>(i_ebs->first,i_ebs->second));
  }
  
  // now work with extrusion blocks
  
  // create new Element Block map
  map<int,ElementBlock> neweblocks;
  
  // loop all existing extrudeblocks
  for (i_ebs = extrudeblocks.begin(); i_ebs != extrudeblocks.end(); ++i_ebs){
    // get existing eblock
    ElementBlock extrudeblock = i_ebs->second;
    int numele = extrudeblock.GetNumEle();
    extrudeblock.Print(cout,true);
    
    // loop through all its elements to create new connectivity
    map<int,vector<int> > newconn;
    for (int iele = 0; iele < numele; ++iele) {
      vector<int> actelenodes = extrudeblock.GetEleNodes(iele);
      vector<int> newelenodes = RiseNodes(highestnode,actelenodes);
      highestnode += actelenodes.size();
      newconn.insert(pair<int,vector<int> >(iele,newelenodes));
      
      // classwork
      PreElement copyele(iele,extrudeblock.GetShape());
      copyele.SetNodeIds(actelenodes.size(), &actelenodes[0]);
      //copyele.Print(cout);
      for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
        RCP<EXODUS::PreNode> actnode = GetNode(actelenodes[inode]);
        int newid = newelenodes[inode];
        //PreNode* newnode = ExtrapolateNode(&actnode, newid);
        //actnode->Print(cout);
      }
    }
    string newname = "anotherblock";
    ElementBlock neweblock(extrudeblock.GetShape(),newconn,newname);
    neweblock.Print(cout, true);
    neweblocks.insert(pair<int,ElementBlock>(i_ebs->first,neweblock));
  }
    

  
  return;
  
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::Soshextrusion::~Soshextrusion()
{
  return;
}

bool EXODUS::Soshextrusion::CheckExtrusion(const ElementBlock eblock)
{
  const ElementBlock::Shape myshape = eblock.GetShape();
  const string myname = eblock.GetName();
  if (myshape == ElementBlock::shell4)
    if (myname.compare(0,7,"extrude") == 0) return true;
  return false;
}

vector<int> EXODUS::Soshextrusion::RiseNodes(const int highestnode, const vector<int> elenodes)
{
  vector<int> newnodes(elenodes.size());
  int newid = highestnode +1;
  for (int i = 0; i < signed(elenodes.size()); ++i){
    newnodes[i] = newid;
    newid ++;
  }
  return newnodes;
}

EXODUS::PreNode* EXODUS::Soshextrusion::ExtrapolateNode(PreNode& actnode, const int newid)
{
  EXODUS::PreNode* newnode = new PreNode(actnode);
  return newnode;
}


#endif
