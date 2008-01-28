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
  
  map<int,ElementBlock> ebs = GetElementBlocks();
  map<int,ElementBlock>::const_iterator i_ebs;
  map<int,ElementBlock> extrudeblocks;
  
  // loop through all EBlocks to check for extrusion blocks
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs ){
    bool toextrude = CheckExtrusion(i_ebs->second);
    if (toextrude){
      ElementBlock extrudeblock = i_ebs->second;
      extrudeblock.Print(cout,true);
      int numele = extrudeblock.GetNumEle();
      //map<int,vector<int> > conn = extrudeblock.GetEleConn();
      for (int iele = 0; iele < numele; ++iele) {
        vector<int> actelenodes = extrudeblock.GetEleNodes(iele);
        for (int inode = 0; inode < signed(actelenodes.size()); ++inode) {
          RCP<PreNode> actnode = GetNode(actelenodes[inode]);
          actnode->Print(cout);
        }
        
      }
    }
    
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


#endif
