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

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
Soshextrusion::Soshextrusion(string exofilename,double thickness,int layers) :
Mesh(exofilename)
{
  for (int i = 0; i < GetNumNodes(); ++i) {
    RCP<PreNode> actnode = GetNode(i);
    actnode->Print(cout);
  }
//  // go through all entities in Mesh
//  for (int i = 0; i < GetNumEntities(); ++i) {
//    RCP<Entity> actEntity = GetEntity(i);
//    //cout << actEntity->GetElementType() << endl;
//    cout << "NumEle: " << actEntity->GetNumEle() << endl;
//    cout << "NumNod: " << actEntity->GetNumNodes() << endl;
//    cout << "NumNodPEle: " << actEntity->GetNumNodpElem() << endl;
//    cout << "EntityID: " << actEntity->GetEntityId() << endl;
//    //cout << "EntityType: " << actEntity->EntityType << endl;
//    switch (actEntity->GetEntityType()){
//    case Entity::elem_blk:{
//      
//      cout << "I am a block" << endl;
//      int blk_id = 1+actEntity->GetEntityId();
//      cout << "entity id in block: " << blk_id << endl;
//      cout << "ExoID: " << GetExoId() << endl;
//      
//      int sizeofconnect = actEntity->GetNumNodpElem()*actEntity->GetNumEle();
//      const int* conn = actEntity->Cont();
//      for (int i = 0; i < sizeofconnect; ++i) {
//        cout << conn[i] << " , ";
//        int inode = conn[i] - 1;
//        RCP<PreNode> actnode = GetNode(inode);
//        actnode->Print(cout);
//        
//      }
//      cout << endl;
//      
//      // create base elements
//      vector<RCP<PreElement> > basels;
//      int numele = actEntity->GetNumEle();
//      basels.resize(numele);
//      for (int i  = 0; i < numele; ++i) {
//        PreElement basele = PreElement(i,PreElement::quad4);
//        int ids[4];
//        ids[0] = conn[i*actEntity->GetNumNodpElem()];
//        ids[1] = conn[i*actEntity->GetNumNodpElem()+1];
//        ids[2] = conn[i*actEntity->GetNumNodpElem()+2];
//        ids[3] = conn[i*actEntity->GetNumNodpElem()+3];
//        basele.SetNodeIds(4,ids);
//        basele.Print(cout);
//        basels[i] = rcp(new PreElement(basele));
//      }
//      for (int i = 0; i < numele; ++i) {
//        RCP<PreElement> safedbasele = basels[i];
//        safedbasele->Print(cout);
//      }
//      break;
//    }
//    default: cout << "I am not a block!" << endl;
//    }
//    
//  }
  
  return;
  
}
/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
Soshextrusion::~Soshextrusion()
{
  return;
}



#endif
