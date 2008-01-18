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
#ifdef EXODUS
#include "pre_exodus_soshextrusion.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
Soshextrusion::Soshextrusion(string exofilename,double thickness,int layers) :
Mesh(exofilename)
{
  // go through all entities in Mesh
  for (int i = 0; i < num_entities_; ++i) {
    RCP<Entity> actEntity = GetEntity(i);
    cout << actEntity->GetElementType() << endl;
    cout << "NumEle: " << actEntity->GetNumEle() << endl;
    cout << "NumNod: " << actEntity->GetNumNodes() << endl;
    cout << "NumNodPEle: " << actEntity->GetNumNodpElem() << endl;
    cout << "EntityID: " << actEntity->GetEntityId() << endl;
    //cout << "EntityType: " << actEntity->EntityType << endl;
    switch (actEntity->elem_blk){
    case Entity::elem_blk:{
      cout << "I am a block" << endl;
      int blk_id = 1+actEntity->GetEntityId();
      cout << "entity id in block: " << blk_id << endl;
      cout << "ExoID: " << exoid_ << endl;
      
      int sizeofconnect = actEntity->GetNumNodpElem()*actEntity->GetNumEle();
      //int* connect;
      int connect[sizeofconnect];
      //int connect[int(actEntity->GetNumEle())][int(actEntity->GetNumNodpElem())];
      int error;
      int elem_blk_ids[20];
      error = ex_get_elem_blk_ids (exoid_, elem_blk_ids);
      if (error != 0) dserror("exo error returned");
      cout << "eleblockids: " << (*elem_blk_ids) << endl;
      error = ex_get_elem_conn(exoid_,blk_id,connect);
      if (error != 0) dserror("exo error returned");
      for (int i = 0; i < sizeofconnect; ++i) {
        cout << connect[i] << "," << endl;
      }
      break;
    }
    default: cout << "I am not a block!" << endl;
    }
    
  }
  
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
