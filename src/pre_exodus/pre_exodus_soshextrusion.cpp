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
  
  cout << num_entities_ << endl;
  
  cout << "hello" << endl;
  
  for (int i = 0; i < num_entities_; ++i) {
    RCP<Entity> actEntity = GetEntity(i);
    actEntity->Print(cout);
    
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
