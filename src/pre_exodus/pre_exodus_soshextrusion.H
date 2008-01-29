/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_soshextrusion.H

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
#ifndef PRE_EXODUS_SOSHEXTRUSION_H_
#define PRE_EXODUS_SOSHEXTRUSION_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;
using namespace Teuchos;

namespace EXODUS
{

/*!
\class Soshextrusion

\brief Soshextrusion handles the extrusion of a solid-shell mesh body based
on a surface mesh


\author maf (frenzel@lnm.mw.tum.de)
*/
class Soshextrusion : public Mesh
{
public:
  //! constructor
  Soshextrusion(string exofilename, double thickness, int layers);
  
  //! destructor
  virtual ~Soshextrusion();
  
  bool CheckExtrusion(const ElementBlock eblock);
  
  vector<int> RiseNodes(const int highestnode, const vector<int> elenodes);
  
  //PreNode* ExtrapolateNode(PreNode& actnode, const int newid);
};

} // of namepsace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

