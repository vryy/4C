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

Mesh SolidShellExtrusion(EXODUS::Mesh basemesh, double thickness, int layers);

bool CheckExtrusion(const EXODUS::ElementBlock eblock);

vector<int> RiseNodeIds(int highestid, map<int,int> pair, const vector<int> elenodes);

vector<double> ExtrudeNodeCoords(const vector<double> basecoords, double distance);

const map<int,set<int> > NodeToEleConn(const EXODUS::ElementBlock eblock);

const map<int,vector<int> > EleNeighbors(const EXODUS::ElementBlock eblock, const map<int,set<int> > node_conn);

bool FindinVec(const int id,vector<int> vec);

int FindEdgeNeighbor(const vector<int> nodes, const int actnode, const int wrong_dir_node);

void PrintMap(ostream& os,const map<int,vector<int> > mymap);
void PrintMap(ostream& os,const map<int,set<int> > mymap);
void PrintVec(ostream& os, const vector<int> actvec);
void PrintSet(ostream& os, const set<int> actset);

} // of namespace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

