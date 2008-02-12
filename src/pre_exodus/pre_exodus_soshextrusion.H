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
bool CheckExtrusion(const EXODUS::SideSet sideset);

vector<double> ExtrudeNodeCoords(const vector<double> basecoords, const double distance,
    const int layer, const int numlayers, const vector<double> normal);

vector<double> NodeToAvgNormal(const int node,const vector<int> elenodes,const vector<double> refnormdir,
    const map<int,set<int> >& nodetoele, const map<int,vector<int> >& ele_conn, const EXODUS::Mesh& basemesh);

vector<double> AverageNormal(const vector<double> normal, const vector<vector<double> > nbr_normals);

vector<double> Normal(int head1, int origin, int head2, const EXODUS::Mesh& basemesh);

void CheckNormDir(vector<double> checkn,const vector<double> refn);

vector<double> MeanVec(const vector<vector<double> > baseVecs);

const map<int,set<int> > NodeToEleConn(const EXODUS::ElementBlock eblock);
const map<int,set<int> > NodeToEleConn(const map<int,vector<int> > ele_conn);

const map<int,vector<int> > EleNeighbors(const EXODUS::ElementBlock eblock, const map<int,set<int> >& node_conn);
const map<int,vector<int> > EleNeighbors(const map<int,vector<int> >ele_conn, const map<int,set<int> >& node_conn);

bool FindinVec(const int id,vector<int> vec);

int FindEdgeNeighbor(const vector<int> nodes, const int actnode, const int wrong_dir_node);

vector<int> FindNodeNeighbors(const vector<int> nodes, const int centernode);

enum ExtrusionType{eblock,sideset};

} // of namespace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

