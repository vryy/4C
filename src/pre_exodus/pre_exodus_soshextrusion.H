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

  Mesh SolidShellExtrusion(EXODUS::Mesh& basemesh, double thickness, int layers, int seedid, int gmsh, int concat2loose);
  
  bool CheckExtrusion(const EXODUS::ElementBlock eblock);
  bool CheckExtrusion(const EXODUS::SideSet sideset);
  bool CheckExtrusion(const EXODUS::NodeSet nodeset);
  
  int RepairTwistedExtrusion(const double thickness,const int numlayers,const int initelesign, int& highestnid,
      const map<int,vector<int> >& encloseconn, const map<int,vector<int> >& layerstack,
      map<int,vector<int> >& newconn, map<int,vector<double> >& newnodes,map<int,set<int> >& ext_node_conn,
      const map<int,vector<double> >& avgnode_normals, map<int,vector<int> >& node_pair,
      const map<int,int>& inv_node_pair);
  
  vector<double> ExtrudeNodeCoords(const vector<double> basecoords, const double distance,
      const int layer, const int numlayers, const vector<double> normal);
  
  vector<double> NodeToAvgNormal(const int node,const vector<int> elenodes,const vector<double> refnormdir,
      const map<int,set<int> >& nodetoele, const map<int,vector<int> >& ele_conn, const EXODUS::Mesh& basemesh,
      bool check_norm_scalarproduct = true);
  
  vector<double> AverageNormal(const vector<double> normal, const vector<vector<double> > nbr_normals);
  
  vector<double> Normal(int head1, int origin, int head2, const EXODUS::Mesh& basemesh);
  vector<double> Normal(int head1, int origin, int head2,const map<int,vector<double> >& coords);
  
  void CheckNormDir(vector<double>& checkn,const vector<double> refn);
  
  vector<double> MeanVec(const vector<vector<double> > baseVecs);
  
  const map<int,set<int> > NodeToEleConn(const map<int,vector<int> > ele_conn);
  
  const map<int,vector<int> > EleNeighbors(const map<int,vector<int> >ele_conn, const map<int,set<int> >& node_conn);
  
  const set<int> FreeEdgeNodes(const map<int,vector<int> >& ele_conn, const map<int,vector<int> >& ele_nbrs);
  
  const set<int> FindExtrudedNodes(const set<int>& freedgenodes, const map<int,vector<int> >& nodepair,const set<int>& ns);
  
  set<int> FreeFaceNodes(const set<int>& freedgenodes, const map<int,vector<int> >& nodepair);
  
  bool FindinVec(const int id,vector<int> vec);
  int  FindPosinVec(const int id,vector<int> vec);
  
  int FindEdgeNeighbor(const vector<int> nodes, const int actnode, const int wrong_dir_node);
  
  vector<int> FindNodeNeighbors(const vector<int> nodes, const int centernode);
  
  enum ExtrusionType{eblock,sideset};
  
  void PlotStartEleGmsh(const int eleid, const vector<int> elenodes, 
      const EXODUS::Mesh& basemesh, const int nodeid, const vector<double> normal);
  
  void PlotEleNbrs(const vector<int> centerele,const vector<int> nbrs, const map<int,vector<int> >& baseconn,
      const EXODUS::Mesh& basemesh, const int nodeid, const vector<double> normal,
      const map<int,set<int> >& node_conn, const map<int,vector<double> >avg_nn);
  
  void PlotEleConnGmsh(const map<int,vector<int> >& conn, const map<int,vector<double> >& nodes);

} // of namespace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

