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
#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"


using namespace Teuchos;

namespace EXODUS
{

  Mesh SolidShellExtrusion(EXODUS::Mesh& basemesh, double thickness, int layers, int seedid, int gmsh,
          int concat2loose, int diveblocks, const string cline, const vector<double> coordcorr);

  bool CheckExtrusion(const EXODUS::ElementBlock eblock);
  bool CheckExtrusion(const EXODUS::SideSet sideset);
  bool CheckExtrusion(const EXODUS::NodeSet nodeset);
  bool CheckFlatEx(const EXODUS::NodeSet nodeset);

  int RepairTwistedExtrusion(const double thickness,
      const map<int,double>& nd_ts, // map of variable thicknesses at node
      const int numlayers,const int initelesign, int& highestnid,
      const map<int,vector<int> >& encloseconn, const map<int,vector<int> >& layerstack,
      std::map<int,vector<int> >& newconn, std::map<int,vector<double> >& newnodes,map<int,std::set<int> >& ext_node_conn,
      const std::map<int,vector<double> >& avgnode_normals, std::map<int,vector<int> >& node_pair,
      const std::map<int,int>& inv_node_pair);

  vector<double> ExtrudeNodeCoords(const vector<double> basecoords, const double distance,
      const int layer, const int numlayers, const vector<double> normal);

  vector<double> NodeToAvgNormal(const int node,const vector<int> elenodes,const vector<double> refnormdir,
      const std::map<int,std::set<int> >& nodetoele, const std::map<int,vector<int> >& ele_conn, const EXODUS::Mesh& basemesh,
      bool check_norm_scalarproduct = true);

  vector<double> AverageNormal(const vector<double> normal, const vector<vector<double> > nbr_normals);

  vector<double> Normal(int head1, int origin, int head2, const EXODUS::Mesh& basemesh);
  vector<double> Normal(int head1, int origin, int head2,const std::map<int,vector<double> >& coords);

  void CheckNormDir(vector<double>& checkn,const vector<double> refn);

  vector<double> MeanVec(const vector<vector<double> > baseVecs);

  const std::map<int,std::set<int> > NodeToEleConn(const std::map<int,vector<int> > ele_conn);

  const std::map<int,vector<int> > EleNeighbors(const std::map<int,vector<int> >ele_conn, const std::map<int,std::set<int> >& node_conn);

  const std::set<int> FreeEdgeNodes(const std::map<int,vector<int> >& ele_conn, const std::map<int,vector<int> >& ele_nbrs);

  const std::set<int> FindExtrudedNodes(const std::set<int>& freedgenodes, const std::map<int,vector<int> >& nodepair,const std::set<int>& ns);

  std::set<int> FreeFaceNodes(const std::set<int>& freedgenodes, const std::map<int,vector<int> >& nodepair);

  bool FindinVec(const int id,vector<int> vec);
  int  FindPosinVec(const int id,vector<int> vec);

  int FindEdgeNeighbor(const vector<int> nodes, const int actnode, const int wrong_dir_node);

  vector<int> FindNodeNeighbors(const vector<int> nodes, const int centernode);

  enum ExtrusionType{eblock,sideset};

  std::map<int,vector<int> > ExtrusionErrorOutput(const int secedgenode,const int todo_counter,
      const std::set<int>& doneles, const std::map<int,vector<int> >& ele_conn, const std::set<int>& todo_eleset);

  void PlotEleGmsh(const vector<int> elenodes, const std::map<int,vector<double> >& nodes, const int id);
  void PlotStartEleGmsh(const int eleid, const vector<int> elenodes,
      const EXODUS::Mesh& basemesh, const int nodeid, const vector<double> normal);

  void PlotEleNbrs(const vector<int> centerele,const vector<int> nbrs, const std::map<int,vector<int> >& baseconn,
      const EXODUS::Mesh& basemesh, const int nodeid, const vector<double> normal,
      const std::map<int,std::set<int> >& node_conn, const std::map<int,vector<double> >avg_nn);

  void PlotEleConnGmsh(const std::map<int,vector<int> >& conn, const std::map<int,vector<double> >& nodes, const int plot_eleid=-1);
  void PlotEleConnGmsh(const std::map<int,vector<int> >& conn, const std::map<int,vector<double> >& nodes, const std::map<int,vector<int> >& leftovers);

} // of namespace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

