/*----------------------------------------------------------------------*/
/*! \file

\brief solid-shell body creation by extruding surface


\level 1

Here everything related with solid-shell body extrusion
*/
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_SOSHEXTRUSION_H_
#define PRE_EXODUS_SOSHEXTRUSION_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"



namespace EXODUS
{
  Mesh SolidShellExtrusion(EXODUS::Mesh& basemesh, double thickness, int layers, int seedid,
      int gmsh, int concat2loose, int diveblocks, const std::string cline,
      const std::vector<double> coordcorr);

  bool CheckExtrusion(const EXODUS::ElementBlock eblock);
  bool CheckExtrusion(const EXODUS::SideSet sideset);
  bool CheckExtrusion(const EXODUS::NodeSet nodeset);
  bool CheckFlatEx(const EXODUS::NodeSet nodeset);

  int RepairTwistedExtrusion(const double thickness,
      const std::map<int, double>& nd_ts,  // map of variable thicknesses at node
      const int numlayers, const int initelesign, int& highestnid,
      const std::map<int, std::vector<int>>& encloseconn,
      const std::map<int, std::vector<int>>& layerstack, std::map<int, std::vector<int>>& newconn,
      std::map<int, std::vector<double>>& newnodes, std::map<int, std::set<int>>& ext_node_conn,
      const std::map<int, std::vector<double>>& avgnode_normals,
      std::map<int, std::vector<int>>& node_pair, const std::map<int, int>& inv_node_pair);

  std::vector<double> ExtrudeNodeCoords(const std::vector<double> basecoords, const double distance,
      const int layer, const int numlayers, const std::vector<double> normal);

  std::vector<double> NodeToAvgNormal(const int node, const std::vector<int> elenodes,
      const std::vector<double> refnormdir, const std::map<int, std::set<int>>& nodetoele,
      const std::map<int, std::vector<int>>& ele_conn, const EXODUS::Mesh& basemesh,
      bool check_norm_scalarproduct = true);

  std::vector<double> AverageNormal(
      const std::vector<double> normal, const std::vector<std::vector<double>> nbr_normals);

  std::vector<double> Normal(int head1, int origin, int head2, const EXODUS::Mesh& basemesh);
  std::vector<double> Normal(
      int head1, int origin, int head2, const std::map<int, std::vector<double>>& coords);

  void CheckNormDir(std::vector<double>& checkn, const std::vector<double> refn);

  std::vector<double> MeanVec(const std::vector<std::vector<double>> baseVecs);

  const std::map<int, std::set<int>> NodeToEleConn(const std::map<int, std::vector<int>> ele_conn);

  const std::map<int, std::vector<int>> EleNeighbors(const std::map<int, std::vector<int>> ele_conn,
      const std::map<int, std::set<int>>& node_conn);

  const std::set<int> FreeEdgeNodes(const std::map<int, std::vector<int>>& ele_conn,
      const std::map<int, std::vector<int>>& ele_nbrs);

  const std::set<int> FindExtrudedNodes(const std::set<int>& freedgenodes,
      const std::map<int, std::vector<int>>& nodepair, const std::set<int>& ns);

  std::set<int> FreeFaceNodes(
      const std::set<int>& freedgenodes, const std::map<int, std::vector<int>>& nodepair);

  bool FindinVec(const int id, std::vector<int> vec);
  int FindPosinVec(const int id, std::vector<int> vec);

  int FindEdgeNeighbor(const std::vector<int> nodes, const int actnode, const int wrong_dir_node);

  std::vector<int> FindNodeNeighbors(const std::vector<int> nodes, const int centernode);

  enum ExtrusionType
  {
    eblock,
    sideset
  };

  std::map<int, std::vector<int>> ExtrusionErrorOutput(const int secedgenode,
      const int todo_counter, const std::set<int>& doneles,
      const std::map<int, std::vector<int>>& ele_conn, const std::set<int>& todo_eleset);

  void PlotEleGmsh(const std::vector<int> elenodes, const std::map<int, std::vector<double>>& nodes,
      const int id);
  void PlotStartEleGmsh(const int eleid, const std::vector<int> elenodes,
      const EXODUS::Mesh& basemesh, const int nodeid, const std::vector<double> normal);

  void PlotEleNbrs(const std::vector<int> centerele, const std::vector<int> nbrs,
      const std::map<int, std::vector<int>>& baseconn, const EXODUS::Mesh& basemesh,
      const int nodeid, const std::vector<double> normal,
      const std::map<int, std::set<int>>& node_conn,
      const std::map<int, std::vector<double>> avg_nn);

  void PlotEleConnGmsh(const std::map<int, std::vector<int>>& conn,
      const std::map<int, std::vector<double>>& nodes, const int plot_eleid = -1);
  void PlotEleConnGmsh(const std::map<int, std::vector<int>>& conn,
      const std::map<int, std::vector<double>>& nodes,
      const std::map<int, std::vector<int>>& leftovers);

}  // namespace EXODUS

#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
