#ifndef CUT_MESH_H
#define CUT_MESH_H

#include <list>
#include <map>
#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

#include "cut_boundingbox.H"
#include "cut_point.H"
#include "cut_facet.H"

namespace GEO
{
namespace CUT
{

class Node;
class Edge;
class Side;
class LevelSetSide;
class Element;

class PointPool;

class Point;
class Line;
class Facet;
class VolumeCell;

class BoundaryCell;
class Tri3BoundaryCell;
class Quad4BoundaryCell;

class IntegrationCell;
class Hex8IntegrationCell;
class Tet4IntegrationCell;
class Wedge6IntegrationCell;
class Pyramid5IntegrationCell;


/// Keeps the whole thing. Volume mesh, cut surface, cuts
class Mesh
{
public:

  Mesh( double norm=1, Teuchos::RCP<PointPool> pp=Teuchos::null, bool cutmesh=false );

  Element * CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  Side * CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  Element * CreateTet4( int eid, const std::vector<int> & nids );

  Element * CreatePyramid5( int eid, const std::vector<int> & nids );

  Element * CreateWedge6( int eid, const std::vector<int> & nids );

  Element * CreateHex8( int eid, const std::vector<int> & nids );

  Side * CreateTri3( int sid, const std::vector<int> & nids );

  Side * CreateTri6( int sid, const std::vector<int> & nids );

  Side * CreateQuad4( int sid, const std::vector<int> & nids );

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side );

  Line* NewLine( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element );

  Facet* NewFacet( const std::vector<Point*> & points, Side * side, bool cutsurface );

  VolumeCell* NewVolumeCell( const std::set<Facet*> & facets,
                             const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
                             Element * element );

  Tri3BoundaryCell* NewTri3Cell( const Epetra_SerialDenseMatrix & xyz, VolumeCell * volume, Facet * facet );

  Quad4BoundaryCell* NewQuad4Cell( const Epetra_SerialDenseMatrix & xyz, VolumeCell * volume, Facet * facet );

  Hex8IntegrationCell* NewHex8Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Tet4IntegrationCell* NewTet4Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Tet4IntegrationCell* NewTet4Cell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, VolumeCell * cell );

  Wedge6IntegrationCell* NewWedge6Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Pyramid5IntegrationCell* NewPyramid5Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

//   void SelfCut();

  void Cut( Mesh & mesh, std::set<Element*> & elements_done );

  void Cut( Side & side, const std::set<Element*> & done, std::set<Element*> & elements_done );

  void Cut( LevelSetSide & side );

  void MakeFacets();

  void MakeVolumeCells();

  void FindNodePositions();

  void FindLSNodePositions();

  void FindNodalDOFSets();

  void CreateIntegrationCells();

#ifdef DEBUG
  void TestElementVolume();
#endif

  void Status();

  void PrintFacets();

  void DumpGmsh( std::string name );

  void DumpGmshIntegrationcells( std::string name );

  Node* GetNode( int nid ) const;

  Node* GetNode( int nid, const double * xyz, double lsv=0.0 );

  Node* GetNode( const std::set<int> & nids, const double * xyz );

  Edge* GetEdge( Node* begin, Node* end );

  const std::vector<GEO::CUT::Side*> & GetSides( int sid );

  Side* GetSide( const std::set<int> & nids );

  Side* GetSide( int sid,
                 const std::vector<int> & nids,
                 const CellTopologyData * top_data );

  Element* GetElement( int eid );

  Element* GetElement( int eid,
                       const std::vector<int> & nids,
                       const CellTopologyData & top_data );

  bool WithinBB( const Epetra_SerialDenseMatrix & xyz );

  bool WithinBB( Element & element );

private:

  void DumpGmsh( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype );

  Edge* GetEdge( const std::set<int> & nids,
                 const std::vector<Node*> & nodes,
                 const CellTopologyData & edge_topology );
  Side* GetSide( int sid,
                 const std::set<int> & nids,
                 const std::vector<Node*> & nodes,
                 const std::vector<Edge*> & edges,
                 const CellTopologyData & side_topology );

  bool setup_;
  double norm_;
  Teuchos::RCP<PointPool> pp_;
  BoundingBox bb_;
  bool cutmesh_;

  //std::list<Teuchos::RCP<Point> > points_;
  std::list<Teuchos::RCP<Line > > lines_;
  std::list<Teuchos::RCP<Facet> > facets_;
  std::list<Teuchos::RCP<VolumeCell> > cells_;
  std::list<Teuchos::RCP<BoundaryCell> > boundarycells_;
  std::list<Teuchos::RCP<IntegrationCell> > integrationcells_;

  std::map<int, Teuchos::RCP<Node> > nodes_;
  std::map<std::set<int>, Teuchos::RCP<Edge> > edges_;
  std::map<std::set<int>, Teuchos::RCP<Side> > sides_;
  std::map<int, std::vector<Side*> > cut_sides_;
  std::map<int, Teuchos::RCP<Element> > elements_;

  std::map<std::set<int>, Node*> shadow_nodes_;
  std::list<Teuchos::RCP<Element> > shadow_elements_;
};


}
}

#endif
