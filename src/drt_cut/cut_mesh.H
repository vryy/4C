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
class Options;

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

  Mesh( const Options & options, double norm=1, Teuchos::RCP<PointPool> pp=Teuchos::null, bool cutmesh=false );

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

  void NewLine( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element, std::vector<Line*> * newlines=NULL );

  bool NewLinesBetween( const std::vector<Point*> & line, Side * cut_side1, Side * cut_side2, Element * cut_element, std::vector<Line*> * newlines=NULL );

  Facet* NewFacet( const std::vector<Point*> & points, Side * side, bool cutsurface );

  VolumeCell* NewVolumeCell( const plain_facet_set & facets,
                             const std::map<std::pair<Point*, Point*>, plain_facet_set> & volume_lines,
                             Element * element );

  Tri3BoundaryCell* NewTri3Cell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points );

  Quad4BoundaryCell* NewQuad4Cell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points );

  Hex8IntegrationCell* NewHex8Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Tet4IntegrationCell* NewTet4Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Tet4IntegrationCell* NewTet4Cell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, VolumeCell * cell );

  Wedge6IntegrationCell* NewWedge6Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  Pyramid5IntegrationCell* NewPyramid5Cell( Point::PointPosition position, const std::vector<Point*> & points, VolumeCell * cell );

  bool DetectSelfCut();

  void Cut( Mesh & mesh, plain_element_set & elements_done, int recursion );

  void Cut( Side & side, const plain_element_set & done, plain_element_set & elements_done, int recursion );

  void Cut( LevelSetSide & side );

  void RectifyCutNumerics();

  void MakeCutLines();

  void MakeFacets();

  void MakeVolumeCells();

  void FindNodePositions();

  void FindLSNodePositions();

  void FindFacetPositions();

  void FindNodalDOFSets( bool include_inner );

  void CreateIntegrationCells( int count, bool levelset );

  void RemoveEmptyVolumeCells();

  void SimplifyIntegrationCells();

  void TestElementVolume( bool fatal );

  void TestElementVolume( DRT::Element::DiscretizationType shape, Element & e, bool fatal );

  void PrintCellStats();

  void Status();

  void PrintFacets();

  void DumpGmsh( std::string name );

  void DumpGmshVolumeCells( std::string name );

  void DumpGmshIntegrationCells( std::string name );

  void NewNodesFromPoints( std::map<Point*, Node*> & nodemap );

  Node* GetNode( int nid ) const;

  Node* GetNode( int nid, const double * xyz, double lsv=0.0 );

#if 0
  Node* GetNode( int nid, Point * p, double lsv=0.0 );
#endif

  Node* GetNode( const plain_int_set & nids, const double * xyz, double lsv=0.0 );

  Edge* GetEdge( Node* begin, Node* end );

  const std::vector<GEO::CUT::Side*> & GetSides( int sid );

  Side* GetSide( const plain_int_set & nids );

  Side* GetSide( int sid,
                 const std::vector<int> & nids,
                 const CellTopologyData * top_data );

#if 0
  Side* GetSide( int sid,
                 const std::vector<Point*> & points,
                 const CellTopologyData * top_data );
#endif

  Side* GetSide( int sid,
                 const std::vector<Node*> & nodes,
                 const CellTopologyData * top_data );

  Element* GetElement( int eid );

  Element* GetElement( int eid,
                       const std::vector<int> & nids,
                       const CellTopologyData & top_data,
                       bool active=true );

#if 0
  Element* GetElement( int eid,
                       const std::vector<Point*> & nids,
                       const CellTopologyData & top_data );
#endif

  Element* GetElement( int eid,
                       const std::vector<Node*> & nids,
                       const CellTopologyData & top_data,
                       bool active=true );

  bool WithinBB( const Epetra_SerialDenseMatrix & xyz );

  bool WithinBB( Element & element );

  Teuchos::RCP<PointPool> Points() { return pp_; }

  const Options & CreateOptions() const { return options_; }

  const std::list<Teuchos::RCP<VolumeCell> > & VolumeCells() const { return cells_; }

  const std::map<plain_int_set, Teuchos::RCP<Edge> > & Edges() const { return edges_; }

  void CreateSideIds();

  void AssignOtherVolumeCells( const Mesh & other );

  void TestVolumeSurface();

  void TestFacetArea();

private:

  void DumpGmsh( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype );

  Edge* GetEdge( const plain_int_set & nids,
                 const std::vector<Node*> & nodes,
                 const CellTopologyData & edge_topology );
  Side* GetSide( int sid,
                 const plain_int_set & nids,
                 const std::vector<Node*> & nodes,
                 const std::vector<Edge*> & edges,
                 const CellTopologyData & side_topology );

  GEO::CUT::Line * NewLineInternal( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element );

  /// setup phase flag
  bool setup_;

  /// options container
  const Options & options_;

  /// mesh dependent point lookup norm
  double norm_;

  /// shared point storage with geometry based access (octtree)
  Teuchos::RCP<PointPool> pp_;

  /// bounding box of this mesh
  BoundingBox bb_;

  /// (output) flag for cut mesh
  bool cutmesh_;

  //! @name Containers that hold all those mesh objects
  /// Plain pointers are used within the library. Memory management is done here.

  std::list<Teuchos::RCP<Line > > lines_;
  std::list<Teuchos::RCP<Facet> > facets_;
  std::list<Teuchos::RCP<VolumeCell> > cells_;
  std::list<Teuchos::RCP<BoundaryCell> > boundarycells_;
  std::list<Teuchos::RCP<IntegrationCell> > integrationcells_;

  /// nodes by unique id
  std::map<int, Teuchos::RCP<Node> > nodes_;

  /// edges by set of nodal ids
  std::map<plain_int_set, Teuchos::RCP<Edge> > edges_;

  /// sides by set of nodal ids
  std::map<plain_int_set, Teuchos::RCP<Side> > sides_;

  /// cut side with (non-unique) side id > -1
  std::map<int, std::vector<Side*> > cut_sides_;

  /// elements by unique id
  std::map<int, Teuchos::RCP<Element> > elements_;

  /// internally generated nodes by nodal ids of element nodes
  std::map<plain_int_set, Node*> shadow_nodes_;

  /// internally generated element (quadratic to linear conversation)
  std::list<Teuchos::RCP<Element> > shadow_elements_;

  //@}
};


}
}

#endif
