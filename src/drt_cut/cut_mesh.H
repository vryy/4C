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

#include "cut_tetgen.H"
#include "cut_boundingbox.H"

namespace GEO
{
namespace CUT
{

class Node;
class Edge;
class Side;
class LinearSide;
class LevelSetSide;
class Element;

class PointPool;

class Point;
class Line;
class Facet;

/// cell generator callback
class CellGenerator
{
public:

  virtual ~CellGenerator() {}

  virtual void NonCut( Element* element ) {}

#ifdef QHULL
  virtual void Generate( Element* element, const tetgenio & out ) = 0;
#endif
};

/// Keeps the whole thing. Volume mesh, cut surface, cuts
class Mesh
{
public:

  Mesh( double norm=1, Teuchos::RCP<PointPool> pp=Teuchos::null, bool cutmesh=false );

#ifdef QHULL

  void AddTetgen( const tetgenio & out );

  void ExtractTetgen( tetgenio & out );

#endif

  void FillComplete();

  Element * CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  Side * CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  Element * CreateTet4( int eid, const std::vector<int> & nids );

  Element * CreateTet10( int eid, const std::vector<int> & nids );

  Element * CreatePyramid5( int eid, const std::vector<int> & nids );

  Element * CreateWedge6( int eid, const std::vector<int> & nids );

  Element * CreateHex8( int eid, const std::vector<int> & nids );

  Element * CreateHex20( int eid, const std::vector<int> & nids );

  Element * CreateHex27( int eid, const std::vector<int> & nids );

  Side * CreateTri3( int sid, const std::vector<int> & nids );

  Side * CreateTri6( int sid, const std::vector<int> & nids );

  Side * CreateQuad4( int sid, const std::vector<int> & nids );

  Side * CreateQuad8( int sid, const std::vector<int> & nids );

  Side * CreateQuad9( int sid, const std::vector<int> & nids );

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side );

  Line* NewLine( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element );

  Facet* NewFacet( const std::vector<Point*> & points, Side * side, bool cutsurface );

  void SelfCut();

  void Cut( Mesh & mesh, std::set<Element*> & elements_done );

  void Cut( LinearSide & side, const std::set<Element*> & done, std::set<Element*> & elements_done );

  void Cut( LevelSetSide & side );

  void MakeFacets();

  void FindNodePositions();

  void FindLSNodePositions();

  void Status();

  void PrintFacets();

  void DumpGmsh( std::string name );

  void GenerateTetgen( CellGenerator * generator=NULL );

  Node* GetNode( int nid, const double * xyz, double lsv=0.0 );

  Node* GetNode( const std::set<int> & nids, const double * xyz );

  Edge* GetEdge( Node* begin, Node* end );

  const std::vector<GEO::CUT::Side*> & GetSides( int sid );

  Side* GetSide( const std::set<int> & nids );

  Side* GetSide( int sid,
                 const std::vector<int> & nids,
                 const CellTopologyData * top_data );

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

  double norm_;
  Teuchos::RCP<PointPool> pp_;
  BoundingBox bb_;
  bool cutmesh_;

  //std::list<Teuchos::RCP<Point> > points_;
  std::list<Teuchos::RCP<Line > > lines_;
  std::list<Teuchos::RCP<Facet> > facets_;

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
