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
//#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"



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
class ArbitraryBoundaryCell;

class IntegrationCell;
class Hex8IntegrationCell;
class Tet4IntegrationCell;
class Wedge6IntegrationCell;
class Pyramid5IntegrationCell;



/// Keeps the whole thing. Volume mesh, cut surface, cuts
/*!
 * There is one background mesh and one mesh for the cut surface. These meshes
 * have different objects but share the same points via the point pool.
 *
 * Mesh does the memory management for the whole thing. Therefore, all
 * creation of cut library objects is done via the mesh.
 */
class Mesh
{
public:


//  /*!
//  \brief This class holds data for all volumecells for that dofsets (dofset numbers) of non-row nodes have
//     to be communicated between processors.
//   */
//  class DofSetData
//  {
//  public:
//
//
//    //! constructor for creating volumecells dofset data for communication between processors (initialize missing dofset numbers with -1)
//    DofSetData(
//        std::vector<Point*> cut_points,
//        int peid,
//        const std::vector<int> nodeIdsToCommunicate
//    ) :
//      peid_(peid)
//    {
//      for(std::vector<Point*>::iterator p=cut_points.begin(); p!=cut_points.end(); ++p)
//      {
//        LINALG::Matrix<3,1> coords(true);
//        const double * x = (*p)->X();
//        for(int i=0; i<3; i++)
//        {
//          coords(i) = x[i];
//        }
//        cut_points_coords_.push_back(coords);
//      }
//
//      // initialize the dofset number map with dofsetnumber=-1 for each node
//      for(int i=0; i< (int)nodeIdsToCommunicate.size(); i++)
//      {
//        node_dofsetnumber_map_.insert(pair<int,int>(nodeIdsToCommunicate[i], -1));
//      }
//
//      if(cut_points.size() == 0) dserror("no cut_points for current volumecell in DofSetData!");
//    }
//
//    //! constructor for creating volumecells dofset data for communication between processors (constructor during the Robin round)
//    DofSetData(
//        std::vector<LINALG::Matrix<3,1> > cut_points_coords, // coordinates of cut_points
//        int peid,
//        const std::map<int,int> node_dofsetnumber_map
//    ) :
//      peid_(peid)
//      //node_dofsetnumber_map_(node_dofsetnumber_map)
//    {
//      std::copy(node_dofsetnumber_map.begin(), node_dofsetnumber_map.end(), std::inserter(node_dofsetnumber_map_, node_dofsetnumber_map_.begin()) );
//      std::copy(cut_points_coords.begin(), cut_points_coords.end(), std::inserter(cut_points_coords_, cut_points_coords_.begin()) );
//    }
//
//
//    //! print the dofset data to screen
//    void print()
//    {
//      // print volumecell information
//      std::cout << "Volumecell-DofSetData: " << endl;
//
//      // print parent element Id
//      std::cout << "\tparent element id: " << peid_ << endl;
//
//      // print node_dofsetnumber_map
//      std::cout << "\tnode_dofsetnumber_map" << endl;
//      for(std::map<int,int>::iterator i=node_dofsetnumber_map_.begin(); i!= node_dofsetnumber_map_.end(); ++i)
//      {
//        std::cout << "\t\tnodeId \t" << i->first << "\t dofsetnumber \t" << i->second << endl;
//      }
//    }
//
//    //! destructor
//    ~DofSetData()
//    {
//      return;
//    }
//
//    std::vector<LINALG::Matrix<3,1> > cut_points_coords_;
//    int peid_;                                // parent element Id for volumecell
//    std::map<int,int> node_dofsetnumber_map_; // node Ids and dofset numbers that for that data has to be communicated
//
//
//  private:
//
//  }; // end class DofSetData


  Mesh( Options & options, double norm=1, Teuchos::RCP<PointPool> pp=Teuchos::null, bool cutmesh=false );

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

  ArbitraryBoundaryCell* NewArbitraryCell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points,
      const DRT::UTILS::GaussIntegration& gaussRule, const LINALG::Matrix<3,1>& normal );

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

  /// check if there are nodes whose position is undecided, return whether undecided node positions available
  bool CheckForUndecidedNodePositions( std::map<int,int> & undecided_node );

  void FindNodalDOFSets( bool include_inner );

//  /// fill parallel DofSetData with information that has to be communicated
//  void FillParallelDofSetData(RCP<std::vector<DofSetData> > parallel_dofSetData_);

  void CreateIntegrationCells( int count, bool levelset );

  void MomentFitGaussWeights(bool include_inner, std::string Bcellgausstype);

  void RemoveEmptyVolumeCells();

  void SimplifyIntegrationCells();

  void TestElementVolume( bool fatal );

  void TestElementVolume( DRT::Element::DiscretizationType shape, Element & e, bool fatal );

  void PrintCellStats();

  void Status();

  void PrintFacets();

  void DumpGmsh( std::string name );

  void DumpGmshVolumeCells( std::string name, bool include_inner );

  void DumpGmshIntegrationCells( std::string name );

  void DumpGmshVolumeCells( std::string name );

  void NewNodesFromPoints( std::map<Point*, Node*> & nodemap );

  void GetNodeMap( std::map<int, Node*> & nodemap);

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
  Element* GetSubElement( int parenteid, int subeid,
                       const std::vector<int> & nids,
                       const CellTopologyData & top_data,
                       bool active=true );

#endif

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

  Options & CreateOptions() { return options_; }

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
  Options & options_;

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


//  RCP<std::vector<DofSetData> > parallel_dofSetData_ ;

  //@}
};


}
}

#endif
