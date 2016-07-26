/*!-----------------------------------------------------------------------------------------------*
\file cut_mesh.cpp

\brief class that holds information about a mesh that is cut or about a cutmesh that cuts another mesh

<pre>
\level 3
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include <boost/bind.hpp>

#include "../drt_lib/drt_discret.H"

#include "cut_mesh.H"

#include "cut_point_impl.H"
#include "cut_levelsetside.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_integrationcell.H"
#include "cut_output.H"
#include "cut_parallel.H"

#include "../drt_geometry/element_volume.H"

// search
#include "../drt_geometry/searchtree.H"

#include "../drt_lib/drt_colors.H"

/*-------------------------------------------------------------------------------------*
 * constructor
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Mesh::Mesh( Options & options, double norm, Teuchos::RCP<PointPool> pp, bool cutmesh, int myrank )
  : options_( options ),
    norm_( norm ),
    pp_( pp ),
    cutmesh_( cutmesh ),
    myrank_( myrank )
{
  if ( pp_ == Teuchos::null )
  {
    pp_ = Teuchos::rcp( new PointPool( norm ) ); // create a new octTree based pointpool that stores all points
  }
}


/*-------------------------------------------------------------------------------------*
 * creates a new element, dependent on distype
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Mesh::CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::hex8:
    return CreateHex8( eid, nids );
  case DRT::Element::tet4:
    return CreateTet4( eid, nids );
  case DRT::Element::pyramid5:
    return CreatePyramid5( eid, nids );
  case DRT::Element::wedge6:
    return CreateWedge6( eid, nids );
  default:
    throw std::runtime_error( "unsupported distype" );
  }
  return NULL;
}


/*-------------------------------------------------------------------------------------*
 * creates a new side, dependent on distype
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Mesh::CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::quad4:
    return CreateQuad4( sid, nids );
  case DRT::Element::tri3:
    return CreateTri3( sid, nids );
  default:
    throw std::runtime_error( "unsupported distype" );
  }
  return NULL;
}


/*-------------------------------------------------------------------------------------*
 * creates a new tet4 element based on element id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Mesh::CreateTet4( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}


/*-------------------------------------------------------------------------------------*
 * creates a new pyramid5 element based on element id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Mesh::CreatePyramid5( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Pyramid<5> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}


/*-------------------------------------------------------------------------------------*
 * creates a new wedge6 element based on element id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Mesh::CreateWedge6( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Wedge<6> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}


/*-------------------------------------------------------------------------------------*
 * creates a new hex8 element based on element id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Mesh::CreateHex8( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}


/*-------------------------------------------------------------------------------------*
 * creates a new tri3 side based on side id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Mesh::CreateTri3( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<3> >();
  return GetSide( sid, nids, top_data );
}


/*-------------------------------------------------------------------------------------*
 * creates a new quad4 side based on side id and node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Mesh::CreateQuad4( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  return GetSide( sid, nids, top_data );
}


/*-------------------------------------------------------------------------------------*
 * creates a new point, optional information about cut-edge and cut-side
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::Mesh::NewPoint( const double * x, Edge * cut_edge, Side * cut_side,double tolerance )
{
  bb_.AddPoint( x ); // add the point to the mesh's bounding box
  //Point* p = pp_->NewPoint( x, cut_edge, cut_side, setup_ ? SETUPNODECATCHTOL : MINIMALTOL );
  Point* p = pp_->NewPoint( x, cut_edge, cut_side, tolerance ); // add the point in the point pool
#if 0
  std::cout << "Mesh::NewPoint: ";
  p->Print();
  std::cout << "\n";
#endif
  return p;
}


/*-------------------------------------------------------------------------------------*
 * creates a new line
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::NewLine( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element,
                              std::vector<Line*> * newlines )
{
  if ( p1==p2 )
    throw std::runtime_error( "no line between same point" );

  // If there is a line between those points already, return it. Otherwise
  // create a new one.
  Line * line = p1->CommonLine( p2 );
  if ( line==NULL )
  {
    plain_edge_set edges;
    p1->CommonEdge( p2, edges );
    for ( plain_edge_set::iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      std::vector<Point*> line_points;
      e->CutPointsIncluding( p1, p2, line_points );

      NewLinesBetween( line_points, cut_side1, cut_side2, cut_element, newlines );
    }

    if ( edges.size()==0 )
      NewLineInternal( p1, p2, cut_side1, cut_side2, cut_element );
  }

  //return line;
}


/*-------------------------------------------------------------------------------------------------*
 * Create new line between the two given cut points that are in given two cut sides
 *-------------------------------------------------------------------------------------------------*/
GEO::CUT::Line * GEO::CUT::Mesh::NewLineInternal( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2,
                                                  Element * cut_element )
{
  if ( p1==p2 )
    throw std::runtime_error( "no line between same point" );

  // If there is a line between those points already, return it. Otherwise
  // create a new one.
  Line * line = p1->CommonLine( p2 );
  if ( line==NULL )
  {
    line = new Line( p1, p2, cut_side1, cut_side2, cut_element );
    lines_.push_back( Teuchos::rcp( line ) );
#if 0
    std::cout << "Mesh::NewLineInternal: ";
    p1->Print();
    std::cout << "--";
    p2->Print();
    std::cout << "\n";
#endif
  }
  else // line already exists. just add cut side details to the line
  {
    if ( cut_side1 )
      line->AddSide( cut_side1 );
    if ( cut_side2 )
      line->AddSide( cut_side2 );
    if ( cut_element )
      line->AddElement( cut_element );
  }
  return line;
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::Mesh::NewLinesBetween( const std::vector<Point*> & line, Side * cut_side1,
                                      Side * cut_side2, Element * cut_element, std::vector<Line*> * newlines )
{
  bool hasnewlines = false;
  std::vector<Point*>::const_iterator i = line.begin();
  if ( i!=line.end() )
  {
    Point * bp = *i;
    for ( ++i; i!=line.end(); ++i )
    {
      Point * ep = *i;
      Line * l = NewLineInternal( bp, ep, cut_side1, cut_side2, cut_element );
      if ( newlines!=NULL )
        newlines->push_back( l );
      bp = ep;
    }
    hasnewlines = true;
  }
  return hasnewlines;
}


/*-------------------------------------------------------------------------------------*
 * creates a new facet, consists of points, additional bool if it is a facet on a cutsurface
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Facet* GEO::CUT::Mesh::NewFacet( const std::vector<Point*> & points, Side * side, bool cutsurface )
{
  if ( points.size()==0 )
    throw std::runtime_error( "empty facet" );

  std::vector<Point*>::const_iterator i=points.begin();
  plain_facet_set facets = ( *i )->Facets();
  for ( ++i; i!=points.end(); ++i )
  {
    Point * p = *i;
    p->Intersection( facets );
    if ( facets.size()==0 )
    {
      break;
    }
  }

  for ( plain_facet_set::iterator j=facets.begin(); j!=facets.end(); ++j )
  {
    Facet * f = *j;
    if ( f->Equals( points ) )
    {
      if ( side->IsCutSide() )
      {
        f->ExchangeSide( side, true );
      }
      return f;
    }
  }

  Facet* f = new Facet( *this, points, side, cutsurface );
  facets_.push_back( Teuchos::rcp( f ) );
#if 0
  std::cout << "Mesh::NewFacet: ";
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    p->Print();
    std::cout << " ";
  }
  std::cout << "\n";
#endif
  return f;
}


/*-------------------------------------------------------------------------------------*
 * creates a new volumecell, consists of facets
 *-------------------------------------------------------------------------------------*/
GEO::CUT::VolumeCell* GEO::CUT::Mesh::NewVolumeCell( const plain_facet_set & facets,
                                                     const std::map<std::pair<Point*, Point*>, plain_facet_set > & volume_lines,
                                                     Element * element )
{
  VolumeCell * c = new VolumeCell( facets, volume_lines, element );
  cells_.push_back( Teuchos::rcp( c ) ); // store the pointer in mesh's cells_
  return c;
}


/*-------------------------------------------------------------------------------------*
 * creates a new tri3 boundary cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Tri3BoundaryCell* GEO::CUT::Mesh::NewTri3Cell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points )
{
  if ( points.size()!=3 )
    throw std::runtime_error( "expect 3 points" );
#ifdef DEBUGCUTLIBRARY
  plain_point_set pointtest;
  pointtest.insert( points.begin(), points.end() );
  if ( points.size()!=pointtest.size() )
  {
    throw std::runtime_error( "point used more than once in boundary cell" );
  }
#endif
  Epetra_SerialDenseMatrix xyz( 3, 3 );
  for ( int i=0; i<3; ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Tri3BoundaryCell * bc = new Tri3BoundaryCell( xyz, facet, points );
  boundarycells_.push_back( Teuchos::rcp( bc ) );

  return bc;
}


/*-------------------------------------------------------------------------------------*
 * creates a new quad4 boundary cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Quad4BoundaryCell* GEO::CUT::Mesh::NewQuad4Cell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points )
{
  if ( points.size()!=4 )
    throw std::runtime_error( "expect 4 points" );
#ifdef DEBUGCUTLIBRARY
  plain_point_set pointtest;
  pointtest.insert( points.begin(), points.end() );
  if ( points.size()!=pointtest.size() )
  {
    throw std::runtime_error( "point used more than once in boundary cell" );
  }
#endif
  Epetra_SerialDenseMatrix xyz( 3, 4 );
  for ( int i=0; i<4; ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Quad4BoundaryCell * bc = new Quad4BoundaryCell( xyz, facet, points );
  boundarycells_.push_back( Teuchos::rcp( bc ) );

  return bc;
}


/*-------------------------------------------------------------------------------------*
 * creates a new ??? boundary cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::ArbitraryBoundaryCell* GEO::CUT::Mesh::NewArbitraryCell( VolumeCell * volume, Facet * facet, const std::vector<Point*> & points,
    const DRT::UTILS::GaussIntegration& gaussRule, const LINALG::Matrix<3,1>& normal )
{
#ifdef DEBUGCUTLIBRARY
  plain_point_set pointtest;
  pointtest.insert( points.begin(), points.end() );
  if ( points.size()!=pointtest.size() )
  {
    throw std::runtime_error( "point used more than once in boundary cell" );
  }
#endif
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  ArbitraryBoundaryCell * bc = new ArbitraryBoundaryCell( xyz, facet, points, gaussRule, normal );
  boundarycells_.push_back( Teuchos::rcp( bc ) );

  return bc;
}


/*-------------------------------------------------------------------------------------*
 * creates a new hex8 integration cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Hex8IntegrationCell* GEO::CUT::Mesh::NewHex8Cell( Point::PointPosition position,
                                                            const std::vector<Point*> & points,
                                                            VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Hex8IntegrationCell * c = new Hex8IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );

  return c;
}


/*-------------------------------------------------------------------------------------*
 * creates a new tet4 integration cell, based on points
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Tet4IntegrationCell* GEO::CUT::Mesh::NewTet4Cell( Point::PointPosition position,
                                                            const std::vector<Point*> & points,
                                                            VolumeCell * cell )
{
  if ( points.size()!=4 )
    throw std::runtime_error( "wrong number of cell points" );
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Tet4IntegrationCell * c = new Tet4IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );

  return c;
}


/*-------------------------------------------------------------------------------------*
 * creates a new tet4 integration cell, based on xyz coordinates
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Tet4IntegrationCell* GEO::CUT::Mesh::NewTet4Cell( Point::PointPosition position,
                                                            const Epetra_SerialDenseMatrix & xyz,
                                                            VolumeCell * cell )
{
  std::vector<Point*> points;   // empty list of points
  Tet4IntegrationCell * c = new Tet4IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}


/*-------------------------------------------------------------------------------------*
 * creates a new wedge6 integration cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Wedge6IntegrationCell* GEO::CUT::Mesh::NewWedge6Cell( Point::PointPosition position,
                                                                const std::vector<Point*> & points,
                                                                VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Wedge6IntegrationCell * c = new Wedge6IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );

  return c;
}


/*-------------------------------------------------------------------------------------*
 * creates a new pyramid5 integration cell
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Pyramid5IntegrationCell* GEO::CUT::Mesh::NewPyramid5Cell( Point::PointPosition position,
                                                                    const std::vector<Point*> & points,
                                                                    VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }

  Pyramid5IntegrationCell * c = new Pyramid5IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );

  return c;
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection in the self cut           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::BuildSelfCutTree()
{

  // constructor for the search tree of the background mesh
  selfcuttree_ = Teuchos::rcp(new GEO::SearchTree(5));   // tree_depth_ 4 is reasonable, 5 is possible

  // extent the bounding volume of the root of the search tree to prevent symmetry issues
  LINALG::Matrix<3,2> boundingvolume = bb_.GetBoundingVolume();
//  boundingvolume(0,1) += 1e-4;
//  boundingvolume(1,1) += 1e-4;
//  boundingvolume(2,1) += 1e-4;

  // initializes the search tree of the background mesh
  selfcuttree_->initializeTree(boundingvolume, GEO::TreeType(GEO::OCTTREE));

  // inserts all linear elements into the search tree of the cutter mesh
  for ( std::map<int, Side* >::iterator i=shadow_sides_.begin();
        i!=shadow_sides_.end();
        ++i )
  {
    int sid = i->first;
    Side* s = i->second;
    selfcuttree_->insertElement(sid);
    selfcutbvs_[sid] = s->GetBoundingVolume().GetBoundingVolume();  // better kdop?
  }

  // builds the static search tree of the background mesh
  if (selfcutbvs_.size() != 0) // *********************************************************** possible in case of parallel computing and using only relevant elements
  {
    selfcuttree_->buildStaticSearchTree(selfcutbvs_);
  }

}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::BuildStaticSearchTree()
{

  // constructor for the search tree of the background mesh
  searchtree_ = Teuchos::rcp(new GEO::SearchTree(5));   // tree_depth_ 4 is reasonable, 5 is possible

  // extent the bounding volume of the root of the search tree to prevent symmetry issues
  LINALG::Matrix<3,2> boundingvolume = bb_.GetBoundingVolume();
  boundingvolume(0,1) += 1e-4;
  boundingvolume(1,1) += 1e-4;
  boundingvolume(2,1) += 1e-4;

  // initializes the search tree of the background mesh
  searchtree_->initializeTree(boundingvolume, GEO::TreeType(GEO::OCTTREE));

  // inserts all linear elements into the search tree of the background mesh
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    int eid = i->first;
    Element* e = &*i->second;
    searchtree_->insertElement(eid);
    boundingvolumes_[eid] = e->GetBoundingVolume().GetBoundingVolume();
  }

  // inserts all quadratic elements into the search tree of the background mesh
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i = shadow_elements_.begin();
      i != shadow_elements_.end(); ++i)
  {
    int eid = i->first;
    Element* e = &*i->second;
    searchtree_->insertElement(eid);
    boundingvolumes_[eid] = e->GetBoundingVolume().GetBoundingVolume();
  }

  // builds the static search tree of the background mesh
  if (boundingvolumes_.size() != 0) // *********************************************************** possible in case of parallel computing and using only relevant elements
  {
    searchtree_->buildStaticSearchTree(boundingvolumes_);
  }

}

/*---------------------------------------------------------------------*
 * Cuts the background elements of the mesh with all the cut sides     *
 * Called by Tetmeshintersection (but not for normal meshintersection) *
 *---------------------------------------------------------------------*/
void GEO::CUT::Mesh::Cut( Mesh & mesh, plain_element_set & elements_done, int recursion )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 6/6 --- Cut_Finalize --- CUT (incl. tetmesh-cut)" );

  plain_element_set my_elements_done;

  // perform the cut for each side of the cut_mesh_
  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    mesh.Cut( side, elements_done, my_elements_done, recursion );
  }

  std::copy( my_elements_done.begin(),
             my_elements_done.end(),
             std::inserter( elements_done, elements_done.begin() ) );
}


/*------------------------------------------------------------------------*
 * Cuts the background elements of the mesh with this considered cut side *
 * Called by Tetmeshintersection (but not for normal meshintersection)    *
 *-----------------------------------------------------------------------*/
void GEO::CUT::Mesh::Cut( Side & side, const plain_element_set & done, plain_element_set & elements_done, int recursion )
{
  BoundingBox sidebox( side ); // define a bounding box around the maybe twisted side to determine a preselection of cutting sides and elements
  plain_element_set elements;

#if(0)
  {
    // REMARK: do not use pp_->CollectElements anymore
    // it can happen that some intersections between elements and sides
    // are not detected, because the octtree bounding boxes can lie
    // within a real background element. Such small bounding boxes
    // do not have to contain any points adjacent to elements,
    // then elements have not been found
    //
    // schott 10/2012
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 6/6 --- Cut_Finalize --- preselection of possible cut between" );

    pp_->CollectElements( sidebox, elements ); // find involved elements (octree-based)
  }

#else
  // use a brute force preselection
  //
  // is there an overlap between the side's bounding box
  // and the elements bounding box ?
  // if yes -> try to find a cut
  //
  // REMARK: if the brute force preselection is to slow, one could think about
  // an additional octtree for elements to get a logartithmic complexity
  //

  {

    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 6/6 --- Cut_Finalize --- preselection of possible cut between" );

    // preselection of possible cut between linear elements and the current side
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end(); i++)
    {
      Element* e= &*(i->second);
      BoundingBox elementbox;
      elementbox.Assign(*e);
      if(elementbox.Within(1.0, sidebox))
      {
        if(elements.count(e) == 0)
        {
          elements.insert(e);
        }
      }
    }
    // preselection of possible cut between shadow elements of quadratic elements and the current side
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        i++)
    {
      Element* e= &*(i->second);
      BoundingBox elementbox;
      elementbox.Assign(*e);
      if(elementbox.Within(1.0, sidebox))
      {
        if(elements.count(e) == 0)
        {
          elements.insert(e);
        }
      }
    }

  }
#endif


  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 6/6 --- Cut_Finalize --- cutting sides with elements" );


  // perform the cut of this side for each involved element
  for ( plain_element_set::iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    if ( done.count( e )==0 )
    {
      if ( e->Cut( *this, side, recursion ) )
      {
        elements_done.insert( e );
      }
    }
  }

  }
}


/*-------------------------------------------------------------------------------------*
 * Cuts the background elements with this levelset side
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::Cut( LevelSetSide & side )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.Cut( *this, side, 0 );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.Cut( *this, side, 0 );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * ???
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::RectifyCutNumerics()
{
  for ( std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator i=edges_.begin();
        i!=edges_.end();
        ++i )
  {
    Edge * e = &*i->second;
    e->RectifyCutNumerics();
  }
}

/*------------------------------------------------------------------------------------------------*
 * detects if a side of the cut mesh possibly collides with an element of the background mesh     *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::SearchCollisions(Mesh & cutmesh)
{

  const std::map<plain_int_set, Teuchos::RCP<Side> > & cutsides = cutmesh.Sides();
  for (std::map<plain_int_set, Teuchos::RCP<Side> >::const_iterator i = cutsides.begin();
      i != cutsides.end(); ++i)
  {
    Side* cutside = &*i->second;
    LINALG::Matrix<3, 2> cutsideBV =
        cutside->GetBoundingVolume().GetBoundingVolume();
    if (boundingvolumes_.size() != 0) // ******************************************************* possible in case of parallel computing and using only relevant elements
    {
      std::set<int> collisions;
      // search collision between cutside and elements in the static search tree
      searchtree_->searchCollisions(boundingvolumes_, cutsideBV, 0, collisions);
      // add the cutside to all the elements which have been found
      for (std::set<int>::iterator ic = collisions.begin();
          ic != collisions.end(); ++ic)
      {
        int collisionid = *ic;
        std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find(
            collisionid);
        if (ie != elements_.end())
        {
          Element & e = *ie->second;
          e.AddCutFace(cutside);
        }
        std::map<int, Teuchos::RCP<Element> >::iterator ise =
            shadow_elements_.find(collisionid);
        if (ise != shadow_elements_.end())
        {
          Element & e = *ise->second;
          e.AddCutFace(cutside);
        }
      }
    }
  }

}

/*-------------------------------------------------------------------------------------*
 * finds intersections between sides and edges                                         *
 *                                                                         wirtz 08/14 *
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindCutPoints(int recursion)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 4/6 --- Cut_MeshIntersection --- FindCutPoints" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.FindCutPoints( *this, recursion);
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.FindCutPoints( *this, recursion);
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }

}

/*-------------------------------------------------------------------------------------*
 * create cut lines based on the point cloud
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeCutLines()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 4/6 --- Cut_MeshIntersection --- MakeCutLines" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeCutLines( *this);
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeCutLines( *this);
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * create facets based on the cut lines
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeFacets()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 4/6 --- Cut_MeshIntersection --- MakeFacets" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeFacets( *this );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeFacets( *this );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * create volumecells based on created facets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeVolumeCells()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 4/6 --- Cut_MeshIntersection --- MakeVolumeCells" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeVolumeCells( *this );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MakeVolumeCells( *this );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * find node positions and propagate the positions to facets, points and volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindNodePositions()
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindNodePositions" );


  // On multiple cuts former outside positions can become inside
  // positions. Thus reset all outside positions.
  pp_->ResetOutsidePoints();

  // get nodal positions from elements
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.FindNodePositions();
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.FindNodePositions();
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  // find undecided nodes
  // * for serial simulations all node positions should be set
  // * for parallel simulations there can be some undecided nodes
//  CheckForUndecidedNodePositions();

}


/*-------------------------------------------------------------------------------------*
 * Find node positions and propagate the positions to points.
 * Might need modification for more difficult cut situations with shadow elements.
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindLSNodePositions()
{
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
        i!=nodes_.end();
        ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    if ( p->Position()==Point::undecided )
    {
      // we have to take into account the tolerance for which a node lies on the levelset-side
      double lsv = n->LSV();
      if ( lsv > 0.0) //REFERENCETOL )
      {
        p->Position( Point::outside );
      }
      else if ( lsv < 0.0)//-REFERENCETOL )
      {
        p->Position( Point::inside );
      }
      else //
      {
        //throw std::runtime_error( "undecided nodal point on levelset
        //surface" );
        std::cout << "UNDECIDED CUT POSITION! SHOULD THIS HAPPEN?!" << " lsv= " << std::setprecision(24) << lsv << std::endl;
        p->Position( Point::oncutsurface );
      }
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * find facet positions for remaining facets, points, volumecells that have not been found
 * using FindNodePositions()
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindFacetPositions()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindFacetPositions" );


  plain_volumecell_set undecided;

  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * c = &**i;
//     if ( c->Empty() )
//       continue;
    if ( c->Position()==Point::undecided )
    {
      const plain_facet_set & facets = c->Facets();

  //    bool haveundecided = false;

      Point::PointPosition position = Point::undecided;
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        Point::PointPosition fp = f->Position();
        switch ( fp )
        {
        case Point::undecided:
  //        haveundecided = true;
          break;
        case Point::oncutsurface:
          break;
        case Point::inside:
        case Point::outside:
          if ( position!=Point::undecided and position!=fp )
          {
            throw std::runtime_error( "mixed facet set" );
          }
          position = fp;
          break;
        }
      }

      // set any undecided facets in a decided volume cell

      if ( position != Point::undecided )
      {
        c->Position( position );
      }
      else
      {
        undecided.insert( c );
      }
    }
  }

  // look to neighbouring volume cells to decide the position of this cell

  // this is a slow algorithm, but there should only be very few undecided
  // volume cells, so it should be fine.

  if ( cells_.size()==undecided.size() and cells_.size() > 0)
    throw std::runtime_error( "all volume cells undecided and volume cells available" );

  while ( undecided.size() > 0 )
  {
    unsigned size = undecided.size();
    for ( plain_volumecell_set::iterator ui=undecided.begin(); ui!=undecided.end(); )
    {
      VolumeCell * c = *ui;
      bool done = false;
      const plain_facet_set & facets = c->Facets();
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        {
          VolumeCell * nc = f->Neighbor( c );
          if ( nc!=NULL )
          {
            Point::PointPosition np = nc->Position();
            switch ( np )
            {
            case Point::undecided:
              if ( undecided.count( nc )==0 )
                throw std::runtime_error( "uncomplete set of undecided volume cells" );
              break;
            case Point::oncutsurface:
              throw std::runtime_error( "illegal volume position" );
              break;
            case Point::inside:
            case Point::outside:
              if ( f->OnCutSide() )
                c->Position( np==Point::inside ? Point::outside : Point::inside );
              else
                c->Position( np );
              done = true;
              break;
            }
            if ( done )
              break;
          }
        }
      }
      if ( done )
      {
        set_erase( undecided, ui );
      }
      else
      {
        ++ui;
      }
    }
    if ( size == undecided.size() )
      throw std::runtime_error( "no progress in volume cell position" );
  }

  // second pass

//   for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin(); i!=cells_.end(); ++i )
//   {
//     VolumeCell * c = &**i;
//     const plain_facet_set & facets = c->Facets();

//     bool haveundecided = false;
//     bool havecutsurface = false;
//     GEO::CUT::Point::PointPosition position = GEO::CUT::Point::undecided;
//     for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
//     {
//       Facet * f = *i;
//       GEO::CUT::Point::PointPosition fp = f->Position();
//       switch ( fp )
//       {
//       case GEO::CUT::Point::undecided:
//         haveundecided = true;
//         break;
//       case GEO::CUT::Point::oncutsurface:
//         havecutsurface = true;
//         break;
//       case GEO::CUT::Point::inside:
//       case GEO::CUT::Point::outside:
//         if ( position!=GEO::CUT::Point::undecided and position!=fp )
//         {
//           throw std::runtime_error( "mixed facet set" );
//         }
//         position = fp;
//       }
//     }

    // this is a bold assumption.

//     if ( position == GEO::CUT::Point::undecided )
//     {
//       if ( havecutsurface )
//         position = GEO::CUT::Point::inside;
//       else
//         position = GEO::CUT::Point::outside;
//     }

    // set any undecided facets in a decided volume cell

//     if ( haveundecided and position != GEO::CUT::Point::undecided )
//     {
//       for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
//       {
//         Facet * f = *i;
//         GEO::CUT::Point::PointPosition fp = f->Position();
//         if ( fp==GEO::CUT::Point::undecided )
//         {
//           f->Position( position );
//         }
//       }
//     }
//   }

#if 0
  // If there are any undecided facets left, those should be outside. This can
  // only happen in test cases where there is no edge cut to determine a
  // genuine facet position.
  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin();
        i!=facets_.end();
        ++i )
  {
    Facet * f = &**i;
    if ( f->Position()==Point::undecided )
    {
      f->Position( Point::outside );
    }
  }
#endif
}


/*-------------------------------------------------------------------------------*
 | fill the vector with nids of nodes with undecided position,       schott 03/12 |
 *-------------------------------------------------------------------------------*/
bool GEO::CUT::Mesh::CheckForUndecidedNodePositions( std::map<int, int> &           undecided_nodes,
                                                     std::map<plain_int_set, int> & undecided_shadow_nodes)
{
  // return whether undecided node positions available
  bool undecided_node_positions = false;

  undecided_nodes.clear();

  // find nodes with undecided node positions
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();

    // check for undecided position
    if ( pos==Point::undecided )
    {

      if(n->Id() < 0)
      {
        // continue for shadow nodes, shadow nodes will be handled afterwards
        continue;
      }
      // insert pair of nid and undecided PointPosition
      undecided_nodes.insert(std::pair<int,int>(n->Id(), pos));

      // set undecided_node_positions to true
      undecided_node_positions = true;
    }

  }

  //----------------------------------
  // do the same for possible shadow nodes
  //----------------------------------

  // find nodes with undecided node positions for shadow nodes of e.g. hex20 elements
  for ( std::map<plain_int_set, Node* >::iterator i=shadow_nodes_.begin(); i!=shadow_nodes_.end(); ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();

    // check for undecided position
    if ( pos==Point::undecided )
    {

      if(n->Id() >= 0)
      {
        throw std::runtime_error( "this cannot be a shadow node, a shadow node should have negative nid" );
      }

      // for hex20 elements, boundary shadow nodes are identified by the eight side-nids of the quad8 side
      // the inner shadow node is identified by the 20 nodes of the hex20 element

      // insert pair of sorted node-Ids of side or element to identify the unique shadow node and undecided PointPosition
      undecided_shadow_nodes.insert(std::pair<plain_int_set, int>(i->first, pos));

      // set undecided_node_positions to true
      undecided_node_positions = true;
    }

  }


  return undecided_node_positions;
}


/*-------------------------------------------------------------------------------------*
 * still used???
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindNodalDOFSets( bool include_inner )
{
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
        i!=nodes_.end();
        ++i )
  {
    Node * n = &*i->second;
    n->FindDOFSets( include_inner );
  }

  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * cell = &**i;
    cell->ConnectNodalDOFSets( include_inner );
  }
}


/*-------------------------------------------------------------------------------------*
 * Execute Tessellation with QHULL for each element to generate integrationcells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::CreateIntegrationCells( int count, bool tetcellsonly )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.CreateIntegrationCells( *this, count+1, tetcellsonly );
    }
    catch ( std::runtime_error & err )
    {
      if( count > 0 )
        //dserror("Error occurred in a recursive call i.e. in a call from TetMeshIntersection.");
        std::cout << "Error occurred in a recursive call i.e. in a call from TetMeshIntersection." << std::endl;

      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.CreateIntegrationCells( *this, count+1, tetcellsonly );
    }
    catch ( std::runtime_error & err )
    {
      if( count > 0 )
        std::cout << "Error occurred in a recursive call i.e. in a call from TetMeshIntersection." << std::endl;

      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * Call the moment fitting method for each element to generate the Gaussian integration rule
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MomentFitGaussWeights(bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype)
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MomentFitGaussWeights( *this, include_inner, Bcellgausstype );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.MomentFitGaussWeights( *this, include_inner, Bcellgausstype );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * Call the DirectDivergence method for each element to generate the Gaussian integration rule
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DirectDivergenceGaussRule(bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype)
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.DirectDivergenceGaussRule( *this, include_inner, Bcellgausstype );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.DirectDivergenceGaussRule( *this, include_inner, Bcellgausstype );
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::RemoveEmptyVolumeCells()
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.RemoveEmptyVolumeCells();
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    try
    {
      e.RemoveEmptyVolumeCells();
    }
    catch ( std::runtime_error & err )
    {
      DebugDump(&e,__FILE__,__LINE__);
      throw;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 * Not Implemented, is called, but is an empty function
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::SimplifyIntegrationCells()
{
  // There are obscure bugs. Do not try to be clever.
#if 0
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * vc = &**i;
    vc->SimplifyIntegrationCells( *this );
  }
#endif
}


/*-------------------------------------------------------------------------------------*
 * test if for all elements the element volume is equal to the volume of all integration cells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestElementVolume( bool fatal, INPAR::CUT::VCellGaussPts VCellGP )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    TestElementVolume( e.Shape(), e, fatal, VCellGP );
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    TestElementVolume( e.Shape(), e, fatal, VCellGP );
  }
}


/*-------------------------------------------------------------------------------------*
 * Find the difference between the volume of background element and the sum of volume
 * of all integration cells. There should be no difference between these two
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestElementVolume( DRT::Element::DiscretizationType shape, Element & e, bool fatal, INPAR::CUT::VCellGaussPts VCellGP )
{
  if ( e.IsCut() )
  {
    const std::vector<Node*> & nodes = e.Nodes();
    Epetra_SerialDenseMatrix xyze( 3, nodes.size() );
    double max_norm = 0.0;
    for ( unsigned i=0; i<nodes.size(); ++i )
    {
      nodes[i]->Coordinates( &xyze( 0, i ) );
      LINALG::Matrix<3,1> vec;
      nodes[i]->Coordinates( &vec(0,0) );
      if(vec.Norm2() > max_norm)
        max_norm = vec.Norm2();
    }

    double ev = GEO::ElementVolume( e.Shape(), xyze );

    int numgp = 0;
    int numic = 0;
    int numbc = 0;
    double cv = 0;
    double ba = 0;
    const plain_volumecell_set & cells = e.VolumeCells();
    if(VCellGP == INPAR::CUT::VCellGaussPts_Tessellation)
    {
      for ( plain_volumecell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        VolumeCell * vc = *i;
        numic += vc->IntegrationCells().size();
        numgp += vc->NumGaussPoints( shape );
        cv += vc->Volume();

        const plain_boundarycell_set & bcells = vc->BoundaryCells();
        numbc += bcells.size();
        for ( plain_boundarycell_set::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
        {
          BoundaryCell * bc = *i;
          ba += bc->Area();
        }
      }
    }
    else if(VCellGP == INPAR::CUT::VCellGaussPts_DirectDivergence)
    {
      for ( plain_volumecell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        VolumeCell * vc = *i;
        //numic += vc->IntegrationCells().size();
//        numgp += vc-> //vc->NumGaussPoints( shape );
        cv += vc->Volume();

        const plain_boundarycell_set & bcells = vc->BoundaryCells();
        numbc += bcells.size();
        for ( plain_boundarycell_set::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
        {
          BoundaryCell * bc = *i;
          ba += bc->Area();
        }
      }
    }

#if 0 //to print the number of gauss points in each volumecell
    std::vector<int> numgpeach(cells.size()), hex(cells.size()), tet(cells.size()),
        pyramid(cells.size()), wedge(cells.size());

    for( unsigned i=0;i<cells.size();i++ )
    {
      hex[i] = 0; tet[i] = 0;
      pyramid[i] = 0; wedge[i] = 0;
    }

    for ( plain_volumecell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      VolumeCell * vc = *i;
      numic += vc->IntegrationCells().size();

      for( unsigned j=0;j<vc->IntegrationCells().size();j++ )
      {
        IntegrationCell * icc = vc->IntegrationCells()[j];
        if( icc->Shape() == DRT::Element::hex8 )
          hex[i-cells.begin()]++;
        if( icc->Shape() == DRT::Element::tet4 )
          tet[i-cells.begin()]++;
        if( icc->Shape() == DRT::Element::pyramid5 )
          pyramid[i-cells.begin()]++;
        if( icc->Shape() == DRT::Element::wedge6 )
          wedge[i-cells.begin()]++;
      }

      numgpeach[i-cells.begin()] = vc->NumGaussPoints( shape );
    }
    std::cout<<"Number of Gauss points in each volume = \n";
    for(unsigned i=0;i<cells.size();i++)
    {
      std::cout<<"vc no = "<<i<<"\t"<<"numgp = "<<numgpeach[i]<<"\n";
      std::cout<<"\tnumber of Hex     = "<<hex[i]<<"\n";
      std::cout<<"\tnumber of Tet     = "<<tet[i]<<"\n";
      std::cout<<"\tnumber of Pyramid = "<<pyramid[i]<<"\n";
      std::cout<<"\tnumber of Wedge   = "<<wedge[i]<<"\n";
    }
#endif

    double volume_error = ( ev-cv )/ev;
//Uncomment for a lot of output.
//#ifdef DEBUGCUTLIBRARY
//    std::cout << "#vc=" << cells.size()
//              << " #ic=" << numic
//              << " #bc=" << numbc
//              << " #gp=" << numgp
//              << " \t-- "
//              << ev << "  "
//              << cv << "  "
//              << ev-cv << "  "
//              << volume_error
//              << " \t-- "
//              << ba
//              << "\n";
//#endif


//    if ( fatal and fabs( volume_error ) > 1e-5 )       //Normal test...
    if(fatal and fabs( ev-cv )/(max_norm) > LINSOLVETOL)
    {
      std::stringstream err;
      err << std::setprecision(10) << " !!!!! volume test failed: !!!!! "
          << "eleID=" << e.Id() << "  "
          << "ve=" << ev << "  "
          << "vc=" << cv << "  "
          << "vd= "<< ev-cv << "  "
          << "err=" << volume_error << "  "
          << "tol=" << LINSOLVETOL*max_norm;
      std::cout << err.str() << std::endl;

#ifdef DEBUGCUTLIBRARY
    if(VCellGP == INPAR::CUT::VCellGaussPts_Tessellation)
    {
    std::cout << "#vc=" << cells.size()
              << " #ic=" << numic
              << " #bc=" << numbc
              << " #gp=" << numgp
              << " \t-- "
              << ev << "  "
              << cv << "  "
              << ev-cv << "  "
              << volume_error
              << " \t-- "
              << ba
              << "\n";
    }
#endif

    //Cut test is written for level-set cases as well.
//    if(LINSOLVETOL < fabs( volume_error ) )
    throw std::runtime_error( err.str() );

    }
    else
    {
#if DEBUGCUTLIBRARY
      std::cout << GREEN << "Passed Element Volume Test!" << END_COLOR << std::endl;
#endif
    }
  }
}

/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::PrintCellStats()
{
  const int vectorlength = 21;
  unsigned cutted = 0;
  std::vector<int> numvc( vectorlength, 0 );
  std::map<DRT::Element::DiscretizationType, std::vector<int> > numcells;
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    if ( e.IsCut() )
    {
      cutted += 1;
      const plain_volumecell_set & volumecells = e.VolumeCells();
      numvc[std::min( static_cast<int>( volumecells.size()-1 ), static_cast<int>( numvc.size()-1 ) )] += 1;
      for ( plain_volumecell_set::const_iterator i=volumecells.begin(); i!=volumecells.end(); ++i )
      {
        VolumeCell * vc = *i;
        std::map<DRT::Element::DiscretizationType, int> cell_count;
        const plain_integrationcell_set & cells = vc->IntegrationCells();
        for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          IntegrationCell * cell = *i;
          cell_count[cell->Shape()] += 1;
        }
        const plain_boundarycell_set & bcells = vc->BoundaryCells();
        for ( plain_boundarycell_set::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
        {
          BoundaryCell * bcell = *i;
          cell_count[bcell->Shape()] += 1;
        }
        for ( std::map<DRT::Element::DiscretizationType, int>::iterator i=cell_count.begin(); i!=cell_count.end(); ++i )
        {
          DRT::Element::DiscretizationType shape = i->first;
          int count = i->second;
          std::map<DRT::Element::DiscretizationType, std::vector<int> >::iterator j = numcells.find( shape );
          if ( j==numcells.end() )
          {
            numcells[shape] = std::vector<int>( vectorlength, 0 );
          }
          numcells[shape][std::min( static_cast<int>( numcells[shape].size()-1 ), count-1 )] += 1;
        }
      }
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = *i->second;
    if ( e.IsCut() )
    {
      cutted += 1;
      const plain_volumecell_set & volumecells = e.VolumeCells();
      numvc[std::min( static_cast<int>( volumecells.size()-1 ), static_cast<int>( numvc.size()-1 ) )] += 1;
      for ( plain_volumecell_set::const_iterator i=volumecells.begin(); i!=volumecells.end(); ++i )
      {
        VolumeCell * vc = *i;
        std::map<DRT::Element::DiscretizationType, int> cell_count;
        const plain_integrationcell_set & cells = vc->IntegrationCells();
        for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          IntegrationCell * cell = *i;
          cell_count[cell->Shape()] += 1;
        }
        const plain_boundarycell_set & bcells = vc->BoundaryCells();
        for ( plain_boundarycell_set::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
        {
          BoundaryCell * bcell = *i;
          cell_count[bcell->Shape()] += 1;
        }
        for ( std::map<DRT::Element::DiscretizationType, int>::iterator i=cell_count.begin(); i!=cell_count.end(); ++i )
        {
          DRT::Element::DiscretizationType shape = i->first;
          int count = i->second;
          std::map<DRT::Element::DiscretizationType, std::vector<int> >::iterator j = numcells.find( shape );
          if ( j==numcells.end() )
          {
            numcells[shape] = std::vector<int>( vectorlength, 0 );
          }
          numcells[shape][std::min( static_cast<int>( numcells[shape].size()-1 ), count-1 )] += 1;
        }
      }
    }
  }

  std::cout << "#elements = " << cutted << " of " << ( elements_.size() + shadow_elements_.size() ) << " total\n";
  std::cout << "volume  cells: ";
  //std::streamsize w = std::cout.width( 4 );
  //std::copy( numvc.begin(), numvc.end(), std::ostream_iterator<int>( std::cout, " " ) );
  for ( std::vector<int>::iterator i=numvc.begin(); i!=numvc.end(); ++i )
  {
    int c = *i;
    if ( c != 0 )
    {
      std::cout << std::setw( 4 ) << c << " ";
    }
    else
    {
      std::cout << "     ";
    }
  }
  std::cout << "\n";

  std::cout << "               ";
  for ( int i=1; i<vectorlength; ++i )
  {
    std::cout << std::setw( 4 ) << i << " ";
  }
  std::cout << "   *\n";

  for ( std::map<DRT::Element::DiscretizationType, std::vector<int> >::iterator i=numcells.begin(); i!=numcells.end(); ++i )
  {
    DRT::Element::DiscretizationType shape = i->first;
    std::vector<int> & nc = i->second;
    std::cout << DRT::DistypeToString( shape ) << "\tcells: ";
    //std::copy( nc.begin(), nc.end(), std::ostream_iterator<int>( std::cout, " " ) );
    for ( std::vector<int>::iterator i=nc.begin(); i!=nc.end(); ++i )
    {
      int c = *i;
      if ( c != 0 )
      {
        std::cout << std::setw( 4 ) << c << " ";
      }
      else
      {
        std::cout << "     ";
      }
    }
    std::cout << "\n";
  }
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::Status()
{
#if 0
  std::cout << "#points:    " << pp_->size() << "\n"
            << "#lines:     " << lines_.size() << "\n"
            << "#facets:    " << facets_.size() << "\n"
            << "#nodes:     " << nodes_.size() << "\n"
            << "#edges:     " << edges_.size() << "\n"
            << "#sides:     " << sides_.size() << "\n"
            << "#elements:  " << elements_.size() << "\n"
            << "#snodes:    " << shadow_nodes_.size() << "\n"
            << "#selements: " << shadow_elements_.size() << "\n"
            << "\n";

//   std::cout << "GetElement: ";
//   std::copy( nids.begin(), nids.end(), std::ostream_iterator<int>( std::cout, " " ) );
//   std::cout << "\n";

  std::cout << "edges: {\n";
  for ( std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator i=edges_.begin();
        i!=edges_.end();
        ++i )
  {
    const plain_int_set & s = i->first;
    std::cout << "  ";
    std::copy( s.begin(), s.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
  }
  std::cout << "}\n";
  pp_->Print( std::cout );
#endif

#ifdef DEBUGCUTLIBRARY
  std::string name;

  if ( cutmesh_ )
  {
    name = "cutmesh";
  }
  else
  {
    name = "mesh";
  }

  std::string linefile = name + "_line.plot";
  std::ofstream lf( linefile.c_str() );

  for ( std::list<Teuchos::RCP<Line > >::iterator i=lines_.begin();
        i!=lines_.end();
        ++i )
  {
    ( *i )->Plot( lf );
  }

  std::string edgefile = name + "_edge.plot";
  std::ofstream ef( edgefile.c_str() );

  for ( std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator i=edges_.begin();
        i!=edges_.end();
        ++i )
  {
    i->second->Plot( ef );
  }
#endif
}


/*-------------------------------------------------------------------------------------*
 * print all facets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::PrintFacets()
{
//   for ( std::list<Teuchos::RCP<Point> >::iterator i=points_.begin();
//         i!=points_.end();
//         ++i )
//   {
//     Point & p = **i;
//     p.Print();
//     std::cout << "\n";
//   }

  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin();
        i!=facets_.end();
        ++i )
  {
    Facet & f = **i;
    f.Print();
  }
}


/*-------------------------------------------------------------------------------------*
 * Write full Gmsh Output
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmsh( std::string name )
{
  std::ofstream file( name.c_str() );
 // file.precision(32); //higher precicion!

  //###############write all elements & shadow elements###############
  if ( elements_.size() > 0 || shadow_elements_.size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"Elements");
    for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end(); ++i )
      GEO::CUT::OUTPUT::GmshElementDump( file, &(*i->second) );

    for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin(); i!=shadow_elements_.end();++i )
      GEO::CUT::OUTPUT::GmshElementDump( file, &(*i->second) );
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all sides###############
  if (sides_.size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"Sides");
    for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin(); i!=sides_.end(); ++i )
      GEO::CUT::OUTPUT::GmshSideDump( file, &(*i->second) );
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all nodes###############
  if (sides_.size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"Nodes");
    for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
      GEO::CUT::OUTPUT::GmshNodeDump( file, &(*i->second) );
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all points in pointpool###############
  if (pp_->GetPoints().size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"PoolPoints");
    const RCPPointSet & points = pp_->GetPoints();
    for ( RCPPointSet::const_iterator i= points.begin() ; i != points.end(); ++i )
      GEO::CUT::OUTPUT::GmshPointDump( file, &(*(*i)),(*i)->Id());
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all edges###############
  if (edges_.size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"Edges");
    for ( std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator i=edges_.begin(); i!=edges_.end(); ++i )
      GEO::CUT::OUTPUT::GmshEdgeDump( file, &(*i->second) );
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all lines###############
  if (lines_.size() > 0)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"Lines");
    for ( std::list<Teuchos::RCP<Line > >::iterator i=lines_.begin(); i!=lines_.end(); ++i )
      GEO::CUT::OUTPUT::GmshLineDump( file, &(**i) );
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //###############write all facets (or bacially the facet points)###############
  if (facets_.size() > 0)
  {
//    GEO::CUT::OUTPUT::GmshNewSection(file,"Facet_Points");
//    for ( std::list<Teuchos::RCP<Facet > >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
//      GEO::CUT::OUTPUT::GmshFacetDump(file,&(**i),"points");
//    GEO::CUT::OUTPUT::GmshEndSection(file);

    //###############write all facets (or bacially the facet lines)###############
    GEO::CUT::OUTPUT::GmshNewSection(file,"Facet_Lines");
    for ( std::list<Teuchos::RCP<Facet > >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
      GEO::CUT::OUTPUT::GmshFacetDump(file,&(**i),"lines");
    GEO::CUT::OUTPUT::GmshEndSection(file);

    //###############write all triangulated facets ###############
    GEO::CUT::OUTPUT::GmshNewSection(file,"Facets");
    for ( std::list<Teuchos::RCP<Facet > >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
      GEO::CUT::OUTPUT::GmshFacetDump(file,&(**i),"sides");
    GEO::CUT::OUTPUT::GmshEndSection(file);

    //###############write all cut facets all ###############
    GEO::CUT::OUTPUT::GmshNewSection(file,"cut_Facets");
    for ( std::list<Teuchos::RCP<Facet > >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      if ((*i)->ParentSide()->IsCutSide())
        GEO::CUT::OUTPUT::GmshFacetDump(file,&(**i),"sides",true);
    }
    GEO::CUT::OUTPUT::GmshEndSection(file);

    //###############write all triangulated facets all ###############
    GEO::CUT::OUTPUT::GmshNewSection(file,"ele_Facets");
    for ( std::list<Teuchos::RCP<Facet > >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      if (!(*i)->ParentSide()->IsCutSide())
        GEO::CUT::OUTPUT::GmshFacetDump(file,&(**i),"sides",true);
    }
    GEO::CUT::OUTPUT::GmshEndSection(file);
  }

  //#############write level set information from cut if level set side exists ################
  bool haslevelsetside = false;
  //Does one element have a level set side?
  for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end();i++)
  {
    Element* ele = &*i->second;
    haslevelsetside = ele->HasLevelSetSide();
    if(haslevelsetside)
      break;
  }

  if(haslevelsetside)
  {
    GEO::CUT::OUTPUT::GmshNewSection(file,"LevelSetValues");
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end();i++)
      GEO::CUT::OUTPUT::GmshLevelSetValueDump(file,&(*i->second));
    GEO::CUT::OUTPUT::GmshEndSection(file);

    GEO::CUT::OUTPUT::GmshNewSection(file,"LevelSetGradient");
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end();i++)
      GEO::CUT::OUTPUT::GmshLevelSetGradientDump(file,&(*i->second));
    GEO::CUT::OUTPUT::GmshEndSection(file);

    GEO::CUT::OUTPUT::GmshNewSection(file,"LevelSetOrientation");
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end();i++)
      GEO::CUT::OUTPUT::GmshLevelSetOrientationDump(file,&(*i->second));
    GEO::CUT::OUTPUT::GmshEndSection(file);

#ifdef DEBUGCUTLIBRARY
    file << "View \"LevelSetZeroShape\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin(); i!=elements_.end();i++)
      GEO::CUT::OUTPUT::GmshLevelSetValueZeroSurfaceDump(file,&(*i->second));
    GEO::CUT::OUTPUT::GmshEndSection(file);
#endif
  }
}

/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmshVolumeCells( std::string name, bool include_inner )
{
  std::ofstream file( name.c_str() );
  int count = 0;

  file.setf(std::ios::scientific,std::ios::floatfield);
  file.precision(16);

  file << "View \"VolumeCells\" {\n";
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * vc = &**i;

//    if ( true  ) // std::cout all volumecells - inside and outside
    if ( include_inner or vc->Position()!=Point::inside )
    {
      const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
      for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
            i!=integrationcells.end();
            ++i )
      {
        IntegrationCell * ic = *i;
        ic->DumpGmsh( file, &count );
      }
    }

    count += 1;
  }
  file << "};\n";

  file << "View \"Elements, NumVcs\" {\n";
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element * e = &*i->second;
    const plain_volumecell_set & volumes = e->VolumeCells();
    count = volumes.size();
    for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
    {
      VolumeCell * vc = *i;
      if ( include_inner or vc->Position()!=Point::inside )
      {
        const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
        for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
              i!=integrationcells.end();
              ++i )
        {
          IntegrationCell * ic = *i;
          ic->DumpGmsh( file, &count );
        }
      }
    }
  }
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element * e = &*i->second;
    const plain_volumecell_set & volumes = e->VolumeCells();
    count = volumes.size();
    for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
    {
      VolumeCell * vc = *i;
      if ( include_inner or vc->Position()!=Point::inside )
      {
        const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
        for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
              i!=integrationcells.end();
              ++i )
        {
          IntegrationCell * ic = *i;
          ic->DumpGmsh( file, &count );
        }
      }
    }
  }
  file << "};\n";

  file << "View \"Node-ID\" {\n";
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    const double * x = p->X();
    file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" << n->Id() << "};\n";
  }
  file << "};\n";

#if(0)
  Teuchos::RCP<PointPool> points = Points();

  file << "View \"Points\" {\n";
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * vc = &**i;

//    if ( true  ) // std::cout all volumecells - inside and outside
    if ( include_inner or vc->Position()!=Point::inside )
    {
      PointSet points;
      vc->GetAllPoints(*this, points);

      for ( PointSet::iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        Point * p = *i;
        const double * x = p->X();

        file << "VP(" << x[0] << "," << x[1] << "," << x[2] << "){" << x[0] << "," << x[1] << "," << x[2] << "};\n";
      }
    }

  }
  file << "};\n";
#endif

  file << "View \"Node-Positions\" {\n";
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    const double * x = p->X();
    file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" <<  p->Position() << "};\n";
  }
  file << "};\n";


  //Does there exist a Level Set cut side?
  bool haslevelsetside = false;
  for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
      i!=elements_.end();i++)
  {
    Element* ele = &*i->second;
    haslevelsetside = ele->HasLevelSetSide();

    if(haslevelsetside)
      break;
  }

  if(haslevelsetside)
  {
    file << "View \"LevelSetValues\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();i++)
    {
      Element* ele = &*i->second;
      GEO::CUT::OUTPUT::GmshLevelSetValueDump(file,ele);
    }
    file << "};\n";

    file << "View \"LevelSetGradient\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();i++)
    {
      Element* ele = &*i->second;
      GEO::CUT::OUTPUT::GmshLevelSetGradientDump(file,ele);
    }
    file << "};\n";

    file << "View \"LevelSetOrientation\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();i++)
    {
      Element* ele = &*i->second;
      GEO::CUT::OUTPUT::GmshLevelSetOrientationDump(file,ele);
    }
    file << "};\n";

#ifdef DEBUGCUTLIBRARY
    file << "View \"LevelSetZeroShape\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();i++)
    {
      Element* ele = &*i->second;
      GEO::CUT::OUTPUT::GmshLevelSetValueZeroSurfaceDump(file,ele);
    }
    file << "};\n";
#endif
  }
}


/*-------------------------------------------------------------------------------------*
 * Write all integration cells and boundarycells intro GMSH output file
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmshIntegrationCells( std::string name )
{

  std::ofstream file( name.c_str() );
  file << "View \"IntegrationCells\" {\n";
  for ( std::list<Teuchos::RCP<IntegrationCell> >::iterator i=integrationcells_.begin();
        i!=integrationcells_.end();
        ++i )
  {
    IntegrationCell * ic = &**i;
    ic->DumpGmsh( file );
  }
  file << "};\n";

  //Write BCs for outside VolumeCell:
  file << "View \"BoundaryCells\" {\n";
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
      i!=cells_.end();
      ++i )
  {
    VolumeCell * volcell = &**i;
    if(volcell->Position()==GEO::CUT::Point::outside)
    {
      plain_boundarycell_set bc_cells = volcell->BoundaryCells();
      for ( plain_boundarycell_set::iterator j=bc_cells.begin();
          j!=bc_cells.end();
          ++j )
      {
        BoundaryCell * bc = *j;
        if(bc->IsValid())
          bc->DumpGmsh(file);
      }
    }
  }
  file << "};\n";

  //Write BCs for inside VolumeCell:
  file << "View \"BoundaryCellsNormal\" {\n";
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
      i!=cells_.end();
      ++i )
  {
    VolumeCell * volcell = &**i;
    if(volcell->Position()==GEO::CUT::Point::outside)
    {
      plain_boundarycell_set bc_cells = volcell->BoundaryCells();
      for ( plain_boundarycell_set::iterator j=bc_cells.begin();
            j!=bc_cells.end();
            ++j )
      {
        BoundaryCell * bc = *j;
        if(bc->IsValid())
          bc->DumpGmshNormal(file);
      }
    }
  }


  file << "};\n";

}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmshVolumeCells( std::string name )
{
  std::ofstream file( name.c_str() );
  file << "View \"VolumeCells\" {\n";
  for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
      i!=elements_.end();i++)
  {
    Element& ele = *i->second;
    const plain_volumecell_set cells = ele.VolumeCells();
    for(plain_volumecell_set::const_iterator j=cells.begin();j!=cells.end();j++)
    {
      VolumeCell *vcc = *j;
      vcc->DumpGmsh(file);

    }
#if 0
    if(cells.size()==1)
    {
      const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
      for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
                    i!=integrationcells.end();
                    ++i )
              {
                IntegrationCell * ic = *i;
                ic->DumpGmsh( file, &count );
              }

    }
/*    else
    {

    }*/
#endif
  }
  file << "};\n";

  file << "View \"BoundaryCells\" {\n";
  bool haslevelsetside = false;
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=boundarycells_.begin();
        i!=boundarycells_.end();
        ++i )
  {
    BoundaryCell * bc = &**i;

    if(bc->GetFacet()->BelongsToLevelSetSide())
      haslevelsetside = true;

     if ( bc->IsValid() )
      bc->DumpGmsh( file );
  }
  file << "};\n";

  if(haslevelsetside)
  {
    file << "View \"LevelSetInfoOnFacet\" {\n";
    for(std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();i++)
    {
      Element* ele = &*i->second;
      const plain_facet_set facets = ele->Facets();
      for(plain_facet_set::const_iterator j=facets.begin();j!=facets.end();j++)
      {
        Facet *facet = *j;

        if(facet->OnCutSide())
        {
          LINALG::Matrix<3,1> facet_triang_midpoint_coord(true);

          if(facet->IsTriangulated())
          {
            std::vector<std::vector<Point*> > facet_triang = facet->Triangulation();
            Point* facet_triang_midpoint = (facet_triang[0])[0];

            //Choose midpoint of facet?
            facet_triang_midpoint->Coordinates(&facet_triang_midpoint_coord(0,0));
          }
          else
          {
            LINALG::Matrix<3,1> cur;
            std::vector<Point*> pts =facet->Points();
            for( std::vector<Point*>::iterator i=pts.begin();i!=pts.end();i++ )
            {
              Point* p1 = *i;
              p1->Coordinates(cur.A());
              facet_triang_midpoint_coord.Update(1.0,cur,1.0);
            }
            facet_triang_midpoint_coord.Scale(1.0/pts.size());
          }

          file << "VT(";
          file << facet_triang_midpoint_coord( 0, 0 ) << ","
              << facet_triang_midpoint_coord( 1, 0 ) << ","
              << facet_triang_midpoint_coord( 2, 0 );

          file << "){";

          std::vector<double> normal = ele->GetLevelSetGradient(facet_triang_midpoint_coord);//facet->GetLevelSetFacetNormal(ele);
          file << normal[0] << "," << normal[1] << "," << normal[2];

          file << "};\n";
        }
      }
    }
  }
  file << "};\n";

}

/*-------------------------------------------------------------------------------------*
 * DebugDump to call before runtime error                                     ager 04/15
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DebugDump(GEO::CUT::Element* ele, std::string file, int line)
{
  if (file != "" || line != -1)
    std::cout << "GEO::CUT::Mesh::DebugDump called in " << file << " in line " << line << std::endl;


  std::stringstream str;
  str << ".full_debug_cut_output." << myrank_ << "_CUTFAIL.pos";
  std::string filename(GEO::CUT::OUTPUT::GenerateGmshOutputFilename(str.str()));
  DumpGmsh(filename);

  ele->DebugDump();
}

/*-------------------------------------------------------------------------------------*
 * ? -> used in cut_tetmeshintersection
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::NewNodesFromPoints( std::map<Point*, Node*> & nodemap )
{
  int nid = nodes_.size();
  for ( std::map<Point*, Node*>::iterator i=nodemap.begin(); i!=nodemap.end(); ++i )
  {
    Point * p = i->first;
    while ( nodes_.count( nid ) > 0 )
      nid += 1;
    Node * n = GetNode( nid, p->X() );
    i->second = n;
    nid += 1;
  }
}


/*-------------------------------------------------------------------------------------*
 * get a map of node id and the pointer to the node
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::GetNodeMap( std::map<int, Node *> & nodemap )
{

    for (std::map<int, Teuchos::RCP<Node> >::iterator i = nodes_.begin();
         i != nodes_.end();
         i++ )
    {
        int nid = i->first;

        nodemap.insert(std::pair< int, Node* >(nid, &*(i->second)) );
    }

}


/*-------------------------------------------------------------------------------------*
    Returns the node with given id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid ) const
{
  std::map<int, Teuchos::RCP<Node> >::const_iterator i = nodes_.find( nid );
  if ( i != nodes_.end() )
  {
    return &*i->second;
  }
  return NULL;
}

/*-------------------------------------------------------------------------------------*
    Returns the unique shadow node
    identified by given nids of a quad8 boundary side or all nodes of hex20 element
    for the inner shadow node
    Remark: Do not identify a shadow node via its Id, since it is not unique over processors
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( const plain_int_set & nids ) const
{
  std::map<plain_int_set, Node*>::const_iterator i=shadow_nodes_.find( nids );
  if ( i != shadow_nodes_.end() )
  {
    return &*i->second;
  }
  return NULL;
}


/*-------------------------------------------------------------------------------------*
    If node with the given id exists return the node, else create a new node with
    given coordinates and levelset value
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid, const double * xyz, double lsv, double tolerance )
{
  std::map<int, Teuchos::RCP<Node> >::iterator i = nodes_.find( nid );
  if ( i != nodes_.end() )
  {
    return &*i->second;
  }
  if ( xyz==NULL )
    throw std::runtime_error( "cannot create node without coordinates" );

//   Point * p = pp_->GetPoint( xyz, NULL, NULL, MINIMALTOL );
//   if ( p!=NULL )
//   {
//     // We already have a point at this location. See if there is a node to
//     // it. If so, that's the one.
//     for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
//           i!=nodes_.end();
//           ++i )
//     {
//       Teuchos::RCP<Node> n = i->second;
//       if ( n->point()==p )
//       {
//         nodes_[nid] = n;
//         return NULL;
//       }
//     }
//   }
  Point * p = NewPoint( xyz, NULL, NULL, tolerance );
  Node * n = new Node( nid, p, lsv );
  nodes_[nid] = Teuchos::rcp( n );
#ifdef DRT_CUT_DUMPCREATION
  if ( cutmesh_ )
    std::cout << "GetCutNode( " << nid << ", ";
  else
    std::cout << "GetNode( " << nid << ", ";
  std::copy( xyz, xyz+3, std::ostream_iterator<double>( std::cout, ", " ) );
  std::cout << lsv << " );\n";
#endif
  return n;
}

#if 0
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid, Point * p, double lsv )
{
  std::map<int, Teuchos::RCP<Node> >::iterator i = nodes_.find( nid );
  if ( i != nodes_.end() )
  {
    Node * n = &*i->second;
    if ( n->point()!=p )
    {
      throw std::runtime_error( "node id with different point exists" );
    }
    return n;
  }
  Node * n = new Node( nid, p, lsv );
  nodes_[nid] = Teuchos::rcp( n );
  return n;
}
#endif


/*-------------------------------------------------------------------------------------*
 * GetNode routine for shadow nodes (e.g. center nodes of hex20 element),
 * find the unique shadow node (if existent) using the eight quad8 nodes of a quad8 side as key
 * in case of a boundary center node, and use the 20 nodes of the hex20 element to identify the
 * inner shadow/center node of the hex20 element
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( const plain_int_set & nids, const double * xyz, double lsv )
{
  std::map<plain_int_set, Node*>::iterator i=shadow_nodes_.find( nids );
  if ( i!=shadow_nodes_.end() )
  {
    return &*i->second;
  }
  int nid = - shadow_nodes_.size() - 1;
  if ( nodes_.find( nid )!=nodes_.end() )
  {
    throw std::runtime_error( "shadow node already exists" );
  }

  // Remark: the nid of a shadow node is not unique over processors, consequently numbered with negative integers on each proc
  Node * n = GetNode( nid, xyz, lsv );
  shadow_nodes_[nids] = n;
  return n;
}


/*-------------------------------------------------------------------------------------*
 * get the edge with begin node and end node
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Edge* GEO::CUT::Mesh::GetEdge( Node* begin, Node* end )
{
  if ( begin->point()==end->point() )
    throw std::runtime_error( "edge between same point" );

  plain_int_set nids;
  nids.insert( begin->Id() );
  nids.insert( end->Id() );

  std::vector<Node*> nodes( 2 );
  nodes[0] = begin;
  nodes[1] = end;

  return GetEdge( nids, nodes, *shards::getCellTopologyData< shards::Line<2> >() );
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Edge* GEO::CUT::Mesh::GetEdge( const plain_int_set & nids, const std::vector<Node*> & nodes, const CellTopologyData & edge_topology )
{
  std::map<plain_int_set, Teuchos::RCP<Edge> >::iterator i = edges_.find( nids );
  if ( i != edges_.end() )
  {
    return &*i->second;
  }

  // Search for edges between those points from an other mesh. This might
  // happen if mesh and cut mesh use the same nodes and thus the same
  // edges. This will happen with double cuts.

  Point * p1 = nodes[0]->point();
  Point * p2 = nodes[1]->point();

  plain_edge_set edges;
  p1->CommonEdge( p2, edges );
  for ( plain_edge_set::iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    Point * ep1 = e->BeginNode()->point();
    Point * ep2 = e->EndNode()  ->point();
    if ( p1==ep1 and p2==ep2 )
    {
      edges_[nids] = Teuchos::rcp( e, false );
      return e;
    }
    if ( p1==ep2 and p2==ep1 )
    {
      edges_[nids] = Teuchos::rcp( e, false );
      return e;
    }
  }

  Edge * e = NULL;
  switch ( edge_topology.key )
  {
  case shards::Line<2>::key :
//     if ( nodes[0]->point()==nodes[1]->point() )
//       throw std::runtime_error( "edge between same point" );
    e = new Edge( nodes );
    break;
  default:
    throw std::runtime_error( "unsupported edge topology" );
  }
  edges_[nids] = Teuchos::rcp( e );
  return e;
}

/*-------------------------------------------------------------------------------------*
 * get the side that contains the nodes with the following node ids
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Mesh::GetSide( std::vector<int>& nids ) const
{
  // create a sorted vector
  plain_int_set node_ids;

  for(unsigned int i=0; i< nids.size(); i++)
  {
    node_ids.insert(nids[i]);
  }

  std::map<plain_int_set, Teuchos::RCP<Side> >::const_iterator i = sides_.find( node_ids );
  if ( i != sides_.end() )
  {
    return &*i->second;
  }
  return NULL;
}


/*-------------------------------------------------------------------------------------*
 * ???
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid, const std::vector<int> & nids, const CellTopologyData * top_data )
{
  unsigned nc = top_data->node_count;
  std::vector<Node*> nodes;
  nodes.reserve( nc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( nids[i], static_cast<double*>( NULL ) ) );
  }

  return GetSide( sid, nodes, top_data );
}

#if 0
GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid,
                                         const std::vector<Point*> & points,
                                         const CellTopologyData * top_data )
{
  unsigned nc = top_data->node_count;
  std::vector<Node*> nodes;
  nodes.reserve( nc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( points[i], static_cast<double*>( NULL ) ) );
  }

  return GetSide( sid, nodes, top_data );
}
#endif


/*-------------------------------------------------------------------------------------*
 * ???
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid,
                                         const std::vector<Node*> & nodes,
                                         const CellTopologyData * top_data )
{
  unsigned ec = top_data->edge_count;
  std::vector<Edge*> edges;
  edges.reserve( ec );

  plain_int_set nidset;

  for ( unsigned i=0; i<ec; ++i )
  {
    const CellTopologyData_Subcell & edge = top_data->edge[i] ;
    const CellTopologyData & edge_topology = *edge.topology;

    std::vector<Node*> edge_nodes;
    plain_int_set edge_nids;
    edge_nodes.reserve( edge_topology.node_count );
    for ( unsigned j=0; j<edge_topology.node_count; ++j )
    {
      edge_nids.insert( nodes[edge.node[j]]->Id() );
      edge_nodes.push_back( nodes[edge.node[j]] );
    }
    edges.push_back( GetEdge( edge_nids, edge_nodes, edge_topology ) );

    std::copy( edge_nids.begin(), edge_nids.end(),
               std::inserter( nidset, nidset.begin() ) );
  }

  Side * s = GetSide( sid, nidset, nodes, edges, *top_data );

  return s;
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid,
                                         const plain_int_set & nids,
                                         const std::vector<Node*> & nodes,
                                         const std::vector<Edge*> & edges,
                                         const CellTopologyData & side_topology )
{
  std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i = sides_.find( nids );
  if ( i != sides_.end() )
  {
    return &*i->second;
  }
  //Facet * f = new Facet;
  //facets_.push_back( Teuchos::rcp( f ) );
  Side * s = NULL;
  switch ( side_topology.key )
  {
  case shards::Triangle<3>::key :
    s = new ConcreteSide<DRT::Element::tri3>( sid, nodes, edges );
    break;
  case shards::Quadrilateral<4>::key :
    s = new ConcreteSide<DRT::Element::quad4>( sid, nodes, edges );
    break;
  default:
    throw std::runtime_error( "unsupported side topology" );
  }
  sides_[nids] = Teuchos::rcp( s );
  int seid = - shadow_sides_.size() - 1;
  // Remark: the seid of a shadow node is not unique over processors, consequently numbered with negative integers on each proc
  shadow_sides_[seid] = s;
  return s;
}


/*-------------------------------------------------------------------------------------*
 * Returns the element with given id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid )
{
  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( eid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }
  return NULL;
}


/*-------------------------------------------------------------------------------------*
 * If element with the given id exists return the element, else create a new element
 * with given node ids. All details of the element are in cell topology data
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid,
                                               const std::vector<int> & nids,
                                               const CellTopologyData & top_data,
                                               bool active )
{
  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( eid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }

  unsigned nc = top_data.node_count;

  /*if( nc != nids.size() )
    dserror("node details does not match with topology data");*/

  std::vector<Node*> nodes;
  nodes.reserve( nc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( nids[i], static_cast<double*>( NULL ) ) );
  }

  return GetElement( eid, nodes, top_data, active );
}

#if 0
GEO::CUT::Element* GEO::CUT::Mesh::GetSubElement( int parenteid, int subeid,
                                               const std::vector<int> & nids,
                                               const CellTopologyData & top_data,
                                               bool active )
{
  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( subeid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }

  unsigned nc = top_data.node_count;
  std::vector<Node*> nodes;
  nodes.reserve( nc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( nids[i], static_cast<double*>( NULL ) ) );
  }

  return GetSubElement( parenteid, subeid, nodes, top_data, active );
}
#endif

#if 0
GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid,
                                               const std::vector<Point*> & points,
                                               const CellTopologyData & top_data )
{
  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( eid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }

  unsigned nc = top_data.node_count;
  std::vector<Node*> nodes;
  nodes.reserve( nc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( points[i], static_cast<double*>( NULL ) ) );
  }

  return GetElement( eid, nodes, top_data );
}
#endif


/*-------------------------------------------------------------------------------------*
    Create a new element with given nodes. All details of the element are in
    cell topology data
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid,
                                               const std::vector<Node*> & nodes,
                                               const CellTopologyData & top_data,
                                               bool active )
{
  static int shards_to_baci[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 26, 20, 25, 24, 22, 21, 23, -1
  };

  unsigned sc = top_data.side_count;

  std::vector<Side*> sides;

  sides.reserve( sc );


  for ( unsigned i=0; i<sc; ++i )
  {
    const CellTopologyData_Subcell & side = top_data.side[i] ;
    const CellTopologyData & side_topology = *side.topology;

    std::vector<Node*> side_nodes;
    std::vector<int> side_nids;
    side_nodes.reserve( side_topology.node_count );
    side_nids .reserve( side_topology.node_count );
    for ( unsigned j=0; j<side_topology.node_count; ++j )
    {
      int nid = shards_to_baci[side.node[j]];
      side_nids .push_back( nodes[nid]->Id() );
      side_nodes.push_back( nodes[nid] );
    }

    std::vector<Edge*> side_edges;
    side_edges.reserve( side_topology.edge_count );
    for ( unsigned j=0; j<side_topology.edge_count; ++j )
    {
      const CellTopologyData_Subcell & edge = side_topology.edge[j] ;
      const CellTopologyData & edge_topology = *edge.topology;

      std::vector<Node*> edge_nodes;
      plain_int_set edge_nids;
      edge_nodes.reserve( edge_topology.node_count );
      for ( unsigned j=0; j<edge_topology.node_count; ++j )
      {
        edge_nids.insert( side_nids[edge.node[j]] );
        edge_nodes.push_back( side_nodes[edge.node[j]] );
      }

      side_edges.push_back( GetEdge( edge_nids, edge_nodes, edge_topology ) );
    }

    plain_int_set side_nidset;
    std::copy( side_nids.begin(), side_nids.end(),
               std::inserter( side_nidset, side_nidset.begin() ) );
    sides.push_back( GetSide( -1, side_nidset, side_nodes, side_edges, side_topology ) );
  }

  Element * e = NULL;
  switch ( top_data.key )
  {
  case shards::Tetrahedron<4>::key :
    e = new ConcreteElement<DRT::Element::tet4>( eid, sides, nodes, active );
    break;

  case shards::Hexahedron<8>::key :
    e = new ConcreteElement<DRT::Element::hex8>( eid, sides, nodes, active );
    break;

  case shards::Pyramid<5>::key :
    e = new ConcreteElement<DRT::Element::pyramid5>( eid, sides, nodes, active );
    break;

  case shards::Wedge<6>::key :
    e = new ConcreteElement<DRT::Element::wedge6>( eid, sides, nodes, active );
    break;

  default:
    throw std::runtime_error( "unsupported element topology" );
  }
  if ( eid > -1 )
  {
    elements_[eid] = Teuchos::rcp( e );
  }
  else
  {
    int seid = - shadow_elements_.size() - 1;
    // Remark: the seid of a shadow node is not unique over processors, consequently numbered with negative integers on each proc
    shadow_elements_[seid] = Teuchos::rcp( e );
  }
  return e;
}
#if 0

GEO::CUT::Element* GEO::CUT::Mesh::GetSubElement( int parenteid,
                                                  int subeid,
                                               const std::vector<Node*> & nodes,
                                               const CellTopologyData & top_data,
                                               bool active )
{
  static int shards_to_baci[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 26, 20, 25, 24, 22, 21, 23, -1
  };

  unsigned sc = top_data.side_count;

  std::vector<Side*> sides;

  sides.reserve( sc );


  for ( unsigned i=0; i<sc; ++i )
  {
    const CellTopologyData_Subcell & side = top_data.side[i] ;
    const CellTopologyData & side_topology = *side.topology;

    std::vector<Node*> side_nodes;
    std::vector<int> side_nids;
    side_nodes.reserve( side_topology.node_count );
    side_nids .reserve( side_topology.node_count );
    for ( unsigned j=0; j<side_topology.node_count; ++j )
    {
      int nid = shards_to_baci[side.node[j]];
      side_nids .push_back( nodes[nid]->Id() );
      side_nodes.push_back( nodes[nid] );
    }

    std::vector<Edge*> side_edges;
    side_edges.reserve( side_topology.edge_count );
    for ( unsigned j=0; j<side_topology.edge_count; ++j )
    {
      const CellTopologyData_Subcell & edge = side_topology.edge[j] ;
      const CellTopologyData & edge_topology = *edge.topology;

      std::vector<Node*> edge_nodes;
      plain_int_set edge_nids;
      edge_nodes.reserve( edge_topology.node_count );
      for ( unsigned j=0; j<edge_topology.node_count; ++j )
      {
        edge_nids.insert( side_nids[edge.node[j]] );
        edge_nodes.push_back( side_nodes[edge.node[j]] );
      }

      side_edges.push_back( GetEdge( edge_nids, edge_nodes, edge_topology ) );
    }

    plain_int_set side_nidset;
    std::copy( side_nids.begin(), side_nids.end(),
               std::inserter( side_nidset, side_nidset.begin() ) );
    sides.push_back( GetSide( -1, side_nidset, side_nodes, side_edges, side_topology ) );
  }

  Element * e = NULL;
  switch ( top_data.key )
  {
  case shards::Tetrahedron<4>::key :
    e = new ConcreteElement<DRT::Element::tet4>( subeid, sides, nodes, active );
    break;

  case shards::Hexahedron<8>::key :
    e = new ConcreteElement<DRT::Element::hex8>( subeid, sides, nodes, active );
    break;

  case shards::Pyramid<5>::key :
    e = new ConcreteElement<DRT::Element::pyramid5>( subeid, sides, nodes, active );
    break;

  case shards::Wedge<6>::key :
    e = new ConcreteElement<DRT::Element::wedge6>( subeid, sides, nodes, active );
    break;

  default:
    throw std::runtime_error( "unsupported element topology" );
  }
  if ( subeid > -1 )
  {
   dserror(" this case should not be called in sub-function");
    elements_[parenteid] = Teuchos::rcp( e );
  }
  else
  {
    shadow_elements_[parenteid].push_back( Teuchos::rcp( e ) );
//    shadow_elements_.push_back( Teuchos::rcp( e ) );
  }
  return e;
}
#endif



/*-------------------------------------------------------------------------------------*
 * check if xyz-coordinates lie within the mesh's bounding box
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::Mesh::WithinBB( const Epetra_SerialDenseMatrix & xyz )
{
  return bb_.Within( norm_, xyz );
}


/*-------------------------------------------------------------------------------------*
 * check if the element lies within the bounding box
 *-------------------------------------------------------------------------------------*/
bool GEO::CUT::Mesh::WithinBB( Element & element )
{
  return bb_.Within( norm_, element );
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::CreateSideIds(int lastid)
{
  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side * s = &*i->second;
    if ( s->Id() < 0 )
    {
      lastid += 1;
      s->SetId( lastid );
    }
  }

  cutmesh_ = true;
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::AssignOtherVolumeCells( const Mesh & other )
{
  const std::list<Teuchos::RCP<VolumeCell> > & other_cells = other.VolumeCells();
  plain_volumecell_set cells;
  for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=other_cells.begin();
        i!=other_cells.end();
        ++i )
  {
    VolumeCell * vc = &**i;

    const plain_facet_set & facets = vc->Facets();
    plain_facet_set cut_facets;
    std::remove_copy_if( facets.begin(), facets.end(),
                         std::inserter( cut_facets, cut_facets.begin() ),
                         not boost::bind( &Facet::OnCutSide, _1 ) );

    plain_side_set cut_sides;
    std::transform( cut_facets.begin(), cut_facets.end(),
                    std::inserter( cut_sides, cut_sides.begin() ),
                    boost::bind( &Facet::ParentSide, _1 ) );

    if ( cut_sides.size() > 1 )
    {
      plain_side_set::iterator i = cut_sides.begin();
      Side * s1 = *i;
      ++i;
      Side * s2 = *i;
      Element * e = s1->CommonElement( s2 );
      if ( e==NULL )
        throw std::runtime_error( "no common element on cut sides" );
      e->AssignOtherVolumeCell( vc );
    }
    else
    {
      cells.insert( vc );
    }
  }

  if ( cells.size() > 0 )
  {
    std::stringstream str;
    str << cells.size() << " volume cells left. Need to handle those.";
    throw std::runtime_error( str.str() );
  }
}


/*-------------------------------------------------------------------------------------*
 * Broken test needs fixing before being run properly.
 * The idea is to test if the boundary of the boundary cells (i.e. lines) created include
 * the cut facet given.
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestVolumeSurface()
{
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = &**i;
    vc->TestSurface();
  }
}


/*-------------------------------------------------------------------------------------*
 * Tests if the area of the facet created by boundary cells is equal on both sides.
 *
 * WARNING: If set sharply (i.e. 10^-11) this will cause some test cases to fail for Combust,
 *          cut is done in local coordinates. Thus this makes it likely that the cut
 *          is sensitive to 10^-11? Why?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestFacetArea(bool istetmeshintersection)
{
  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = &**i;

    // This is a crude test. We do not demand so much here...
    double tolerance = 1e-7;
    //double tolerance = 1e-11; // sharper tolerance
    f->TestFacetArea( tolerance, istetmeshintersection );
  }
}
