/*!-----------------------------------------------------------------------------------------------*
\file cut_mesh.cpp

\brief class that holds information about a mesh that is cut or about a cutmesh that cuts another mesh

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include <boost/bind.hpp>

#include "../drt_lib/drt_discret.H"

#include "cut_mesh.H"

#include "cut_creator.H"
#include "cut_point_impl.H"
#include "cut_levelsetside.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_integrationcell.H"

#include "cut_parallel.H"

#include "../drt_geometry/element_volume.H"


/*-------------------------------------------------------------------------------------*
 * constructor
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Mesh::Mesh( Options & options, double norm, Teuchos::RCP<PointPool> pp, bool cutmesh )
  : setup_( true ),
    options_( options ),
    norm_( norm ),
    pp_( pp ),
    cutmesh_( cutmesh )
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
GEO::CUT::Point* GEO::CUT::Mesh::NewPoint( const double * x, Edge * cut_edge, Side * cut_side )
{
  bb_.AddPoint( x ); // add the point to the mesh's bounding box
  //Point* p = pp_->NewPoint( x, cut_edge, cut_side, setup_ ? SETUPNODECATCHTOL : MINIMALTOL );
  Point* p = pp_->NewPoint( x, cut_edge, cut_side, MINIMALTOL ); // add the point in the point pool
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


/*---------------------------------------------------------------*
 * check whether any one of the cut sides cuts another cut side
 * Our algorithm fails in such cases
 *---------------------------------------------------------------*/
bool GEO::CUT::Mesh::DetectSelfCut()
{
  std::vector<Side*> mysides;
  mysides.reserve( sides_.size() );

  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side * side = &*i->second;
    mysides.push_back( side );
  }

  std::sort( mysides.begin(), mysides.end() ); // sorting to use std::binary_search

  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    {
      BoundingBox sidebox( side );
      plain_side_set sides;
      pp_->CollectSides( sidebox, sides );  // all sides that fall within BB of considered side
      sides.erase( &side );
      for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); )
      {
        Side * s = *i;
        if ( not std::binary_search( mysides.begin(), mysides.end(), s ) or
             side.HaveCommonNode( *s ) )  // if common nodes exist -- these sides are connected (always??)
        {
          set_erase( sides, i );
        }
        else
        {
          ++i;
        }
      }
      for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        {
          bool normal_cut  = side.FindCutPoints( *this, NULL, *s, 0 );
          bool reverse_cut = s->FindCutPoints( *this, NULL, side, 0 );
          if ( normal_cut or reverse_cut )
          {
            //std::cout << side << "\n" << ( *s ) << "\n";
            return true;
          }
        }
      }
    }
  }
  return false;
}

#if 0
void GEO::CUT::Mesh::SelfCut()
{
  plain_facet_set facets;
  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    {
      BoundingBox sidebox( side );
      plain_side_set sides;
      pp_->CollectSides( sidebox, sides );
      sides.erase( &side );
      for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); )
      {
        Side * s = *i;
        if ( side.HaveCommonEdge( *s ) )
        {
          sides.erase( i++ );
        }
        else
        {
          ++i;
        }
      }
      for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        {
          side.FindCutPoints( *this, NULL, *s );
          s->FindCutPoints( *this, NULL, side );
        }
      }
      bool cut = false;
      for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        {
          bool normal_cut  = side.FindCutLines( *this, NULL, *s );
          bool reverse_cut = s->FindCutLines( *this, NULL, side );
          if ( normal_cut or reverse_cut )
            cut = true;
        }
      }
      if ( cut )
      {
        for ( plain_side_set::iterator i=sides.begin(); i!=sides.end(); ++i )
        {
          Side * s = *i;
          {
            SideSideCutFilter filter( &side, s );
            side.MakeOwnedSideFacets( *this, filter, facets );
          }
        }
        side.MakeSideCutFacets( *this, NULL, facets );
      }
    }
  }
  for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    f->CreateLinearElements( *this );
  }
}
#endif


/*-----------------------------------------------------------------*
 * Cuts the background elements of the mesh with all the cut sides *
 *-----------------------------------------------------------------*/
void GEO::CUT::Mesh::Cut( Mesh & mesh, plain_element_set & elements_done, int recursion )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- CUT (incl. tetmesh-cut)" );

  if ( DetectSelfCut() )
    throw std::runtime_error( "cut surface with self cut not supported" );

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
 *------------------------------------------------------------------------*/
void GEO::CUT::Mesh::Cut( Side & side, const plain_element_set & done, plain_element_set & elements_done, int recursion )
{
  BoundingBox sidebox( side ); // define a bounding box around the maybe twisted side to determine a preselection of cutting sides and elements
  plain_element_set elements;
  
#if(0)
  // REMARK: do not use pp_->CollectElements anymore
  // it can happen that some intersections between elements and sides
  // are not detected, because the octtree bounding boxes can lie
  // within a real background element. Such small bounding boxes
  // do not have to contain any points adjacent to elements,
  // then elements have not been found
  //
  // schott 10/2012
  
  //pp_->CollectElements( sidebox, elements ); // find involved elements (octree-based)

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
  for(std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
      i!=shadow_elements_.end();
      i++)
  {
    Element* e= &**i;
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

#endif


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
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.Cut( *this, side, 0 );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.Cut( *this, side, 0 );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
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


/*-------------------------------------------------------------------------------------*
 * create cut lines based on the point cloud
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeCutLines()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- MakeCutLines" );

  Creator creator;
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeCutLines( *this, creator );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeCutLines( *this, creator );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }

  creator.Execute( *this );
}


/*-------------------------------------------------------------------------------------*
 * create facets based on the cut lines
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeFacets()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- MakeFacets" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeFacets( *this );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeFacets( *this );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
}


/*-------------------------------------------------------------------------------------*
 * create volumecells based on created facets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MakeVolumeCells()
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- MakeVolumeCells" );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeVolumeCells( *this );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MakeVolumeCells( *this );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
}


/*-------------------------------------------------------------------------------------*
 * find node positions and propagate the positions to facets, points and volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::FindNodePositions()
{
  // On multiple cuts former outside positions can become inside
  // positions. Thus reset all outside positions.
  pp_->ResetOutsidePoints();

  // get nodal positions from elements
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.FindNodePositions();
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.FindNodePositions();
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  // find undecided nodes
  // * for serial simulations all node positions should be set
  // * for parallel simulations there can be some undecided nodes
//  CheckForUndecidedNodePositions();

}


/*-------------------------------------------------------------------------------------*
 * ?
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
      double lsv = n->LSV();
      if ( lsv > 0 )
      {
        p->Position( Point::outside );
      }
      else if ( lsv < 0 )
      {
        p->Position( Point::inside );
      }
      else
      {
        //throw std::runtime_error( "undecided nodal point on levelset
        //surface" );
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
bool GEO::CUT::Mesh::CheckForUndecidedNodePositions( std::map<int, int> & undecided_nodes )
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
        throw std::runtime_error( "node with node-Id <0 found" );
      }
      // insert pair of nid and undecided PointPosition
      undecided_nodes.insert(std::pair<int,int>(n->Id(), pos));

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
void GEO::CUT::Mesh::CreateIntegrationCells( int count, bool levelset )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.CreateIntegrationCells( *this, count+1, levelset );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.CreateIntegrationCells( *this, count+1, levelset );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
}


/*-------------------------------------------------------------------------------------*
 * Call the moment fitting method for each element to generate the Gaussian integration rule
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::MomentFitGaussWeights(bool include_inner, std::string Bcellgausstype)
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MomentFitGaussWeights( *this, include_inner, Bcellgausstype );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.MomentFitGaussWeights( *this, include_inner, Bcellgausstype );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
}


/*-------------------------------------------------------------------------------------*
 * Call the DirectDivergence method for each element to generate the Gaussian integration rule
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DirectDivergenceGaussRule(bool include_inner, std::string Bcellgausstype)
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.DirectDivergenceGaussRule( *this, include_inner, Bcellgausstype );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.DirectDivergenceGaussRule( *this, include_inner, Bcellgausstype );
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
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
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.RemoveEmptyVolumeCells();
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
#ifndef DEBUGCUTLIBRARY
    try
    {
#endif
      e.RemoveEmptyVolumeCells();
#ifndef DEBUGCUTLIBRARY
    }
    catch ( std::runtime_error & err )
    {
      e.DebugDump();
      throw;
    }
#endif
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
void GEO::CUT::Mesh::TestElementVolume( bool fatal )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    TestElementVolume( e.Shape(), e, fatal );
  }
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
    TestElementVolume( e.Shape(), e, fatal );
  }
}


/*-------------------------------------------------------------------------------------*
 * Find the difference between the volume of background element and the sum of volume
 * of all integration cells. There should be no difference between these two
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestElementVolume( DRT::Element::DiscretizationType shape, Element & e, bool fatal )
{
  if ( e.IsCut() )
  {
    const std::vector<Node*> & nodes = e.Nodes();
    Epetra_SerialDenseMatrix xyze( 3, nodes.size() );
    for ( unsigned i=0; i<nodes.size(); ++i )
    {
      nodes[i]->Coordinates( &xyze( 0, i ) );
    }

    double ev = GEO::ElementVolume( e.Shape(), xyze );

    int numgp = 0;
    int numic = 0;
    int numbc = 0;
    double cv = 0;
    double ba = 0;
    const plain_volumecell_set & cells = e.VolumeCells();
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

#ifdef DEBUGCUTLIBRARY
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
#endif

    if ( fatal and fabs( volume_error ) > 1e-5 )
    {
      std::stringstream err;
      err << " !!!!!!!!!!! volume test failed: !!!!!!!!!!!!!!!"
          << "eleID=" << e.Id() << "  "
          << "ve=" << ev << "  "
          << "vc=" << cv << "  "
          << "vd= "<< ev-cv << "  "
          << "err=" << volume_error;
      cout << err.str() << endl;
//      throw std::runtime_error( err.str() );
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
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element & e = **i;
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
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmsh( std::string name )
{
  std::ofstream file( name.c_str() );
  file << "View \"" << name << "\" {\n";
  if ( elements_.size() > 0 )
  {
    for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
          i!=elements_.end();
          ++i )
    {
      Element & e = *i->second;
      //e.DumpGmsh( file );
      {
        const std::vector<Node*> & nodes = e.Nodes();
        char elementtype;
        switch ( nodes.size() )
        {
        case 8:
          elementtype = 'H';
          break;
        case 4:
          elementtype = 'S';
          break;
        case 6:
          elementtype = 'I';
          break;
        default:
          throw std::runtime_error( "unknown element type" );
        }
        DumpGmsh( file, nodes, elementtype );
      }
    }
    for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
          i!=shadow_elements_.end();
          ++i )
    {
      Element & e = **i;
      const std::vector<Node*> & nodes = e.Nodes();
      char elementtype;
      switch ( nodes.size() )
      {
      case 8:
        elementtype = 'H';
        break;
      case 4:
        elementtype = 'S';
        break;
      case 6:
        elementtype = 'I';
        break;
      default:
        throw std::runtime_error( "unknown element type" );
      }
      DumpGmsh( file, nodes, elementtype );
    }
  }
  else
  {
    for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
          i!=sides_.end();
          ++i )
    {
      Side & s = *i->second;
      {
        const std::vector<Node*> & nodes = s.Nodes();
        char elementtype;
        switch ( nodes.size() )
        {
        case 3:
          elementtype = 'T';
          break;
        case 4:
          elementtype = 'Q';
          break;
        default:
          throw std::runtime_error( "unknown element type" );
        }
        DumpGmsh( file, nodes, elementtype );
      }
    }
  }
  file << "};\n";
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmsh( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype )
{
  file << "S" << elementtype
       << "(";
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    double x[3];
    n->Coordinates( x );
    if ( i!=nodes.begin() )
      file << ",";
    file << x[0] << "," << x[1] << "," << x[2];
  }
  file << "){";
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    if ( i!=nodes.begin() )
      file << ",";
    file << p->Position();
    //file << n->DofSets().size();
  }
  file << "};\n";
}


/*-------------------------------------------------------------------------------------*
 * ?
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::DumpGmshVolumeCells( std::string name, bool include_inner )
{
  std::ofstream file( name.c_str() );
  int count = 0;

  file << "View \"VolumeCells\" {\n";
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * vc = &**i;

//    if ( true  ) // cout all volumecells - inside and outside
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
  for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
        i!=shadow_elements_.end();
        ++i )
  {
    Element * e = &**i;
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


  file << "View \"Node-Positions\" {\n";
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    const double * x = p->X();
    file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" <<  p->Position() << "};\n";
  }
  file << "};\n";
}


/*-------------------------------------------------------------------------------------*
 * ?
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

  file << "View \"BoundaryCells\" {\n";
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=boundarycells_.begin();
        i!=boundarycells_.end();
        ++i )
  {
    BoundaryCell * bc = &**i;
    if ( bc->IsValid() )
      bc->DumpGmsh( file );
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
/*	  else
	  {

	  }*/
#endif
  }
  file << "};\n";

  file << "View \"BoundaryCells\" {\n";
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=boundarycells_.begin();
        i!=boundarycells_.end();
        ++i )
  {
    BoundaryCell * bc = &**i;
    if ( bc->IsValid() )
      bc->DumpGmsh( file );
  }
  file << "};\n";
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
    If node with the given id exists return the node, else create a new node with
    given coordinates and levelset value
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid, const double * xyz, double lsv )
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
  Point * p = NewPoint( xyz, NULL, NULL );
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
 * ?
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
 * get all sides (subsides) that belong to the parent side with sid
 *-------------------------------------------------------------------------------------*/
const std::vector<GEO::CUT::Side*> & GEO::CUT::Mesh::GetSides( int sid )
{
  std::map<int, std::vector<Side*> >::iterator i=cut_sides_.find( sid );
  if ( i!=cut_sides_.end() )
  {
    return i->second;
  }
  throw std::runtime_error( "no side with given id" );
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
  if ( sid > -1 )
  {
    cut_sides_[sid].push_back( s );
  }
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
    shadow_elements_.push_back( Teuchos::rcp( e ) );
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
void GEO::CUT::Mesh::CreateSideIds()
{
#if 0
  int localmaxid = cut_sides_.rbegin()->first;
  int globalmaxid;

  int err = MPI_Allreduce( &localmaxid, &globalmaxid, 1, MPI_INT, MPI_MAX, comm );
  if ( err!=0 )
    throw std::runtime_error( "mpi error" );

  int myrank;
  int numproc;
  MPI_Comm_rank( comm, &myrank );
  MPI_Comm_size( comm, &numproc );

  //std::vector<int>
#endif

  int lastid = 0;
  if ( cut_sides_.size() > 0 )
  {
    lastid = cut_sides_.rbegin()->first;
  }

  //int numelementsides = sides_.size() - cut_sides_.size();

  for ( std::map<plain_int_set, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side * s = &*i->second;
    if ( s->Id() < 0 )
    {
      lastid += 1;
      s->SetId( lastid );
      cut_sides_[lastid].push_back( s );
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
 * ???
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
 * ???
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::Mesh::TestFacetArea()
{
  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = &**i;

    // This is a crude test. We do not demand so much here...
    f->TestFacetArea( 1e-7 );
  }
}


