
#include <stack>

#include "cut_tetmesh.H"
#include "cut_levelsetside.H"

#ifdef QHULL
extern "C" {
#include <qhull/qhull_a.h>
}
#endif

namespace GEO
{
  namespace CUT
  {
    namespace
    {
      class NullFile
      {
      public:
        NullFile()
          : f_( NULL )
        {
        }

        ~NullFile()
        {
          if ( f_ != NULL )
            fclose( f_ );
        }

        operator FILE* ()
        {
          if ( f_ == NULL )
            f_ = fopen( "/dev/null", "w" );
          return f_;
        }

      private:
        FILE * f_;
      };
    }
  }
}


GEO::CUT::TetMesh::TetMesh( const std::vector<Point*> & points,
                            const std::set<Facet*> & facets )
  : points_( points ),
    facets_( facets )
{
  std::vector<std::vector<int> > original_tets;

  CallQHull( points, original_tets );

  tets_.reserve( original_tets.size() );
  for ( std::vector<std::vector<int> >::iterator i=original_tets.begin(); i!=original_tets.end(); ++i )
  {
    std::vector<int> & t = *i;
    std::vector<Point*> tet;
    tet.reserve( t.size() );
    for ( std::vector<int>::iterator i=t.begin(); i!=t.end(); ++i )
    {
      tet.push_back( points[*i] );
    }
    if ( IsValidTet( tet ) )
    {
      //tets_.push_back( std::vector<int>() );
      //std::swap( tets_.back(), t );
      tets_.push_back( t );
    }
  }

  TestUsedPoints( tets_ );

  Init();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteCells();
#endif

  // find the cut surface

  FindTriFacets();
  FindCutSurfaceTris();

  //RotateNonMatchingCutTets();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteSurfaceCells();
#endif

  DeactivateCutSurfaceTets();
  ActivateCutSurfaceTris();
  DeactivateOverlappingTris();
  TestCutSurface();

  // find the active tets

  ActivateSingleTetsOnCutSurface();
  ActivateDoubleTetsOnCutSurface();
  ActivateUsedCutSurfaceTris();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteActiveCells();
  GmshWriteActiveSurfaceCells();
#endif

  ClearExternalTets();

  FixBrokenTets();
  FillCutSides( sides_xyz_ );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteSurfaceTris();
#endif

  TestUsedPoints( tets_ );

  // Clear all tets that have been deactivated because they are too small. We
  // might lose some points in rare cases. That is alright.
  ClearExternalTets();
}

void GEO::CUT::TetMesh::Init()
{
  unsigned numtets = tets_.size();

  for ( unsigned i=0; i<numtets; ++i )
  {
    const std::vector<int> & t = tets_[i];
    tet_entities_.push_back( Entity<4>( i, Handle<4>( &t[0] ) ) );
  }

  for ( std::vector<Entity<4> >::iterator i=tet_entities_.begin();
        i!=tet_entities_.end();
        ++i )
  {
    Entity<4> & e = *i;
    e.CreateChildren( tet_surfaces_ );
  }

  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & e = i->second;
    e.CreateChildren( tet_lines_ );
  }

  accept_tets_.resize( tets_.size() );
  std::fill( accept_tets_.begin(), accept_tets_.end(), true );

  accept_tris_.resize( tet_surfaces_.size() );
  std::fill( accept_tris_.begin(), accept_tris_.end(), true );

  active_surface_tris_.resize( tet_surfaces_.size() );
  std::fill( active_surface_tris_.begin(), active_surface_tris_.end(), false );
}

void GEO::CUT::TetMesh::CallQHull( const std::vector<Point*> & points,
                                   std::vector<std::vector<int> > & tets )
{
#ifdef QHULL

  const int dim = 3;
  const int n = points.size();

  if ( n < 4 )
  {
    throw std::runtime_error( "illegal element topology" );
  }
  if ( n == 4 )
  {
    throw std::runtime_error( "no need to triangulate" );
  }

  std::vector<double> coordinates( dim*n );

  for ( int i=0; i<n; ++i )
  {
    Point * p = points[i];
    p->Coordinates( &coordinates[dim*i] );
  }

  boolT ismalloc = false;

  // a set of option we try to process the input with
  // Qz seems to be required for rotational symmetric input
  std::vector<std::string> options;
  options.push_back( "qhull d Qt Qbb Qc Pp" );
  options.push_back( "qhull d Qt Qbb Qc Qz Pp" );

  // If you want some debugging information replace the 0 pointer
  // with stdout or some other file open for writing.

  FILE * outfile = 0;
//#define QHULL_DEBUG_OUTPUT
#ifdef QHULL_DEBUG_OUTPUT
  FILE * errfile = stderr;
#else
  static NullFile errfile;
#endif

#if 0
  static FILE * debug_f;
  if ( debug_f==NULL )
    debug_f = fopen( "qhull.debug", "w" );
  for ( int i=0; i<n; ++i )
  {
    Point * p = points[i];
    const double * x = p->X();
    fprintf( debug_f, "% .20f % .20f % .20f ", x[0], x[1], x[2] );
  }
  fprintf( debug_f, "\n");
  fflush( debug_f );
#endif

  for ( std::vector<std::string>::iterator i=options.begin(); i!=options.end(); ++i )
  {
    std::string & ostr = *i;
    if ( not qh_new_qhull( dim, n, &coordinates[0], ismalloc, const_cast<char*>( ostr.c_str() ), outfile, errfile ) )
    {
      // triangulate non-simplicial facets
      qh_triangulate ();

      facetT *facet;
      vertexT *vertex, **vertexp;
      int nf = 0;

      FORALLfacets
      {
        if ( not facet->upperdelaunay )
          nf++;

        // Double check
        if ( not facet->simplicial )
        {
          throw std::runtime_error("Qhull returned non-simplicial facets -- try delaunayn with different options");
        }
      }

      tets.reserve( nf );

      FORALLfacets
      {
        if ( not facet->upperdelaunay )
        {
          std::vector<int> ids;
          ids.reserve( dim+1 );

          FOREACHvertex_(facet->vertices)
          {
            // if delaunayn crashes, enable this check
#if 0
            if (j > dim)
            {
              std::runtime_error("internal error. Qhull returned non-tetsicial facets");
            }
#endif

            int p = qh_pointid(vertex->point);
            if ( p >= n )
            {
              throw std::runtime_error( "new node in delaunay" );
            }
            ids.push_back( p );
          }

          tets.push_back( ids );
        }
      }
    }

    // free long memory
    qh_freeqhull( not qh_ALL );

    // free short memory and memory allocator
    int curlong, totlong;
    qh_memfreeshort( &curlong, &totlong );

    if ( curlong or totlong )
    {
      std::stringstream str;
      str << "did not free " << totlong << " bytes of long memory (" << curlong << " pieces)";
      throw std::runtime_error( str.str() );
    }

    if ( tets.size() > 0 )
    {
      std::set<int> used_points;
      for ( std::vector<std::vector<int> >::iterator i=tets.begin(); i!=tets.end(); ++i )
      {
        std::vector<int> & t = *i;
        std::copy( t.begin(), t.end(), std::inserter( used_points, used_points.begin() ) );
      }
      if ( used_points.size() == points.size() )
      {
        return;
      }
      //throw std::runtime_error( "failed to triangulate all points" );

      // failed! start a new iteration.
      tets.clear();
    }
  }

#if 1
  // debug output to be read by qhull_test programm
  FILE * debug_f = fopen( "qhull.debug", "w" );
  for ( int i=0; i<n; ++i )
  {
    Point * p = points[i];
    const double * x = p->X();
    fprintf( debug_f, "% .20f % .20f % .20f ", x[0], x[1], x[2] );
  }
  fprintf( debug_f, "\n");
  fclose( debug_f );
#endif

  throw std::runtime_error( "qhull failed" );

#else
  throw std::runtime_error( "qhull not available" );
#endif
}

bool GEO::CUT::TetMesh::IsValidTet( const std::vector<Point*> & t )
{
  std::set<Side*> sides;
  FindCommonSides( t, sides );
  if ( sides.size()==0 )
  {
    std::set<Facet*> facets;
    FindCommonFacets( t, facets );
    if ( facets.size()==0 )
    {
      return true;
    }
    for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->IsTriangulated() )
      {
        return false;
      }
    }
    return true;
  }
  if ( sides.size()==1 )
  {
    if ( dynamic_cast<LevelSetSide*>( *sides.begin() ) != NULL )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::TetMesh::TestUsedPoints( const std::vector<std::vector<int> > & tets )
{
  std::set<int> used_points;
  for ( std::vector<std::vector<int> >::const_iterator i=tets.begin(); i!=tets.end(); ++i )
  {
    const std::vector<int> & t = *i;
    std::copy( t.begin(), t.end(), std::inserter( used_points, used_points.begin() ) );
  }
  if ( used_points.size() != points_.size() )
  {
    throw std::runtime_error( "failed to triangulate all points" );
  }
}

void GEO::CUT::TetMesh::FixBrokenTets()
{
  for ( std::vector<std::vector<int> >::iterator i=tets_.begin();
        i!=tets_.end();
        ++i )
  {
    std::vector<int> & t = *i;

    // create planes consisting of 3 nodes each
    LINALG::Matrix<3,1> p0( points_[t[0]]->X() );
    LINALG::Matrix<3,1> p1( points_[t[1]]->X() );
    LINALG::Matrix<3,1> p2( points_[t[2]]->X() );
    LINALG::Matrix<3,1> p3( points_[t[3]]->X() );

    LINALG::Matrix<3,1> v01;
    LINALG::Matrix<3,1> v02;
    LINALG::Matrix<3,1> v03;

    v01.Update( 1, p1, -1, p0, 0 );
    v02.Update( 1, p2, -1, p0, 0 );
    v03.Update( 1, p3, -1, p0, 0 );

    // create 4 normal vectors to each tet surface plane
    LINALG::Matrix<3,1> nplane012;

    // cross product
    nplane012(0) = v01(1)*v02(2) - v01(2)*v02(1);
    nplane012(1) = v01(2)*v02(0) - v01(0)*v02(2);
    nplane012(2) = v01(0)*v02(1) - v01(1)*v02(0);

    // compute norm (area) of plane
    double norm012 = nplane012.Norm2();

    // compute normal distance of point to plane of the three remaining points
    double distance = nplane012.Dot( v03 );

    double vol_tet = distance / 6.0;

    // Deactivate all tets that are too small. We might still need the tet to
    // create a cut surface tri. Afterwards we will discard it.

    if ( fabs( vol_tet ) < 1e-10 )
    {
      accept_tets_[i - tets_.begin()] = false;
    }
    else if ( fabs( distance / norm012 ) < 1e-7 )
    {
      accept_tets_[i - tets_.begin()] = false;
    }

    // tet numbering wrong exchange 1 with 3
    if ( distance < 0 )
    {
      std::swap( t[1], t[3] );
    }
  }
}

void GEO::CUT::TetMesh::FindTriFacets()
{
  std::map<Handle<3>, std::set<Facet*> > double_tris;

  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & tri = i->second;

    std::vector<Point*> tri_points;
    tri.PointVector( points_, tri_points );

    std::set<Facet*> facets;
    FindCommonFacets( tri_points, facets );

#if 0
    // remove all facets that do not belong to the volume cell
    std::set<Facet*> intersection;
    std::set_intersection( facets.begin(), facets.end(),
                           facets_.begin(), facets_.end(),
                           std::inserter( intersection, intersection.begin() ) );
    std::swap( facets, intersection );
#endif

    if ( facets.size() == 1 )
    {
      tri_facets_[&tri] = *facets.begin();
    }
    else if ( facets.size() > 1 )
    {
      std::swap( double_tris[tri.GetHandle()], facets );
    }
  }

  if ( double_tris.size() == 0 )
    return;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteFacetTris();
#endif

  for ( std::map<Handle<3>, std::set<Facet*> >::iterator itris=double_tris.begin();
        itris!=double_tris.end();
    )
  {
    Entity<3> & tri = tet_surfaces_[itris->first];
    std::set<Facet*> & facets = itris->second;

    std::vector<Point*> tri_points;
    tri.PointVector( points_, tri_points );

    bool found = false;
    for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->IsTriangle( tri_points ) )
      {
        tri_facets_[&tri] = f;
        found = true;
        break;
      }
    }

    if ( found )
    {
      double_tris.erase( itris++ );
    }
    else
    {
      ++itris;
    }
  }

  if ( double_tris.size() == 0 )
    return;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteFacetTris();
#endif

  for ( std::map<Entity<3>*, Facet*>::iterator i=tri_facets_.begin(); i!=tri_facets_.end(); ++i )
  {
    Entity<3> & tri = *i->first;
    Facet * f = i->second;

    const std::vector<Entity<2>*> & lines = tri.Children();
    for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
    {
      Entity<2> & line = **i;
      bool isfacetline = f->IsLine( points_[line[0]], points_[line[1]] );
      if ( not isfacetline )
      {
        const std::vector<Entity<3>*> & tris = line.Parents();
        for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
        {
          Entity<3> * other_tri = *i;
          if ( tri_facets_.count( other_tri )==0 )
          {
            std::map<Handle<3>, std::set<Facet*> >::iterator j =
              double_tris.find( other_tri->GetHandle() );
            if ( j != double_tris.end() )
            {
              if ( j->second.count( f )==0 )
              {
                throw std::runtime_error( "assigned facet not in facet list" );
              }
              tri_facets_[other_tri] = f;
              double_tris.erase( j );
            }
          }
        }
      }
    }
  }

  if ( double_tris.size() == 0 )
    return;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteFacetTris();
#endif

  while ( double_tris.size() > 0 )
  {
    std::set<Entity<3>*> patch;
    std::stack<Entity<3>*> work;
    std::set<Entity<2>*> patch_lines;

    work.push( &tet_surfaces_[double_tris.begin()->first] );
    std::set<Facet*> & facets = double_tris.begin()->second;

    //Facet * f = *facets.begin();

    while ( not work.empty() )
    {
      Entity<3> * tri = work.top();
      work.pop();

      patch.insert( tri );

      const std::vector<Entity<2>*> & lines = tri->Children();
      for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        Entity<2> & line = **i;
        bool isfacetline = false;
        for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
        {
          Facet * f = *i;
          if ( f->IsLine( points_[line[0]], points_[line[1]] ) )
          {
            isfacetline = true;
            break;
          }
        }
        if ( not isfacetline )
        {
          const std::vector<Entity<3>*> & tris = line.Parents();
          for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
          {
            Entity<3> * other_tri = *i;
            if ( tri_facets_.count( other_tri )==0 and patch.count( other_tri )==0 )
            {
              std::map<Handle<3>, std::set<Facet*> >::iterator j =
                double_tris.find( other_tri->GetHandle() );
              if ( j != double_tris.end() )
              {
                work.push( other_tri );
                std::set<Facet*> & other_facets = j->second;
                std::set<Facet*> intersection;
                std::set_intersection( facets.begin(), facets.end(),
                                       other_facets.begin(), other_facets.end(),
                                       std::inserter( intersection, intersection.begin() ) );
                std::swap( facets, intersection );
                if ( facets.size()==0 )
                  throw std::runtime_error( "no common facet in patch" );
              }
            }
          }
        }
        else
        {
          patch_lines.insert( &line );
        }
      }
    }

    std::set<Point*> patch_points;
    for ( std::set<Entity<2>*>::iterator i=patch_lines.begin(); i!=patch_lines.end(); ++i )
    {
      Entity<2> & line = **i;
      patch_points.insert( points_[line[0]] );
      patch_points.insert( points_[line[1]] );
    }

    std::vector<Point*> final_points;
    final_points.reserve( patch_points.size() );
    final_points.assign( patch_points.begin(), patch_points.end() );

    for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->NumPoints()==patch_points.size() and f->Contains( final_points ) )
      {
        for ( std::set<Entity<3>*>::iterator i=patch.begin(); i!=patch.end(); ++i )
        {
          Entity<3> * tri = *i;
          tri_facets_[tri] = f;
          double_tris.erase( tri->GetHandle() );
        }
        patch.clear();
        break;
      }
    }

    // These tris are done, even if we did not find a facet. There are broken
    // tris...
    for ( std::set<Entity<3>*>::iterator i=patch.begin(); i!=patch.end(); ++i )
    {
      Entity<3> * tri = *i;
      double_tris.erase( tri->GetHandle() );
    }
  }

  if ( double_tris.size() > 0 )
    throw std::runtime_error( "failed to assign triangle facets" );
}

/// find all tris at a cut surface
void GEO::CUT::TetMesh::FindCutSurfaceTris()
{
  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & tri = i->second;
    tri.FindCutEntities( points_, facets_, cutsurface_ );
  }
}

#if 0
void GEO::CUT::TetMesh::RotateNonMatchingCutTets()
{
  std::map<Facet*, std::vector<Entity<3>*> >::iterator ic = cutsurface_.find( NULL );
  if ( ic!=cutsurface_.end() )
  {
    std::vector<Entity<3>*> & tris = ic->second;

    std::set<int> tri_ids;
    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;
      tri_ids.insert( tri->Id() );
    }

    for ( std::vector<Entity<3>*>::iterator it=tris.begin(); it!=tris.end(); ++it )
    {
      Entity<3> * my_tri = *it;

      if ( tri_ids.count( my_tri->Id() )==0 )
        continue;

      // find my neighbor that needs to be rotated with me
      std::set<Entity<3>*> neighbors;
      my_tri->Neighbors( neighbors );

      for ( std::set<Entity<3>*>::iterator i=neighbors.begin();
            i!=neighbors.end();
        )
      {
        Entity<3> * other = *i;
        if ( tri_ids.count( other->Id() )==0 )
        {
          neighbors.erase( i++ );
        }
        else
        {
          ++i;
        }
      }

      Entity<3> * other_tri = NULL;
      if ( neighbors.size()==1 )
      {
        other_tri = *neighbors.begin();
      }
      else
      {
        throw std::runtime_error( "only exact matches for triangle rotation supported" );
      }

      // remove done tris (before we are done, since the pointer changes!)
      tri_ids.erase( other_tri->Id() );
      tri_ids.erase( my_tri->Id() );

      // test the findings
      const std::vector<Entity<4>*> & my_parents = my_tri->Parents();
      const std::vector<Entity<4>*> & other_parents = other_tri->Parents();

      if ( my_parents.size()!=1 or other_parents.size()!=1 )
      {
        throw std::runtime_error( "expect one tet for each surface tri" );
      }

      Entity<4> * my_tet = my_parents[0];
      Entity<4> * other_tet = other_parents[0];

      int p1 = my_tet->OtherPoint( my_tri );
      int p2 = other_tet->OtherPoint( other_tri );
      Point * p = points_[p1];

      if ( p1!=p2 or p->Position()==Point::oncutsurface )
      {
        throw std::runtime_error( "tets must share a non cut surface point" );
      }

      // rotate

      Entity<2> * middle_line = my_tri->CommonChild( other_tri );
      Entity<3> * middle_tri = my_tet->CommonChild( other_tet );

      int points[2];
      points[0] = my_tri->OtherPoint( middle_line );
      points[1] = other_tri->OtherPoint( middle_line );
      std::sort( points, points+2 );
      Handle<2> new_line_handle( points );

      Handle<3> new_tri_handle = new_line_handle + p1;
      Handle<3> my_tri_handle = new_line_handle + ( *middle_line )[0];
      Handle<3> other_tri_handle = new_line_handle + ( *middle_line )[1];

      Handle<4> my_tet_handle = my_tri_handle + p1;
      Handle<4> other_tet_handle = other_tri_handle + p1;

      middle_line = SwapHandle( middle_line, new_line_handle, tet_lines_ );

      middle_tri = SwapHandle( middle_tri, new_tri_handle, tet_surfaces_ );
      my_tri = SwapHandle( my_tri, my_tri_handle, tet_surfaces_ );
      other_tri = SwapHandle( other_tri, other_tri_handle, tet_surfaces_ );

      SwapTetHandle( my_tet, my_tet_handle );
      SwapTetHandle( other_tet, other_tet_handle );

      my_tet->CreateChildren( tet_surfaces_ );
      other_tet->CreateChildren( tet_surfaces_ );

      my_tri->FindCutEntities( points_, facets_, cutsurface_ );
      other_tri->FindCutEntities( points_, facets_, cutsurface_ );
    }

#if 0
    std::vector<Entity<3>*>::iterator i =
      std::remove_if( tris.begin(), tris.end(),
                      std::bind1st( std::equal_to<Entity<3>*>(),
                                    static_cast<Entity<3>*>( NULL ) ) );
    tris.erase( i, tris.end() );

    if ( tris.size() > 0 )
      throw std::runtime_error( "unhandled tris left" );
#endif
    if ( tri_ids.size() > 0 )
      throw std::runtime_error( "unhandled tris left" );

    cutsurface_.erase( ic );
  }
}
#endif

/// deactivate all tets that are on a cut surface with all four points
void GEO::CUT::TetMesh::DeactivateCutSurfaceTets()
{
  // All tets that are not at a cutsurface are included --- remove all
  // the others.

  // There might be external tets that are connected to cut sides but
  // fail the cutsurface test. Deactivate any such tet.

  unsigned numtets = tets_.size();
  for ( unsigned i=0; i<numtets; ++i )
  {
    if ( tet_entities_[i].AllOnCutSurface( points_ ) )
    {
      accept_tets_[i] = false;
    }
  }
}

void GEO::CUT::TetMesh::ActivateCutSurfaceTris()
{
  for ( std::map<Facet*, std::vector<Entity<3>*> >::iterator i=cutsurface_.begin();
        i!=cutsurface_.end();
        ++i )
  {
    std::vector<Entity<3>*> & tris = i->second;

    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;

      active_surface_tris_[tri->Id()] = true;
    }
  }
}

void GEO::CUT::TetMesh::DeactivateOverlappingTris()
{
  // Test if the set of tris we have is unique. It might not be, if
  // qhull returns non matching interfaces.
  //
  // The list is not guaranteed to find the right pairs. But it finds
  // the right number of subsets.
  for ( std::map<Facet*, std::vector<Entity<3>*> >::iterator i=cutsurface_.begin();
        i!=cutsurface_.end();
        ++i )
  {
    std::vector<Entity<3>*> & tris = i->second;

    std::vector<std::set<int> > unique_tris;

    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> & tri = **i;

      bool found = false;
      for ( std::vector<std::set<int> >::iterator i=unique_tris.begin();
            i!=unique_tris.end();
            ++i )
      {
        std::set<int> & triset = *i;
        unsigned size = triset.size();
        std::copy( tri(), tri()+3, std::inserter( triset, triset.begin() ) );
        if ( size!=triset.size() )
        {
          found = true;
          break;
        }
      }
      if ( not found )
      {
        std::set<int> newset;
        unique_tris.push_back( newset );
        std::copy( tri(), tri()+3, std::inserter( unique_tris.back(),
                                                  unique_tris.back().begin() ) );
      }
    }
    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> & tri = **i;

      if ( unique_tris.size()==0 )
      {
        throw std::runtime_error( "no tris here" );
      }
      if ( unique_tris.size()>1 )
      {
        accept_tris_[tri.Id()] = false;
      }
    }
  }
}

/// if there is just one tet at an cut surface it is active
void GEO::CUT::TetMesh::ActivateSingleTetsOnCutSurface()
{
  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & tri = i->second;

    if ( accept_tris_[tri.Id()] )
    {
      const std::vector<Entity<4>*> & tets = tri.Parents();
      if ( tets.size()==1 )
      {
        if ( active_surface_tris_[tri.Id()] )
        {
          accept_tets_[tets[0]->Id()] = true;
        }
        //else if ( tri.AtElementSide( points_, facets_ ) )
        else
        {
          std::map<Entity<3>*, Facet*>::iterator j = tri_facets_.find( &tri );
          if ( j != tri_facets_.end() )
          {
            Facet * f = j->second;
            if ( facets_.count( f ) > 0 and f->SideId() < 0 )
            {
              accept_tets_[tets[0]->Id()] = true;
            }
          }
        }
      }
    }
  }
}

/// activate the right tets at the cut surface
void GEO::CUT::TetMesh::ActivateDoubleTetsOnCutSurface()
{
  for ( std::vector<int>::iterator i = std::find( accept_tets_.begin(), accept_tets_.end(), 1 );
        i != accept_tets_.end();
        i = std::find( ++i, accept_tets_.end(), 1 ) )
  {
    unsigned pos = i - accept_tets_.begin();
    Entity<4> & e = tet_entities_[pos];
    ActivateSurroundingTets( &e );
  }
}

void GEO::CUT::TetMesh::ActivateSurroundingTets( Entity<4> * tet )
{
  const std::vector<Entity<3>*> & tris = tet->Children();
  for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> * tri = *i;
    if ( not active_surface_tris_[tri->Id()] and accept_tris_[tri->Id()] )
    {
      const std::vector<Entity<4>*> & tets = tri->Parents();
      if ( tets.size()==2 )
      {
        Entity<4> * other = NULL;
        if ( tets[0]==tet )
        {
          other = tets[1];
        }
        else if ( tets[1]==tet )
        {
          other = tets[0];
        }
        else
        {
          throw std::runtime_error( "expected parent tet not in list" );
        }
        if ( not accept_tets_[other->Id()] )
        {
          accept_tets_[other->Id()] = true;
          ActivateSurroundingTets( other );
        }
      }
      else if ( tets.size()>2 )
      {
        throw std::runtime_error( "too many tets at tri surface" );
      }
    }
  }
}

void GEO::CUT::TetMesh::ActivateUsedCutSurfaceTris()
{
  for ( std::vector<int>::iterator i = std::find( accept_tets_.begin(), accept_tets_.end(), 1 );
        i != accept_tets_.end();
        i = std::find( ++i, accept_tets_.end(), 1 ) )
  {
    unsigned pos = i - accept_tets_.begin();
    Entity<4> & e = tet_entities_[pos];
    const std::vector<Entity<3>*> & tris = e.Children();
    for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;
      accept_tris_[tri->Id()] = true;
    }
  }
}

/// test for illegal (cutted) cut surfaces
void GEO::CUT::TetMesh::TestCutSurface()
{
  std::set<Entity<2>*> inner_lines;

  for ( std::map<Facet*, std::vector<Entity<3>*> >::iterator i=cutsurface_.begin();
        i!=cutsurface_.end();
        ++i )
  {
    Facet * f = i->first;
    std::vector<Entity<3>*> & tris = i->second;

    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;

      // test for gaps that show illegal tets
      //
      // Here we recognize triangulation errors that need to be fixed
      // manually. There are no fixes yet.

      const std::vector<Entity<2>*> & lines = tri->Children();
      for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        Entity<2> & line = **i;
        int count = line.CountMarkedParents( active_surface_tris_, accept_tris_ );
        bool isline = f->IsLine( points_[line[0]], points_[line[1]] );
#if 1
        if ( count > 2 )
        {
          throw std::runtime_error( "fork in the road" );
        }
#endif
        if ( count == 1 and not isline )
        {
          //throw std::runtime_error( "gap in cut facet" );
          inner_lines.insert( &line );
        }
#if 0
        if ( count == 2 and isline )
        {
          throw std::runtime_error( "two surface elements at line" );
        }
#endif
      }
    }
  }

  if ( inner_lines.size()>0 )
  {
    // So we have lines within facets that are not connected to two
    // tris. Select all tris these lines are connected with, provided
    //
    // - the tri is a cut surface tri
    // - the tri is (not yet) activated at the surface
    // - the tri has at least one parent with a non cut surface point

    std::set<Entity<3>*> inner_tris;
    for ( std::set<Entity<2>*>::iterator i=inner_lines.begin(); i!=inner_lines.end(); ++i )
    {
      Entity<2> * line = *i;
      const std::vector<Entity<3>*> & tris = line->Parents();
      for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
      {
        Entity<3> * tri = *i;
        if ( not active_surface_tris_[tri->Id()] and tri->AllOnCutSurface( points_ ) )
        {
          const std::vector<Entity<4>*> & inner_tets = tri->Parents();
          for ( std::vector<Entity<4>*>::const_iterator i=inner_tets.begin();
                i!=inner_tets.end();
                ++i )
          {
            Entity<4> * tet = *i;
            Point * p = points_[tet->OtherPoint( tri )];
            if ( p->Position()!=Point::oncutsurface )
            {
              inner_tris.insert( tri );
              break;
            }
          }
        }
      }
    }

    for ( std::set<Entity<3>*>::iterator i=inner_tris.begin(); i!=inner_tris.end(); ++i )
    {
      Entity<3> * tri = *i;

      // activate all tris found in that way
      active_surface_tris_[tri->Id()] = true;

#if 1
      std::cout << "WARNING: Activate non-facet tri " << tri->Id() << "\n";
#endif

      const std::vector<Entity<2>*> & lines = tri->Children();
      for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        Entity<2> * line = *i;
        inner_lines.erase( line );
      }
    }

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
    GmshWriteBrokenSurfaceTris( inner_tris );
#endif

    if ( inner_lines.size()>0 )
    {
      throw std::runtime_error( "remaining gap in cut facet" );
    }
  }
}

/// remove deactivated tets
void GEO::CUT::TetMesh::ClearExternalTets()
{
  for ( std::vector<int>::iterator i = std::find( accept_tets_.begin(), accept_tets_.end(), 0 );
        i != accept_tets_.end();
        i = std::find( ++i, accept_tets_.end(), 0 ) )
  {
    unsigned pos = i - accept_tets_.begin();
    tets_[pos].clear();
  }
}

/// Get tri coordinates at the cut facets in the proper order
void GEO::CUT::TetMesh::FillCutSides( std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  for ( std::map<Facet*, std::vector<Entity<3>*> >::iterator i=cutsurface_.begin();
        i!=cutsurface_.end();
        ++i )
  {
    Facet * f = i->first;
    std::vector<Entity<3>*> & tris = i->second;

    std::vector<Epetra_SerialDenseMatrix> & side_coords = sides_xyz[f];

    for ( std::vector<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;

      if ( active_surface_tris_[tri->Id()] and
           accept_tris_[tri->Id()] )
      {
        bool done = false;
        const std::vector<Entity<4>*> & tets = tri->Parents();
        for ( std::vector<Entity<4>*>::const_iterator i=tets.begin();
              i!=tets.end();
              ++i )
        {
          Entity<4> * tet = *i;
          std::vector<int> & original_tet = tets_[tet->Id()];

          //if ( accept_tets_[tet->Id()] )
          if ( original_tet.size() > 0 )
          {
            if ( done )
            {
              throw std::runtime_error( "double tets at cut surface" );
            }

            done = true;

            bool found = false;
            for ( int i=0; i<4; ++i )
            {
              std::vector<int> side( 3 );
              for ( int j=0; j<3; ++j )
              {
                side[j] = original_tet[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
              }
              if ( tri->Equals( side ) )
              {
                found = true;

                Epetra_SerialDenseMatrix xyz( 3, 3 );
                points_[side[0]]->Coordinates( &xyz( 0, 0 ) );
                points_[side[1]]->Coordinates( &xyz( 0, 1 ) );
                points_[side[2]]->Coordinates( &xyz( 0, 2 ) );
                side_coords.push_back( xyz );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
                surface_tris_.push_back( side );
#endif
                break;
              }
            }
            if ( not found )
            {
              throw std::runtime_error( "failed to find side" );
            }
          }
        }
        if ( not done )
        {
          throw std::runtime_error( "failed to find tet" );
        }
      }
    }
  }
}

#ifdef TETMESH_GMSH_DEBUG_OUTPUT

void GEO::CUT::TetMesh::GmshWriteCells()
{
  std::ofstream file( "delaunaycells.pos" );
  file << "View \"delaunaycells\" {\n";
  for ( std::vector<std::vector<int> >::iterator i=tets_.begin(); i<tets_.end(); ++i )
  {
    std::vector<int> & t = *i;
    GmshWriteTet( file, i-tets_.begin(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteActiveCells()
{
  std::ofstream file( "activecells.pos" );
  file << "View \"activecells\" {\n";
  for ( std::vector<std::vector<int> >::iterator i=tets_.begin(); i<tets_.end(); ++i )
  {
    unsigned pos = i - tets_.begin();
    if ( accept_tets_[pos] )
    {
      std::vector<int> & t = *i;
      GmshWriteTet( file, i-tets_.begin(), t );
    }
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteBrokenSurfaceTris( const std::set<Entity<3>*> & tris )
{
  std::ofstream file( "brokentris.pos" );
  file << "View \"brokentris\" {\n";
  for ( std::set<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> & tri = **i;
    std::vector<int> t( tri(), tri()+3 );
    GmshWriteTri( file, tri.Id(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteFacetTris()
{
  std::map<Facet*, int> my_facets;
  for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    my_facets[*i] = my_facets.size();
  }

  std::ofstream file( "facettris.pos" );
  file << "View \"facettris\" {\n";
  for ( std::map<Entity<3>*, Facet*>::iterator i=tri_facets_.begin(); i!=tri_facets_.end(); ++i )
  {
    Entity<3> & tri = *i->first;
    Facet * f = i->second;
    std::vector<int> t( tri(), tri()+3 );
    GmshWriteConnect( file, "ST", t );
    for ( std::vector<int>::const_iterator j=t.begin(); j!=t.end(); ++j )
    {
      //Point * p = points_[*j];
      if ( j != t.begin() )
        file << ",";
      file << my_facets[f];
    }
    file << "};  // " << tri.Id() << "\n";
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteSurfaceCells()
{
  std::ofstream file( "surfacecells.pos" );
  file << "View \"surfacecells\" {\n";
  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & tri = i->second;
    std::vector<int> t( tri(), tri()+3 );
    GmshWriteTri( file, tri.Id(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteActiveSurfaceCells()
{
  std::ofstream file( "activesurfacecells.pos" );
  file << "View \"activesurfacecells\" {\n";
  for ( std::map<Handle<3>, Entity<3> >::iterator i=tet_surfaces_.begin();
        i!=tet_surfaces_.end();
        ++i )
  {
    Entity<3> & tri = i->second;
    if ( active_surface_tris_[tri.Id()] and
         accept_tris_[tri.Id()] )
    {
      std::vector<int> t( tri(), tri()+3 );
      GmshWriteTri( file, tri.Id(), t );
    }
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteSurfaceTris()
{
  std::ofstream file( "surfacetris.pos" );
  file << "View \"surfacetris\" {\n";
  for ( std::vector<std::vector<int> >::iterator i=surface_tris_.begin(); i<surface_tris_.end(); ++i )
  {
    std::vector<int> & t = *i;
    GmshWriteTri( file, i-surface_tris_.begin(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTri( std::ostream & file, int eid, const std::vector<int> & t )
{
  GmshWriteConnect( file, "ST", t );
  GmshWritePosition( file, eid, t );
}

void GEO::CUT::TetMesh::GmshWriteTet( std::ostream & file, int eid, const std::vector<int> & t )
{
  GmshWriteConnect( file, "SS", t );
  GmshWritePosition( file, eid, t );
}

void GEO::CUT::TetMesh::GmshWriteConnect( std::ostream & file, std::string name, const std::vector<int> & t )
{
  file << name << "(";
  for ( std::vector<int>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = points_[*j];
    if ( j != t.begin() )
      file << ",";
    file << p->X()[0] << ","
         << p->X()[1] << ","
         << p->X()[2];
  }
  file << "){";
}

void GEO::CUT::TetMesh::GmshWritePosition( std::ostream & file, int eid, const std::vector<int> & t )
{
  for ( std::vector<int>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = points_[*j];
    if ( j != t.begin() )
      file << ",";
    file << p->Position();
  }
  file << "};  // " << eid << "\n";
}

#endif
