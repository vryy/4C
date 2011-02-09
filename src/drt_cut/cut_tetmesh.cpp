
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

  CallQHull( points_, original_tets );

  tets_.reserve( original_tets.size() );
  for ( std::vector<std::vector<int> >::iterator i=original_tets.begin(); i!=original_tets.end(); ++i )
  {
    std::vector<int> & t = *i;
    std::vector<Point*> tet;
    tet.reserve( t.size() );
    for ( std::vector<int>::iterator i=t.begin(); i!=t.end(); ++i )
    {
      tet.push_back( points_[*i] );
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
  FindOverlappingTriFacets();
  RemoveSimpleOverlappingTriFacets();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteSurfaceCells();
#endif

  ActivateCutSurfaceTris();

  // find the active tets

  DeactivateCutSurfaceTets();
  ActivateSingleTetsOnFacets();
  ActivateDoubleTetsOnCutSurface();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteActiveCells();
  //GmshWriteActiveSurfaceCells();
#endif

  RemoveFullOverlappingTriFacets();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  //GmshWriteActiveCells();
  GmshWriteActiveSurfaceCells();
#endif

  TestCutSurface();

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
  options.push_back( "qhull d Qt Qbb Qc QJ Pp" );

  // If you want some debugging information replace the 0 pointer
  // with stdout or some other file open for writing.

  FILE * outfile = 0;
//#define QHULL_DEBUG_OUTPUT
#ifdef QHULL_DEBUG_OUTPUT
#if 0
  static FILE * errfile;
  if ( errfile==NULL )
    errfile = fopen( "qhull_error.log", "w" );
#else
  FILE * errfile = stderr;
#endif
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
          if ( facet->vertices )
          {
            std::vector<int> ids;
            ids.reserve( dim+1 );
            vertexT *  vertex = NULL;

            //FOREACHvertex_(facet->vertices)
            for ( void ** vertexp = & facet->vertices->e[0].p;
                  ( vertex = static_cast<vertexT*>( * vertexp++ ) );
              )
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

#ifdef QHULL_DEBUG_OUTPUT
  fflush( errfile );
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

    // remove all facets that do not belong to the volume cell
    std::set<Facet*> intersection;
    std::set_intersection( facets.begin(), facets.end(),
                           facets_.begin(), facets_.end(),
                           std::inserter( intersection, intersection.begin() ) );
    std::swap( facets, intersection );

    if ( facets.size() == 1 )
    {
      Facet * f = *facets.begin();
      tri_facets_[&tri] = f;
      facet_info_[f].tris_.insert( &tri );
    }
    else if ( facets.size() > 1 )
    {
      std::swap( double_tris[tri.GetHandle()], facets );
    }
  }

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteFacetTris();
#endif

  if ( double_tris.size() == 0 )
    return;

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
          facet_info_[f].tris_.insert( tri );
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

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteFacetTris();
#endif

  if ( double_tris.size() > 0 )
    throw std::runtime_error( "failed to assign triangle facets" );
}

void GEO::CUT::TetMesh::FindOverlappingTriFacets()
{
  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    Facet * f = i->first;
    FacetInfo & fi = i->second;
    std::set<Entity<3>*> & tris = fi.tris_;

    // find all lines surrounding the facet

    std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
    f->GetLines( lines );

    std::set<Entity<2>*> trace_lines;
    for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=lines.begin();
          i!=lines.end();
          ++i )
    {
      const std::pair<Point*, Point*> & l = i->first;
      Handle<2> h = MakeHandle( l.first, l.second );
      std::map<Handle<2>, Entity<2> >::iterator j = tet_lines_.find( h );
      if ( j!=tet_lines_.end() )
      {
        Entity<2> & line = j->second;
        trace_lines.insert( &line );
      }
    }

    // find number of tris at all lines

    std::map<Entity<2>*, std::set<Entity<3>*> > line_usage;
    for ( std::set<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Entity<3> * tri = *i;
      const std::vector<Entity<2>*> & lines = tri->Children();
      for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        Entity<2> * line = *i;
        line_usage[line].insert( tri );
      }
    }

    std::list<OverlappingTriSets> overlappingsets;

    // Find all possible forks. This creates more ots than needed in
    // complicated cases. But at least there are all cases covered. The merge
    // is done later on.

    for ( std::map<Entity<2>*, std::set<Entity<3>*> >::iterator i=line_usage.begin();
          i!=line_usage.end();
          ++i )
    {
      Entity<2> * line = i->first;
      std::set<Entity<3>*> & line_tris = i->second;
      if ( line_tris.size() > 2 )
      {
        if ( trace_lines.count( line )==0 )
        {
          // find all combinations at this line
          for ( std::set<Entity<3>*>::iterator i=line_tris.begin(); i!=line_tris.end(); ++i )
          {
            Entity<3> * tri1 = *i;
            std::set<Entity<3>*>::iterator j = i;
            for ( ++j; j!=line_tris.end(); ++j )
            {
              Entity<3> * tri2 = *j;
              overlappingsets.push_back( OverlappingTriSets( trace_lines, line, tri1, tri2 ) );
            }
          }
        }
        else
        {
          // this is a trace line, thus there is just one tri
          for ( std::set<Entity<3>*>::iterator i=line_tris.begin(); i!=line_tris.end(); ++i )
          {
            Entity<3> * tri = *i;
            overlappingsets.push_back( OverlappingTriSets( trace_lines, line, tri ) );
          }
        }
      }
    }

    // There are no illegal lines, just use the normal trace.
    if ( overlappingsets.size()==0 )
    {
      overlappingsets.push_back( OverlappingTriSets( trace_lines ) );

      // There are no options. The order of lines in input does not matter in
      // any way. Ignore it.
      line_usage.clear();
    }

    for ( std::list<OverlappingTriSets>::iterator i=overlappingsets.begin();
          i!=overlappingsets.end();
          ++i )
    {
      OverlappingTriSets & ots = *i;

      while ( not ots.Done() )
      {
        Entity<2> * line = ots.NextLine( line_usage );

        std::set<Entity<3>*> new_tris;
        ots.NewLineTris( tris, line, new_tris );
        if ( new_tris.size()==0 )
        {
          // Empty line.
          ots.EmptyTriSet( line );
        }
        else
        {
          // Simple new triangle or multiple choices.

          // make a new ots for each triangle that is available at the line
          // at hand
          std::map<Entity<3>*, OverlappingTriSets*> new_ots;
          std::set<Entity<3>*>::iterator i = new_tris.begin();
          new_ots[*i] = &ots;
          for ( ++i; i!=new_tris.end(); ++i )
          {
            overlappingsets.push_back( ots );
            new_ots[*i] = &overlappingsets.back();
          }
          for ( std::set<Entity<3>*>::iterator i=new_tris.begin(); i!=new_tris.end(); ++i )
          {
            Entity<3> * tri = *i;
            if ( not new_ots[tri]->InsertTri( tri ) )
            {
              new_ots[tri]->EmptyTriSet( line );
            }
          }
        }
      }
    }

    // remove duplicates
    fi.SetOverlappingTriSets( overlappingsets );
  }
}

void GEO::CUT::TetMesh::RemoveSimpleOverlappingTriFacets()
{
//     unsigned matches = std::count_if( fullset.begin(), fullset.end(),
//                                       std::bind1st( std::equal_to<int>(), true ) );

  std::vector<FacetInfo*> reconstruct_sets;

  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    Facet * f = i->first;
    FacetInfo & fi = i->second;
    std::set<Entity<3>*> & tris = fi.tris_;
    std::vector<OverlappingTriSets> & overlappingsets = fi.overlappingsets_;
    std::vector<OverlappingTriSets*> matches;
    fi.FindFullOTS( matches );
    if ( matches.size() == 1 )
    {
      OverlappingTriSets * ots = matches[0];
      std::swap( tris, ots->tset_ );
      fi.done_ = true;
    }
    else if ( matches.size() == 0 )
    {
      std::cout << "WARNING: " << overlappingsets.size()
                << " non-full matches of facet " << f->SideId() << ":\n";
      for ( std::vector<OverlappingTriSets>::iterator i=overlappingsets.begin();
            i!=overlappingsets.end();
            ++i )
      {
        OverlappingTriSets & ots = *i;
        std::cout << "\t(";
        for ( std::set<Entity<3>*>::iterator i=ots.tset_.begin(); i!=ots.tset_.end(); ++i )
        {
          Entity<3> & tri = **i;
          std::cout << tri;
        }
        std::cout << ")\n";
      }

      //fi.done_ = true;
      reconstruct_sets.push_back( &fi );
    }
    else if ( matches.size() > 1 )
    {
      if ( f->IsTriangulated() )
      {
        std::vector<Point*> points;
        for ( std::vector<OverlappingTriSets*>::iterator i=matches.begin(); i!=matches.end(); ++i )
        {
          OverlappingTriSets * ots = *i;
          bool all_triangulated = true;
          for ( std::set<Entity<3>*>::iterator i=ots->tset_.begin(); i!=ots->tset_.end(); ++i )
          {
            Entity<3> * tri = *i;
            tri->PointVector( points_, points );
            if ( not f->IsTriangulatedSide( points ) )
            {
              all_triangulated = false;
              break;
            }
          }
          if ( all_triangulated )
          {
            std::swap( tris, ots->tset_ );
            fi.done_ = true;
            break;
          }
        }
      }
    }
  }

  if ( reconstruct_sets.size() > 0 )
  {
    std::vector<std::pair<std::vector<OverlappingTriSets>::iterator,
      std::vector<OverlappingTriSets>::iterator> > reconstruct_range;
    std::vector<std::vector<OverlappingTriSets>::iterator> reconstruct_counter;
    reconstruct_range  .reserve( reconstruct_sets.size() );
    reconstruct_counter.reserve( reconstruct_sets.size() );
    for ( std::vector<FacetInfo*>::iterator i=reconstruct_sets.begin();
          i!=reconstruct_sets.end();
          ++i )
    {
      FacetInfo * fi = *i;
      std::vector<OverlappingTriSets> & otsl = fi->overlappingsets_;
      reconstruct_range  .push_back( std::make_pair( otsl.begin(), otsl.end() ) );
      reconstruct_counter.push_back( otsl.begin() );
    }

    std::vector<std::set<Entity<3>*> > surface_tris;
    std::vector<std::vector<std::vector<OverlappingTriSets>::iterator> > counter;
    while ( reconstruct_counter.back() != reconstruct_range.back().second )
    {
      std::set<Entity<2>*> open_surface_lines;

      for ( std::vector<std::vector<OverlappingTriSets>::iterator>::iterator i=reconstruct_counter.begin();
            i!=reconstruct_counter.end();
            ++i )
      {
        OverlappingTriSets & ots = **i;
        std::copy( ots.openlines_.begin(), ots.openlines_.end(),
                   std::inserter( open_surface_lines, open_surface_lines.begin() ) );
      }

      std::set<Entity<3>*> st;
      if ( FillSurfaceGaps( open_surface_lines, st ) )
      {
        counter.push_back( reconstruct_counter );
        surface_tris.push_back( st );
      }

      // next combination
      for ( unsigned i=0; i<reconstruct_sets.size(); ++i )
      {
        ++reconstruct_counter[i];
        if ( reconstruct_counter[i] == reconstruct_range[i].second )
        {
          if ( i+1 < reconstruct_sets.size() )
          {
            reconstruct_counter[i] = reconstruct_range[i].first;
          }
        }
        else
        {
          break;
        }
      }
    }

    if ( counter.size()==0 )
    {
      throw std::runtime_error( "not able to fill surface gaps" );
    }

    // We are undecided which ots to use. Postphone it. Treat them like the
    // full non-unique ots.

    for ( unsigned i=0; i<counter.size(); ++i )
    {
      std::set<Entity<3>*> & st = surface_tris[i];
      std::vector<std::vector<OverlappingTriSets>::iterator> & c = counter[i];

      // collect and activate broken tris
      std::vector<Entity<3>*> local_broken_tris;
      local_broken_tris.reserve( st.size() );
      for ( std::set<Entity<3>*>::iterator i=st.begin(); i!=st.end(); ++i )
      {
        Entity<3> * tri = *i;
        if ( tri->AllOnCutSurface( points_ ) )
        {
          //std::cout << "WARNING: Activate non-facet tri " << tri->Id() << "\n";
          active_surface_tris_[tri->Id()] = true;
          local_broken_tris.push_back( tri );
        }
        else
        {
          if ( tri->HasOnCutSurface( points_ ) )
          {
            throw std::runtime_error( "broken tri partly on cut surface -- rotation needed" );
          }
        }
      }

      if ( counter.size() == 1 )
      {
        // there is a unique gap fill-in in this facet

        // take the combination that worked for all gaps
        for ( unsigned i=0; i!=reconstruct_sets.size(); ++i )
        {
          FacetInfo * fi = reconstruct_sets[i];

          OverlappingTriSets & ots = *c[i];
          fi->ExchangeTriSet( ots.tset_, active_surface_tris_ );
          fi->done_ = true;
        }

        // fill all gaps
        std::copy( local_broken_tris.begin(), local_broken_tris.end(),
                   std::inserter( broken_tris_, broken_tris_.begin() ) );
      }
      else
      {
        // there are more options

        std::set<OverlappingTriSets*> select;

        // set marker to process later
        for ( unsigned i=0; i!=reconstruct_sets.size(); ++i )
        {
          FacetInfo * fi = reconstruct_sets[i];
          OverlappingTriSets & ots = *c[i];

          if ( fi->IsUnique() )
          {
            fi->ExchangeTriSet( ots.tset_, active_surface_tris_ );
            fi->done_ = true;
          }
          else
          {
            fi->done_ = false;

            // mark the ots that it needs to be looked at again
            ots.full_ = true;

            // mark all tris as (possible) cut surface tris
            fi->SetMarker( active_surface_tris_, true );

            select.insert( &ots );
          }
        }

        // remember the set of broken tris that applies to this particular set
        // of overlapping tris
        if ( select.size() > 0 )
        {
          std::swap( broken_tri_selection_[select], local_broken_tris );
        }
        else
        {
          // fill all gaps
          std::copy( local_broken_tris.begin(), local_broken_tris.end(),
                     std::inserter( broken_tris_, broken_tris_.begin() ) );
        }
      }
    }

//     for ( std::vector<FacetInfo*>::iterator i=reconstruct_sets.begin();
//           i!=reconstruct_sets.end();
//           ++i )
//     {
//       FacetInfo * fi = *i;
//       std::vector<OverlappingTriSets*> matches;
//       fi->FindFullOTS( matches );
//       if ( matches.size()==1 )
//       {
//         OverlappingTriSets * ots = matches[0];
//         fi->ExchangeTriSet( ots->tset_, active_surface_tris_ );
//         fi->done_ = true;
//       }
//     }
  }

#ifdef DEBUGCUTLIBRARY
  GmshWriteTriSet( "brokentris1", broken_tris_ );
#endif
}

void GEO::CUT::TetMesh::RemoveFullOverlappingTriFacets()
{
  // At this point the tets are assigned and can be used to test which version
  // is to be used.

  std::set<OverlappingTriSets*> select;

  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    //Facet * f = i->first;
    FacetInfo & fi = i->second;

    if ( fi.done_ )
      continue;

    //std::set<Entity<3>*> & tris = fi.tris_;
    //std::vector<OverlappingTriSets> & overlappingsets = fi.overlappingsets_;
    std::vector<OverlappingTriSets*> matches;
    fi.FindFullOTS( matches );
//     if ( matches.size() > 1 )
    {
#ifdef DEBUGCUTLIBRARY
      for ( unsigned i=0; i!=matches.size(); ++i )
      {
        OverlappingTriSets * ots = matches[i];
        std::stringstream str;
        str << "ots" << i;
        GmshWriteTriSet( str.str(), ots->tset_ );
      }
#endif
      std::vector<OverlappingTriSets*> cutsurface_ots;
      for ( std::vector<OverlappingTriSets*>::iterator i=matches.begin(); i!=matches.end(); ++i )
      {
        OverlappingTriSets * ots = *i;
        bool all_accepted = true;
        for ( std::set<Entity<3>*>::iterator i=ots->tset_.begin(); i!=ots->tset_.end(); ++i )
        {
          Entity<3> * tri = *i;
          if ( not HasAcceptedTet( *tri ) )
          {
            all_accepted = false;
            break;
          }
        }
        if ( all_accepted )
        {
          cutsurface_ots.push_back( ots );
        }
      }
      if ( cutsurface_ots.size()==0 )
      {
        throw std::runtime_error( "no set with active tet parents" );
      }
      if ( cutsurface_ots.size() > 1 )
      {
        throw std::runtime_error( "ambiguous sets" );
      }

      OverlappingTriSets * ots = cutsurface_ots[0];
      fi.ExchangeTriSet( ots->tset_, active_surface_tris_ );

      select.insert( ots );
    }
  }

  //bool found_state = false;
  for ( std::map<std::set<OverlappingTriSets*>, std::vector<Entity<3>*> >::iterator ibts=broken_tri_selection_.begin();
        ibts!=broken_tri_selection_.end();
        ++ibts )
  {
    const std::set<OverlappingTriSets*> & select_state = ibts->first;
    std::vector<Entity<3>*> & broken = ibts->second;

    bool found = true;
    for ( std::set<OverlappingTriSets*>::const_iterator i=select_state.begin();
          i!=select_state.end();
          ++i )
    {
      OverlappingTriSets * ots = *i;
      if ( select.count( ots )==0 )
      {
        found = false;
        break;
      }
    }
    if ( found )
    {
      std::copy( broken.begin(), broken.end(),
                 std::inserter( broken_tris_, broken_tris_.begin() ) );

      // remove used broken tris
      broken_tri_selection_.erase( ibts );

      //found_state = true;
      break;
    }
  }
//   if ( not found_state )
//   {
//     throw std::runtime_error( "" );
//   }

  // deactivate unused broken tris
  for ( std::map<std::set<OverlappingTriSets*>, std::vector<Entity<3>*> >::iterator i=broken_tri_selection_.begin();
        i!=broken_tri_selection_.end();
        ++i )
  {
    std::vector<Entity<3>*> & broken = i->second;
    for ( std::vector<Entity<3>*>::iterator i=broken.begin(); i!=broken.end(); ++i )
    {
      Entity<3> * tri = *i;
      active_surface_tris_[tri->Id()] = false;
    }
  }

  // Activate used broken tris again, in case there is some overlap with the
  // unused ones. (There is!)
  for ( std::set<Entity<3>*>::iterator i=broken_tris_.begin(); i!=broken_tris_.end(); ++i )
  {
    Entity<3> * tri = *i;
    active_surface_tris_[tri->Id()] = true;
  }

// #ifdef DEBUGCUTLIBRARY
//   {
//     int count = 0;
//     for ( std::set<OverlappingTriSets*>::iterator i=select.begin();
//           i!=select.end();
//           ++i )
//     {
//       OverlappingTriSets * ots = *i;
//       std::stringstream str;
//       str << "select" << count++;
//       GmshWriteTriSet( str.str(), ots->tset_ );
//     }
//   }
// #endif

#ifdef DEBUGCUTLIBRARY
  GmshWriteTriSet( "brokentris2", broken_tris_ );
#endif
}

bool GEO::CUT::TetMesh::FillSurfaceGaps( std::set<Entity<2>*> & surface_lines,
                                         std::set<Entity<3>*> & surface_tris )
{
  while ( surface_lines.size() > 0 )
  {
    bool line_found = false;
    unsigned numtris = surface_tris.size();

    std::map<Handle<1>, Entity<1> > tet_points;

    for ( std::set<Entity<2>*>::iterator i=surface_lines.begin(); i!=surface_lines.end(); ++i )
    {
      Entity<2> * line = *i;
      line->CreateChildren( tet_points );
    }

    for ( std::map<Handle<1>, Entity<1> >::iterator i=tet_points.begin(); i!=tet_points.end(); ++i )
    {
      Entity<1> & point = i->second;
      const std::vector<Entity<2>*> & lines = point.Parents();
      if ( lines.size()==2 )
      {
        Entity<2> * l1 = lines[0];
        Entity<2> * l2 = lines[1];

        std::vector<Entity<3>*> tri1 = l1->Parents();
        std::vector<Entity<3>*> tri2 = l2->Parents();
        std::sort( tri1.begin(), tri1.end() );
        std::sort( tri2.begin(), tri2.end() );
        std::vector<Entity<3>*> intersection;
        std::set_intersection( tri1.begin(), tri1.end(),
                               tri2.begin(), tri2.end(),
                               std::back_inserter( intersection ) );
        if ( intersection.size()==1 )
        {
          Entity<3> * tri = intersection[0];

          if ( //tri_facets_.count( tri )==0 and
               surface_tris.count( tri )==0 )
          {
            surface_tris.insert( tri );

            const std::vector<Entity<2>*> & lines = tri->Children();
            for ( std::vector<Entity<2>*>::const_iterator i=lines.begin();
                  i!=lines.end();
                  ++i )
            {
              Entity<2> * line = *i;
              if ( surface_lines.count( line ) > 0 )
              {
                surface_lines.erase( line );
              }
              else
              {
                surface_lines.insert( line );
              }
            }
            line_found = true;
          }
        }
        else if ( intersection.size() > 1 )
        {
          throw std::runtime_error( "cannot have more than one common tri between lines" );
        }
        else
        {
          // there are more that one tris at this point
        }
        if ( line_found )
        {
          // after each change we need to rebuild the tet_points connection
          break;
        }
      }
    }

    for ( std::map<Handle<1>, Entity<1> >::iterator i=tet_points.begin(); i!=tet_points.end(); ++i )
    {
      Entity<1> & point = i->second;
      point.Disconnect();
    }

    if ( not line_found or numtris == surface_tris.size() )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::TetMesh::ActivateCutSurfaceTris()
{
  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    Facet * f = i->first;
    FacetInfo & fi = i->second;
    if ( f->SideId() > -1 )
    {
      std::set<Entity<3>*> & tris = fi.tris_;
      for ( std::set<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
      {
        Entity<3> * tri = *i;

        active_surface_tris_[tri->Id()] = true;
      }
    }
  }
}

/// test for illegal (cutted) cut surfaces
void GEO::CUT::TetMesh::TestCutSurface()
{
  std::set<Entity<2>*> inner_lines;

  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    Facet * f = i->first;
    FacetInfo & fi = i->second;
    if ( f->SideId() > -1 )
    {
      std::set<Entity<3>*> & tris = fi.tris_;

      for ( std::set<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
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
          int count = line.CountMarkedParents( active_surface_tris_ );
          bool isline = f->IsLine( points_[line[0]], points_[line[1]] );
#if 1
          if ( count > 2 )
          {
            throw std::runtime_error( "fork in the road" );
          }
#endif
          if ( count == 1 and not isline )
          {
            throw std::runtime_error( "gap in cut facet" );
            //inner_lines.insert( &line );
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
      broken_tris_.insert( tri );

#if 1
      std::cout << "WARNING: Last minute activation of non-facet tri " << tri->Id() << "\n";
#endif

      const std::vector<Entity<2>*> & lines = tri->Children();
      for ( std::vector<Entity<2>*>::const_iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        Entity<2> * line = *i;
        inner_lines.erase( line );
      }
    }
#ifdef TETMESH_GMSH_DEBUG_OUTPUT
    GmshWriteTriSet( "brokentris3", broken_tris_ );
#endif

    if ( inner_lines.size()>0 )
    {
      throw std::runtime_error( "remaining gap in cut facet" );
    }
  }
}

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

/// if there is just one tet at an cut surface it is active
void GEO::CUT::TetMesh::ActivateSingleTetsOnFacets()
{
  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    //Facet * f = i->first;
    FacetInfo & fi = i->second;
    if ( fi.IsUnique() )
    {
      std::set<Entity<3>*> & tris = fi.tris_;
      for ( std::set<Entity<3>*>::iterator i=tris.begin(); i!=tris.end(); ++i )
      {
        Entity<3> & tri = **i;
        const std::vector<Entity<4>*> & tets = tri.Parents();
        if ( tets.size()==1 )
        {
          accept_tets_[tets[0]->Id()] = true;
        }
      }
    }
  }
}

/// activate the right tets at the cut surface
void GEO::CUT::TetMesh::ActivateDoubleTetsOnCutSurface()
{
  if ( std::find( accept_tets_.begin(), accept_tets_.end(), true )==accept_tets_.end() )
  {
    throw std::runtime_error( "no accepted tets to begin with" );
  }
  for ( std::vector<int>::iterator i = std::find( accept_tets_.begin(), accept_tets_.end(), true );
        i != accept_tets_.end();
        i = std::find( ++i, accept_tets_.end(), true ) )
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
    if ( not active_surface_tris_[tri->Id()] )
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

/// remove deactivated tets
void GEO::CUT::TetMesh::ClearExternalTets()
{
  for ( std::vector<int>::iterator i = std::find( accept_tets_.begin(), accept_tets_.end(), false );
        i != accept_tets_.end();
        i = std::find( ++i, accept_tets_.end(), false ) )
  {
    unsigned pos = i - accept_tets_.begin();
    tets_[pos].clear();
  }
}

/// Get tri coordinates at the cut facets in the proper order
void GEO::CUT::TetMesh::FillCutSides( std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  // rotate broken tris and find unique facet

  while ( broken_tris_.size() > 0 )
  {
    std::set<Entity<3>*> patch;
    patch.insert( *broken_tris_.begin() );
    broken_tris_.erase( broken_tris_.begin() );
    std::set<Entity<2>*> common_lines;
    std::set<Entity<3>*>::iterator ib = broken_tris_.begin();
    while ( ib!=broken_tris_.end() )
    {
      Entity<3> * tri = *ib;
      bool found = false;
      for ( std::set<Entity<3>*>::iterator i=patch.begin(); i!=patch.end(); ++i )
      {
        Entity<3> * patch_tri = *i;
        Entity<2> * line = tri->CommonChild( patch_tri );
        if ( line != NULL )
        {
          common_lines.insert( line );
          patch.insert( tri );
          broken_tris_.erase( tri );
          ib = broken_tris_.begin();
          found = true;
          break;
        }
      }
      if ( not found )
      {
        ++ib;
      }
    }

    std::set<int> all_points;
    std::set<int> surface_points;
    std::set<Entity<2>*> surface_lines;
    for ( std::set<Entity<3>*>::iterator i=patch.begin(); i!=patch.end(); ++i )
    {
      Entity<3> & tri = **i;
      const std::vector<Entity<2>*> lines = tri.Children();
      std::copy( lines.begin(), lines.end(), std::inserter( surface_lines, surface_lines.begin() ) );
      std::copy( tri(), tri()+3, std::inserter( all_points, all_points.begin() ) );
    }
    for ( std::set<Entity<2>*>::iterator i=common_lines.begin(); i!=common_lines.end(); ++i )
    {
      Entity<2> * line = *i;
      surface_lines.erase( line );
    }
    for ( std::set<Entity<2>*>::iterator i=surface_lines.begin(); i!=surface_lines.end(); ++i )
    {
      Entity<2> & line = **i;
      std::copy( line(), line()+2, std::inserter( surface_points, surface_points.begin() ) );
    }

    if ( all_points.size()==surface_points.size() )
    {
      std::vector<int> points;
      points.reserve( surface_points.size() );
      Entity<2> * line = *surface_lines.begin();
      points.push_back( ( *line )[0] );
      points.push_back( ( *line )[1] );
      surface_lines.erase( line );
      while ( surface_lines.size() > 0 )
      {
        unsigned numlines = surface_lines.size();
        for ( std::set<Entity<2>*>::iterator i=surface_lines.begin(); i!=surface_lines.end(); ++i )
        {
          Entity<2> & line = **i;
          if ( line[0]==points.back() )
          {
            points.push_back( line[1] );
            surface_lines.erase( i );
            break;
          }
          else if ( line[1]==points.back() )
          {
            points.push_back( line[0] );
            surface_lines.erase( i );
            break;
          }
        }
        if ( numlines == surface_lines.size() )
        {
          throw std::runtime_error( "unsupported patch topology" );
        }
      }

      if ( points.front()!=points.back() )
      {
        throw std::runtime_error( "patch not closed" );
      }
      points.pop_back();

      std::vector<std::vector<int> > sides;
      FindProperSides( patch, sides );

      if ( sides.size()!=patch.size() )
        throw std::runtime_error( "wrong number of sides" );

      std::map<int, int> index;
      unsigned numpoints = points.size();
      for ( unsigned i=0; i<numpoints; ++i )
      {
        index[points[i]] = i;
      }

      bool match = true;
      for ( unsigned step=0; step<patch.size()-1; ++step )
      {
        // rotate
        for ( std::vector<std::vector<int> >::iterator i=sides.begin();
              i!=sides.end();
              ++i )
        {
          std::vector<int> & side = *i;
          for ( std::vector<int>::iterator i=side.begin(); i!=side.end(); ++i )
          {
            int & p = *i;
            p = points[( index[p]+1 ) % numpoints];
          }
        }

        match = true;
        std::map<Facet*, std::vector<std::vector<int> > > patch_facets;
        std::vector<Point*> tri_points( 3 );
        for ( std::vector<std::vector<int> >::iterator i=sides.begin();
              i!=sides.end();
              ++i )
        {
          std::vector<int> & side = *i;
          for ( int i=0; i<3; ++i )
          {
            tri_points[i] = points_[side[i]];
          }
          std::set<Facet*> facets;
          FindCommonFacets( tri_points, facets );
          if ( facets.size()==1 )
          {
            Facet * f = *facets.begin();
            patch_facets[ f ].push_back( side );
          }
          else if ( facets.size() > 1 )
          {
            throw std::runtime_error( "facet not unique" );
          }
          else
          {
            match = false;
            break;
          }
        }
        if ( match )
        {
          {
            std::cout << "WARNING: Insert cut surface patch with " << patch.size() << " tris: { ";
            for ( std::set<Entity<3>*>::iterator i=patch.begin(); i!=patch.end(); ++i )
            {
              Entity<3> & tri = **i;
              if ( i!=patch.begin() )
                std::cout << ", ";
              std::cout << tri.Id();
            }
            std::cout << " }. Expect volume degeneration.\n";
          }
          for ( std::map<Facet*, std::vector<std::vector<int> > >::iterator i=patch_facets.begin();
                i!=patch_facets.end();
                ++i )
          {
            Facet * f = i->first;
            std::vector<std::vector<int> > & sides = i->second;

            std::vector<Epetra_SerialDenseMatrix> & side_coords = sides_xyz[f];
            CollectCoordinates( sides, side_coords );
          }
          break;
        }
      }
      if ( not match )
      {
        throw std::runtime_error( "unambiguous facet needed" );
      }
    }
    else
    {
      throw std::runtime_error( "unsupported patch topology" );
    }
  }

  for ( std::map<Facet*, FacetInfo>::iterator i=facet_info_.begin();
        i!=facet_info_.end();
        ++i )
  {
    Facet * f = i->first;
    FacetInfo & fi = i->second;
    if ( f->SideId() > -1 )
    {
      std::set<Entity<3>*> & tris = fi.tris_;

      std::vector<Epetra_SerialDenseMatrix> & side_coords = sides_xyz[f];

      std::vector<std::vector<int> > sides;
      FindProperSides( tris, sides );
      CollectCoordinates( sides, side_coords );
    }
  }
}

void GEO::CUT::TetMesh::FindProperSides( const std::set<Entity<3>*> & tris, std::vector<std::vector<int> > & sides )
{
  sides.reserve( tris.size() );
  for ( std::set<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> * tri = *i;

    bool done = false;
    const std::vector<Entity<4>*> & tets = tri->Parents();
    for ( std::vector<Entity<4>*>::const_iterator i=tets.begin();
          i!=tets.end();
          ++i )
    {
      Entity<4> * tet = *i;
      std::vector<int> & original_tet = tets_[tet->Id()];

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
            sides.push_back( std::vector<int>() );
            std::swap( sides.back(), side );
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

void GEO::CUT::TetMesh::CollectCoordinates( const std::vector<std::vector<int> > & sides,
                                            std::vector<Epetra_SerialDenseMatrix> & side_coords )
{
  for ( std::vector<std::vector<int> >::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    const std::vector<int> & side = *i;
    Epetra_SerialDenseMatrix xyz( 3, 3 );
    points_[side[0]]->Coordinates( &xyz( 0, 0 ) );
    points_[side[1]]->Coordinates( &xyz( 0, 1 ) );
    points_[side[2]]->Coordinates( &xyz( 0, 2 ) );
    side_coords.push_back( xyz );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
    surface_tris_.push_back( side );
#endif
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
    file << "};  // " << tri.Id() << "; "
         << tri[0] << ", "
         << tri[1] << ", "
         << tri[2]
         << "\n";
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
    if ( active_surface_tris_[tri.Id()] )
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

void GEO::CUT::TetMesh::GmshWriteTriSet( const std::string & name, const std::set<Entity<3>*> & tris )
{
  std::string filename = name + ".pos";
  std::ofstream file( filename.c_str() );
  file << "View \"" << name << "\" {\n";
  for ( std::set<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> & tri = **i;
    std::vector<int> t( tri(), tri()+3 );
    GmshWriteTri( file, tri.Id(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTetSet( const std::string & name, const std::set<Entity<4>*> & tets )
{
  std::string filename = name + ".pos";
  std::ofstream file( filename.c_str() );
  file << "View \"" << name << "\" {\n";
  for ( std::set<Entity<4>*>::const_iterator i=tets.begin(); i!=tets.end(); ++i )
  {
    Entity<4> & tet = **i;
    std::vector<int> t( tet(), tet()+4 );
    GmshWriteTet( file, tet.Id(), t );
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
