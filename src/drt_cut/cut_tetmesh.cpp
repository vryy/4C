
#include <stack>

#include "cut_tetmesh.H"
#include "cut_tetmeshintersection.H"
#include "cut_levelsetside.H"
#include "cut_volumecell.H"

extern "C" {
#include <qhull/qhull_a.h>
}

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
                            const plain_facet_set & facets,
                            bool project )
  : points_( points ),
    facets_( facets )
{
  std::vector<std::vector<int> > original_tets;

  CallQHull( points_, original_tets, project );

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

#if 0
  // if this is not the first cut, it might be fine not to have all points
  TestUsedPoints( tets_ );
#endif

  Init();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteCells();
#endif
}

bool GEO::CUT::TetMesh::FillFacetMesh()
{
  for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( not FillFacet( f ) )
      return false;
    if ( f->HasHoles() )
    {
      const plain_facet_set & holes = f->Holes();
      for ( plain_facet_set::const_iterator i=holes.begin(); i!=holes.end(); ++i )
      {
        Facet * h = *i;
        if ( not FillFacet( h ) )
          return false;
      }
    }
  }
  return true;
}

void GEO::CUT::TetMesh::CreateElementTets( Mesh & mesh,
                                           Element * element,
                                           const plain_volumecell_set & cells,
                                           const plain_side_set & cut_sides,
                                           int count,
                                           bool tetcellsonly)
{
  FixBrokenTets();

  if ( FillFacetMesh() )
  {
    for ( plain_volumecell_set::const_iterator i=cells.begin();
          i!=cells.end();
          ++i )
    {
      VolumeCell * vc = *i;

      Domain<4> cell_domain;
      PlainEntitySet<4> & cell_members = cell_domain.Members();
      PlainEntitySet<3> & cell_border  = cell_domain.Border();

#ifdef DEBUGCUTLIBRARY
      std::vector<Side*> facet_sides;
#endif

      const plain_facet_set & facets = vc->Facets();
      for ( plain_facet_set::const_iterator i=facets.begin();
            i!=facets.end();
            ++i )
      {
        Facet * f = *i;
        FacetMesh & fm = facet_mesh_[f];
        const PlainEntitySet<3> & tris = fm.SurfaceTris();

#ifdef DEBUGCUTLIBRARY
        facet_sides.push_back( f->ParentSide() );
#endif

        std::copy( tris.begin(), tris.end(), std::inserter( cell_border, cell_border.begin() ) );
      }

      for ( plain_facet_set::const_iterator i=facets.begin();
            i!=facets.end();
            ++i )
      {
        Facet * f = *i;
        SeedDomain( cell_domain, f );
      }

      if ( cell_domain.Empty() )
      {
        // Emergency. If this is the only volume cell within the element (an
        // element is surrounded by cut surfaces), try and force any available
        // tet into the volume cell. (This is a special case that turns up in
        // debug situations only.)

        if ( cells.size()==1 )
        {
          for ( plain_facet_set::const_iterator i=facets.begin();
                i!=facets.end();
                ++i )
          {
            Facet * f = *i;
            SeedDomain( cell_domain, f, true );
          }
        }

        if ( cell_domain.Empty() )
        {
          // Assume the volume cell in question is degenerated and does not
          // contain any tets.
          continue;
        }
      }

      cell_domain.Fill();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      GmshWriteTriSet( "cell_border" , cell_border  );
      GmshWriteTetSet( "cell_members", cell_members );
#endif

      std::vector<std::vector<Point*> > tets;
      tets.reserve( cell_members.size() );

      for ( PlainEntitySet<4>::iterator i=cell_members.begin(); i!=cell_members.end(); ++i )
      {
        Entity<4> & t = **i;
        if ( accept_tets_[t.Id()] )
        {
          std::vector<int> & fixedtet = tets_[t.Id()];
          if ( fixedtet.size()!=4 )
            throw std::runtime_error( "confused" );
          tets.push_back( std::vector<Point*>( 4 ) );
          std::vector<Point*> & tet = tets.back();
          for ( int i=0; i<4; ++i )
          {
            tet[i] = points_[fixedtet[i]];
          }
        }
      }

      // all facets (whose facet-mesh are just tri3s here) which are on cut surface obtain the coordinates
      // to create boundary integration cells then
      std::map<Facet*, std::vector<Point*> > sides_xyz;

      for ( plain_facet_set::const_iterator i=facets.begin();
            i!=facets.end();
            ++i )
      {
        Facet * f = *i;
        if ( f->OnCutSide() )
        {
          FacetMesh & fm = facet_mesh_[f];
          const PlainEntitySet<3> & tris = fm.SurfaceTris();
          std::vector<Point*> & side_coords = sides_xyz[f]; // create entry for facet and get a reference to the facets side coordinates
          std::vector<std::vector<int> > sides;
          FindProperSides( tris, sides, &cell_members );
          CollectCoordinates( sides, side_coords ); // fill the side coordinates, if all the side's coordinates are on cut surface
        }
      }

      vc->CreateTet4IntegrationCells( mesh, tets, sides_xyz );

      if ( vc->Empty() )
      {
        throw std::runtime_error( "empty volume cell detected" );
      }
    }
  }
  else
  {
    //if ( count <= 3 )
    {
      TetMeshIntersection intersection( mesh.CreateOptions(), element, tets_, accept_tets_, points_, cut_sides );
      intersection.Cut( mesh, element, cells, count ,tetcellsonly);
    }
  }
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

//   active_surface_tris_.resize( tet_surfaces_.size() );
//   std::fill( active_surface_tris_.begin(), active_surface_tris_.end(), false );
}

void GEO::CUT::TetMesh::CallQHull( const std::vector<Point*> & points,
                                   std::vector<std::vector<int> > & tets,
                                   bool project )
{
  const int dim = 3;
  const int n = points.size();

  if ( n < 4 )
  {
    throw std::runtime_error( "illegal element topology" );
  }
//   if ( n == 4 )
//   {
//     throw std::runtime_error( "no need to triangulate" );
//   }

  std::vector<double> coordinates( dim*n );

  if ( project )
  {
    LINALG::Matrix<3,1> m;
    m = 0;
    double scale = 1./n;
    for ( int i=0; i<n; ++i )
    {
      Point * p = points[i];
      LINALG::Matrix<3,1> x( p->X() );
      m.Update( scale, x, 1 );
    }
    double length = 0;
    LINALG::Matrix<3,1> l;
    for ( int i=0; i<n; ++i )
    {
      Point * p = points[i];
      LINALG::Matrix<3,1> x( p->X() );
      l = m;
      l.Update( 1, x, -1 );
      double n = l.Norm2();
      length = std::max( n, length );
    }
#ifdef DEBUGCUTLIBRARY
    std::ofstream pointfile( "points.plot" );
    pointfile << m( 0 ) << " " << m( 1 ) << " " << m( 2 ) << "\n";
#endif
    for ( int i=0; i<n; ++i )
    {
      Point * p = points[i];
      LINALG::Matrix<3,1> x( p->X() );
      l = m;
      l.Update( 1, x, -1 );
      double n = l.Norm2();
      l.Scale( length/n );
      l.Update( 1, m, 1 );
      std::copy( l.A(), l.A()+3, &coordinates[dim*i] );
#ifdef DEBUGCUTLIBRARY
      pointfile << l( 0 ) << " " << l( 1 ) << " " << l( 2 ) << "\n";
#endif
    }
  }
  else
  {
    for ( int i=0; i<n; ++i )
    {
      Point * p = points[i];
      p->Coordinates( &coordinates[dim*i] );
    }
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
      plain_int_set used_points;
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

  throw std::runtime_error( "qhull failed: Maybe the wrong version is used. Check your installation." );
}

bool GEO::CUT::TetMesh::IsValidTet( const std::vector<Point*> & t )
{
  plain_side_set sides;
  FindCommonSides( t, sides );
  if ( sides.size()==0 )
  {
    plain_facet_set facets;
    FindCommonFacets( t, facets );
    if ( facets.size()==0 )
    {
      return true;
    }
    for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); ++i )
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

  plain_facet_set facets;
  FindCommonFacets( t, facets );
  if ( facets.size()==0 )
  {
    return true;
  }
  for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->IsTriangulated() )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::TetMesh::TestUsedPoints( const std::vector<std::vector<int> > & tets )
{
  plain_int_set used_points;
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

    // compute normal distance of point to plane of the three remaining points
    double distance = nplane012.Dot( v03 );

    // compute norm (area) of plane
    //double norm012 = nplane012.Norm2();

    double vol_tet = distance / 6.0;

    // Deactivate all tets that are too small. We might still need the tet to
    // create a cut surface tri. Afterwards we will discard it.

    if ( fabs( vol_tet ) < VOLUMETOL )
    {
      accept_tets_[i - tets_.begin()] = false;
    }
//     else if ( fabs( distance / norm012 ) < 1e-7 )
//     {
//       accept_tets_[i - tets_.begin()] = false;
//     }

    // tet numbering wrong exchange 1 with 3
    if ( distance < 0 )
    {
      std::swap( t[1], t[3] );
    }
  }
}

void GEO::CUT::TetMesh::FindProperSides( const PlainEntitySet<3> & tris,
                                         std::vector<std::vector<int> > & sides,
                                         const PlainEntitySet<4> * members )
{
  sides.reserve( tris.size() );
  for ( PlainEntitySet<3>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> * tri = *i;

    bool done = false;
    const std::vector<Entity<4>*> & tets = tri->Parents();
    for ( std::vector<Entity<4>*>::const_iterator i=tets.begin();
          i!=tets.end();
          ++i )
    {
      Entity<4> * tet = *i;

      if ( members!=NULL and members->count( tet )==0 )
        continue;

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

/// collects the coordinates for the tri3 sides of the facet if all its points are on cut surface
void GEO::CUT::TetMesh::CollectCoordinates( const std::vector<std::vector<int> > & sides,
                                            std::vector<Point*> & side_coords )
{
  for ( std::vector<std::vector<int> >::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    const std::vector<int> & side = *i;

    Point * p1 = points_[side[0]];
    Point * p2 = points_[side[1]];
    Point * p3 = points_[side[2]];

    if ( p1->Position()==Point::oncutsurface and
         p2->Position()==Point::oncutsurface and
         p3->Position()==Point::oncutsurface )
    {
      side_coords.push_back( p1 );
      side_coords.push_back( p2 );
      side_coords.push_back( p3 );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      surface_tris_.push_back( side );
#endif
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

void GEO::CUT::TetMesh::GmshWriteTriSet( const std::string & name, const PlainEntitySet<3> & tris )
{
  std::string filename = name + ".pos";
  std::ofstream file( filename.c_str() );
  file << "View \"" << name << "\" {\n";
  for ( PlainEntitySet<3>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
  {
    Entity<3> & tri = **i;
    std::vector<int> t( tri(), tri()+3 );
    GmshWriteTri( file, tri.Id(), t );
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTetSet( const std::string & name, const PlainEntitySet<4> & tets )
{
  std::string filename = name + ".pos";
  std::ofstream file( filename.c_str() );
  file << "View \"" << name << "\" {\n";
  for ( PlainEntitySet<4>::const_iterator i=tets.begin(); i!=tets.end(); ++i )
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
