
#include <iterator>

#include <Epetra_SerialDenseMatrix.h>

#include "cut_element.H"
#include "cut_volumecell.H"
#include "cut_facet.H"

#ifdef QHULL
extern "C" {
#include <qhull/qhull_a.h>
}
#endif

namespace GEO
{
  namespace CUT
  {

struct PendingTet
{
  std::vector<VolumeCell*> owners_;
  std::vector<std::vector<Point*> > sides_;
};

struct VolumeCellTet
{
  void AddTet( std::vector<std::vector<Point*> >::iterator it, VolumeCell * owner )
  {
    cell_tets_[owner].push_back( it );

    std::vector<Point*> & t = *it;

    for ( int i=0; i<4; ++i )
    {
      std::vector<Point*> side( t );
      side.erase( side.begin() + i );

      if ( side.size()!=3 )
        throw std::runtime_error( "illegal side size" );

      bool all_oncutsurface = true;
      for ( std::vector<Point*>::iterator i=side.begin(); i!=side.end(); ++i )
      {
        Point * p = *i;
        if ( p->Position() != Point::oncutsurface )
        {
          all_oncutsurface = false;
          break;
        }
      }

      if ( all_oncutsurface )
      {
        if ( OnSameCutSurface( side ) )
        {
          cut_sides_[it].push_back( side );

          std::sort( side.begin(), side.end() );
          onsurface_tet_sides_[side] = owner;
        }
        else
        {
          std::sort( side.begin(), side.end() );
          internal_tet_sides_[side] = owner;
        }
      }
    }
  }

  bool OnSameCutSurface( const std::vector<Point*> & side )
  {
    std::set<Side*> cut_sides = side[0]->CutSides();
    side[1]->Intersection( cut_sides );
    side[2]->Intersection( cut_sides );

    for ( std::set<Side*>::iterator i=cut_sides.begin(); i!=cut_sides.end(); ++i )
    {
      Side * s = *i;
      if ( s->Id() > -1 )
      {
        return true;
      }
    }
    return false;
  }

  void AddPendingTet( std::vector<std::vector<Point*> >::iterator it,
                      const std::vector<VolumeCell*> & owners )
  {
    std::vector<Point*> & tet = *it;
    pending_tets_[it].owners_ = owners;
    std::vector<std::vector<Point*> > & sides = pending_tets_[it].sides_;

    for ( int i=0; i<4; ++i )
    {
      std::vector<Point*> side( tet );
      side.erase( side.begin() + i );

      //if ( not OnSameCutSurface( side ) )
      {
        std::sort( side.begin(), side.end() );
        sides.push_back( side );
      }
    }
  }

  void FillPendingTets( const std::set<VolumeCell*> & cells )
  {
    // look for cut surface tris that are known on one side but not on the
    // other

    for ( std::map<std::vector<std::vector<Point*> >::iterator, PendingTet >::iterator i=pending_tets_.begin();
          i!=pending_tets_.end();
      )
    {
      std::vector<std::vector<Point*> >::iterator j = i->first;
      //std::vector<Point*> & tet = *j;
      std::vector<std::vector<Point*> > & sides = i->second.sides_;
      std::vector<VolumeCell*> & owners = i->second.owners_;

      bool found = false;

      if ( owners.size()==2 )
      {
        for ( std::vector<std::vector<Point*> >::iterator is=sides.begin(); not found and is!=sides.end(); ++is )
        {
          std::vector<Point*> & side = *is;

          if ( side.size()!=3 )
          {
            throw std::runtime_error( "wrong number of points" );
          }

          std::vector<Point*> sorted_side( side );
          std::sort( sorted_side.begin(), sorted_side.end() );
          std::map<std::vector<Point*>, VolumeCell*>::iterator l = onsurface_tet_sides_.find( sorted_side );
          if ( l!=onsurface_tet_sides_.end() )
          {
            VolumeCell * other_owner = l->second;
            VolumeCell * owner;
            if ( other_owner==owners[0] )
            {
              owner = owners[1];
            }
            else if ( other_owner==owners[1] )
            {
              owner = owners[0];
            }
            else
            {
              throw std::runtime_error( "confused" );
            }

            AddTet( j, owner );
            pending_tets_.erase( i++ );
            found = true;
            break;
          }
        }
      }
      if ( not found )
      {
        ++i;
      }
    }

#if 0
    // Get tet volume cell ownership from element side facets.

    for ( std::map<std::vector<std::vector<Point*> >::iterator, PendingTet >::iterator i=pending_tets_.begin();
          i!=pending_tets_.end();
      )
    {
      std::vector<std::vector<Point*> >::iterator j = i->first;
      //std::vector<Point*> & tet = *j;
      std::vector<std::vector<Point*> > & sides = i->second.sides_;

      bool found = false;

      for ( std::vector<std::vector<Point*> >::iterator is=sides.begin(); not found and is!=sides.end(); ++is )
      {
        std::vector<Point*> & side = *is;

        if ( side.size()!=3 )
        {
          throw std::runtime_error( "wrong number of points" );
        }

        std::set<Facet*> facets = side[0]->Facets();
        side[1]->Intersection( facets );
        side[2]->Intersection( facets );

        //for ( std::set<Facet*>::iterator l=facets.begin(); not found and l!=facets.end(); ++l )
        if ( facets.size()==1 )
        {
          //Facet * f = *l;
          Facet * f = *facets.begin();
          if ( f->SideId() < 0 )
          {
            std::vector<VolumeCell*> owners;
            for ( std::set<VolumeCell*>::iterator ic=cells.begin(); ic!=cells.end(); ++ic )
            {
              VolumeCell * cell = *ic;
              if ( cell->OwnsFacet( f ) )
              {
                owners.push_back( cell );
              }
            }
            if ( owners.size()==0 )
            {
              //throw std::runtime_error( "facet without volume" );
            }
            else if ( owners.size()==1 )
            {
              AddTet( j, owners[0] );
              pending_tets_.erase( i++ );
              found = true;
              break;
            }
          }
        }
      }
      if ( not found )
      {
        ++i;
      }
    }
#endif

    // Get tet volume cell ownership from neighboring tets.

    for ( std::size_t numpending = pending_tets_.size(); numpending > 0; numpending = pending_tets_.size() )
    {
      for ( std::map<std::vector<std::vector<Point*> >::iterator, PendingTet >::iterator i=pending_tets_.begin();
            i!=pending_tets_.end();
        )
      {
        // Tets with all points on the cut surface.

        // Find the tet sides that are not on a cut surface. There have to be
        // some. Find another tet that shares the same surface. Use the same volume.

        std::vector<std::vector<Point*> >::iterator j = i->first;
        //std::vector<Point*> & tet = *j;
        std::vector<std::vector<Point*> > & sides = i->second.sides_;

        bool found = false;

        for ( std::vector<std::vector<Point*> >::iterator l=sides.begin(); l!=sides.end(); ++l )
        {
          std::vector<Point*> & side = *l;
          std::map<std::vector<Point*>, VolumeCell*>::iterator k = internal_tet_sides_.find( side );
          if ( k!=internal_tet_sides_.end() )
          {
            VolumeCell * owner = k->second;
            AddTet( j, owner );
            pending_tets_.erase( i++ );
            found = true;
            break;
          }
        }

        if ( not found )
        {
          ++i;
        }
      }
      if ( numpending == pending_tets_.size() )
      {
        throw std::runtime_error( "no progress" );
      }
    }
  }

  std::map<VolumeCell*, std::vector<std::vector<std::vector<Point*> >::iterator> > cell_tets_;

  std::map<std::vector<Point*>, VolumeCell*> internal_tet_sides_;

  std::map<std::vector<Point*>, VolumeCell*> onsurface_tet_sides_;

  std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > > cut_sides_;

  //std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > > pending_tets_;
  std::map<std::vector<std::vector<Point*> >::iterator, PendingTet > pending_tets_;
};

#if 1
void GmshWriteTri( std::ofstream & file, const std::vector<Point*> & t )
{
  file << "ST(";
  for ( std::vector<Point*>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = *j;
    if ( j != t.begin() )
      file << ",";
    file << p->X()[0] << ","
         << p->X()[1] << ","
         << p->X()[2];
  }
  file << "){";
  for ( std::vector<Point*>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = *j;
    if ( j != t.begin() )
      file << ",";
    file << p->Position();
  }
  file << "};\n";
}
void GmshWriteTet( std::ofstream & file, const std::vector<Point*> & t )
{
  file << "SS(";
  for ( std::vector<Point*>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = *j;
    if ( j != t.begin() )
      file << ",";
    file << p->X()[0] << ","
         << p->X()[1] << ","
         << p->X()[2];
  }
  file << "){";
  for ( std::vector<Point*>::const_iterator j=t.begin(); j!=t.end(); ++j )
  {
    Point * p = *j;
    if ( j != t.begin() )
      file << ",";
    file << p->Position();
  }
  file << "};\n";
}
#endif

  }
}

void GEO::CUT::Element::Delaunay( Mesh & mesh, std::vector<std::vector<Point*> > & tets )
{
#ifdef QHULL

  std::set<Point*> cut_points;
  GetAllPoints( cut_points );

  std::vector<Point*> points;
  points.reserve( cut_points.size() );
  points.assign( cut_points.begin(), cut_points.end() );

  const int dim = 3;
  const int n = cut_points.size();

  if ( n < 4 )
  {
    throw std::runtime_error( "illegal element topology" );
  }

  Epetra_SerialDenseMatrix coordinates( dim, n );

  for ( int i=0; i<n; ++i )
  {
    Point * p = points[i];
    p->Coordinates( &coordinates( 0, i ) );
  }

  boolT ismalloc = false;
  const char * options = "qhull d Qt Qbb Qc Qz Pp";
  //const char * options = "qhull d Qt Qbb Qc Qx";

  // If you want some debugging information replace the 0 pointer
  // with stdout or some other file open for writing.

  FILE * outfile = 0;
  FILE * errfile = stderr;

  if ( not qh_new_qhull( dim, n, &coordinates( 0, 0 ), ismalloc, const_cast<char*>( options ), outfile, errfile ) )
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
        std::vector<Point*> t;
        t.reserve( dim+1 );

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
          t.push_back( points[p] );
        }

        // make sure we do not get tets at any side due to rounding errors

        std::vector<Point*>::iterator is=t.begin();
        std::set<Side*> sides = ( *is )->CutSides();
        for ( ++is; is!=t.end(); ++is )
        {
          Point * p = *is;
          p->Intersection( sides );
          if ( sides.size()==0 )
          {
            break;
          }
        }

        if ( sides.size()==0 )
        {
          tets.push_back( t );
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

#if 0

//     std::copy( tets.begin(), tets.end(), std::ostream_iterator<int>( std::cout, " " ) );
//     std::cout << "\n";

    {
      std::ofstream file( "delaunaycells.pos" );
      file << "View \"delaunaycells\" {\n";
      for ( std::vector<std::vector<Point*> >::iterator i=tets.begin(); i<tets.end(); ++i )
      {
        std::vector<Point*> & t = *i;
        GmshWriteTet( file, t );
      }
      file << "};\n";
    }

#endif

    VolumeCellTet volume_tets;

    std::vector<std::vector<Point*> > broken_tets;

#define DELAUNAY_RETHROW
#ifdef DELAUNAY_RETHROW
    try
    {
#endif

    // repair delaunay tets if the cut surface is penetrated

    // This might cause problems if there are multiple cuts. The assumption is
    // that each cut separates inside and outside.

    for ( std::vector<std::vector<Point*> >::iterator i=tets.begin(); i<tets.end(); ++i )
    {
      std::vector<Point*> & t = *i;
      bool is_oncutsurface = false;
      bool is_inside = false;
      bool is_outside = false;
      for ( std::vector<Point*>::iterator j=t.begin(); j!=t.end(); ++j )
      {
        Point * p = *j;
        switch ( p->Position() )
        {
        case Point::undecided:
          throw std::runtime_error( "must not have undecided nodal positions at this point" );
        case Point::oncutsurface:
          is_oncutsurface = true;
          break;
        case Point::inside:
          is_inside = true;
          break;
        case Point::outside:
          is_outside = true;
          break;
        }
      }

      if ( is_inside and is_outside )
      {
        if ( not is_oncutsurface )
        {
          throw std::runtime_error( "penetration that does not touch cut surface" );
        }

//         std::vector<Point*>::iterator is=t.begin();
//         std::set<Side*> sides = ( *is )->CutSides();
//         for ( ++is; is!=t.end(); ++is )
//         {
//           Point * p = *is;
//           p->Intersection( sides );
//           if ( sides.size()==0 )
//           {
//             break;
//           }
//         }

//         if ( sides.size() > 0 )
//         {
//           // All tet points on the same side! This is most certainly a
//           // delaunay error!

//           t.clear();
//           continue;
//         }

        std::vector<Point*> oncutsurface;
        std::vector<Point*> inside;
        std::vector<Point*> outside;

        for ( std::vector<Point*>::iterator j=t.begin(); j!=t.end(); ++j )
        {
          Point * p = *j;
          switch ( p->Position() )
          {
          case Point::undecided:
            throw std::runtime_error( "must not have undecided nodal positions at this point" );
          case Point::oncutsurface:
            oncutsurface.push_back( p );
            break;
          case Point::inside:
            inside.push_back( p );
            break;
          case Point::outside:
            outside.push_back( p );
            break;
          }
        }

        if ( outside.size()==1 and inside.size()==1 )
        {
          Point * op = outside[0];
          Point * ip = inside[0];

          Edge * e = op->CommonEdge( ip );
          if ( e!=NULL )
          {
            bool point_on_edge = false;
            for ( std::vector<Point*>::iterator j=oncutsurface.begin(); j!=oncutsurface.end(); ++j )
            {
              Point * cp = *j;
              if ( cp->IsCut( e ) )
              {
                point_on_edge = true;
                break;
              }
            }
            if ( point_on_edge )
            {
              // A flat tet. Ignore it!
              t.clear();
            }
            else
            {
              broken_tets.push_back( t );
            }
          }
          else
          {
            broken_tets.push_back( t );
          }
        }
        else
        {
          broken_tets.push_back( t );
          //throw std::runtime_error( "delaunay tet repair needed" );
        }

        // flip sides here
      }
    }

    if ( broken_tets.size()>0 )
    {
#if 1
      std::ofstream file( "brokentets.pos" );
      file << "View \"brokentets\" {\n";
      for ( std::vector<std::vector<Point*> >::iterator i=broken_tets.begin();
            i!=broken_tets.end();
            ++i )
      {
        std::vector<Point*> & t = *i;
        GmshWriteTet( file, t );
      }
      file << "};\n";
      file.close();
#endif
      throw std::runtime_error( "delaunay tet repair needed" );
    }

    for ( std::vector<std::vector<Point*> >::iterator it=tets.begin(); it<tets.end(); ++it )
    {
      std::vector<Point*> & t = *it;
      if ( t.size()>0 )
      {
        std::vector<VolumeCell*> owners;
        for ( std::set<VolumeCell*>::iterator i=cells_.begin(); i!=cells_.end(); ++i )
        {
          VolumeCell * cell = *i;
          if ( cell->OwnsTet( t ) )
          {
            owners.push_back( cell );
          }
        }
        switch ( owners.size() )
        {
        case 0:
          throw std::runtime_error( "tet without owner" );
        case 1:
        {
          volume_tets.AddTet( it, owners[0] );
          break;
        }
        default:
        {
          volume_tets.AddPendingTet( it, owners );
        }
        }
      }
    }

    volume_tets.FillPendingTets( cells_ );

    for ( std::map<VolumeCell*, std::vector<std::vector<std::vector<Point*> >::iterator> >::iterator i=volume_tets.cell_tets_.begin();
          i!=volume_tets.cell_tets_.end();
          ++i )
    {
      VolumeCell * cell = i->first;
      std::vector<std::vector<std::vector<Point*> >::iterator> & tets = i->second;

      cell->SetIntegrationCells( mesh, tets, volume_tets.cut_sides_ );
    }

#ifdef DELAUNAY_RETHROW
    }
    catch ( std::runtime_error & err )
    {
#endif

#if 1
    {
      int count = 0;
      for ( std::map<VolumeCell*, std::vector<std::vector<std::vector<Point*> >::iterator> >::iterator i=volume_tets.cell_tets_.begin();
            i!=volume_tets.cell_tets_.end();
            ++i )
      {
        //VolumeCell * cell = i->first;
        std::vector<std::vector<std::vector<Point*> >::iterator> & tets = i->second;

        std::stringstream str;
        str << "volumecelltets" << count;
        std::ofstream file( ( str.str() + ".pos" ).c_str() );
        file << "View \"" << str.str() << "\" {\n";

        for ( std::vector<std::vector<std::vector<Point*> >::iterator>::iterator i=tets.begin();
              i!=tets.end();
              ++i )
        {
          std::vector<Point*> & t = **i;
          GmshWriteTet( file, t );
        }
        file << "};\n";

        str << "boundary";

        std::ofstream bfile( ( str.str() + ".pos" ).c_str() );
        bfile << "View \"" << str.str() << "\" {\n";
        for ( std::vector<std::vector<std::vector<Point*> >::iterator>::iterator i=tets.begin();
              i!=tets.end();
              ++i )
        {
          std::vector<std::vector<Point*> >::iterator it = *i;
          std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > >::iterator j = volume_tets.cut_sides_.find( it );
          if ( j!=volume_tets.cut_sides_.end() )
          {
            std::vector<std::vector<Point*> > & sides = j->second;
            for ( std::vector<std::vector<Point*> >::iterator i=sides.begin(); i!=sides.end(); ++i )
            {
              std::vector<Point*> & t = *i;
              GmshWriteTri( bfile, t );
            }
          }
        }
        bfile << "};\n";

        count += 1;
      }
    }

    {
      std::ofstream file( "pendingtets.pos" );
      file << "View \"pendingtets\" {\n";
      for ( std::map<std::vector<std::vector<Point*> >::iterator, PendingTet >::iterator i=volume_tets.pending_tets_.begin();
            i!=volume_tets.pending_tets_.end();
            ++i )
      {
        std::vector<std::vector<Point*> >::iterator j = i->first;
        std::vector<Point*> & t = *j;
        GmshWriteTet( file, t );
      }
      file << "};\n";
    }

#ifdef DELAUNAY_RETHROW
    throw;
    }
#endif

#endif

    if ( volume_tets.cell_tets_.size() != cells_.size() )
    {
      throw std::runtime_error( "empty volume cells" );
    }
  }
  else
    throw std::runtime_error( "qhull failed" );

#endif
}
