#ifndef CUT_TETMESH_H
#define CUT_TETMESH_H

#include <queue>

#include "cut_facet.H"
#include "cut_side.H"

#ifdef DEBUGCUTLIBRARY
#define TETMESH_GMSH_DEBUG_OUTPUT
#endif

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Facet;
    class VolumeCell;

    /*!
    \brief Mesh of tet elements that is used to triangulate a VolumeCell

      We use qhull for triangulation. qhull assumes convex volumes, thus it
      might create more tets than we need. Furthermore, it is not guaranteed
      to reconstruct the cut surface. We need to do post-processing of the
      qhull triangulation to get proper integration cells.

      \author u.kue
      \date 01/11
     */
    class TetMesh
    {

      /// Unique description of an entity in a tet mesh.
      /*!
        The node numbers of an entity, sorted.
       */
      template <int length>
      class Handle
      {
      public:

        Handle()
        {
        }

        explicit Handle( const int * points )
        {
          std::copy( points, points+length, points_ );
          std::sort( points_, points_+length );
        }

        Handle( const int * points, int skip )
        {
          int count = 0;
          for ( int i=0; i<length+1; ++i )
          {
            if ( i!=skip )
            {
              points_[count] = points[i];
              count += 1;
            }
          }
        }

        bool operator<( const Handle<length> & other ) const
        {
          for ( int i=0; i<length; ++i )
          {
            if ( points_[i] < other.points_[i] )
              return true;
            else if ( points_[i] > other.points_[i] )
              return false;
          }
          return false;
        }

        bool operator==( const Handle<length> & other ) const
        {
          bool fine = true;
          for ( int i=0; fine and i<length; ++i )
          {
            fine = fine and points_[i] == other.points_[i];
          }
          return fine;
        }

        bool operator!=( const Handle<length> & other ) const
        {
          return not ( *this ) == other;
        }

        int operator[]( int i ) const { return points_[i]; }

        const int * operator()() const { return points_; }

        bool Equals( const std::vector<int> & points )
        {
          for ( std::vector<int>::const_iterator i=points.begin();
                i!=points.end();
                ++i )
          {
            int p = *i;
            if ( std::find( points_, points_+length, p )==points_+length )
              return false;
          }
          return true;
        }

        bool Contains( int p )
        {
          return std::find( points_, points_+length, p )!=points_+length;
        }

        Handle<length+1> operator+( int p ) const
        {
          int points[length+1];
          std::copy( points_, points_+length, points );
          points[length] = p;
          std::sort( points, points+( length+1 ) );
          return Handle<length+1>( points );
        }

        friend std::ostream & operator<<( std::ostream & stream, const TetMesh::Handle<length> & h )
        {
          stream << "{";
          for ( int i=0; i<length; ++i )
          {
            if ( i > 0 )
              stream << ",";
            stream << h[i];
          }
          return stream << "}";
        }

      private:
        int points_[length];
      };

      /// An entity in a tet mesh. Might be tet, tri, line, point.
      template <int length>
      class Entity
      {
      public:

        Entity()
        {
        }

        Entity( int id, const Handle<length> & handle )
          : id_( id ),
            handle_( handle )
        {
        }

        int operator[]( int i ) const { return handle_[i]; }

        const int * operator()() const { return handle_(); }

        bool Equals( const std::vector<int> & points ) { return handle_.Equals( points ); }

        bool Contains( int p ) { return handle_.Contains( p ); }

        int Id() const { return id_; }

        void SetId( int i )
        {
          id_ = i;
        }

        const Handle<length> & GetHandle() const
        {
          return handle_;
        }

        void SetHandle( const Handle<length> & handle )
        {
          handle_ = handle;
        }

        void AddParent( Entity<length+1> * parent )
        {
          parents_.push_back( parent );
          parent->AddChild( this );
        }

        void RemoveParent( Entity<length+1> * parent )
        {
          typename std::vector<Entity<length+1>*>::iterator i =
            std::find( parents_.begin(), parents_.end(), parent );
          if ( i!=parents_.end() )
          {
            parents_.erase( i );
          }
        }

        void AddChild( Entity<length-1> * child )
        {
          children_.push_back( child );
        }

        void RemoveChild( Entity<length-1> * child )
        {
          typename std::vector<Entity<length-1>*>::iterator i =
            std::find( children_.begin(), children_.end(), child );
          if ( i!=children_.end() )
          {
            children_.erase( i );
          }
          //child->RemoveParent( this );
        }

        void Disconnect()
        {
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=children_.begin();
                i!=children_.end();
                ++i )
          {
            Entity<length-1> * child = *i;
            child->RemoveParent( this );
          }
          children_.clear();
          for ( typename std::vector<Entity<length+1>*>::const_iterator i=parents_.begin();
                i!=parents_.end();
                ++i )
          {
            Entity<length+1> * parent = *i;
            parent->RemoveChild( this );
          }
          parents_.clear();
        }

        /// create all my children
        void CreateChildren( std::map<Handle<length-1>, Entity<length-1> > & entities )
        {
          for ( int i=0; i<length; ++i )
          {
            Handle<length-1> h( handle_(), i );
            typename std::map<Handle<length-1>, Entity<length-1> >::iterator j = entities.find( h );
            if ( j!=entities.end() )
            {
              j->second.AddParent( this );
            }
            else
            {
              int id = entities.size();
              Entity<length-1> & e = entities[h];
              e.SetId( id );
              e.SetHandle( h );
              e.AddParent( this );
            }
          }
        }

        bool IsCut( const std::vector<Point*> & points, Facet * facet )
        {
          for ( int i=0; i<length; ++i )
          {
            Point * p = points[handle_[i]];
            if ( not p->IsCut( facet ) )
              return false;
          }
          return true;
        }

        const std::vector<Entity<length+1>*> & Parents() const { return parents_; }

        const std::vector<Entity<length-1>*> & Children() const { return children_; }

      private:

        int id_;
        Handle<length> handle_;
        std::vector<Entity<length+1>*> parents_;
        std::vector<Entity<length-1>*> children_;
      };

#ifdef CUT_USE_SORTED_VECTOR
      template <int length>
      class PlainEntitySet : public sorted_vector<Entity<length>*>
      {
      };
#else
      template <int length>
      class PlainEntitySet : public std::set<Entity<length>*>
      {
      };
#endif

      template <int length>
      class Domain
      {
      public:

        bool Empty() const { return members_.size()==0; }

        PlainEntitySet<length> & Members() { return members_; }

        PlainEntitySet<length-1> & Border() { return border_; }

        bool Contains( Entity<length> * m ) { return members_.count( m ) > 0; }

        void Add( Entity<length> * m, bool check_done=true )
        {
          const typename std::vector<Entity<length-1>*> & bs = m->Children();
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=bs.begin(); i!=bs.end(); ++i )
          {
            Entity<length-1> * b = *i;
            if ( check_done and done_border_.count( b ) > 0 )
            {
              throw std::runtime_error( "border entity has been visited before" );
            }
            if ( border_.count( b ) > 0 )
            {
              border_.erase( b );
              done_border_.insert( b );
            }
            else
            {
              border_.insert( b );
            }
          }
          members_.insert( m );
        }

        void Fill()
        {
          if ( members_.size()==0 )
          {
            throw std::runtime_error( "cannot fill without some members" );
          }

          // mark the external border as done since this will be the domains
          // boundary
          for ( typename PlainEntitySet<length-1>::iterator i=border_.begin(); i!=border_.end(); ++i )
          {
            Entity<length-1> * b = *i;
            bool skip = false;
            const std::vector<Entity<length>*> & ms = b->Parents();
            for ( typename std::vector<Entity<length>*>::const_iterator i=ms.begin(); i!=ms.end(); ++i )
            {
              Entity<length> * m = *i;
              if ( Contains( m ) )
              {
                skip = true;
                break;
              }
            }
            if ( not skip )
            {
              done_border_.insert( b );
            }
          }

          PlainEntitySet<length> stack;

          for ( typename PlainEntitySet<length>::iterator i=members_.begin();
                i!=members_.end();
                ++i )
          {
            Entity<length> * m = *i;
            PushNewNeighbors( stack, m );
          }

          while ( stack.size() > 0 )
          {
            Entity<length> * m = *stack.begin();
            stack.erase( m );

            Add( m, false );

            PushNewNeighbors( stack, m );
          }

          if ( border_.size()!=0 )
            throw std::runtime_error( "failed to fill domain" );
        }

      private:

        void PushNewNeighbors( PlainEntitySet<length> & stack, Entity<length> * m )
        {
          const std::vector<Entity<length-1>*> & bs = m->Children();
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=bs.begin();
                i!=bs.end();
                ++i )
          {
            Entity<length-1> * b = *i;
            if ( done_border_.count( b )==0 )
            {
              const std::vector<Entity<length>*> & ms = b->Parents();
              for ( typename std::vector<Entity<length>*>::const_iterator i=ms.begin(); i!=ms.end(); ++i )
              {
                Entity<length> * m = *i;
                if ( not Contains( m ) )
                {
                  stack.insert( m );
                }
              }
            }
          }
        }

        PlainEntitySet<length>   members_;
        PlainEntitySet<length-1> border_;
        PlainEntitySet<length-1> done_border_;
      };

      /// tracking of cut facet to find its tris
      class FacetMesh
      {
      public:

        bool Fill( TetMesh * tm, Facet * facet )
        {
          Domain<3> domain;
          PlainEntitySet<2> & lines = domain.Border();
          if ( FindTrace( tm, facet, lines ) )
          {
            while ( lines.size() > 0 )
            {
              bool match = false;
              for ( PlainEntitySet<2>::iterator i=lines.begin(); i!=lines.end(); ++i )
              {
                Entity<2> * l = *i;
                Entity<3> * tri = FindUniqueTri( tm, facet, domain, l );

                if ( tri!=NULL )
                {
                  match = true;
                  domain.Add( tri );
                  break;
                }
              }
              if ( not match )
              {
                // No way to fill the facet. Sorry.
                return false;
              }
            }

            std::swap( domain.Members(), tris_ );

            return true;
          }
          return false;
        }

        bool FindTrace( TetMesh * tm, Facet * facet, PlainEntitySet<2> & trace_lines )
        {
          std::map<std::pair<Point*, Point*>, plain_facet_set > lines;
          facet->GetLines( lines );

          for ( std::map<std::pair<Point*, Point*>, plain_facet_set >::iterator i=lines.begin();
                i!=lines.end();
                ++i )
          {
            const std::pair<Point*, Point*> & l = i->first;
            Handle<2> h = tm->MakeHandle( l.first, l.second );
            std::map<Handle<2>, Entity<2> >::iterator j = tm->tet_lines_.find( h );
            if ( j!=tm->tet_lines_.end() )
            {
              Entity<2> & line = j->second;
              trace_lines.insert( &line );
            }
            else
            {
              return false;
            }
          }
          return true;
        }

        const PlainEntitySet<3> & SurfaceTris() const { return tris_; }

      private:

        /// Find the tri at the given line that belongs to the facet. If there
        /// is no such tri (or it is not unique), return NULL.
        Entity<3> * FindUniqueTri( TetMesh * tm, Facet * facet, Domain<3> & domain, Entity<2> * l )
        {
          Entity<3> * tri = NULL;
          const std::vector<Entity<3>*> & tris = l->Parents();
          for ( std::vector<Entity<3>*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
          {
            Entity<3> * t = *i;
            if ( not domain.Contains( t ) )
            {
              if ( t->IsCut( tm->points_, facet ) )
              {
                if ( tri==NULL )
                {
                  tri = t;
                }
                else
                {
                  //throw std::runtime_error( "double tri on cut facet" );
                  return NULL;
                }
              }
            }
          }
          return tri;
        }

        PlainEntitySet<3> tris_;
      };

    public:

      /// Construct a tet mesh by calling qhull and selecting all the tets and
      /// tris that are needed.
      /*!
        This is where the work is done.
       */
      TetMesh( const std::vector<Point*> & points,
               const plain_facet_set & facets,
               bool project );

      void CreateElementTets( Mesh & mesh,
                              Element * element,
                              const plain_volumecell_set & cells,
                              const plain_side_set & cut_sides,
                              int count,
                              bool levelset,
                              bool tetcellsonly=false);

      void CreateVolumeCellTets();

      /// Return the list of tets.
      const std::vector<std::vector<int> > & Tets() const { return tets_; }

//       /// Return the tris for each facet
//       const std::map<Facet*, std::vector<Point*> > & SidesXYZ() const { return sides_xyz_; }

    private:

      void FixBrokenTets();

      void FindProperSides( const PlainEntitySet<3> & tris,
                            std::vector<std::vector<int> > & sides,
                            const PlainEntitySet<4> * members=NULL );

      void CollectCoordinates( const std::vector<std::vector<int> > & sides,
                               std::vector<Point*> & side_coords );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT

      void GmshWriteCells();

      void GmshWriteActiveCells();

      void GmshWriteSurfaceCells();

      void GmshWriteSurfaceTris();

      void GmshWriteTriSet( const std::string & name, const PlainEntitySet<3> & tris );

      void GmshWriteTetSet( const std::string & name, const PlainEntitySet<4> & tets );

      void GmshWriteTri( std::ostream & file, int eid, const std::vector<int> & t );

      void GmshWriteTet( std::ostream & file, int eid, const std::vector<int> & t );

      void GmshWriteConnect( std::ostream & file, std::string name, const std::vector<int> & t );

      void GmshWritePosition( std::ostream & file, int eid, const std::vector<int> & t );

#endif

      void Init();

      void CallQHull( const std::vector<Point*> & points,
                      std::vector<std::vector<int> > & tets,
                      bool project );

      bool IsValidTet( const std::vector<Point*> & t );

      void TestUsedPoints( const std::vector<std::vector<int> > & tets );

      bool FillFacetMesh();

      void SwapTetHandle( Entity<4> * tet, const Handle<4> & handle )
      {
        tet->Disconnect();
        tet->SetHandle( handle );

        std::vector<int> & t = tets_[tet->Id()];
        t.clear();
        std::copy( handle(), handle()+4, std::back_inserter( t ) );
      }

      template <int length>
        Entity<length> * SwapHandle( Entity<length> * e,
                                     const Handle<length> & handle,
                                     std::map<Handle<length>, Entity<length> > & entities )
      {
        Entity<length> & new_e = entities[handle];
        new_e.SetId( e->Id() );
        new_e.SetHandle( handle );
        e->Disconnect();
        entities.erase( e->GetHandle() );
        return &new_e;
      }

      Handle<2> MakeHandle( Point * p1, Point * p2 )
      {
        std::vector<Point*>::const_iterator i1 = std::find( points_.begin(), points_.end(), p1 );
        std::vector<Point*>::const_iterator i2 = std::find( points_.begin(), points_.end(), p2 );
        if ( i1==points_.end() or i2==points_.end() )
        {
          throw std::runtime_error( "point not in list" );
        }
        int points[2];
        points[0] = i1 - points_.begin();
        points[1] = i2 - points_.begin();
        std::sort( points, points+2 );
        return Handle<2>( points );
      }

      Handle<3> MakeHandle( Point * p1, Point * p2, Point * p3 )
      {
        std::vector<Point*>::const_iterator i1 = std::find( points_.begin(), points_.end(), p1 );
        std::vector<Point*>::const_iterator i2 = std::find( points_.begin(), points_.end(), p2 );
        std::vector<Point*>::const_iterator i3 = std::find( points_.begin(), points_.end(), p3 );
        if ( i1==points_.end() or i2==points_.end() or i3==points_.end() )
        {
          throw std::runtime_error( "point not in list" );
        }
        int points[3];
        points[0] = i1 - points_.begin();
        points[1] = i2 - points_.begin();
        points[2] = i3 - points_.begin();
        std::sort( points, points+3 );
        return Handle<3>( points );
      }

      bool FillFacet( Facet * f )
      {
        FacetMesh & cf = facet_mesh_[f];
        if ( not cf.Fill( this, f ) )
        {
          facet_mesh_.clear();
          return false;
        }
        return true;
      }

      void SeedDomain( Domain<4> & cell_domain, Facet * f, bool force=false )
      {
        if ( force or not f->OnCutSide() )
        {
          FacetMesh & fm = facet_mesh_[f];
          const PlainEntitySet<3> & tris = fm.SurfaceTris();
          for ( PlainEntitySet<3>::const_iterator i=tris.begin();
                i!=tris.end();
                ++i )
          {
            Entity<3> * t = *i;
            const std::vector<Entity<4>*> & tets = t->Parents();
            if ( tets.size()==1 )
            {
              if ( not cell_domain.Contains( tets[0] ) )
                cell_domain.Add( tets[0] );
            }
          }
        }
      }

      const std::vector<Point*> & points_;
      const plain_facet_set & facets_;
      std::vector<std::vector<int> > tets_;

      std::vector<int> accept_tets_;

      std::vector<Entity<4> > tet_entities_;
      std::map<Handle<3>, Entity<3> > tet_surfaces_;
      std::map<Handle<2>, Entity<2> > tet_lines_;
      //std::map<Handle<1>, Entity<1> > tet_points_;

      std::map<Facet*, FacetMesh> facet_mesh_;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      std::vector<std::vector<int> > surface_tris_;
#endif
    };
  }
}

#endif
