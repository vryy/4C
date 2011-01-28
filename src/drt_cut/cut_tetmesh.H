#ifndef CUT_TETMESH_H
#define CUT_TETMESH_H

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "cut_point.H"
#include "cut_facet.H"
#include "cut_side.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#ifdef DEBUGCUTLIBRARY
#define TETMESH_GMSH_DEBUG_OUTPUT
#endif

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Facet;


    class TetMesh
    {
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

      private:
        int points_[length];
      };

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

        const Handle<length> & GetHandle()
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

        /// find all the cut facets this entity belongs to
        void FindCutEntities( const std::vector<Point*> & points,
                              const std::set<Facet*> & facets,
                              std::map<Facet*, std::vector<Entity<3>*> > & cutsurface )
        {
          if ( AllOnCutSurface( points ) )
          {
            std::vector<Point*> mypoints;
            this->PointVector( points, mypoints );
            std::set<Facet*> cut_facets;
            FindCommonFacets( mypoints, cut_facets );
            if ( cut_facets.size() > 0 )
            {
              for ( std::set<Facet*>::iterator i=cut_facets.begin();
                    i!=cut_facets.end();
                    ++i )
              {
                Facet * f = *i;
                if ( facets.count( f ) > 0 and f->SideId() >= 0 )
                {
                  if ( f->IsTriangulated() )
                  {
                    if ( f->IsTriangulatedSide( mypoints ) )
                    {
                      cutsurface[f].push_back( this );
                    }
                  }
                  else
                  {
                    cutsurface[f].push_back( this );
                  }
                }
              }
            }
#if 0
            else if ( parents_.size()==1 )
            {
              // Collect all "surfaces" that are at the cut surface but do not
              // belong to any cut facet. If there is one parent (we are at a
              // surface) and its other point does not belong to the cut
              // surface, we have a match.
              //
              // These are the "surfaces" that need to be rotated. Along with
              // its parents.

              Entity<length+1> * parent = parents_[0];
              Point * p = points[parent->OtherPoint( this )];
              if ( p->Position()!=Point::oncutsurface )
              {
                cutsurface[NULL].push_back( this );
              }
            }
#endif
          }
        }

        /// test if this entity belongs entirely to a cut surface
        bool AllOnCutSurface( const std::vector<Point*> & points )
        {
          for ( int i=0; i<length; ++i )
          {
            Point * p = points[handle_[i]];
            if ( p->Position()!=Point::oncutsurface )
              return false;
          }
          return true;
        }

#if 0
        bool AtElementSide( const std::vector<Point*> & points,
                            const std::set<Facet*> & facets )
        {
          std::vector<Point*> mypoints;
          this->PointVector( points, mypoints );
          std::set<Side*> cut_sides;
          FindCommonSides( mypoints, cut_sides );
          if ( cut_sides.size()==0 )
          {
            return false;
          }
          for ( std::set<Side*>::iterator i=cut_sides.begin(); i!=cut_sides.end(); ++i )
          {
            Side * side = *i;
            if ( side->Id() < 0 )
            {
              // this needs to be unique!
              std::set<Facet*> cut_facets;
              FindCommonFacets( mypoints, cut_facets );
              for ( std::set<Facet*>::iterator i=cut_facets.begin();
                    i!=cut_facets.end();
                    ++i )
              {
                Facet * f = *i;
                if ( facets.count( f ) > 0 )
                {
                  return true;
                }
              }
            }
          }
          return false;
        }
#endif

        int CountMarkedParents( const std::vector<int> & marker1,
                                const std::vector<int> & marker2 )
        {
          int count = 0;
          for ( typename std::vector<Entity<length+1>*>::iterator i=parents_.begin();
                i!=parents_.end();
                ++i )
          {
            Entity<length+1> * pt = *i;
            if ( marker1[pt->Id()] and marker2[pt->Id()] )
            {
              count += 1;
            }
          }
          return count;
        }

        int OtherPoint( Entity<length-1> * child )
        {
          for ( int i=0; i<length; ++i )
          {
            if ( not child->Contains( handle_[i] ) )
              return handle_[i];
          }
          throw std::runtime_error( "malformed entity" );
        }

        Entity<length-1> * CommonChild( Entity<length> * other )
        {
          const std::vector<Entity<length-1>*> & other_children = other->Children();
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=children_.begin();
                i!=children_.end();
                ++i )
          {
            Entity<length-1> * child = *i;
            if ( std::find( other_children.begin(), other_children.end(), child )!=other_children.end() )
              return child;
          }
          throw std::runtime_error( "exactly one common child required" );
        }

        const std::vector<Entity<length+1>*> & Parents() const { return parents_; }

        const std::vector<Entity<length-1>*> & Children() const { return children_; }

        void Neighbors( std::set<Entity<length>*> & neighbors )
        {
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=children_.begin();
                i!=children_.end();
                ++i )
          {
            Entity<length-1> * child = *i;
            const std::vector<Entity<length>*> & n = child->Parents();
            std::copy( n.begin(), n.end(), std::inserter( neighbors, neighbors.begin() ) );
          }
          neighbors.erase( this );
        }

        void PointVector( const std::vector<Point*> & points, std::vector<Point*> & mypoints )
        {
          mypoints.reserve( length );
          for ( int i=0; i<length; ++i )
          {
            Point * p = points[handle_[i]];
            mypoints.push_back( p );
          }
        }

      private:

        int id_;
        Handle<length> handle_;
        std::vector<Entity<length+1>*> parents_;
        std::vector<Entity<length-1>*> children_;
      };

    public:

      TetMesh( const std::vector<Point*> & points,
               const std::set<Facet*> & facets );

      const std::vector<std::vector<int> > & Tets() const { return tets_; }

      const std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & SidesXYZ() const { return sides_xyz_; }

    private:

      void FindTriFacets();

      /// find all tris at a cut surface
      void FindCutSurfaceTris();

#if 0
      void RotateNonMatchingCutTets();
#endif

      /// deactivate all tets that are on a cut surface with all four points
      void DeactivateCutSurfaceTets();

      void ActivateCutSurfaceTris();

      void DeactivateOverlappingTris();

      /// test for illegal (cutted) cut surfaces
      void TestCutSurface();

      /// if there is just one tet at an cut surface it is active
      void ActivateSingleTetsOnCutSurface();

      /// activate the right tets at the cut surface
      void ActivateDoubleTetsOnCutSurface();

      void ActivateSurroundingTets( Entity<4> * tet );

      void ActivateUsedCutSurfaceTris();

      /// remove deactivated tets
      void ClearExternalTets();

      void FixBrokenTets();

      /// Get tri coordinates at the cut facets in the proper order
      void FillCutSides( std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT

      void GmshWriteCells();

      void GmshWriteActiveCells();

      void GmshWriteBrokenSurfaceTris( const std::set<Entity<3>*> & tris );

      void GmshWriteFacetTris();

      void GmshWriteSurfaceCells();

      void GmshWriteActiveSurfaceCells();

      void GmshWriteSurfaceTris();

      void GmshWriteTri( std::ostream & file, int eid, const std::vector<int> & t );

      void GmshWriteTet( std::ostream & file, int eid, const std::vector<int> & t );

      void GmshWriteConnect( std::ostream & file, std::string name, const std::vector<int> & t );

      void GmshWritePosition( std::ostream & file, int eid, const std::vector<int> & t );

#endif

      void Init();

      void CallQHull( const std::vector<Point*> & points,
                      std::vector<std::vector<int> > & tets );

      bool IsValidTet( const std::vector<Point*> & t );

      void TestUsedPoints( const std::vector<std::vector<int> > & tets );

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

      const std::vector<Point*> & points_;
      const std::set<Facet*> & facets_;
      std::vector<std::vector<int> > tets_;

      std::vector<int> accept_tets_;
      std::vector<int> accept_tris_;
      std::vector<int> active_surface_tris_;

      std::vector<Entity<4> > tet_entities_;
      std::map<Handle<3>, Entity<3> > tet_surfaces_;
      std::map<Handle<2>, Entity<2> > tet_lines_;
      //std::map<Handle<1>, Entity<1> > tet_points_;

      std::map<Entity<3>*, Facet*> tri_facets_;

      std::map<Facet*, std::vector<Entity<3>*> > cutsurface_;
      std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > sides_xyz_;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      std::vector<std::vector<int> > surface_tris_;
#endif
    };

  }
}

#endif
