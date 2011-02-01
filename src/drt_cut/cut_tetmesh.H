#ifndef CUT_TETMESH_H
#define CUT_TETMESH_H

#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <iostream>

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

    /// Mesh of tet elements that is used to triangulate a VolumeCell
    /*!
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

        bool HasOnCutSurface( const std::vector<Point*> & points )
        {
          for ( int i=0; i<length; ++i )
          {
            Point * p = points[handle_[i]];
            if ( p->Position()==Point::oncutsurface )
              return true;
          }
          return false;
        }

        int CountMarkedParents( const std::vector<int> & marker )
        {
          int count = 0;
          for ( typename std::vector<Entity<length+1>*>::iterator i=parents_.begin();
                i!=parents_.end();
                ++i )
          {
            Entity<length+1> * pt = *i;
            if ( marker[pt->Id()] )
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
          return NULL;
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
          mypoints.clear();
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

      /// Representation of one possible triangulation of a facet
      /*!
        We find the triangulation of each facet from the qhull tet
        mesh. Sometimes, this triangulation is ambiguous. Thus we need to
        collect all possible interpretations of a facet and select the right
        one later on.
       */
      struct OverlappingTriSets
      {
        std::set<Entity<2>*> lset_;
        std::set<Entity<3>*> tset_;
        std::set<Entity<2>*> dlset_;
        std::set<Entity<2>*> openlines_;
        bool full_;

        OverlappingTriSets()
          : full_( true )
        {}

        bool Done()
        {
          return lset_.size()==0;
        }

        Entity<2> * NextLine() { return *lset_.begin(); }

        void NewLineTris( const std::set<Entity<3>*> & tris,
                          Entity<2> * line,
                          std::set<Entity<3>*> & new_tris )
        {
          const std::vector<Entity<3>*> & line_tris = line->Parents();
          for ( std::vector<Entity<3>*>::const_iterator i=line_tris.begin();
                i!=line_tris.end();
                ++i )
          {
            Entity<3> * tri = *i;
            if ( tris.count( tri ) > 0 and tset_.count( tri )==0 )
            {
              new_tris.insert( tri );
            }
          }
        }

        void EmptyTriSet( Entity<2> * line )
        {
          full_ = false;
          lset_.erase( line );
          dlset_.insert( line );
          openlines_.insert( line );
        }

        void InsertTri( Entity<3> * tri )
        {
          const std::vector<Entity<2>*> tri_lines = tri->Children();
          bool reject = false;
          for ( std::vector<Entity<2>*>::const_iterator i=tri_lines.begin();
                i!=tri_lines.end();
                ++i )
          {
            Entity<2> * line = *i;
            if ( lset_.count( line ) == 0 and dlset_.count( line ) > 0 )
            {
              // reject this tri, if it would overlap with any other tri
              full_ = false;
              reject = true;
              break;
            }
          }
          for ( std::vector<Entity<2>*>::const_iterator i=tri_lines.begin();
                i!=tri_lines.end();
                ++i )
          {
            Entity<2> * line = *i;
            dlset_.insert( line );
            if ( lset_.count( line ) > 0 )
            {
              lset_.erase( line );
            }
            else if ( not reject )
            {
              lset_.insert( line );
            }
          }
          if ( not reject )
          {
            tset_.insert( tri );
          }
        }
      };

      /// Collection of facet triangulations.
      struct FacetInfo
      {
        std::set<Entity<3>*> tris_;
        std::list<OverlappingTriSets> overlappingsets_;
        bool done_;

        FacetInfo()
          : done_( false )
        {}

        bool IsUnique()
        {
          return overlappingsets_.size()==1;
        }

        bool IsFull()
        {
          for ( std::list<OverlappingTriSets>::iterator i=overlappingsets_.begin();
                i!=overlappingsets_.end();
                ++i )
          {
            OverlappingTriSets & ots = *i;
            if ( ots.full_ )
              return true;
          }
          return false;
        }

        void SetMarker( std::vector<int> & marker, bool value )
        {
          for ( std::set<Entity<3>*>::iterator i=tris_.begin(); i!=tris_.end(); ++i )
          {
            Entity<3> * tri = *i;
            marker[tri->Id()] = value;
          }
        }

        void ExchangeTriSet( std::set<Entity<3>*> & tris, std::vector<int> & marker )
        {
          SetMarker( marker, false );
          std::swap( tris_, tris );
          SetMarker( marker, true );
        }
      };

    public:

      /// Construct a tet mesh by calling qhull and selecting all the tets and
      /// tris that are needed.
      /*!
        This is where the work is done.
       */
      TetMesh( const std::vector<Point*> & points,
               const std::set<Facet*> & facets );

      /// Return the list of tets.
      const std::vector<std::vector<int> > & Tets() const { return tets_; }

      /// Return the tris for each facet
      const std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & SidesXYZ() const { return sides_xyz_; }

    private:

      void FindTriFacets();

      void FindOverlappingTriFacets();

      void RemoveSimpleOverlappingTriFacets();

      bool FillSurfaceGaps( std::set<Entity<2>*> & surface_lines,
                            std::set<Entity<3>*> & surface_tris );

      void RemoveFullOverlappingTriFacets();

#if 0
      void RotateNonMatchingCutTets();
#endif

      /// deactivate all tets that are on a cut surface with all four points
      void DeactivateCutSurfaceTets();

      void ActivateCutSurfaceTris();

      /// test for illegal (cutted) cut surfaces
      void TestCutSurface();

      /// if there is just one tet at an cut surface it is active
      void ActivateSingleTetsOnFacets();

      /// activate the right tets at the cut surface
      void ActivateDoubleTetsOnCutSurface();

      void ActivateSurroundingTets( Entity<4> * tet );

      //void ActivateUsedCutSurfaceTris();

      /// remove deactivated tets
      void ClearExternalTets();

      void FixBrokenTets();

      /// Get tri coordinates at the cut facets in the proper order
      void FillCutSides( std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

      void FindProperSides( const std::set<Entity<3>*> & tris, std::vector<std::vector<int> > & sides );

      void CollectCoordinates( const std::vector<std::vector<int> > & sides,
                               std::vector<Epetra_SerialDenseMatrix> & side_coords );

#ifdef TETMESH_GMSH_DEBUG_OUTPUT

      void GmshWriteCells();

      void GmshWriteActiveCells();

      void GmshWriteFacetTris();

      void GmshWriteSurfaceCells();

      void GmshWriteActiveSurfaceCells();

      void GmshWriteSurfaceTris();

      void GmshWriteTriSet( const std::string & name, const std::set<Entity<3>*> & tris );

      void GmshWriteTetSet( const std::string & name, const std::set<Entity<4>*> & tets );

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

      const std::vector<Point*> & points_;
      const std::set<Facet*> & facets_;
      std::vector<std::vector<int> > tets_;

      std::vector<int> accept_tets_;
      std::vector<int> active_surface_tris_;

      std::vector<Entity<4> > tet_entities_;
      std::map<Handle<3>, Entity<3> > tet_surfaces_;
      std::map<Handle<2>, Entity<2> > tet_lines_;
      //std::map<Handle<1>, Entity<1> > tet_points_;

      std::map<Entity<3>*, Facet*> tri_facets_;
      std::map<Facet*, FacetInfo> facet_info_;

      std::set<Entity<3>*> broken_tris_;

      //std::map<Facet*, std::vector<Entity<3>*> > cutsurface_;
      std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > sides_xyz_;

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      std::vector<std::vector<int> > surface_tris_;
#endif
    };

    template <int length>
    std::ostream & operator<<( std::ostream & stream, const TetMesh::Handle<length> & h )
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

    template <int length>
    std::ostream & operator<<( std::ostream & stream, const TetMesh::Entity<length> & e )
    {
      stream << "( id=" << e.Id() << ", " << e.GetHandle() << " )";
      return stream;
    }

  }
}

#endif
