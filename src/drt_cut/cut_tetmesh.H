/*---------------------------------------------------------------------*/
/*!
\file cut_tetmesh.H

\brief Tessellate the element by means of a Delauney triangulation

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_TETMESH_H
#define CUT_TETMESH_H

#include <queue>

#include "cut_facet.H"
#include "cut_side.H"

#ifdef DEBUGCUTLIBRARY
#define TETMESH_GMSH_DEBUG_OUTPUT
#endif

// FLAGS FOR EXTENDED OUTPUT
//#define TETMESH_EXTENDED_DEBUG_OUTPUT
//#define QHULL_EXTENDED_DEBUG_OUTPUT

// FLAGS FOR NEW IMPLEMENTATION TESTING
// Probably don't want to use:
//#define REMOVE_ALL_TETS_ON_CUTSIDE

// Probably want to use:
#define DIFF_QHULL_CALL
#define NEW_SEED_DOMAIN
#define NEW_POSTOL_TET

//Extra safety test: (leads to error if REMOVE_ALL_TETS_ON_CUTSIDE is defined in cut_test hex8_quad4_alex38)
//#define TEST_EMPTY_VC_IS_COMPLETELY_ON_CUTSIDE


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
          return not(( *this ) == other);
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
        /// This creats from, for instance a tet, 4 surfaces and also stores handle to access "parent" element.
        /// Same thing is done for a surface where lines are created as children.
        /// Or at least as far as my understanding goes.... //m.w
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

        /* Do all the points in this tri cut the facet?
         * I.e. make sure all points are on the facet. Otherwise return false!
         */
        bool IsCut( const std::vector<Point*> & points, Facet * facet )
        {
          for ( int i=0; i<length; ++i )
          {
            Point * p = points[handle_[i]];
            if ( not p->IsCut( facet ) )
            {
              return false;
            }
          }
          return true;
        }

        /*!
        \brief Stores connectivity information i.e. surface information for a line (what tet-surfaces a line is contributing to).
         */
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

        /// Is the domain empty?
        bool Empty() const { return members_.size()==0; }

        /// Members contains for a TET-domain the TETs and for a TRI-domain the TRIs
        PlainEntitySet<length> & Members() { return members_; }

        /// Border contains for a TET-domain the surfaces and for a TRI-domain the lines
        PlainEntitySet<length-1> & Border() { return border_; }

        /// Does this member already exist in the domain?
        bool Contains( Entity<length> * m ) { return members_.count( m ) > 0; }

        /* Add a tri or a tet to the domain (either a mesh containing TETs for a cell or a mesh containing TRIs for a surface).
         * Check the border of the added element:
         *    If it is on the border of the added domain i.e. for a facet its defining lines. Add the border to done_border_.
         *    If it has an entity which is in the interior of the domain (i.e. triangulization of a facet), add this to the
         *       border_.
         *    If an entity (tet or tri) is added but one of its children (tri or line) have already been added to done_border_
         *       throw an error as this can not happen.
         */
        void Add( Entity<length> * m, bool check_done=true )
        {
          const typename std::vector<Entity<length-1>*> & bs = m->Children();
          for ( typename std::vector<Entity<length-1>*>::const_iterator i=bs.begin(); i!=bs.end(); ++i )
          {
            Entity<length-1> * b = *i;
            //Each interior border_ should be visited twice (either surface for tets or lines for facet mesh),
            // and on the boundary (of the (facet or tet)/domain) only once.
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

        /*!
          \brief Fill the inner domain of the volume cell for tets.
              SeedDomain() has already been called, i.e. tets associated with boundary-tris are already covered.
         */
        void Fill()
        {
          if ( members_.size()==0 )
          {
            throw std::runtime_error( "cannot fill without some members" );
          }

          // mark the external border as done since this will be the domains
          // boundary
          // Check if a border_ (i.e. tri) has been forgotten on the boundary. And if so add it to the done_border_.
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

          //Adds tets not in cell_domain already. These are the interior tets not associated with
          //  an element boundary.
          for ( typename PlainEntitySet<length>::iterator i=members_.begin();
                i!=members_.end();
                ++i )
          {
            Entity<length> * m = *i;
            PushNewNeighbors( stack, m );
          }

          //Loop over the stack of "interior tets" associated with border_ entities (i.e. surfaces)
          //  and add these to the members_ and done_border_
          while ( stack.size() > 0 )
          {
            Entity<length> * m = *stack.begin();
            stack.erase( m );

            Add( m, false );

            PushNewNeighbors( stack, m );
          }

          //Throw error if the domain is not filled.
          //  All border_ should be converted into border_done_!!!
          if ( border_.size()!=0 )
            throw std::runtime_error( "failed to fill domain" );
        }

      private:

        /*!
          \brief Find if a border of a tet (i.e. tri-surface) is not contained in done_border_,
               if not then check if any of the tets associated with the border (i.e. tri-surface) is not
               contained in the Domain. If so then insert the tet to the stack.
               QUESTION: Is this testing a bit too much?
         */
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

        /* Create a facet mesh for a facet in a TetMesh (if possible)
         *   Make sure that the boundaries (i.e. lines) of the facet exists in the TetMesh.
         *    If NO ->  return false
         *    If Yes -> Make sure that this line is associated with an unique tri on this facet
         *       If NO ->  return false
         *       If Yes -> Add the created tris.
         *                 return true
         */
        bool Fill( TetMesh * tm, Facet * facet )
        {
          Domain<3> domain;
          PlainEntitySet<2> & lines = domain.Border();
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
          facet->PrintPointIds();
#endif
          if ( FindTrace( tm, facet, lines ) )
          {
            if(lines.size()<3)
            {
              std::cout << "lines.size(): " << lines.size() << std::endl;
              dserror("How can a facet contain less than 3 lines?");
            }
            //Loop through the lines. The lines are connected to the domain. I.e. when a tri is added it's "boundary lines"
            //  of the domain are removed from the boundary, but its interior lines are added.
            //Example: For a facet with 4 lines: An unique tri is found. It has 2 outer boundaries and one inner boundary
            //          This leads to that lines.size()=3 after wards 2 outer removed and 1 inner added.
            while ( lines.size() > 0 )
            {
              bool match = false;
              for ( PlainEntitySet<2>::iterator i=lines.begin(); i!=lines.end(); ++i )
              {
                Entity<2> * l = *i;
                Entity<3> * tri = FindUniqueTri( tm, facet, domain, l );

                if ( tri!=NULL )
                {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
                  std::cout << "ADDED tri!" <<std::endl;
                  std::cout << "tri->GetHandle()[0]: " << tri->GetHandle()[0]<< ", tri->GetHandle()[1]: " << tri->GetHandle()[1] << ", tri->GetHandle()[2]: " << tri->GetHandle()[2] << std::endl;
#endif
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

        /* Get the lines of the (un-triangulated?) facets.
           Create a handle for the TetMesh of the line. Check that the line exists in the previously created (TetInit(..)) tet_lines_.
           If it exist          -  add lines of facet to trace_lines (i.e. border_) and continue.
           If it DOES NOT exist -  return false. This would mean that a line of one facet does not exist in the
                                   lines created by the tet mesh.
                                   Probably the tetmesh runs through a cut-side, as such all facet lines will not be
                                   included.
         */
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
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
              std::cout << "DID NOT FIND HANDLE FROM tet_lines_ map!!!!" << std::endl;
#endif
              return false;
            }
          }
          return true;
        }

        const PlainEntitySet<3> & SurfaceTris() const { return tris_; }

      private:


        /* See if there exists an unique tri for this line. Taking into consideration not to add tris double.
         * And also that
         *  Get all the tris connected to this line. Now test each of them:
         *    Check if this tri is already added in the domain?
         *      If Yes -> Continue to next tri.
         *      If NO  -> Check if the tri lies on the Facet (IsCut(..))
         *          If NO  -> Continue to next tri.
         *          If Yes -> Check if a tri has already been assigned for this line.
         *              If Yes -> return NULL (i.e. no unique tri exists)
         *              If NO  -> Add tri and continue loop over tris.
         */
        Entity<3> * FindUniqueTri( TetMesh * tm, Facet * facet, Domain<3> & domain, Entity<2> * l )
        {
          Entity<3> * tri = NULL;
          const std::vector<Entity<3>*> & tris = l->Parents(); //TRIs connected to this line.
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
                  //Might lead to level set cut entering here for a cut with more than one cut-nodes!
                  //It is handled in TetMeshIntersection. Could possibly be handled more elegantly if information about
                  // the connectivity of the facet is provided here.
                  //  I.e. a tri has to be on the facet. In this case the triangulated point has to be added as a
                  //  cut-point information, and not as is now without any info.

//                  facet->PrintPointIds();
//
//                  std::cout << "l: " << std::endl;
//                  std::cout << "l->GetHandle()[0]: " << l->GetHandle()[0] << ", l->GetHandle()[1]: " << l->GetHandle()[1] << std::endl;
//
//                  std::cout << "t: " << std::endl;
//                  std::cout << "t->GetHandle()[0]: " << t->GetHandle()[0]<< ", t->GetHandle()[1]: " << t->GetHandle()[1] << ", t->GetHandle()[2]: " << t->GetHandle()[2] << std::endl;
//
//                  std::cout << "tri: " << std::endl;
//                  std::cout << "tri->GetHandle()[0]: " << tri->GetHandle()[0]<< ", tri->GetHandle()[1]: " << tri->GetHandle()[1] << ", tri->GetHandle()[2]: " << tri->GetHandle()[2] << std::endl;

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
                              bool tetcellsonly=false);

      void CreateVolumeCellTets();

      /// Return the list of tets.
      const std::vector<std::vector<int> > & Tets() const { return tets_; }

//       /// Return the tris for each facet
//       const std::map<Facet*, std::vector<Point*> > & SidesXYZ() const { return sides_xyz_; }

    private:

      /*!
      \brief Calculate the volume of the created delauney-cells and remove cells smaller than VOLUMETOL. Also switches
             the order in which nodes are stored if the tets are numbered the wrong way (deduced from taking the cross product).
       */
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

      /// Initialize valid tets and create children for the tet-cells (i.e. surfaces) and tet-surfaces (i.e. lines)
      ///   and the connectivity of these.
      void Init();

      void CallQHull( const std::vector<Point*> & points,
                      std::vector<std::vector<int> > & tets,
                      bool project );

      /// Somehow checks if tet is valid. Uses intersection for FindCommonSides(t,sides).
      bool IsValidTet( const std::vector<Point*> & t );

#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      double CalcVolumeOfTet(const std::vector<Point*> & t);
#endif

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

      /*!
        \brief Seed Domain for tets. This is done by finding the tets associated with the boundaries
                 of the element (i.e. the "parent" of the facet).
               Facets on the Cut surface are not seeded (unless forced).
               Tets which do not have a tri on a facet (tet in "interior") is not added in this function.
                 This is handled in Fill() function.

               1) Kuttler-way: Get the tris from the boundary-facet. Find it's associated tet.
                        If the tri is associated with only 1 tet -> add this tet to the cell_domain.

               2) Magnus-test: As problems appeared with disappearing integration cells in volume-cells,
                        it was noted that some tets were not added in cell_domain due to small perturbations
                        of the provided nodes. Thus the criterion that as long as the tri is associated to
                        ONE accepted tet it will be added in the cell_domain.
                        Furthermore, if a tri is not associated to any accepted tet but is still associated
                        to an UNIQUE tet (however small) it is added.
       */
      void SeedDomain( Domain<4> & cell_domain, Facet * f, bool force=false )
      {
        if ( force or not f->OnCutSide() )
        {
          FacetMesh & fm = facet_mesh_[f];
          const PlainEntitySet<3> & tris = fm.SurfaceTris();
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
          int tri_counter = 0;
#endif
          for ( PlainEntitySet<3>::const_iterator i=tris.begin();
                i!=tris.end();
                ++i )
          {
            Entity<3> * t = *i;
            const std::vector<Entity<4>*> & tets = t->Parents();
#ifdef NEW_SEED_DOMAIN
            Entity<4>* tet_accepted=tets[0];
            int possible_tets = tets.size();
            for(unsigned j=0; j<tets.size(); j++)
            {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
              std::cout << "tri_counter: " << tri_counter << std::endl;
              std::cout << "tets["   << j << "]->Id(): " << tets[j]->Id() << std::endl;
              std::cout << "tets["   << j << "](0):" << tets[j]->GetHandle()[0];
              std::cout << ", tets[" << j << "](1):" << tets[j]->GetHandle()[1];
              std::cout << ", tets[" << j << "](2):" << tets[j]->GetHandle()[2];
              std::cout << ", tets[" << j << "](3):" << tets[j]->GetHandle()[3];
              std::cout << std::endl;

              std::cout << "accept_tets_[tets[j]->Id()]: " << accept_tets_[tets[j]->Id()] << std::endl;
#endif
              if(not accept_tets_[tets[j]->Id()])
                possible_tets--;
              else
                tet_accepted=tets[j];
            }
#endif

#ifdef NEW_SEED_DOMAIN
            //There exists one ACCEPTED tet corresponding to the tri and, perhaps, other DISREGARDED
            //  tets which are not considered here.
            if(possible_tets==1)
            {
              if ( not cell_domain.Contains( tet_accepted ) )
              {
                cell_domain.Add( tet_accepted );
              }
            }
            //In the case a tri is created on the boundary of the element but is so small it
            //  is not an accepted tet.
            if(possible_tets==0)
            {
              if(tets.size()==1)
              {
                if ( not cell_domain.Contains( tets[0] ) )
                {
                  cell_domain.Add( tets[0] );
                }
              }
            }
#else
//            Kuettler-way
            if ( tets.size()==1 )
            {
              if ( not cell_domain.Contains( tets[0] ) )
                cell_domain.Add( tets[0] );
            }
#endif

//            //Test this to see what happens...
//            if(not tri_of_facet_given_VC)
//            {
//              std::cout << "tets.size(): " << tets.size() << std::endl;
//              for(unsigned j=0; j<tets.size(); j++)
//              {
//                std::cout << "tets["   << j << "]->Id(): " << tets[j]->Id() << std::endl;
//                std::cout << "tets["   << j << "](0):" << tets[j]->GetHandle()[0];
//                std::cout << ", tets[" << j << "](1):" << tets[j]->GetHandle()[1];
//                std::cout << ", tets[" << j << "](2):" << tets[j]->GetHandle()[2];
//                std::cout << ", tets[" << j << "](3):" << tets[j]->GetHandle()[3];
//                std::cout << std::endl;
//
//                std::cout << "accept_tets_[tets[j]->Id()]: " << accept_tets_[tets[j]->Id()] << std::endl;
//              }
//              if (tets.size() > 0 and not cell_domain.Contains( tets[0] ) )
//              {
//                cell_domain.Add( tets[0], false );
//              }
//              //throw std::runtime_error("A boundary tri is not assigned a tet-cell");
//            }

#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
            tri_counter++;
#endif

          }
#ifdef DEBUGCUTLIBRARY
          if(tris.size()==0)
          {
            std::cout << "WARNING: No tris found for this facet!" << std::endl;
          }
//          else if(not atleast_one_tet_in_VC)
//          {
//            std::cout << "WARNING: No tets were accepted from this facet!!! Possibly, Volume Cell too small?" << std::endl;
//          }
#endif
        }
      }

      //The cut points (i.e. nodes of the facets)
      const std::vector<Point*> & points_;
      //Facets of this element and cut-surface
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
