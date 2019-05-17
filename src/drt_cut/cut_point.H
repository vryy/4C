/*---------------------------------------------------------------------*/
/*!
\file cut_point.H

\brief Cut Point

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_POINT_H
#define CUT_POINT_H

#include <fstream>

#include <map>

#include <Teuchos_RCP.hpp>

#include "cut_clnwrapper.H"
#include "../linalg/linalg_fixedsizematrix.H"

#include "cut_utils.H"
#include "cut_tolerance.H"

namespace GEO
{
  namespace CUT
  {
    namespace IMPL
    {
      class PointLineFilter
      {
       public:
        virtual ~PointLineFilter() {}

        virtual bool operator()(Line* line) const = 0;
      };

    }  // namespace IMPL

    /// Points in a mesh. Both nodal points and cut points.
    /*!
     * A point knows its position with respect to the cut interface.
     */
    class Point
    {
     public:
      enum PointPosition
      {
        undecided = 0,
        oncutsurface = -1,
        inside = -2,
        outside = -3
      };

      /// translate PointPosition enumerator to string
      static inline const std::string PointPosition2String(const enum PointPosition& pos)
      {
        switch (pos)
        {
          case undecided:
            return "undecided";
          case oncutsurface:
            return "oncutsurface";
          case inside:
            return "inside";
          case outside:
            return "outside";
          default:
            break;
        }
        return "Unknown PointPosition enum";
      };

      static Point* NewPoint(Mesh& mesh, const double* x, double t, Edge* cut_edge, Side* cut_side,
          double tolerance = 0.0);

      static Point* InsertCut(Edge* cut_edge, Side* cut_side, Node* n);

      /// constructor
      Point(unsigned pid, double tolerance);

      /// destructor
      virtual ~Point(){};

      int Id() const { return pid_; }

      const double* X() const { return x_; }

      PointPosition Position() const { return position_; }

      void Position(PointPosition p);

      /*! \brief Add this edge to the list of edges that are cut by this point */
      void AddEdge(Edge* e);

      /*! \brief Add this side to the list of sides that are cut by this point */
      void AddSide(Side* s);

      /*! \brief Add this element to the list of elements that are cut by this point */
      void AddElement(Element* e)
      {
        if (e != NULL)
        {
          cut_elements_.insert(e);
        }
      }

      /// removing element from the list of elements cut by this point (can be used when position
      /// was determined wrongly)
      void RemoveElement(Element* e) { cut_elements_.erase(e); }

      /// adding pair of side-edge to container of intersected, does not insert if such pair already
      /// exists
      void AddPair(Side* side, Edge* edge, const std::pair<Side*, Edge*>& or_pair);

      /// default call
      void AddPair(Side* side, Edge* edge);

      /// called during merging
      void AddPair(const std::pair<Side*, Edge*>& pair,
          const std::pair<Side*, Edge*>& original_pair, Point* source);

      /// Add all topological connections of intersection of this edge and other edge ( all
      /// necessery pairs, etc)
      void AddEdgeIntersection(Edge* first, Edge* second,
          const std::pair<Side*, Edge*>& original_cut_pair, const std::string& extra_msg = "");

      void AddEdgeIntersection(Edge* first, Edge* second, Side* original_side, Edge* original_edge,
          const std::string& extra_msg = "");

      /// returns true if this point is created by this pair of side and edge intersection
      bool IsCut(Side* side, Edge* edge);

      /*! \brief Returns true if the edge is cut by this point */
      bool IsCut(Edge* e) { return cut_edges_.count(e) > 0; }

      /*! \brief Returns true if the facet is cut by this point */
      bool IsCut(Facet* f) { return facets_.count(f) > 0; }

      /*! \brief Returns true if the side is cut by this point */
      bool IsCut(Side* s) { return cut_sides_.count(s) > 0; }

      /*! \brief Returns true if this point cuts the element */
      bool IsCut(Element* s) { return cut_elements_.count(s) > 0; }

      /// Checks is this cut point is the result of s1 x s2 by checking cut_pairs
      bool IsCut(Side* s1, Side* s2);

      /*
       * Return if this point associated with a facet which is either a cut-side or a marked side
       *  Needed for the tessellation.
       *  TODO: Look into this
       *  winter 02/2017
       */
      bool HasAssociatedBoundaryCellFacet();

      /*! \brief Get the coordinates of this point */
      virtual void Coordinates(double* x) const = 0;

      void Register(Line* line) { lines_.insert(line); }

      void Register(Facet* facet) { facets_.insert(facet); }

      /*!
      \brief Identifies the edges that are cut by considered point and given point
      */
      void CommonEdge(Point* other, plain_edge_set& edges);

      /*!
      \brief Identifies the sides that are cut by considered point and given point
      */
      void CommonSide(Point* other, plain_side_set& sides);

      Line* CommonLine(Point* other);

      Line* CutLine(Side* side, bool unique = true);

      Line* CutLine(Line* line, Side* side, Element* element);

      Line* CutLine(const IMPL::PointLineFilter& filter, bool unique = true);

      Line* CutLine(Line* line, const IMPL::PointLineFilter& filter, bool unique = true);

      void CutLines(const IMPL::PointLineFilter& filter, plain_line_set& cut_lines);

      void CutLines(Side* side, plain_line_set& cut_lines);

      // std::vector<Edge*> CutEdges( Point * other );

      const plain_edge_set& CutEdges() { return cut_edges_; }

      const plain_side_set& CutSides() { return cut_sides_; }

      const std::set<std::pair<Side*, Edge*>>& CutPairs() { return cut_pairs_; }

      Side* CutSide(Side* side, Point* other);

      void CutEdge(Side* side, Line* other_line, std::vector<Edge*>& matches);

      void Print(std::ostream& stream = std::cout) const
      {
        stream << "(" << pid_ << "; " << std::setprecision(16) << x_[0] << ","
               << std::setprecision(16) << x_[1] << "," << std::setprecision(16) << x_[2] << ")";
      }

      void Plot(std::ostream& f, int nid = -50) const
      {
        f << std::setprecision(25) << x_[0] << " " << std::setprecision(25) << x_[1] << " "
          << std::setprecision(25) << x_[2] << " # " << pid_ << " " << nid;
        DumpDoubles(f, X(), 3);
        f << "\n";
      }

      /** \brief check if the given point lies on the edge and return the
       *  corresponding local coordinate of this point on the edge
       *
       *  Return the parameter space coordinate \f$t\f$ of this point \f$x\f$ on
       *  the given \c edge with begin node \f$\bar{x}_{0}\f$ and the end point
       *  \f$\bar{x}_{1}\f$:
       *
       *  x = \sum\limits_{i=0}^{1} N_{i}(t) \bar{x}_{i}
       *
       *  \param edge (in) : edge which we want to check */
      double t(Edge* edge);

      /// Same as above.
      /// But it can happen that origanl intersection point is merged into different location.
      /// However, we want to obtain position of the cut_point
      /// on the edge based on the real intersection( before mergin),to be consistest with rest of
      /// the code
      double t(Edge* edge, const LINALG::Matrix<3, 1>& coord);

      // void t( Edge* edge, double pos ) { t_[edge] = pos; }

      unsigned Pid() const { return pid_; }

      bool NodalPoint(const std::vector<Node*>& nodes) const;

      // checking if the cut_point is close enough to the nodal point
      bool AlmostNodalPoint(const std::vector<Node*>& nodes, double tolerance = 1e-10) const;

      Node* CutNode();

      void Intersection(plain_edge_set& edges);

      void Intersection(plain_side_set& sides);

      void Intersection(plain_facet_set& facets);

      void Intersection(plain_element_set& elements);

      const plain_element_set& Elements() const { return cut_elements_; }

      const plain_facet_set& Facets() const { return facets_; }

      Edge* CommonCutEdge(Side* side);

      /// erase all cut pairs containing this side
      void ErasedContainingCutPairs(Side* cutside);

      /// erased all cut pairs containing this edge
      void ErasedContainingCutPairs(Edge* cutsideedge);
      /// Erase the cutside from this point because it is deleted in the selfcut
      void EraseCutSide(Side* cutside) { cut_sides_.erase(cutside); }

      /// Erase the cutsideedge from this point because it is deleted in the selfcut
      void EraseCutSideEdge(Edge* cutsideedge) { cut_edges_.erase(cutsideedge); }

      /// Return actual point tolerance
      double Tolerance() { return tol_; }

      /// Set new Tolerance (in case points are merged this is necessary)
      void EnlargeTolerance(double newtol)
      {
        if (newtol > tol_)
        {
          tol_ = newtol;
        }
      }

      /// Removes information that this point was created by cut from everywhere (e.g. sides,
      /// edges, element)
      void RemoveConnectivityInfo();

      /// Carefully replace this points by another at the stage before/during creation of cut_lines
      void Replace(Point* p);

      void RemoveEdge(Edge* edge);

      void RemoveSide(Side* side);

      void RemovePair(std::pair<Side*, Edge*>& pair) { cut_pairs_.erase(pair); };

      // move point slighlty (used for merging)
      virtual void MovePoint(const double* new_coord) = 0;

      // Merge point into anothe r
      void Merge(Point* dest);

      /// Dump all information of how this point was created, edges x sides intersections. Useful to
      /// debug reasons crashing
      void DumpConnectivityInfo();

      /// Get the Merged Points
      const std::vector<Point*> GetMergedPoints() { return merged_points_; };

#if CUT_CREATION_INFO
      // Add merged pair. If the coord = NULL, it means it was inserted without actual computation
      // of the coordinate, solely on topology e.g. in the InsertCut
      void AddMergedPair(const std::pair<Side*, Edge*>& inter, LINALG::Matrix<3, 1>* coord);

      const std::map<std::pair<Side*, Edge*>, std::pair<LINALG::Matrix<3, 1>, bool>>&
      GetMergedPointsInfo()
      {
        return merged_points_info_;
      };

      // Add info about based on which side and edge this point was created (with optional message)
      void AddCreationInfo(const std::pair<Side*, Edge*>& cut_pair, const std::string& info);

      // Add creation info, that does not depended on particular cut pair
      void AddAdditionalCreationInfo(const std::string& info);

      const std::string& GetAdditionalCreationInfo();

      // returns empty string if there is no info for this pair
      const std::string& GetCreationInfo(const std::pair<Side*, Edge*>& pair)
      {
        return creation_info_[pair];
      };


      // for particular cut pair gives info if it is really cut or added later on during topology
      // identification (e.g. because of neigboring edges )
      bool IsReallyCut(const std::pair<Side*, Edge*>& pair)
      {
        std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>>::iterator it =
            cut_pairs_info_.find(pair);
        if (it != cut_pairs_info_.end())
        {
          std::pair<Side*, Edge*> found_pair = it->second.first;
          return (found_pair == it->first);
        }
        else
        {
          dserror("Following cut_pair does not exist!");
        }
        return false;
      }

      const std::pair<Side*, Edge*>& AddedFrom(const std::pair<Side*, Edge*>& resultant_pair)
      {
        std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>>::iterator it =
            cut_pairs_info_.find(resultant_pair);
        if (it != cut_pairs_info_.end())
        {
          return (it->second.first);
        }
        else
        {
          dserror("Following cut pair does not exists!");
        }
        // should not reach here, just to make compiler happy
        return cut_pairs_info_.begin()->second.first;
      }

#endif

     protected:
      // internal stored point coordinates (always of dimension 3!)
      double x_[3];

     private:
      // find the unique line that matches the filter criteria
      template <class Filter>
      Line* Find(Filter& filter, bool unique = true);

      template <class Filter>
      void Find(Filter& filter, plain_line_set& cut_lines);

      unsigned pid_;
      PointPosition position_;

      std::map<Edge*, double> t_;

      // precicion of the point
      double tol_;

      plain_edge_set cut_edges_;
      plain_side_set cut_sides_;  // a set with all sides
      plain_element_set cut_elements_;

      plain_line_set lines_;
      plain_facet_set facets_;

      // stores cut pairs of side and edge intersections for this point in associative containers
      // for easy lookup and insert
      std::set<std::pair<Side*, Edge*>> cut_pairs_;

      /// points that were merged here during replacement
      std::vector<Point*> merged_points_;

#if CUT_CREATION_INFO

      // Map that indicates whether cut pairs whether originally (during point creation) or later
      // (during neighbors identification), key is the pair that was added, value is the original
      // pair that result in this add. if they are same it was added directly in intersection, if
      // different - during topology identification Point represent origin source point, that was
      // merged into this one
      std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>> cut_pairs_info_;

      // Map that indicates what intersections (keys of the map)  merged coordinates in the the
      // first part of value pair to the Point* of second part of value parit.Indicates what
      // coordinates were merged in this particular point, due to closeness
      std::map<std::pair<Side*, Edge*>, std::pair<LINALG::Matrix<3, 1>, bool>> merged_points_info_;

      // Displays why particular cut pair is intersected in this point
      std::map<std::pair<Side*, Edge*>, std::string> creation_info_;

      // additional creation info that might be cut pair independent
      std::string additional_creation_info_;

      // Indicates actuall real points, that were merged here (due to mergin in pointpool)
      std::vector<Point*> real_merged_points_;

      // point that this point got merged into
      Point* merged_to_;

#endif
    };  // class Point

    /*--------------------------------------------------------------------------*/
    /* \class ConcretePoint
     *
     * Necessary due to the different problem dimensions. Actually only the
     * constructor and the Coordinates() function change, since the stored
     * point coordinate stays 3. If probDim is smaller than three, the last
     * coordinates are set to zero.
     *
     * \author hiermeier
     * \date 11/16 */
    template <unsigned probDim>
    class ConcretePoint : public Point
    {
     public:
      /// constructor
      ConcretePoint(unsigned pid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
          : Point(pid, tolerance)
      {
        // copy the given coordinates into the base class member variable
        std::copy(x, x + probDim, this->x_);

        // set the remaining entries to zero
        std::fill(this->x_ + probDim, this->x_ + 3, 0.0);

        // The x_ coordinate must be set, before these methods are called!
        if (cut_edge != NULL)
        {
          AddEdge(cut_edge);
        }
        if (cut_side != NULL)
        {
          AddSide(cut_side);
        }

        if ((cut_side != NULL) && (cut_edge != NULL)) this->AddPair(cut_side, cut_edge);
      };

      /// destructor
      virtual ~ConcretePoint(){};

      /** \brief Get the coordinates of this point */
      virtual void Coordinates(double* x) const { std::copy(this->x_, this->x_ + probDim, x); };

      virtual void MovePoint(const double* new_coord);

    };  // class ConcretePoint

    /*--------------------------------------------------------------------------*/
    /** \class PointFactory
     *
     *  This class creates a ConcretePoint object with correct problem dimension.
     *  Don't call this class directly! Use the non-member function instead.
     *
     *  \date 11/16
     *  \author hiermeier */
    class PointFactory
    {
     public:
      PointFactory(){};

      // non-member function to create a concrete point of desired dimension
      Teuchos::RCP<GEO::CUT::Point> CreatePoint(unsigned pid, const double* x, Edge* cut_edge,
          Side* cut_side, double tolerance, int probdim) const
      {
        Teuchos::RCP<GEO::CUT::Point> point = Teuchos::null;

        switch (probdim)
        {
          case 2:
            point = Teuchos::rcp(new ConcretePoint<2>(pid, x, cut_edge, cut_side, tolerance));
            break;
          case 3:
            point = Teuchos::rcp(new ConcretePoint<3>(pid, x, cut_edge, cut_side, tolerance));
            break;
          default:
            dserror("Unsupported problem dimension! (probdim=%d)", probdim);
            break;
        }
        return point;
      };
    };  // class PointFactory

    // non-member function to create a concrete point of desired dimension
    Teuchos::RCP<GEO::CUT::Point> CreatePoint(
        unsigned pid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance);

    inline int EntityId(const Point& p) { return p.Pid(); }

    /// id based comparison for sorting and searching
    template <class T>
    class EntityIdLess
    {
     public:
      bool operator()(const T& p1, const T& p2) const { return EntityId(p1) < EntityId(p2); }

      bool operator()(const T* p1, const T* p2) const { return EntityId(*p1) < EntityId(*p2); }

      bool operator()(const Teuchos::RCP<T> p1, const Teuchos::RCP<T> p2) const
      {
        return EntityId(*p1) < EntityId(*p2);
      }
    };

    typedef EntityIdLess<Point> PointPidLess;

    /// position on edge based comparison for sorting and searching
    class PointPositionLess
    {
     public:
      explicit PointPositionLess(Edge* edge) : edge_(edge) {}

      bool operator()(Point& p1, Point& p2) const
      {
        if (not p1.IsCut(edge_) or not p2.IsCut(edge_))
          throw std::runtime_error("point position compare only on cut edges");
        return p1.t(edge_) < p2.t(edge_);
      }

      bool operator()(Point* p1, Point* p2) const
      {
        if (not p1->IsCut(edge_) or not p2->IsCut(edge_))
          throw std::runtime_error("point position compare only on cut edges");
        return p1->t(edge_) < p2->t(edge_);
      }

     private:
      Edge* edge_;
    };

#ifdef CUT_USE_SORTED_VECTOR
    typedef sorted_vector<Point*, true, PointPidLess> PointSet;
    typedef sorted_vector<Point*, true, PointPositionLess> PointPositionSet;

    typedef sorted_vector<Teuchos::RCP<Point>, true, PointPidLess> RCPPointSet;

#else
    typedef std::set<Point*, PointPidLess> PointSet;
    typedef std::set<Point*, PointPositionLess> PointPositionSet;

    typedef std::set<Teuchos::RCP<Point>, PointPidLess> RCPPointSet;

#endif

    namespace IMPL
    {
      class ExcludeLineFilter : public PointLineFilter
      {
       public:
        ExcludeLineFilter(Line* line, const PointLineFilter& filter) : line_(line), filter_(filter)
        {
        }

        virtual bool operator()(Line* line) const { return line != line_ and filter_(line); }

       private:
        Line* line_;
        const PointLineFilter& filter_;
      };

      class SideCutFilter : public PointLineFilter
      {
       public:
        SideCutFilter(Side* side) : side_(side) {}

        virtual bool operator()(Line* line) const;

       private:
        Side* side_;
      };

      class SideElementCutFilter : public PointLineFilter
      {
       public:
        SideElementCutFilter(Side* side, Element* element) : side_(side), element_(element) {}

        virtual bool operator()(Line* line) const;

       private:
        Side* side_;
        Element* element_;
      };

      class SideSideCutFilter : public PointLineFilter
      {
       public:
        SideSideCutFilter(Side* side1, Side* side2) : side1_(side1), side2_(side2) {}

        virtual bool operator()(Line* line) const;

       private:
        Side* side1_;
        Side* side2_;
      };

    }  // namespace IMPL

    /** \brief Find if the points in element share a common element,
     *  (i.e. do all points lie on the same element?)
     *
     *  A point knows which cut_elements it is associated to.
     *
     *  \author hiermeier \date 11/16 */
    void FindCommonElements(const std::vector<Point*>& element, plain_element_set& elements);

    /// Find if the points in side (which is actually a tet...) share a common side,
    ///   (i.e. do all points lie on the same surface?)
    /// A point knows which cut_sides it is associated to.
    inline void FindCommonSides(const std::vector<Point*>& side, plain_side_set& sides)
    {
      std::vector<Point*>::const_iterator is = side.begin();
      // Get the sides this point cuts
      // A cut-side which is a LevelSet-side is empty but is still present.
      sides = (*is)->CutSides();
      for (++is; is != side.end(); ++is)
      {
        Point* p = *is;
        p->Intersection(sides);
        if (sides.size() == 0)
        {
          break;
        }
      }
    }

    inline void FindCommonSides(Point* p1, Point* p2, Point* p3, plain_side_set& sides)
    {
      sides = p1->CutSides();
      p2->Intersection(sides);
      p3->Intersection(sides);
    }

    inline void FindCommonSides(Point* p1, Point* p2, Point* p3, Point* p4, plain_side_set& sides)
    {
      sides = p1->CutSides();
      p2->Intersection(sides);
      p3->Intersection(sides);
      p4->Intersection(sides);
    }

    /// Find distance between points
    template <unsigned int probDim>
    double DistanceBetweenPoints(
        const LINALG::Matrix<probDim, 1>& coord_a, const LINALG::Matrix<probDim, 1>& coord_b)
    {
      LINALG::Matrix<probDim, 1> diff;
      diff.Update(1, coord_a, -1, coord_b);
      return diff.Norm2();
    }

    /// Find distance between points
    double DistanceBetweenPoints(Point* p1, Point* p2);

    /// Find distance between points
    double DistanceBetweenPoints(Point* p1, const LINALG::Matrix<3, 1>& coord_b);

  }  // namespace CUT
}  // namespace GEO


std::ostream& operator<<(std::ostream& stream, GEO::CUT::Point& point);

std::ostream& operator<<(std::ostream& stream, GEO::CUT::Point* point);

#endif
