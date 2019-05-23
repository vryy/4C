/*!-----------------------------------------------------------------------------------------------*

\brief class representing a geometrical edge

\level 2

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_EDGE_H
#define CUT_EDGE_H

#include "../drt_inpar/inpar_cut.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "cut_node.H"
#include "cut_tolerance.H"

namespace GEO
{
  namespace CUT
  {
    class Mesh;
    class Side;
    class IntersectionBase;

    /*--------------------------------------------------------------------------*/
    /*! \brief Linear edge between two nodes. The edge nodes are always cut points.
     *  There can be further cut points on the edge. */
    class Edge
    {
     public:
      /** \brief Create a new concrete edge object
       *
       *  \param edgetype (in) : element type of the edge
       *  \param nodes    (in) : vector of nodes defining the new edge
       *
       *  \author hiermeier \date 08/16 */
      static Teuchos::RCP<Edge> Create(
          DRT::Element::DiscretizationType edgetype, const std::vector<Node*>& nodes);

      /** \brief Create a new concrete edge object
       *
       *  \param shardskey (in) : unique key equivalent to a element type ( see TRILINOS library )
       *  \param nodes     (in) : vector of nodes defining the new edge
       *
       *  \author hiermeier \date 08/16 */
      static Teuchos::RCP<Edge> Create(unsigned shardskey, const std::vector<Node*>& nodes);

      /// constructor
      explicit Edge(const std::vector<Node*>& nodes)
          : nodes_(nodes), cut_points_(PointPositionLess(this))
      {
        for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
        {
          Node* n = *i;
          n->Register(this);
          n->point()->AddEdge(this);
        }
        selfcutposition_ = Point::undecided;
#if CUT_CREATION_INFO
        std::stringstream id;
        id << BeginNode()->point()->Id() << "->" << EndNode()->point()->Id();
        id_ = id.str();
#endif
      }

      virtual ~Edge() {}

      virtual unsigned ProbDim() const = 0;

      virtual unsigned Dim() const = 0;

      virtual unsigned NumNodes() const = 0;

      virtual DRT::Element::DiscretizationType Shape() const = 0;

      /*! \brief Add the side to the list of sides cut by this edge */
      void Register(Side* side)
      {
        sides_.insert(side);
        for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
        {
          Node* n = *i;
          n->Register(side);
        }
      }

      /*! \brief Check whether the edge is part of this side */
      bool AtSide(Side* side) { return sides_.count(side) > 0; }

      /*! \brief Check whether the edges defined by these two points coincide with
       *  this edge */
      bool Matches(Point* begin, Point* end)
      {
        return ((BeginNode()->point() == begin and EndNode()->point() == end) or
                (BeginNode()->point() == end and EndNode()->point() == begin));
      }

      // virtual bool HasPoint( Point * p ) = 0;

      /*! \brief Get all the sides on which this edge is a part of */
      const plain_side_set& Sides() { return sides_; }

      virtual void GetTouchingPoints(const std::vector<Node*>& nodes, std::vector<Node*>& points,
          INPAR::CUT::CUT_Floattype floattype = INPAR::CUT::floattype_double) = 0;

      /*! \brief Get the intersection points of this edge with the given side and
       * store the cut points in cuts */
      virtual bool Cut(Mesh& mesh, Side& side, PointSet& cuts) = 0;

      /*! Tries to compute side x edge intersection, if edge was parallel to side (using only
       * ComputeDistance and Edge-edge intersection), no real intersection */
      virtual bool JustParallelCut(Mesh& mesh, Side& side, PointSet& cuts, int skip_id = -1) = 0;

      /*! \brief Get the intersection points of this edge with the level-set side
       *  and store the cut points in cuts */
      virtual bool LevelSetCut(Mesh& mesh, Side& side, PointSet& cuts) = 0;

      /*! \brief Add cut_point to the list of points at which this edge cuts some
       *  other elements */
      void AddPoint(Point* cut_point);

      /*! \brief Get the coordinates of the nodes of the edge */
      virtual void Coordinates(double* xyze) = 0;

      template <class T>
      void Coordinates(T& xyze_lineElement)
      {
        if (static_cast<unsigned>(xyze_lineElement.M()) != ProbDim())
          dserror("xyze_lineElement has the wrong number of rows! (probdim = %d)", ProbDim());
        if (static_cast<unsigned>(xyze_lineElement.N()) != NumNodes())
          dserror("xyze_lineElement has the wrong number of columns! (dim = %d)", NumNodes());

        Coordinates(xyze_lineElement.A());
      }

      /*! \brief Print the coordinates of the nodes on screen */
      void Print()
      {
        nodes_[0]->Print();
        for (unsigned i = 1; i < nodes_.size(); ++i)
        {
          std::cout << "--";
          nodes_[i]->Print();
        }
      }

      void Plot(std::ofstream& f)
      {
        f << "# edge\n";
        BeginNode()->Plot(f);
        if (nodes_.size() == 3) MiddleNode()->Plot(f);
        EndNode()->Plot(f);
        f << "\n\n";
      }

      /*! \brief Get the cut points on the edge defined by the edge_start and
       *  edge_end nodes */
      void CutPoint(Node* edge_start, Node* edge_end, std::vector<Point*>& edge_points);

      /*! \brief Unused */
      void CutPoints(Side* side, PointSet& cut_points);

      /*! \brief Unused */
      void CutPointsBetween(Point* begin, Point* end, std::vector<Point*>& line);

      void CutPointsIncluding(Point* begin, Point* end, std::vector<Point*>& line);

      void CutPointsInside(Element* element, std::vector<Point*>& line);

      bool IsCut(Side* side);

      /** \brief Compute the intersection between THIS edge and a \c other
       *  edge.
       *
       *  \param mesh       (in)  : cut mesh, necessary to add cut points to the point pool
       *  \param other      (in)  : pointer to the other edge object
       *  \param side       (in)  : pointer to the other side object(can be NULL)
       *  \param cut_points (out) : if the cut was successful, this set will hold the new cut points
       *  \param tolerance  (out) : tolerance used for the intersection ( defined by the cut kernel
       * )
       *
       *  This routine returns TRUE, if the computation was successful AND the cut
       *  point is within the edge limits.
       *
       *  \author hiermeier \date 12/16 */
      virtual bool ComputeCut(
          Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance) = 0;

      Node* BeginNode() { return nodes_.front(); }

      Node* MiddleNode()
      {
        if (nodes_.size() != 3) throw std::runtime_error("middle node in line3 only");
        return nodes_[2];
      }

      Node* EndNode() { return nodes_[1]; }

      Point* NodeInElement(Element* element, Point* other);

      const std::vector<Node*>& Nodes() const { return nodes_; }

      /// Find common points (excluding cut_points points) between two edges
      void CommonNodalPoints(Edge* edge, std::vector<Point*>& nodes);

      bool FindCutPoints(Mesh& mesh, Element* element, Side& side, Side& other);

      /*!
      \brief Computes the points at which both the sides intersect
       */
      bool FindCutPointsMeshCut(
          Mesh& mesh, Element* element, Side& side, Side& other, PointSet* cutpoints = NULL);

      /*!
      \brief Simplified version of the FindCutPoints as used by the function LevelSetCut
       */
      bool FindCutPointsLevelSet(Mesh& mesh, Element* element, Side& side, Side& other);

      /*!
      \brief Cut points falling on this edge that are common to the two given sides are extracted
       */
      void GetCutPoints(Element* element, Side& side, Side& other, PointSet& cuts);

      /*!
      \brief Cut points falling on this edge that are common to the given edge are extracted
       */
      void GetCutPoints(Edge* other, PointSet& cuts);

      const PointPositionSet& CutPoints() const { return cut_points_; }
      // const std::vector<Point*> & CutPoints() const { return cut_points_; }

      void RectifyCutNumerics();

      /// Returns the selfcutposition of this edge
      Point::PointPosition SelfCutPosition() { return selfcutposition_; }

      /// Gives this edge a selfcutposition and spreads the positional information
      void SelfCutPosition(Point::PointPosition p);

      /// Changes the selfcutposition of this edge and spreads the positional information
      void ChangeSelfCutPosition(Point::PointPosition p);

      /// Erase the cutside from this edge because it is deleted in the selfcut
      void EraseCutSide(Side* cutside) { sides_.erase(cutside); }

      /*!
      \brief Replace the node "nod" of this edge with given node "replwith"
       */
      void replaceNode(Node* nod, Node* replwith);

      // remove information about this point
      void RemovePoint(Point* p) { cut_points_.erase(p); };

      /*!
       \brief Add all topological connections of interesection of this edge and other edge ( all
       necessery pairs, etc)
       */
      void AddConnections(Edge* other, const std::pair<Side*, Edge*>& original_cut_pair);

      void AddConnections(Edge* other, Side* original_side, Edge* original_edge);

#if CUT_CREATION_INFO
      std::string Id() { return id_; };
#endif

     private:
#if CUT_CREATION_INFO
      std::string id_;
#endif

      std::vector<Node*> nodes_;

      plain_side_set sides_;

      //! sorted vector contains all points (end points and cut points) on this edge
      PointPositionSet cut_points_;

      //! every cutsideedge knows its selfcutposition
      Point::PointPosition selfcutposition_;

    };  // class Edge

    /*--------------------------------------------------------------------------*/
    template <unsigned probDim, DRT::Element::DiscretizationType edgeType,
        unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
        unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement>
    class ConcreteEdge : public Edge
    {
     public:
      /// constructor
      ConcreteEdge(const std::vector<Node*>& nodes) : Edge(nodes) {}

      /// destructor
      virtual ~ConcreteEdge(){};

      /// get the element dimension of this edge
      virtual unsigned Dim() const { return dimEdge; }

      /// get the number of nodes of this edge
      virtual unsigned NumNodes() const { return numNodesEdge; }

      /// get the problem dimension
      virtual unsigned ProbDim() const { return probDim; }

      /// get the shape of this edge element
      virtual DRT::Element::DiscretizationType Shape() const { return edgeType; }

      /*! \brief Get the coordinates of the nodes of the side */
      virtual void Coordinates(double* xyze)
      {
        double* x = xyze;
        for (std::vector<Node*>::const_iterator i = Nodes().begin(); i != Nodes().end(); ++i)
        {
          const Node& n = **i;
          n.Coordinates(x);
          x += probDim;
        }
      }

      virtual void GetTouchingPoints(const std::vector<Node*>& nodes, std::vector<Node*>& points,
          INPAR::CUT::CUT_Floattype floattype = INPAR::CUT::floattype_double);

      /*! \brief Handles intersection of two edges that are close to each other */
      virtual bool HandleParallelCut(Edge* other, Side* side, PointSet* cut_points,
          INPAR::CUT::CUT_Floattype floattype = INPAR::CUT::floattype_double);

      virtual bool JustParallelCut(Mesh& mesh, Side& side, PointSet& cuts, int skip_id = -1);

      /*! \brief Get the intersection points of this edge with the given side and
       * store the cut points in cuts */
      virtual bool Cut(Mesh& mesh, Side& side, PointSet& cuts);

      /*! \brief Get the intersection point of THIS edge with the given \c other edge and
       *  store the global cut point in x and the local coordinate of THIS edge in \c pos.
       *
       *  \param mesh       (in)  : cut mesh, necessary to add cut points to the point pool
       *  \param other      (in)  : pointer to the other edge object
       *  \param side       (in)  : pointer to the other side object(can be NULL)
       *  \param cut_points (out) : if the cut was successful, this set will hold the new cut points
       *  \param tolerance  (out) : tolerance used for the intersection ( defined by the cut kernel
       * )
       *
       *  This routine returns TRUE, if the computation was successful AND the cut
       *  point is within the edge limits.
       *
       *  \author hiermeier \date 12/16 */
      virtual bool ComputeCut(
          Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance);

      /*! \brief Get the intersection points of this edge with the level-set side
       *  and store the cut points in cuts */
      virtual bool LevelSetCut(Mesh& mesh, Side& side, PointSet& cuts)
      {
        double blsv = BeginNode()->LSV();
        double elsv = EndNode()->LSV();

        bool cutfound = false;
        // version for single element cuts, here we need to watch for tolerances on
        // nodal cuts
        if (std::abs(blsv) <= REFERENCETOL)
        {
          cuts.insert(Point::InsertCut(this, &side, BeginNode()));
          cutfound = true;
        }
        if (std::abs(elsv) <= REFERENCETOL)
        {
          cuts.insert(Point::InsertCut(this, &side, EndNode()));
          cutfound = true;
        }

        if (not cutfound)
        {
          if ((blsv < 0.0 and elsv > 0.0) or (blsv > 0.0 and elsv < 0.0))
          {
            cutfound = true;
            double z = blsv / (blsv - elsv);

            LINALG::Matrix<probDim, 1> x1;
            LINALG::Matrix<probDim, 1> x2;
            BeginNode()->Coordinates(x1.A());
            EndNode()->Coordinates(x2.A());

            LINALG::Matrix<probDim, 1> x;
            x.Update(-1., x1, 1., x2, 0.);
            x.Update(1., x1, z);
            Point* p = Point::NewPoint(mesh, x.A(), 2. * z - 1., this, &side, 0.0);
            cuts.insert(p);
          }
        }

        //      std::cout <<"LS cut found? --> " << (cutfound ? "TRUE" : "FALSE") << std::endl;
        //      dserror("GEO::CUT::Edge::LevelSetCut -- STOP -- hiermeier 08/16");

        return cutfound;
      }

     private:
      Teuchos::RCP<GEO::CUT::IntersectionBase> IntersectionPtr(
          const DRT::Element::DiscretizationType& sidetype) const;

    };  // class ConcreteEdge

    /*--------------------------------------------------------------------------*/
    class EdgeFactory
    {
     public:
      /// constructor
      EdgeFactory(){};

      /// destructor
      virtual ~EdgeFactory(){};

      /// create the concrete edge object
      Teuchos::RCP<Edge> CreateEdge(
          const DRT::Element::DiscretizationType& edgetype, const std::vector<Node*>& nodes) const;

     private:
      template <DRT::Element::DiscretizationType edgetype>
      Edge* CreateConcreteEdge(const std::vector<Node*>& nodes, int probdim) const
      {
        Edge* e = NULL;
        switch (probdim)
        {
          case 2:
            e = new ConcreteEdge<2, edgetype>(nodes);
            break;
          case 3:
            e = new ConcreteEdge<3, edgetype>(nodes);
            break;
          default:
            dserror("Unsupported problem dimension! (probdim=%d)", probdim);
            break;
        }
        return e;
      }
    };  // class EdgeFactory
  }     // namespace CUT
}  // namespace GEO

#endif /* CUT_EDGE_H */
