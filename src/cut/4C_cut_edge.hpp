/*----------------------------------------------------------------------*/
/*! \file

\brief class representing a geometrical edge

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_EDGE_HPP
#define FOUR_C_CUT_EDGE_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"
#include "4C_cut_node.hpp"
#include "4C_cut_tolerance.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
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
      static Teuchos::RCP<Edge> create(
          Core::FE::CellType edgetype, const std::vector<Node*>& nodes);

      /** \brief Create a new concrete edge object
       *
       *  \param shardskey (in) : unique key equivalent to a element type ( see TRILINOS library )
       *  \param nodes     (in) : vector of nodes defining the new edge
       *
       *  \author hiermeier \date 08/16 */
      static Teuchos::RCP<Edge> create(unsigned shardskey, const std::vector<Node*>& nodes);

      /// constructor
      explicit Edge(const std::vector<Node*>& nodes)
          : nodes_(nodes), cut_points_(PointPositionLess(this))
      {
        for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
        {
          Node* n = *i;
          n->Register(this);
          n->point()->add_edge(this);
        }
        selfcutposition_ = Point::undecided;
#if CUT_CREATION_INFO
        std::stringstream id;
        id << BeginNode()->point()->Id() << "->" << EndNode()->point()->Id();
        id_ = id.str();
#endif
      }

      virtual ~Edge() = default;
      virtual unsigned n_prob_dim() const = 0;

      virtual unsigned n_dim() const = 0;

      virtual unsigned num_nodes() const = 0;

      virtual Core::FE::CellType shape() const = 0;

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
      bool at_side(Side* side) { return sides_.count(side) > 0; }

      /*! \brief Check whether the edges defined by these two points coincide with
       *  this edge */
      bool matches(Point* begin, Point* end)
      {
        return ((begin_node()->point() == begin and end_node()->point() == end) or
                (begin_node()->point() == end and end_node()->point() == begin));
      }

      // virtual bool HasPoint( Point * p ) = 0;

      /*! \brief Get all the sides on which this edge is a part of */
      const plain_side_set& sides() { return sides_; }

      virtual void get_touching_points(const std::vector<Node*>& nodes, std::vector<Node*>& points,
          CutFloatType floattype = floattype_double) = 0;

      /*! \brief Get the intersection points of this edge with the given side and
       * store the cut points in cuts */
      virtual bool cut(Mesh& mesh, Side& side, PointSet& cuts) = 0;

      /*! Tries to compute side x edge intersection, if edge was parallel to side (using only
       * compute_distance and Edge-edge intersection), no real intersection */
      virtual bool just_parallel_cut(Mesh& mesh, Side& side, PointSet& cuts, int skip_id = -1) = 0;

      /*! \brief Get the intersection points of this edge with the level-set side
       *  and store the cut points in cuts */
      virtual bool level_set_cut(Mesh& mesh, Side& side, PointSet& cuts) = 0;

      /*! \brief Add cut_point to the list of points at which this edge cuts some
       *  other elements */
      void add_point(Point* cut_point);

      /*! \brief Get the coordinates of the nodes of the edge */
      virtual void coordinates(double* xyze) = 0;

      template <class T>
      void coordinates(T& xyze_lineElement)
      {
        if (static_cast<unsigned>(xyze_lineElement.num_rows()) != n_prob_dim())
          FOUR_C_THROW(
              "xyze_lineElement has the wrong number of rows! (probdim = %d)", n_prob_dim());
        if (static_cast<unsigned>(xyze_lineElement.num_cols()) != num_nodes())
          FOUR_C_THROW("xyze_lineElement has the wrong number of columns! (dim = %d)", num_nodes());

        coordinates(xyze_lineElement.values());
      }

      /*! \brief Print the coordinates of the nodes on screen */
      void print()
      {
        nodes_[0]->print();
        for (unsigned i = 1; i < nodes_.size(); ++i)
        {
          std::cout << "--";
          nodes_[i]->print();
        }
      }

      void plot(std::ofstream& f)
      {
        f << "# edge\n";
        begin_node()->plot(f);
        if (nodes_.size() == 3) middle_node()->plot(f);
        end_node()->plot(f);
        f << "\n\n";
      }

      /*! \brief Get the cut points on the edge defined by the edge_start and
       *  edge_end nodes */
      void cut_point(Node* edge_start, Node* edge_end, std::vector<Point*>& edge_points);

      /*! \brief Unused */
      void cut_points(Side* side, PointSet& cut_points);

      /*! \brief Unused */
      void cut_points_between(Point* begin, Point* end, std::vector<Point*>& line);

      void cut_points_including(Point* begin, Point* end, std::vector<Point*>& line);

      void cut_points_inside(Element* element, std::vector<Point*>& line);

      bool is_cut(Side* side);

      /** \brief Compute the intersection between THIS edge and a \c other
       *  edge.
       *
       *  \param mesh       (in)  : cut mesh, necessary to add cut points to the point pool
       *  \param other      (in)  : pointer to the other edge object
       *  \param side       (in)  : pointer to the other side object(can be nullptr)
       *  \param cut_points (out) : if the cut was successful, this set will hold the new cut points
       *  \param tolerance  (out) : tolerance used for the intersection ( defined by the cut kernel
       * )
       *
       *  This routine returns TRUE, if the computation was successful AND the cut
       *  point is within the edge limits.
       *
       *  \author hiermeier \date 12/16 */
      virtual bool compute_cut(
          Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance) = 0;

      Node* begin_node() { return nodes_.front(); }

      Node* middle_node()
      {
        if (nodes_.size() != 3) FOUR_C_THROW("middle node in line3 only");
        return nodes_[2];
      }

      Node* end_node() { return nodes_[1]; }

      Point* node_in_element(Element* element, Point* other);

      const std::vector<Node*>& nodes() const { return nodes_; }

      /// Find common points (excluding cut_points points) between two edges
      void common_nodal_points(Edge* edge, std::vector<Point*>& common);

      bool find_cut_points(Mesh& mesh, Element* element, Side& side, Side& other);

      /*!
      \brief Computes the points at which both the sides intersect
       */
      bool find_cut_points_mesh_cut(
          Mesh& mesh, Element* element, Side& side, Side& other, PointSet* cutpoints = nullptr);

      /*!
      \brief Simplified version of the find_cut_points as used by the function LevelSetCut
       */
      bool find_cut_points_level_set(Mesh& mesh, Element* element, Side& side, Side& other);

      /*!
      \brief Cut points falling on this edge that are common to the two given sides are extracted
       */
      void get_cut_points(Element* element, Side& side, Side& other, PointSet& cuts);

      /*!
      \brief Cut points falling on this edge that are common to the given edge are extracted
       */
      void get_cut_points(Edge* other, PointSet& cuts);

      const PointPositionSet& cut_points() const { return cut_points_; }
      // const std::vector<Point*> & CutPoints() const { return cut_points_; }

      void rectify_cut_numerics();

      /// Returns the selfcutposition of this edge
      Point::PointPosition self_cut_position() { return selfcutposition_; }

      /// Gives this edge a selfcutposition and spreads the positional information
      void self_cut_position(Point::PointPosition p);

      /// Changes the selfcutposition of this edge and spreads the positional information
      void change_self_cut_position(Point::PointPosition p);

      /// Erase the cutside from this edge because it is deleted in the selfcut
      void erase_cut_side(Side* cutside) { sides_.erase(cutside); }

      /*!
      \brief Replace the node "nod" of this edge with given node "replwith"
       */
      void replace_node(Node* nod, Node* replwith);

      // remove information about this point
      void remove_point(Point* p) { cut_points_.erase(p); };

      /*!
       \brief Add all topological connections of interesection of this edge and other edge ( all
       necessery pairs, etc)
       */
      void add_connections(Edge* other, const std::pair<Side*, Edge*>& original_cut_pair);

      void add_connections(Edge* other, Side* original_side, Edge* original_edge);

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
    template <unsigned prob_dim, Core::FE::CellType edge_type,
        unsigned dim_edge = Core::FE::dim<edge_type>,
        unsigned num_nodes_edge = Core::FE::num_nodes<edge_type>>
    class ConcreteEdge : public Edge
    {
     public:
      /// constructor
      ConcreteEdge(const std::vector<Node*>& nodes) : Edge(nodes) {}

      /// get the element dimension of this edge
      unsigned n_dim() const override { return dim_edge; }

      /// get the number of nodes of this edge
      unsigned num_nodes() const override { return num_nodes_edge; }

      /// get the problem dimension
      unsigned n_prob_dim() const override { return prob_dim; }

      /// get the shape of this edge element
      Core::FE::CellType shape() const override { return edge_type; }

      /*! \brief Get the coordinates of the nodes of the side */
      void coordinates(double* xyze) override
      {
        double* x = xyze;
        for (std::vector<Node*>::const_iterator i = nodes().begin(); i != nodes().end(); ++i)
        {
          const Node& n = **i;
          n.coordinates(x);
          x += prob_dim;
        }
      }

      void get_touching_points(const std::vector<Node*>& nodes, std::vector<Node*>& touch_nodes,
          CutFloatType floattype = floattype_double) override;

      /*! \brief Handles intersection of two edges that are close to each other */
      virtual bool handle_parallel_cut(
          Edge* other, Side* side, PointSet* cut_points, CutFloatType floattype = floattype_double);

      bool just_parallel_cut(Mesh& mesh, Side& side, PointSet& cuts, int skip_id = -1) override;

      /*! \brief Get the intersection points of this edge with the given side and
       * store the cut points in cuts */
      bool cut(Mesh& mesh, Side& side, PointSet& cuts) override;

      /*! \brief Get the intersection point of THIS edge with the given \c other edge and
       *  store the global cut point in x and the local coordinate of THIS edge in \c pos.
       *
       *  \param mesh       (in)  : cut mesh, necessary to add cut points to the point pool
       *  \param other      (in)  : pointer to the other edge object
       *  \param side       (in)  : pointer to the other side object(can be nullptr)
       *  \param cut_points (out) : if the cut was successful, this set will hold the new cut points
       *  \param tolerance  (out) : tolerance used for the intersection ( defined by the cut kernel
       * )
       *
       *  This routine returns TRUE, if the computation was successful AND the cut
       *  point is within the edge limits.
       *
       *  \author hiermeier \date 12/16 */
      bool compute_cut(
          Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance) override;

      /*! \brief Get the intersection points of this edge with the level-set side
       *  and store the cut points in cuts */
      bool level_set_cut(Mesh& mesh, Side& side, PointSet& cuts) override
      {
        double blsv = begin_node()->lsv();
        double elsv = end_node()->lsv();

        bool cutfound = false;
        // version for single element cuts, here we need to watch for tolerances on
        // nodal cuts
        if (std::abs(blsv) <= REFERENCETOL)
        {
          cuts.insert(Point::insert_cut(this, &side, begin_node()));
          cutfound = true;
        }
        if (std::abs(elsv) <= REFERENCETOL)
        {
          cuts.insert(Point::insert_cut(this, &side, end_node()));
          cutfound = true;
        }

        if (not cutfound)
        {
          if ((blsv < 0.0 and elsv > 0.0) or (blsv > 0.0 and elsv < 0.0))
          {
            cutfound = true;
            double z = blsv / (blsv - elsv);

            Core::LinAlg::Matrix<prob_dim, 1> x1;
            Core::LinAlg::Matrix<prob_dim, 1> x2;
            begin_node()->coordinates(x1.data());
            end_node()->coordinates(x2.data());

            Core::LinAlg::Matrix<prob_dim, 1> x;
            x.update(-1., x1, 1., x2, 0.);
            x.update(1., x1, z);
            Point* p = Point::new_point(mesh, x.data(), 2. * z - 1., this, &side, 0.0);
            cuts.insert(p);
          }
        }

        //      std::cout <<"LS cut found? --> " << (cutfound ? "TRUE" : "FALSE") << std::endl;
        //      FOUR_C_THROW("Core::Geo::Cut::Edge::LevelSetCut -- STOP -- hiermeier 08/16");

        return cutfound;
      }

     private:
      Teuchos::RCP<Core::Geo::Cut::IntersectionBase> intersection_ptr(
          const Core::FE::CellType& sidetype) const;

    };  // class ConcreteEdge

    /*--------------------------------------------------------------------------*/
    class EdgeFactory
    {
     public:
      /// constructor
      EdgeFactory(){};

      /// destructor
      virtual ~EdgeFactory() = default;

      /// create the concrete edge object
      Teuchos::RCP<Edge> create_edge(
          const Core::FE::CellType& edgetype, const std::vector<Node*>& nodes) const;

     private:
      template <Core::FE::CellType edgetype>
      Edge* create_concrete_edge(const std::vector<Node*>& nodes, int probdim) const
      {
        Edge* e = nullptr;
        switch (probdim)
        {
          case 2:
            e = new ConcreteEdge<2, edgetype>(nodes);
            break;
          case 3:
            e = new ConcreteEdge<3, edgetype>(nodes);
            break;
          default:
            FOUR_C_THROW("Unsupported problem dimension! (probdim=%d)", probdim);
            break;
        }
        return e;
      }
    };  // class EdgeFactory
  }     // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
