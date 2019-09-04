/*---------------------------------------------------------------------*/
/*! \file

\brief PointGraph, Graph Algorithm to create Facets from lines and edges

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

// necessary due to usage of graph_t, vertex_t etc.
#include "cut_find_cycles.H"

// forward declarations
namespace GEO
{
  namespace CUT
  {
    class Point;
    class Line;
    class Side;
    class Element;
    class Mesh;

    namespace IMPL
    {
      /// a planar graph that is used to create facets out of lines and points
      class PointGraph
      {
       public:
        /** specifies to which side the pointgraph and the related facets will
         *  belong. See the Side::MakeInternalFacets and Side::MakeOwnedSideFacets
         *  routines for more information. */
        enum Location
        {
          element_side,
          cut_side
        };

        enum Strategy
        {
          all_lines,
          own_lines  // used in levelset cut to create internal facets
        };

        /** \brief create a point graph object
         *
         *  Call this method to get the pointgraph which fits to your problem dimension.
         *
         *  \author hiermeier \date 11/16 */
        static PointGraph *Create(Mesh &mesh, Element *element, Side *side,
            PointGraph::Location location, PointGraph::Strategy strategy);

       protected:
        class Graph
        {
         public:
          /// constructor
          Graph()
          {
            // empty
          }

          virtual ~Graph(){};

          void AddEdge(int row, int col);

          void AddEdge(Point *p1, Point *p2);

          Point *GetPoint(int i);

          void Print(std::ostream &stream = std::cout);

          void PlotAllPoints(std::ostream &stream = std::cout);

          void PlotPoints(Element *element);

          /** Creates maincycles (outer polygons) and holecycles (inner polygons = holes)
           *  of the selfcut graph */
          void FindCycles(Side *side, Cycle &cycle);

          virtual void FindCycles(
              Element *element, Side *side, Cycle &cycle, Location location, Strategy strategy);

          /*!
          \brief Any edge with single point in the graph is deleted
           */
          void FixSinglePoints(Cycle &cycle);

          virtual bool HasSinglePoints(Location location);

          virtual bool HasTouchingEdge(Element *element, Side *side);

          // Simplify connection if the single point lies close to the nodal point
          // and touches the same edges as the nodal point
          // however because of the way pointgraph is constructed only connection from nodal
          // to this point is created
          // In this case we remove connection to nodal point, and treat cut point and "new nodal
          // point"
          virtual bool SimplifyConnections(Element *element, Side *side);

          void GnuplotDumpCycles(const std::string &filename, const std::vector<Cycle> &cycles);

          std::map<int, plain_int_set> graph_;
          std::map<int, Point *> all_points_;
          std::vector<Cycle> main_cycles_;
          std::vector<std::vector<Cycle>> hole_cycles_;

        };  // struct Graph

        /** empty constructor (for derived classes only)
         *
         *  \author hiermeier \date 11/16 */
        PointGraph(unsigned dim) : graph_(CreateGraph(dim)){/* intentionally left blank */};

       public:
        typedef std::vector<Cycle>::iterator facet_iterator;
        typedef std::vector<std::vector<Cycle>>::iterator hole_iterator;

        PointGraph(Mesh &mesh, Element *element, Side *side, Location location, Strategy strategy);

        /// Constructor for the selfcut
        PointGraph(Side *side);

        /// destructor
        virtual ~PointGraph(){};

        facet_iterator fbegin() { return GetGraph().main_cycles_.begin(); }

        facet_iterator fend() { return GetGraph().main_cycles_.end(); }

        hole_iterator hbegin() { return GetGraph().hole_cycles_.begin(); }

        hole_iterator hend() { return GetGraph().hole_cycles_.end(); }

        void Print() { GetGraph().Print(); }

       protected:
        /*! \brief Graph is filled with all edges that are created due to additional
         *  cut points and cut lines */
        void FillGraph(Element *element, Side *side, Cycle &cycle, Strategy strategy);

        /** Graph is filled with all edges of the selfcut: uncutted edges,
         *  selfcutedges and new splitted edges; but no the cutted edges */
        void FillGraph(Side *side, Cycle &cycle);

        /** \brief add cut lines to graph
         *
         *  no need to add any more point to cycle because cut lines just join already
         *  existing points on the edge. making cut lines do not introduce additional
         *  points */
        virtual void AddCutLinesToGraph(Element *element, Side *side, Strategy strategy);

        virtual void BuildCycle(const std::vector<Point *> &edge_points, Cycle &cycle) const;

        /// access the graph of the most derived class
        virtual inline PointGraph::Graph &GetGraph() { return *graph_; }

        virtual inline const Teuchos::RCP<PointGraph::Graph> &GraphPtr() { return graph_; }

       private:
        Teuchos::RCP<GEO::CUT::IMPL::PointGraph::Graph> CreateGraph(unsigned dim);

        Teuchos::RCP<Graph> graph_;
      };  // class PointGraph

      // non-member function
      bool FindCycles(graph_t &g, Cycle &cycle, std::map<vertex_t, LINALG::Matrix<3, 1>> &local,
          std::vector<Cycle> &cycles);

      Teuchos::RCP<PointGraph> CreatePointGraph(Mesh &mesh, Element *element, Side *side,
          const PointGraph::Location &location, const PointGraph::Strategy &strategy);
    }  // namespace IMPL
  }    // namespace CUT
}  // namespace GEO


#endif
