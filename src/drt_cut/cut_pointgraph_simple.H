/*----------------------------------------------------------------------------*/
/**
\file cut_pointgraph_simple.H

\brief Create a simple point graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )

\maintainer Michael Hiermeier

\date Nov 12, 2016

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CUT_CUT_POINTGRAPH_SIMPLE_H_
#define SRC_DRT_CUT_CUT_POINTGRAPH_SIMPLE_H_

#include "cut_pointgraph.H"

namespace GEO {
namespace CUT {
namespace IMPL {
  /*--------------------------------------------------------------------------*/
  /** \brief Simplified version of the base class for 1-D elements
   *
   *  \author hiermeier \date 11/16 */
  class SimplePointGraph_1D : public PointGraph
  {
  public:
    class Graph : public PointGraph::Graph
    {
    public:
      /// constructor
      Graph() : PointGraph::Graph()
      {
        // empty
      }

      /// destructor
      virtual ~Graph() {};

      /** \brief Simplified version of the base class implementation
       *
       *  \author hiermeier \date 11/16 */
      virtual void FindCycles( Element * element, Side * side, Cycle & cycle,
          Location location, Strategy strategy );

    };  // struct Graph

    /// constructor
    SimplePointGraph_1D( Mesh & mesh, Element * element,
        Side * side, Location location, Strategy strategy );

    /// destructor
    virtual ~SimplePointGraph_1D() {};

  protected:
    /** \brief Build the cycle of edge points [derived]
     *
     *  \param edge_points (in) : points of the current edge
     *  \param cycle       (out): created point cycle
     *
     *  \author hiermeier \date 11/16 */
    virtual void BuildCycle(const std::vector<Point*> & edge_points,
        Cycle & cycle) const;

    /** \brief Add cut lines to the graph object [derived]
     *
     *  There are no cut lines in 1-D, thus we jump to AddCutPointsToCycle,
     *  thus it is at least a little bit less confusing.
     *
     *  \author hiermeier \date 11/16 */
    virtual void AddCutLinesToGraph(Element * element, Side * side, Strategy strategy,
        Cycle & cycle);

    /** \brief Add points instead of lines in 1-D
     *
     *  \param element (in) : background element
     *  \param side    (in) : element or cut side (see PointGraph::Location enumerator)
     *  \param cycle   (out): point cycle to be filled
     *
     *  This routine becomes for example important for the LevelSetCut, since a
     *  LevelSet side has no edges nor nodes. But we need the cut point for the boundary
     *  integration cell construction. */
    void AddCutPointsToCycle( Element * element, Side * side, Cycle & cycle );

  }; // class SimplePointGraph_1D

  /*--------------------------------------------------------------------------*/
  /** \brief Modified version of the base class for 2-D elements
   *
   *  \author hiermeier \date 11/16 */
  class SimplePointGraph_2D : public PointGraph
  {
  public:
    class Graph : public PointGraph::Graph
    {
    public:
      /// constructor
      Graph() : PointGraph::Graph(), correct_rotation_direction_( false )
      {
        // empty
      };

      /// destructor
      virtual ~Graph() {};

      /** \brief Overloaded version of the base class implementation
       *
       *  The base class implementation is called first. Afterwards
       *  a split of the created surface facet cycles into line cycles
       *  is performed.
       *
       *  \author hiermeier \date 11/16 */
      virtual void FindCycles( Element * element, Side * side, Cycle & cycle,
          Location location, Strategy strategy );

      inline const std::vector<Cycle> & SurfaceMainCycles() const
      {
        return surface_main_cycles_;
      }

      virtual bool HasSinglePoints( Location location );

      void SetCorrectRotationDirection( bool correct_rotation )
      {
        correct_rotation_direction_ = correct_rotation;
      }

    private:

      /**  \brief split the main cycles which are created by the base class
       *   implementation into line cycles
       *
       *   We store the generated surface facet main cycles as well, since we
       *   can reuse the implementation for the volume cell generation in 2-D.
       *
       *   \author hiermeier \date 01/17 */
      void SplitMainCyclesIntoLineCycles();

    private:
      std::vector<Cycle> surface_main_cycles_;

      bool correct_rotation_direction_;
    };  // struct Graph

    typedef std::vector<Cycle>::const_iterator surface_const_iterator;

    inline surface_const_iterator sbegin() const
    {
      return graph_2d_->SurfaceMainCycles().begin();
    }

    inline surface_const_iterator send() const
    {
      return graph_2d_->SurfaceMainCycles().end();
    }

    inline unsigned NumSurfaces() const
    {
      return graph_2d_->SurfaceMainCycles().size();
    }

    /// constructor
    SimplePointGraph_2D( Mesh & mesh, Element * element,
        Side * side, Location location, Strategy strategy );

    // empty constructor
    SimplePointGraph_2D();

    /// destructor
    virtual ~SimplePointGraph_2D() {};

    void FindLineFacetCycles( const plain_facet_set & line_facets,
        Element * parent_element );

    /** \brief Test and correct the rotation direction of the given cycles such
     *  that the direction fits the rotation direction of the underlying side
     *
     *  The test is necessary, since the created cycles from the boost library seem
     *  to have no preferred rotation direction. First the normal on the parent side
     *  is calculated, afterwards the normal of the given cycle is calculated. If the
     *  inner-product of these normals is positive, the rotation direction is correct,
     *  otherwise the cycle direction is inverted.
     *
     *  \param side   (in)     : parent side ( for 2-D elements this is equivalent to the element )
     *  \param cycles (in/out) : cycles which have to be checked and inverted if necessary
     *
     *  \author hiermeier \date 01/17 */
    static void CorrectRotationDirection( const Side * side, std::vector<Cycle> & cycles );

  private:
    void FillGraphAndCycleWithLineFacets( const plain_facet_set & line_facets, Cycle & cycle );

    void FindCycles( Element * element, Cycle & cycle );

  private:
    Teuchos::RCP<SimplePointGraph_2D::Graph> graph_2d_;
  }; // class SimplePointGraph_2D

} // namespace IMPL
} // namespace CUT
} // namespace GEO


#endif /* SRC_DRT_CUT_CUT_POINTGRAPH_SIMPLE_H_ */
