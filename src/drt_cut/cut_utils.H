#ifndef CUT_UTILS_H
#define CUT_UTILS_H

#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <cstring>

namespace GEO
{
  namespace CUT
  {
    class Node;
    class Edge;
    class Side;
    class Element;

    class Point;
    class Line;
    class Facet;
    class VolumeCell;

    class BoundaryCell;
    class IntegrationCell;

    class Mesh;

    typedef std::set<int>         plain_int_set;

    typedef std::set<Node*>       plain_node_set;
    typedef std::set<Edge*>       plain_edge_set;
    typedef std::set<Side*>       plain_side_set;
    typedef std::set<Element*>    plain_element_set;

    typedef std::set<Point*>      plain_point_set;
    typedef std::set<Line*>       plain_line_set;
    typedef std::set<Facet*>      plain_facet_set;
    typedef std::set<VolumeCell*> plain_volumecell_set;

    typedef std::set<BoundaryCell*>    plain_boundarycell_set;
    typedef std::set<IntegrationCell*> plain_integrationcell_set;

    template <class set>
    void set_erase( set & s, typename set::iterator & i )
    {
      s.erase( i++ );
    }

    inline void DumpDoubles( std::ostream & stream, const double * data, int length )
    {
      std::vector<int> c( length*sizeof( double ) / sizeof( int ) );
      std::memcpy( &c[0], data, c.size()*sizeof( int ) );
      stream << "  {";
      std::copy( c.begin(), c.end(), std::ostream_iterator<int>( stream, "," ) );
      stream << "};";
    }

  }
}

#endif
