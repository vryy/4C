/*---------------------------------------------------------------------*/
/*!
\file cut_utils.H

\brief Utils Header defining ours sets of geometrical objects

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_UTILS_H
#define CUT_UTILS_H

#include <iostream>
#include <iterator>
#include <set>
#include <cstring>

//#define DRT_CUT_DUMPCREATION    // dump object creation to generate test case
#define CUT_USE_SORTED_VECTOR

#ifdef CUT_USE_SORTED_VECTOR
#include "sorted_vector.H"
#else
#include <vector>
#endif

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

#ifdef CUT_USE_SORTED_VECTOR

    typedef sorted_vector<int>         plain_int_set;

    typedef sorted_vector<Node*>       plain_node_set;
    typedef sorted_vector<Edge*>       plain_edge_set;
    typedef sorted_vector<Side*>       plain_side_set;
    typedef sorted_vector<Element*>    plain_element_set;

    typedef sorted_vector<Point*>      plain_point_set;
    typedef sorted_vector<Line*>       plain_line_set;
    typedef sorted_vector<Facet*>      plain_facet_set;
    typedef sorted_vector<VolumeCell*> plain_volumecell_set;

    typedef sorted_vector<std::pair<Point*, Point*> > point_line_set;

    typedef sorted_vector<BoundaryCell*>    plain_boundarycell_set;
    typedef sorted_vector<IntegrationCell*> plain_integrationcell_set;

    template <class set>
    void set_erase( set & s, typename set::iterator & i )
    {
      s.ierase( i );
    }

#else

    typedef std::set<int>         plain_int_set;

    typedef std::set<Node*>       plain_node_set;
    typedef std::set<Edge*>       plain_edge_set;
    typedef std::set<Side*>       plain_side_set;
    typedef std::set<Element*>    plain_element_set;

    typedef std::set<Point*>      plain_point_set;
    typedef std::set<Line*>       plain_line_set;
    typedef std::set<Facet*>      plain_facet_set;
    typedef std::set<VolumeCell*> plain_volumecell_set;

    typedef std::set<std::pair<Point*, Point*> > point_line_set;

    typedef std::set<BoundaryCell*>    plain_boundarycell_set;
    typedef std::set<IntegrationCell*> plain_integrationcell_set;

    template <class set>
    void set_erase( set & s, typename set::iterator & i )
    {
      s.erase( i++ );
    }

#endif

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
