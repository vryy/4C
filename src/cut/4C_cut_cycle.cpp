/*---------------------------------------------------------------------*/
/*! \file

\brief a cylcle of points (basic to create facets)

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_cycle.hpp"

#include "4C_cut_edge.hpp"
#include "4C_cut_output.hpp"

FOUR_C_NAMESPACE_OPEN

bool Core::Geo::Cut::Cycle::is_valid() const
{
  if (points_.size() < 3) return false;

  // ignore cycles with all points on one and the same edge
  {
    plain_edge_set edges;
    common_edges(edges);
    if (edges.size() > 0)
    {
      return false;
    }
  }

  return true;
}

bool Core::Geo::Cut::Cycle::is_cut(Element* element) const
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    if (not p->is_cut(element))
    {
      return false;
    }
  }
  return true;
}

void Core::Geo::Cut::Cycle::add(point_line_set& lines) const
{
  for (unsigned i = 0; i != points_.size(); ++i)
  {
    Point* p1 = points_[i];
    Point* p2 = points_[(i + 1) % points_.size()];

    std::pair<Point*, Point*> line = std::make_pair(p1, p2);
    if (lines.count(line) > 0)
    {
      lines.erase(line);
    }
    else
    {
      std::pair<Point*, Point*> reverse_line = std::make_pair(p2, p1);
      if (lines.count(reverse_line) > 0)
      {
        lines.erase(reverse_line);
      }
      else
      {
        lines.insert(line);
      }
    }
  }
}

void Core::Geo::Cut::Cycle::common_edges(plain_edge_set& edges) const
{
  std::vector<Point*>::const_iterator i = points_.begin();
  if (i != points_.end())
  {
    edges = (*i)->cut_edges();
    for (++i; i != points_.end(); ++i)
    {
      Point* p = *i;
      p->intersection(edges);
      if (edges.size() == 0)
      {
        break;
      }
    }
  }
  else
  {
    edges.clear();
  }
}

void Core::Geo::Cut::Cycle::intersection(plain_side_set& sides) const
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    p->intersection(sides);
    if (sides.size() == 0) break;
  }
}

bool Core::Geo::Cut::Cycle::equals(const Cycle& other)
{
  if (size() != other.size())
  {
    return false;
  }

  for (std::vector<Point*>::const_iterator i = other.points_.begin(); i != other.points_.end(); ++i)
  {
    Point* p = *i;
    // if ( not std::binary_search( sorted.begin(), sorted.end(), p ) )

    if (std::count(points_.begin(), points_.end(), p) != 1)
    {
      return false;
    }
  }
  return true;
}

void Core::Geo::Cut::Cycle::drop_point(Point* p)
{
  std::vector<Point*>::iterator j = std::find(points_.begin(), points_.end(), p);
  if (j != points_.end())
  {
    std::vector<Point*>::iterator prev = j == points_.begin() ? points_.end() : j;
    std::advance(prev, -1);
    std::vector<Point*>::iterator next = j;
    std::advance(next, 1);
    if (next == points_.end()) next = points_.begin();
    if (*prev == *next)
    {
      if (next > j)
      {
        points_.erase(next);
        points_.erase(j);
      }
      else
      {
        points_.erase(j);
        points_.erase(next);
      }
    }
    else
    {
      points_.erase(j);
    }
  }
}

void Core::Geo::Cut::Cycle::test_unique()
{
  PointSet c_copy;
  c_copy.insert(points_.begin(), points_.end());
  if (points_.size() != c_copy.size())
  {
    for (std::vector<Point*>::iterator it = points_.begin(); it != points_.end(); ++it)
    {
      int num_el_erased = c_copy.erase(*it);
      // if element was duplicated we cannot erase it again
      if (num_el_erased == 0)
      {
        int num_occ = std::count(points_.begin(), points_.end(), (*it));
        print();
        std::stringstream str;
        str << "Multiple( " << num_occ << " ) occcurence of point " << (*it)->id()
            << " in the cycle" << std::endl;
        FOUR_C_THROW(str.str());
      }
    }
  }
}

void Core::Geo::Cut::Cycle::gnuplot_dump(std::ostream& stream) const
{
  for (unsigned i = 0; i != points_.size(); ++i)
  {
    Point* p1 = points_[i];
    Point* p2 = points_[(i + 1) % points_.size()];

    p1->plot(stream);
    p2->plot(stream);
    stream << "\n\n";
  }
}


void Core::Geo::Cut::Cycle::gmsh_dump(std::ofstream& file) const
{
  for (unsigned i = 0; i != points_.size(); ++i)
  {
    Point* p1 = points_[i];
    p1->dump_connectivity_info();
    Point* p2 = points_[(i + 1) % points_.size()];
    std::stringstream section_name;
    section_name << "Line" << i;
    Core::Geo::Cut::Output::GmshNewSection(file, section_name.str());
    Core::Geo::Cut::Output::GmshLineDump(file, p1, p2, p1->id(), p2->id(), false, nullptr);
    Core::Geo::Cut::Output::GmshEndSection(file, false);
  }
}


void Core::Geo::Cut::Cycle::reverse() { std::reverse(points_.begin(), points_.end()); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Cycle::print() const
{
  std::cout << "--- Cycle ---" << std::endl;
  for (std::vector<Point*>::const_iterator cit = points_.begin(); cit != points_.end(); ++cit)
    std::cout << "Point " << (*cit)->id() << "\n";
  std::cout << std::endl;
}

FOUR_C_NAMESPACE_CLOSE
