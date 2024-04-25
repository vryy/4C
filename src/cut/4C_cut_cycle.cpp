/*---------------------------------------------------------------------*/
/*! \file

\brief a cylcle of points (basic to create facets)

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_cycle.hpp"

#include "4C_cut_edge.hpp"
#include "4C_cut_output.hpp"

FOUR_C_NAMESPACE_OPEN

bool CORE::GEO::CUT::Cycle::IsValid() const
{
  if (points_.size() < 3) return false;

  // ignore cycles with all points on one and the same edge
  {
    plain_edge_set edges;
    CommonEdges(edges);
    if (edges.size() > 0)
    {
      return false;
    }
  }

  return true;
}

bool CORE::GEO::CUT::Cycle::IsCut(Element* element) const
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    if (not p->IsCut(element))
    {
      return false;
    }
  }
  return true;
}

void CORE::GEO::CUT::Cycle::Add(point_line_set& lines) const
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

void CORE::GEO::CUT::Cycle::CommonEdges(plain_edge_set& edges) const
{
  std::vector<Point*>::const_iterator i = points_.begin();
  if (i != points_.end())
  {
    edges = (*i)->CutEdges();
    for (++i; i != points_.end(); ++i)
    {
      Point* p = *i;
      p->Intersection(edges);
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

void CORE::GEO::CUT::Cycle::Intersection(plain_side_set& sides) const
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    p->Intersection(sides);
    if (sides.size() == 0) break;
  }
}

bool CORE::GEO::CUT::Cycle::Equals(const Cycle& other)
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

void CORE::GEO::CUT::Cycle::DropPoint(Point* p)
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

void CORE::GEO::CUT::Cycle::TestUnique()
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
        Print();
        std::stringstream str;
        str << "Multiple( " << num_occ << " ) occcurence of point " << (*it)->Id()
            << " in the cycle" << std::endl;
        FOUR_C_THROW(str.str());
      }
    }
  }
}

void CORE::GEO::CUT::Cycle::GnuplotDump(std::ostream& stream) const
{
  for (unsigned i = 0; i != points_.size(); ++i)
  {
    Point* p1 = points_[i];
    Point* p2 = points_[(i + 1) % points_.size()];

    p1->Plot(stream);
    p2->Plot(stream);
    stream << "\n\n";
  }
}


void CORE::GEO::CUT::Cycle::GmshDump(std::ofstream& file) const
{
  for (unsigned i = 0; i != points_.size(); ++i)
  {
    Point* p1 = points_[i];
    p1->DumpConnectivityInfo();
    Point* p2 = points_[(i + 1) % points_.size()];
    std::stringstream section_name;
    section_name << "Line" << i;
    CORE::GEO::CUT::OUTPUT::GmshNewSection(file, section_name.str());
    CORE::GEO::CUT::OUTPUT::GmshLineDump(file, p1, p2, p1->Id(), p2->Id(), false, nullptr);
    CORE::GEO::CUT::OUTPUT::GmshEndSection(file, false);
  }
}


void CORE::GEO::CUT::Cycle::reverse() { std::reverse(points_.begin(), points_.end()); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::Cycle::Print() const
{
  std::cout << "--- Cycle ---" << std::endl;
  for (std::vector<Point*>::const_iterator cit = points_.begin(); cit != points_.end(); ++cit)
    std::cout << "Point " << (*cit)->Id() << "\n";
  std::cout << std::endl;
}

FOUR_C_NAMESPACE_CLOSE
