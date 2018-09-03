/*---------------------------------------------------------------------*/
/*!
\file cut_cycle.cpp

\brief a cylcle of points (basic to create facets)

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_cycle.H"
#include "cut_edge.H"


std::ostream& operator<<(std::ostream& stream, const GEO::CUT::Cycle& cycle)
{
  std::copy(cycle.points_.begin(), cycle.points_.end(),
      std::ostream_iterator<GEO::CUT::Point*>(stream, " "));
  return stream << "\n";
}

bool GEO::CUT::Cycle::MakeCycle(const point_line_set& lines, Cycle& cycle)
{
  cycle.clear();
  std::vector<Point*> frompoints;
  std::vector<Point*> topoints;
  frompoints.reserve(lines.size());
  topoints.reserve(lines.size());
  for (point_line_set::const_iterator i = lines.begin(); i != lines.end(); ++i)
  {
    frompoints.push_back(i->first);
    topoints.push_back(i->second);
  }

  cycle.reserve(lines.size());

  std::vector<int> done(lines.size(), false);

  unsigned pos = 0;
  while (not done[pos])
  {
    cycle.push_back(frompoints[pos]);
    done[pos] = true;
    std::vector<Point*>::iterator i =
        std::find(frompoints.begin(), frompoints.end(), topoints[pos]);
    if (i == frompoints.end())
    {
      run_time_error("no cycle: \"to point\" not in \"from list\"");
    }
    pos = std::distance(frompoints.begin(), i);
  }

  if (cycle.size() != lines.size())
  {
    cycle.clear();
    return false;
  }
  return true;
}

bool GEO::CUT::Cycle::IsValid() const
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

bool GEO::CUT::Cycle::IsCut(Element* element) const
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

void GEO::CUT::Cycle::Add(point_line_set& lines) const
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

void GEO::CUT::Cycle::CommonEdges(plain_edge_set& edges) const
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

void GEO::CUT::Cycle::CommonSides(plain_side_set& sides) const
{
  std::vector<Point*>::const_iterator i = points_.begin();
  if (i != points_.end())
  {
    sides = (*i)->CutSides();
    for (++i; i != points_.end(); ++i)
    {
      Point* p = *i;
      p->Intersection(sides);
      if (sides.size() == 0)
      {
        break;
      }
    }
  }
  else
  {
    sides.clear();
  }
}

void GEO::CUT::Cycle::Intersection(plain_side_set& sides) const
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    p->Intersection(sides);
    if (sides.size() == 0) break;
  }
}

bool GEO::CUT::Cycle::Equals(const Cycle& other)
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

void GEO::CUT::Cycle::DropPoint(Point* p)
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

void GEO::CUT::Cycle::TestUnique()
{
  PointSet c_copy;
  c_copy.insert(points_.begin(), points_.end());
  if (points_.size() != c_copy.size()) throw std::runtime_error("double point in cycle");
}

void GEO::CUT::Cycle::GnuplotDump(std::ostream& stream) const
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

void GEO::CUT::Cycle::reverse() { std::reverse(points_.begin(), points_.end()); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Cycle::Print() const
{
  std::cout << "--- Cycle ---" << std::endl;
  for (std::vector<Point*>::const_iterator cit = points_.begin(); cit != points_.end(); ++cit)
    std::cout << "Point " << (*cit) << "\n";
  std::cout << std::endl;
}
