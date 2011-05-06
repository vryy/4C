#ifndef CUT_LINESEGMENT_H
#define CUT_LINESEGMENT_H

#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>

namespace GEO
{
namespace CUT
{
  class Point;
  class Line;
  class Side;
  class Element;
  class Mesh;

class LineSegment
{
public:

  LineSegment( Mesh & mesh, Element * element, Side * side, std::set<Line*> & cut_lines, bool inner );

  LineSegment( Mesh & mesh, Element * element, Side * side, const std::vector<Point*> & points );

  bool Combine( Mesh & mesh, Element * element, Side * side, LineSegment & other );

  bool IsClosed() const { return closed_; }

  Point* BeginPoint() { return facet_points_.front(); }

  Point* EndPoint() { return facet_points_.back(); }

  Line* BeginLine() { return facet_lines_.front(); }

  Line* EndLine() { return facet_lines_.back(); }

  const std::vector<Line*> & Lines() { return facet_lines_; }

  const std::vector<Point*> & Points() { return facet_points_; }

  Side* OnSide( Element * element );

private:

  bool ClosedOnEdge( Mesh & mesh, Element * element, Side * side );

  bool ClosedOnEdge( Mesh & mesh, Element * element, Side * side,
                     Point * begin, Point * end,
                     Line * begin_line, Line * end_line );

  void AppendLines( Mesh & mesh, Element * element, Side * side, const std::vector<Point*> & line );

  void AppendLinesReversed( Mesh & mesh, Element * element, Side * side, const std::vector<Point*> & line );

  void InsertLines( Mesh & mesh, Element * element, Side * side, const std::vector<Point*> & line );

  void InsertLinesReversed( Mesh & mesh, Element * element, Side * side, const std::vector<Point*> & line );

  void CloseGap( Mesh & mesh, Element * element, Side * side, Point * end );

  std::vector<Line*> facet_lines_;
  std::vector<Point*> facet_points_;

  bool closed_;
};

class LineSegmentList
{
public:

  void Create( Mesh & mesh, Element * element, Side * side, bool inner );

  void Create( Mesh & mesh, Element * element, Side * side, std::set<Line*> & lines, bool inner );

  const std::vector<Teuchos::RCP<LineSegment> > & Segments() const { return segments_; }

private:
  std::vector<Teuchos::RCP<LineSegment> > segments_;
};

}
}

#endif
