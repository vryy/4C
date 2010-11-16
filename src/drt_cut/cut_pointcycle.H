#ifndef CUT_POINTCYCLE_H
#define CUT_POINTCYCLE_H

#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>

#include "cut_point.H"

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Side;
    class Element;
    class Facet;
    class Mesh;

class InnerCutPoints
{
public:

  InnerCutPoints( const std::map<Point*, std::vector<Point*>::const_iterator*> & point_positions,
                  const std::vector<Point*> & points );

  bool operator==( const std::vector<Point*> & other ) const;
  bool operator!=( const std::vector<Point*> & other ) const
  {
    return not ( *this == other );
  }

  bool operator==( const InnerCutPoints & other ) const
  {
    return *this == other.Points();
  }

  bool operator!=( const InnerCutPoints & other ) const
  {
    return not ( *this == other );
  }

  const std::vector<Point*> & Points() const { return points_; }

  Point* Front() const { return points_.front(); }
  Point* Back()  const { return points_.back(); }

private:

  std::vector<Point*> points_;
};


class InnerCutPointsLess
{
public:

  InnerCutPointsLess( const std::map<Point*, std::vector<Point*>::const_iterator*> & point_positions )
    : point_positions_( point_positions )
  {
  }

  bool operator()( const InnerCutPoints & l1, const InnerCutPoints & l2 ) const
  {
    Point * p1 = l1.Front();
    Point * p2 = l2.Front();

    if ( p1!=p2 )
    {
      std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i1 = point_positions_.find( p1 );
      std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i2 = point_positions_.find( p2 );

      if ( i1==point_positions_.end() or i2==point_positions_.end() )
      {
        throw std::runtime_error( "inner cut list start point not in cut point list" );
      }

      return *i1->second < *i2->second;
    }
    else
    {
      p1 = l1.Back();
      p2 = l2.Back();

      if ( p1==p2 )
      {
        throw std::runtime_error( "same start and end point in different inner cuts" );
      }

      std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i1 = point_positions_.find( p1 );
      std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i2 = point_positions_.find( p2 );

      if ( i1==point_positions_.end() or i2==point_positions_.end() )
      {
        throw std::runtime_error( "inner cut list end point not in cut point list" );
      }

      return *i1->second < *i2->second;
    }
  }

private:
  const std::map<Point*, std::vector<Point*>::const_iterator*> & point_positions_;
};


class PointCycle
{
public:

  explicit PointCycle( const std::vector<Point*> & facet_points )
    : facet_points_( facet_points )
  {
  }

  bool Split( const std::vector<Point*> & cut_points );

  bool Split( const InnerCutPoints & cut_points )
  {
    return Split( cut_points.Points() );
  }

  bool NeedsSplit( const std::vector<Point*> & cut_points );

  bool MyCut( const std::vector<Point*> & cut_points );

  void CreateFacets( Mesh & mesh, Side * side, std::vector<Facet*> & facets );

  const std::vector<Point*> & FacetPoints() const { return facet_points_; }

private:
  std::vector<Point*> facet_points_;

  Teuchos::RCP<PointCycle> lhs_;
  Teuchos::RCP<PointCycle> rhs_;
};


class PointCycleList
{
public:

  PointCycleList( const PointLineFilter & filter,
                  Side * side,
                  const std::vector<Point*> & facet_points,
                  const std::set<Point*, PointPidLess> & cut_points );

  void CreateFacets( Mesh & mesh, Side * side, std::vector<Facet*> & facets );

private:

  PointCycle pcl_;
};


  }
}

#endif
