#ifndef CUT_FACET_H
#define CUT_FACET_H

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <vector>
#include <stdexcept>

#include "../linalg/linalg_fixedsizematrix.H"

#include "cut_tetgen.H"
#include "cut_point.H"
#include "cut_line.H"

namespace GEO
{
namespace CUT
{
class LinearElement;
class Mesh;

class Facet
{
public:

  Facet( Mesh & mesh, const std::vector<Point*> & points, Side * side, bool cutsurface )
    : points_( points ),
      parentside_( side ),
      planar_( IsPlanar( mesh, points ) ),
      position_( cutsurface ? Point::oncutsurface : Point::undecided )
  {
    for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
    {
      Point * p = *i;
      if ( cutsurface )
      {
        p->Position( Point::oncutsurface );
      }
      p->Register( this );
    }
  }

  void Print()
  {
    for ( std::vector<Point*>::iterator i=points_.begin();
          i!=points_.end();
          ++i )
    {
      Point & p = **i;
      std::cout << p.Pid() << " ";
    }
    std::cout << "\n";
  }

  int SideId();

  void GetPoints( std::set<Point*, PointPidLess> & points );

#ifdef QHULL
  void GenerateTetgen( LinearElement * element, tetgenio::facet & f, int num, int & marker, const std::vector<Point*> & pointlist );
#endif

  void AddHole( const std::vector<Point*> & hole )
  {
    holes_.insert( hole );
  }

  int NumTetgenFacets()
  {
    if ( planar_ )
    {
      return 1;
    }
    else
    {
      if ( holes_.size()>0 )
        throw std::runtime_error( "cannot have holes in a non-planar facet" );
      return triangulation_.size();
    }
  }

  void ExchangeSide( Side * side, bool cutsurface )
  {
    parentside_ = side;
    if ( cutsurface )
    {
      Position( Point::oncutsurface );
      for ( std::vector<Point*>::const_iterator i=points_.begin(); i!=points_.end(); ++i )
      {
        Point * p = *i;
        p->Position( Point::oncutsurface );
      }
    }
  }

  bool Equals( const std::vector<Point*> & facet_points );

  bool IsCutSide( Side * side );

  Point::PointPosition Position() const { return position_; }

  void Position( Point::PointPosition p );

private:

  Facet( const Facet & );
  Facet & operator=( const Facet & );

#ifdef QHULL
  void GeneratePolygon( tetgenio::polygon & p,
                        const std::vector<Point*> & points,
                        const std::vector<Point*> & pointlist,
                        LINALG::Matrix<3,1> & mid );
#endif

  bool IsPlanar( Mesh & mesh, const std::vector<Point*> & points );

  std::vector<Point*> points_;

  std::set<std::vector<Point*> > holes_;

  std::vector<std::vector<Point*> > triangulation_;

  Side* parentside_;

  bool planar_;

  Point::PointPosition position_;
};


}
}

#endif
